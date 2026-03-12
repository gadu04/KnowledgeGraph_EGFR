"""Knowledge Graph Builder for Neo4j.

Pipeline: KG = f(SMILES_Data, Domain_Config)

Pass a *config_path* to switch targets (EGFR → HDAC, etc.) without
touching any Python source.  All domain-knowledge logic is delegated to
ChemistryAnalyzer which reads from the JSON config.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from neo4j import GraphDatabase
from rdkit import Chem

from src.config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, BATCH_SIZE
from src.utils.chemistry import ChemistryAnalyzer, get_scaffold

# Default config bundled with the package
_DEFAULT_CONFIG = Path(__file__).parent.parent / "config" / "domain_config.json"


class KnowledgeGraphBuilder:
    """Config-driven KANO-style Knowledge Graph Builder.

    Parameters
    ----------
    config_path:
        Path to a ``domain_config.json`` file.  Falls back to the bundled
        EGFR config when *None*.
    uri, user, password:
        Neo4j connection parameters; fall back to environment/settings values.
    """

    def __init__(
        self,
        config_path: Optional[str | Path] = None,
        uri: Optional[str] = None,
        user: Optional[str] = None,
        password: Optional[str] = None,
    ) -> None:
        self.uri = uri or NEO4J_URI
        self.user = user or NEO4J_USER
        self.password = password or NEO4J_PASSWORD
        self.driver = GraphDatabase.driver(self.uri, auth=(self.user, self.password))

        resolved = Path(config_path) if config_path else _DEFAULT_CONFIG
        self.analyzer = ChemistryAnalyzer.from_json(resolved)
        print(f"✅ Loaded domain config: {resolved.name}")

    def close(self) -> None:
        self.driver.close()

    def __enter__(self) -> "KnowledgeGraphBuilder":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    # ------------------------------------------------------------------
    # Schema management
    # ------------------------------------------------------------------

    def nuke_and_prepare_db(self) -> None:
        """Clear all data and rebuild indexes/constraints."""
        with self.driver.session() as session:
            session.run("MATCH (n) DETACH DELETE n")
            print("💥 Cleared old data")

            for c in session.run("SHOW CONSTRAINTS YIELD name").data():
                try:
                    session.run(f"DROP CONSTRAINT {c['name']}")
                except Exception:
                    pass

            for i in session.run(
                "SHOW INDEXES YIELD name, type WHERE type <> 'LOOKUP'"
            ).data():
                try:
                    session.run(f"DROP INDEX {i['name']}")
                except Exception:
                    pass

            print("🧹 Cleaned old schema")

            session.run("CREATE CONSTRAINT FOR (m:Molecule) REQUIRE m.smiles IS UNIQUE")
            session.run("CREATE INDEX FOR (s:Scaffold) ON (s.smiles)")
            session.run("CREATE INDEX FOR (t:Target) ON (t.name)")
            session.run("CREATE INDEX FOR (ig:Interaction_Group) ON (ig.name)")
            session.run("CREATE INDEX FOR (moa:MoA) ON (moa.name)")
            session.run("CREATE INDEX FOR (fp:FunctionalGroup) ON (fp.name)")
            print("✅ Created new indexes")

    # ------------------------------------------------------------------
    # Batch import
    # ------------------------------------------------------------------

    def import_batch(self, batch: List[Dict]) -> None:
        """Write one batch of pre-processed molecule records to Neo4j.

        Each record must contain:
          smiles, is_virtual, source, scaffold,
          interaction_groups (list of {name, moa}),
          functional_prompts (list of str),
          target (str),
          docking_affinity / ligand_id  (virtual only)
        """
        query = """
        UNWIND $batch AS row

        // 1. Molecule node
        MERGE (m:Molecule {smiles: row.smiles})
        SET m.is_virtual = row.is_virtual,
            m.source     = row.source

        FOREACH (_ IN CASE WHEN row.is_virtual THEN [1] ELSE [] END |
            SET m.docking_affinity = row.docking_affinity,
                m.ligand_id        = row.ligand_id
        )

        // 2. Scaffold
        MERGE (s:Scaffold {smiles: row.scaffold})
        MERGE (m)-[:HAS_SCAFFOLD]->(s)

        // 3. Functional groups
        FOREACH (fp_name IN row.functional_prompts |
            MERGE (fp:FunctionalGroup {name: fp_name})
            MERGE (m)-[:HAS_FUNCTIONAL_GROUP]->(fp)
        )

        // 4. Interaction groups  →  each carries its own MoA
        FOREACH (ig IN row.interaction_groups |
            MERGE (igNode:Interaction_Group {name: ig.name})
            MERGE (moa:MoA {name: ig.moa})
            MERGE (igNode)-[:ACTS_VIA]->(moa)
            MERGE (m)-[:HAS_INTERACTION_GROUP]->(igNode)
        )

        // 5. Target
        MERGE (t:Target {name: row.target})
        MERGE (m)-[:TESTED_AGAINST]->(t)
        """
        with self.driver.session() as session:
            session.run(query, batch=batch)

    # ------------------------------------------------------------------
    # Data processing
    # ------------------------------------------------------------------

    def process_experimental_molecules(self, df: pd.DataFrame) -> None:
        """Build KG records from experimental data (data_end.csv)."""
        print(f"📊 Processing {len(df)} EXPERIMENTAL molecules...")

        batch_data: List[Dict] = []
        for idx, row in df.iterrows():
            smiles = row["SMILES"]
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            batch_data.append({
                "smiles": smiles,
                "is_virtual": False,
                "source": "Experimental",
                "docking_affinity": None,
                "ligand_id": None,
                "scaffold": get_scaffold(mol),
                "interaction_groups": self.analyzer.get_interaction_groups(mol),
                "functional_prompts": self.analyzer.get_functional_prompts(mol),
                "target": self.analyzer.assign_target(smiles),
            })

            if len(batch_data) >= BATCH_SIZE:
                self.import_batch(batch_data)
                print(f"   ✅ Imported experimental: {idx + 1}/{len(df)}")
                batch_data = []

        if batch_data:
            self.import_batch(batch_data)

        print(f"✅ Completed {len(df)} experimental molecules!")

    def process_denovo_molecules(self, df: pd.DataFrame) -> None:
        """Build KG records from de-novo virtual molecules (DeNovo_Molecule.csv)."""
        print(f"🧪 Processing {len(df)} DE NOVO (Virtual) molecules...")

        batch_data: List[Dict] = []
        for idx, row in df.iterrows():
            smiles = row["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            batch_data.append({
                "smiles": smiles,
                "is_virtual": True,
                "source": "DiffSBDD",
                "docking_affinity": float(row["affinity"]),
                "ligand_id": int(row["ligand_id"]),
                "scaffold": get_scaffold(mol),
                "interaction_groups": self.analyzer.get_interaction_groups(mol),
                "functional_prompts": self.analyzer.get_functional_prompts(mol),
                "target": self.analyzer.assign_target(smiles),
            })

            if len(batch_data) >= BATCH_SIZE:
                self.import_batch(batch_data)
                print(f"   ✅ Imported de novo: {idx + 1}/{len(df)}")
                batch_data = []

        if batch_data:
            self.import_batch(batch_data)

        print(f"✅ Completed {len(df)} de novo molecules!")
