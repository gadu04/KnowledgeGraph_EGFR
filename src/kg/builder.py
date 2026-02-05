"""Knowledge Graph Builder for Neo4j"""
from typing import List, Dict
import pandas as pd
from neo4j import GraphDatabase
from rdkit import Chem

from src.config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD, BATCH_SIZE
from src.utils.chemistry import (
    get_scaffold, predict_target, get_functional_prompts,
    identify_warhead_and_moa
)


class KnowledgeGraphBuilder:
    """KANO-style Knowledge Graph Builder"""
    
    def __init__(self, uri: str = None, user: str = None, password: str = None):
        self.uri = uri or NEO4J_URI
        self.user = user or NEO4J_USER
        self.password = password or NEO4J_PASSWORD
        self.driver = GraphDatabase.driver(self.uri, auth=(self.user, self.password))
    
    def close(self):
        """Close Neo4j connection"""
        self.driver.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def nuke_and_prepare_db(self):
        """Clear database and create indexes"""
        with self.driver.session() as session:
            # Delete all nodes
            session.run("MATCH (n) DETACH DELETE n")
            print("💥 Cleared old data")
            
            # Drop constraints
            constraints = session.run("SHOW CONSTRAINTS YIELD name").data()
            for c in constraints:
                try:
                    session.run(f"DROP CONSTRAINT {c['name']}")
                except:
                    pass
            
            # Drop indexes
            indexes = session.run("SHOW INDEXES YIELD name, type WHERE type <> 'LOOKUP'").data()
            for i in indexes:
                try:
                    session.run(f"DROP INDEX {i['name']}")
                except:
                    pass
            
            print("🧹 Cleaned old schema")
            
            # Create new indexes
            session.run("CREATE CONSTRAINT FOR (m:Molecule) REQUIRE m.smiles IS UNIQUE")
            session.run("CREATE INDEX FOR (s:Scaffold) ON (s.smiles)")
            session.run("CREATE INDEX FOR (t:Target) ON (t.name)")
            session.run("CREATE INDEX FOR (w:Warhead) ON (w.name)")
            session.run("CREATE INDEX FOR (moa:MoA) ON (moa.name)")
            session.run("CREATE INDEX FOR (fp:FunctionalGroup) ON (fp.name)")
            print("✅ Created new indexes")
    
    def import_batch(self, batch: List[Dict]):
        """Import batch of molecules into Neo4j"""
        query = """
        UNWIND $batch AS row
        
        // 1. Create Molecule (experimental or virtual)
        MERGE (m:Molecule {smiles: row.smiles})
        SET m.is_virtual = row.is_virtual,
            m.source = row.source
        
        // Experimental properties
        FOREACH (_ IN CASE WHEN NOT row.is_virtual THEN [1] ELSE [] END |
            SET m.activity = row.activity,
                m.ic50 = row.ic50
        )
        
        // Virtual properties
        FOREACH (_ IN CASE WHEN row.is_virtual THEN [1] ELSE [] END |
            SET m.docking_affinity = row.docking_affinity,
                m.ligand_id = row.ligand_id
        )
        
        // 2. Scaffold
        MERGE (s:Scaffold {smiles: row.scaffold})
        MERGE (m)-[:HAS_SCAFFOLD]->(s)
        
        // 3. Functional Groups
        FOREACH (fp_name IN row.functional_prompts |
            MERGE (fp:FunctionalGroup {name: fp_name})
            MERGE (m)-[:HAS_FUNCTIONAL_GROUP]->(fp)
        )
        
        // 4. Warheads
        FOREACH (w_name IN row.warheads |
            MERGE (w:Warhead {name: w_name})
            MERGE (m)-[:CONTAINS_WARHEAD]->(w)
        )
        
        // 5. MoA
        MERGE (moa:MoA {name: row.moa})
        MERGE (m)-[:ACTS_VIA]->(moa)
        
        // 6. Target
        MERGE (t:Target {name: row.target})
        MERGE (m)-[:TESTED_AGAINST]->(t)
        
        // ❌ XÓA PHẦN NÀY:
        // FOREACH (_ IN CASE 
        //     WHEN NOT row.is_virtual AND row.activity = 1 
        //     THEN [1] 
        //     ELSE [] 
        // END |
        //     MERGE (m)-[:POTENT_AGAINST]->(t)
        // )
        """
        with self.driver.session() as session:
            session.run(query, batch=batch)
    
    def process_experimental_molecules(self, df: pd.DataFrame):
        """Process experimental molecules from data_end.csv"""
        print(f"📊 Processing {len(df)} EXPERIMENTAL molecules...")
        
        batch_data = []
        for idx, row in df.iterrows():
            smiles = row['SMILES']
            mol = Chem.MolFromSmiles(smiles)
            
            if not mol:
                continue
            
            # Chemistry analysis
            scaffold = get_scaffold(mol)
            warheads, moa = identify_warhead_and_moa(mol)
            functional_prompts = get_functional_prompts(mol)
            target = predict_target(mol)
            
            # Experimental properties
            label = str(row['Final_Label']).lower()
            activity = 1 if label == 'active' else 0
            ic50 = float(row.get('IC50 value(nM)', 0))
            
            item = {
                'smiles': smiles,
                'is_virtual': False,
                'source': 'Experimental',
                'activity': activity,
                'ic50': ic50,
                'docking_affinity': None,
                'ligand_id': None,
                'scaffold': scaffold,
                'warheads': warheads,
                'moa': moa,
                'functional_prompts': functional_prompts,
                'target': target
            }
            
            batch_data.append(item)
            
            if len(batch_data) >= BATCH_SIZE:
                self.import_batch(batch_data)
                print(f"   ✅ Imported experimental: {idx + 1}/{len(df)}")
                batch_data = []
        
        # Import remaining
        if batch_data:
            self.import_batch(batch_data)
        
        print(f"✅ Completed {len(df)} experimental molecules!")
    
    def process_denovo_molecules(self, df: pd.DataFrame):
        """Process de novo molecules from DeNovo_Molecule.csv"""
        print(f"🧪 Processing {len(df)} DE NOVO (Virtual) molecules...")
        
        batch_data = []
        for idx, row in df.iterrows():
            smiles = row['smiles']
            mol = Chem.MolFromSmiles(smiles)
            
            if not mol:
                continue
            
            # Chemistry analysis
            scaffold = get_scaffold(mol)
            warheads, moa = identify_warhead_and_moa(mol)
            functional_prompts = get_functional_prompts(mol)
            target = predict_target(mol)
            
            # Virtual properties
            ligand_id = int(row['ligand_id'])
            docking_affinity = float(row['affinity'])
            
            item = {
                'smiles': smiles,
                'is_virtual': True,
                'source': 'DiffSBDD',
                'activity': 0,
                'ic50': None,
                'docking_affinity': docking_affinity,
                'ligand_id': ligand_id,
                'scaffold': scaffold,
                'warheads': warheads,
                'moa': moa,
                'functional_prompts': functional_prompts,
                'target': target
            }
            
            batch_data.append(item)
            
            if len(batch_data) >= BATCH_SIZE:
                self.import_batch(batch_data)
                print(f"   ✅ Imported de novo: {idx + 1}/{len(df)}")
                batch_data = []
        
        # Import remaining
        if batch_data:
            self.import_batch(batch_data)
        
        print(f"✅ Completed {len(df)} de novo molecules!")
