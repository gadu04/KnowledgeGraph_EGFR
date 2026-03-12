"""Chemistry utilities using RDKit.

Module-level helpers (canonicalize_smiles, get_scaffold, get_ecfp4) remain
stateless.  All domain-knowledge-dependent logic lives in ChemistryAnalyzer,
which is initialised from a domain_config dict loaded by KnowledgeGraphBuilder.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs, Fragments
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

# Shared fingerprint generator (stateless, reusable across all instances)
_morgan_gen = GetMorganGenerator(radius=2, fpSize=1024)


# ---------------------------------------------------------------------------
# Stateless helpers
# ---------------------------------------------------------------------------

def canonicalize_smiles(smiles: str) -> Optional[str]:
    """Return canonical SMILES, or None if the string is invalid."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        pass
    return None


def get_scaffold(mol: Chem.Mol) -> Optional[str]:
    """Return the Murcko scaffold SMILES for *mol*."""
    if not mol:
        return None
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core) if core else Chem.MolToSmiles(mol)
    except Exception:
        return Chem.MolToSmiles(mol)


def get_ecfp4(smiles_list: List[str], n_bits: int = 1024) -> np.ndarray:
    """Return an (N, n_bits) array of ECFP4 fingerprints."""
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fps.append(np.array(_morgan_gen.GetFingerprint(mol)))
        else:
            fps.append(np.zeros(n_bits))
    return np.array(fps)


# ---------------------------------------------------------------------------
# Config-driven analyser
# ---------------------------------------------------------------------------

class ChemistryAnalyzer:
    """Domain-knowledge-aware molecule analyser.

    All domain constants (reference drugs, interaction-group SMARTS,
    target-specific rules) are read from a *domain_config* dict so that
    switching targets (EGFR → HDAC, etc.) requires only a config swap.

    Parameters
    ----------
    config:
        Parsed domain config dict, e.g. loaded from ``domain_config.json``.
    """

    def __init__(self, config: Dict) -> None:
        target_info = config["Target_Info"]
        self._generic_label: str = target_info["generic_label"]
        self._similarity_threshold: float = float(target_info["similarity_threshold"])

        # Pre-compute reference fingerprints
        self._ref_fps: List[Dict] = []
        for name, info in config["Reference_Drugs"].items():
            mol = Chem.MolFromSmiles(info["smiles"])
            if mol:
                self._ref_fps.append({
                    "name": name,
                    "fp": _morgan_gen.GetFingerprint(mol),
                    "target": info["target"],
                })

        # Pre-compile Interaction_Group patterns
        self._interaction_groups: List[Dict] = []
        for entry in config["Interaction_Groups"]:
            pattern = Chem.MolFromSmarts(entry["smarts"])
            if pattern:
                self._interaction_groups.append({
                    "name": entry["name"],
                    "pattern": pattern,
                    "moa": entry["moa"],
                })

        # Pre-compile Target_Specific_Rules patterns
        self._target_rules: List[Tuple[str, Chem.Mol]] = []
        for name, smarts in config["Target_Specific_Rules"].items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                self._target_rules.append((name, pattern))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def assign_target(self, smiles: str) -> str:
        """Infer the target label for *smiles* by Tanimoto similarity.

        Returns the target of the most-similar reference drug when the
        similarity exceeds *similarity_threshold*, otherwise the generic label.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return self._generic_label

        query_fp = _morgan_gen.GetFingerprint(mol)
        best_target = self._generic_label
        max_sim = 0.0

        for ref in self._ref_fps:
            sim = DataStructs.TanimotoSimilarity(query_fp, ref["fp"])
            if sim > max_sim:
                max_sim = sim
                best_target = ref["target"]

        return best_target if max_sim >= self._similarity_threshold else self._generic_label

    def get_interaction_groups(self, mol: Chem.Mol) -> List[Dict[str, str]]:
        """Return matched interaction groups for *mol*.

        Each element is ``{"name": <group_name>, "moa": <moa_label>}``.
        If no group matches, returns a single entry with moa
        ``"Reversible_Inhibitor"`` (safe default).
        """
        if not mol:
            return [{"name": "Unknown", "moa": "Unknown"}]

        matched = [
            {"name": ig["name"], "moa": ig["moa"]}
            for ig in self._interaction_groups
            if mol.HasSubstructMatch(ig["pattern"])
        ]
        return matched if matched else [{"name": "Non_Covalent", "moa": "Reversible_Inhibitor"}]

    def get_functional_prompts(self, mol: Chem.Mol) -> List[str]:
        """Extract functional-group labels for *mol* (KANO-style).

        Step 1 – RDKit built-in fragment counters (fr_*).
        Step 2 – Target-specific SMARTS from config.
        """
        if not mol:
            return []

        prompts: List[str] = []

        # Step 1: RDKit Fragments
        for func_name in (f for f in dir(Fragments) if f.startswith("fr_")):
            func = getattr(Fragments, func_name)
            try:
                if func(mol) > 0:
                    prompts.append(func_name.replace("fr_", ""))
            except Exception:
                pass

        # Step 2: target-specific rules from config
        for name, pattern in self._target_rules:
            if mol.HasSubstructMatch(pattern):
                prompts.append(name)

        return list(set(prompts))

    # ------------------------------------------------------------------
    # Convenience factory
    # ------------------------------------------------------------------

    @classmethod
    def from_json(cls, config_path: str | Path) -> "ChemistryAnalyzer":
        """Load a *domain_config.json* file and return a ready analyser."""
        with open(config_path, "r", encoding="utf-8") as fh:
            config = json.load(fh)
        return cls(config)
