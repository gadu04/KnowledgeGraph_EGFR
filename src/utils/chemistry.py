"""Chemistry utilities using RDKit"""
from typing import List, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Fragments, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
import numpy as np

# Domain Knowledge
REF_DRUGS = {
    'Gefitinib': {'smiles': 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1', 'target': 'EGFR_WT'},
    'Erlotinib': {'smiles': 'COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC', 'target': 'EGFR_WT'},
    'Osimertinib': {'smiles': 'C=CC(=O)Nc1cc(Nc2ncnc3cc(N(C)CCN(C)C)c(OC)cc23)c(OC)cc1', 'target': 'EGFR_T790M'},
    'Rociletinib': {'smiles': 'CN(C)CCN(C)c1cc(Nc2ncc(C(F)(F)F)c(Nc3cc(C(=O)N)ccc3)n2)cc(N)n1', 'target': 'EGFR_T790M'},
    'Afatinib': {'smiles': 'CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1O[C@H]3CCOC3', 'target': 'EGFR_WT'}
}

WARHEAD_SMARTS = {
    'Acrylamide': 'C=CC(=O)N',
    'Propynamide': 'C#CC(=O)N',
    'Chloroacetamide': 'ClCC(=O)N',
    'Vinyl_sulfonamide': 'C=CS(=O)(=O)N',
    'Epoxide': 'C1OC1',
    'Michael_Acceptor': 'C=CC(=O)',
}

EGFR_SPECIFIC_SMARTS = {
    'Quinazoline_Core': 'c1cc2ncnc(Nc3ccccc3)c2cc1',
    'Pyrimidine_Core': 'c1cncnc1',
    'Aniline_Group': 'Nc1ccccc1',
    'Fluorine': '[F]',
    'Chlorine': '[Cl]',
    'Bromine': '[Br]',
    'Ether_Link': 'COC'
}

# Initialize reference fingerprints
ref_fps = []
for name, info in REF_DRUGS.items():
    mol = Chem.MolFromSmiles(info['smiles'])
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        ref_fps.append({'name': name, 'fp': fp, 'target': info['target']})


def canonicalize_smiles(smiles: str) -> Optional[str]:
    """Canonicalize SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except:
        pass
    return None


def get_scaffold(mol: Chem.Mol) -> Optional[str]:
    """Get Murcko scaffold"""
    if not mol:
        return None
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core) if core else Chem.MolToSmiles(mol)
    except:
        return Chem.MolToSmiles(mol)


def predict_target(mol: Chem.Mol, threshold: float = 0.35) -> str:
    """Predict target based on similarity to reference drugs"""
    if not mol:
        return "EGFR_Generic"
    
    target_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    max_sim = 0
    best_target = "EGFR_Generic"
    
    for ref in ref_fps:
        sim = DataStructs.TanimotoSimilarity(target_fp, ref['fp'])
        if sim > max_sim:
            max_sim = sim
            best_target = ref['target']
    
    return best_target if max_sim > threshold else "EGFR_Generic"


def get_functional_prompts(mol: Chem.Mol) -> List[str]:
    """Extract functional prompts (KANO style)"""
    if not mol:
        return []
    
    prompts = []
    
    # RDKit Fragments
    fr_funcs = [m for m in dir(Fragments) if m.startswith('fr_')]
    for func_name in fr_funcs:
        func = getattr(Fragments, func_name)
        try:
            if func(mol) > 0:
                clean_name = func_name.replace("fr_", "")
                prompts.append(clean_name)
        except:
            pass
    
    # EGFR-specific SMARTS
    for name, smarts in EGFR_SPECIFIC_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            prompts.append(name)
    
    return list(set(prompts))


def identify_warhead_and_moa(mol: Chem.Mol) -> Tuple[List[str], str]:
    """Identify warheads and mechanism of action"""
    if not mol:
        return [], "Unknown"
    
    found_warheads = []
    for name, smarts in WARHEAD_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_warheads.append(name)
    
    covalent_markers = ['Acrylamide', 'Propynamide', 'Vinyl_sulfonamide', 'Epoxide']
    if any(w in covalent_markers for w in found_warheads):
        moa = "Covalent_Inhibitor"
    else:
        moa = "Reversible_Inhibitor"
    
    return found_warheads, moa


def get_ecfp4(smiles_list: List[str], n_bits: int = 1024) -> np.ndarray:
    """Generate ECFP4 fingerprints"""
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
            fps.append(np.array(fp))
        else:
            fps.append(np.zeros(n_bits))
    return np.array(fps)
