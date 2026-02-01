import pandas as pd
import numpy as np
from neo4j import GraphDatabase
from rdkit import Chem
from rdkit.Chem import AllChem, Fragments
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
import warnings

# Tắt cảnh báo RDKit
warnings.filterwarnings('ignore')

# =============================================================================
# 1. CẤU HÌNH HỆ THỐNG
# =============================================================================
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = "12345678"

# --- MODIFIED: Support both experimental and de novo data ---
EXPERIMENTAL_CSV_PATH = "Data/data_end.csv"
DENOVO_CSV_PATH = "Data/DeNovo_Molecule.csv"

# =============================================================================
# 2. KIẾN THỨC CHUYÊN GIA (DOMAIN KNOWLEDGE) - UNCHANGED
# =============================================================================

# A. Danh sách thuốc chuẩn để suy luận Target (Target Inference)
REF_DRUGS = {
    'Gefitinib': {'smiles': 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1', 'target': 'EGFR_WT'},
    'Erlotinib': {'smiles': 'COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC', 'target': 'EGFR_WT'},
    'Osimertinib': {'smiles': 'C=CC(=O)Nc1cc(Nc2ncnc3cc(N(C)CCN(C)C)c(OC)cc23)c(OC)cc1', 'target': 'EGFR_T790M'},
    'Rociletinib': {'smiles': 'CN(C)CCN(C)c1cc(Nc2ncc(C(F)(F)F)c(Nc3cc(C(=O)N)ccc3)n2)cc(N)n1', 'target': 'EGFR_T790M'},
    'Afatinib': {'smiles': 'CN(C)C/C=C/C(=O)Nc1cc2c(Nc3ccc(F)c(Cl)c3)ncnc2cc1O[C@H]3CCOC3', 'target': 'EGFR_WT'}
}

# B. Danh sách Warhead (Vũ khí) để suy luận MoA
WARHEAD_SMARTS = {
    'Acrylamide': 'C=CC(=O)N',           # Đặc trưng thuốc thế hệ 3
    'Propynamide': 'C#CC(=O)N',
    'Chloroacetamide': 'ClCC(=O)N',
    'Vinyl_sulfonamide': 'C=CS(=O)(=O)N',
    'Epoxide': 'C1OC1',
    'Michael_Acceptor': 'C=CC(=O)',
}

# C. Danh sách Functional Prompts bổ sung (KANO Style)
EGFR_SPECIFIC_SMARTS = {
    'Quinazoline_Core': 'c1cc2ncnc(Nc3ccccc3)c2cc1', 
    'Pyrimidine_Core': 'c1cncnc1',
    'Aniline_Group': 'Nc1ccccc1',
    'Fluorine': '[F]',
    'Chlorine': '[Cl]',
    'Bromine': '[Br]',
    'Ether_Link': 'COC'
}

# =============================================================================
# 3. CÁC HÀM XỬ LÝ HÓA HỌC (RDKIT LOGIC) - UNCHANGED
# =============================================================================

# --- Tạo thư viện Fingerprint cho thuốc chuẩn ---
ref_fps = []
for name, info in REF_DRUGS.items():
    m = Chem.MolFromSmiles(info['smiles'])
    if m:
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
        ref_fps.append({'name': name, 'fp': fp, 'target': info['target']})

def predict_target(mol):
    """Suy luận Target dựa trên độ tương đồng với thuốc chuẩn"""
    if not mol: return "EGFR_Generic"
    
    target_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    max_sim = 0
    best_target = "EGFR_Generic"
    
    for ref in ref_fps:
        sim = DataStructs.TanimotoSimilarity(target_fp, ref['fp'])
        if sim > max_sim:
            max_sim = sim
            best_target = ref['target']
            
    return best_target if max_sim > 0.35 else "EGFR_Generic"

def get_functional_prompts(mol):
    """Trích xuất Functional Prompts (Ý tưởng cốt lõi của bài báo KANO)"""
    if not mol: return []
    prompts = []
    
    # 1. Dùng RDKit Fragments
    fr_funcs = [m for m in dir(Fragments) if m.startswith('fr_')]
    for func_name in fr_funcs:
        func = getattr(Fragments, func_name)
        try:
            if func(mol) > 0:
                clean_name = func_name.replace("fr_", "")
                prompts.append(clean_name)
        except: pass
            
    # 2. Dùng SMARTS thủ công cho EGFR
    for name, smarts in EGFR_SPECIFIC_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            prompts.append(name)
            
    return list(set(prompts))

def identify_warhead_and_moa(mol):
    """Xác định Warhead và Cơ chế tác động"""
    if not mol: return [], "Unknown"
    
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

def get_scaffold(mol):
    """Lấy khung sườn Murcko"""
    if not mol: return None
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core) if core else Chem.MolToSmiles(mol)
    except: return Chem.MolToSmiles(mol)

# =============================================================================
# 4. CLASS XÂY DỰNG KG (NEO4J IMPORTER)
# =============================================================================
class KANO_KG_Builder:
    def __init__(self, uri, user, pwd):
        self.driver = GraphDatabase.driver(uri, auth=(user, pwd))
        
    def close(self):
        self.driver.close()
        
    def nuke_and_prepare_db(self):
        """Xóa sạch DB và tạo Index mới"""
        with self.driver.session() as session:
            session.run("MATCH (n) DETACH DELETE n")
            print("💥 Đã xóa sạch dữ liệu cũ.")
            
            constraints = session.run("SHOW CONSTRAINTS YIELD name").data()
            for c in constraints:
                try: session.run(f"DROP CONSTRAINT {c['name']}")
                except: pass
                
            indexes = session.run("SHOW INDEXES YIELD name, type WHERE type <> 'LOOKUP'").data()
            for i in indexes:
                try: session.run(f"DROP INDEX {i['name']}")
                except: pass
                
            print("🧹 Đã dọn sạch Schema cũ.")
            
            session.run("CREATE CONSTRAINT FOR (m:Molecule) REQUIRE m.smiles IS UNIQUE")
            session.run("CREATE INDEX FOR (s:Scaffold) ON (s.smiles)")
            session.run("CREATE INDEX FOR (fp:FunctionalGroup) ON (fp.name)")
            print("✅ Đã tạo Index mới.")

    def import_batch(self, batch):
        # --- MODIFIED: Support both experimental and virtual molecules ---
        query = """
        UNWIND $batch AS row
        
        // 1. Tạo Molecule với phân biệt experimental vs virtual
        MERGE (m:Molecule {smiles: row.smiles})
        SET m.is_virtual = row.is_virtual,
            m.source = row.source
        
        // --- DIVERGENCE POINT: Experimental vs Virtual ---
        // Experimental molecules: set activity, ic50
        FOREACH (_ IN CASE WHEN NOT row.is_virtual THEN [1] ELSE [] END |
            SET m.activity = row.activity,
                m.ic50 = row.ic50
        )
        
        // Virtual molecules: set docking_affinity, ligand_id
        FOREACH (_ IN CASE WHEN row.is_virtual THEN [1] ELSE [] END |
            SET m.docking_affinity = row.docking_affinity,
                m.ligand_id = row.ligand_id
        )
        
        // 2. Tạo Scaffold (Khung sườn) - SAME FOR BOTH
        MERGE (s:Scaffold {smiles: row.scaffold})
        MERGE (m)-[:HAS_SCAFFOLD]->(s)
        
        // 3. Tạo Functional Prompts (KANO) - SAME FOR BOTH
        FOREACH (fp_name IN row.functional_prompts |
            MERGE (fp:FunctionalGroup {name: fp_name})
            MERGE (m)-[:HAS_FUNCTIONAL_GROUP]->(fp)
        )
        
        // 4. Tạo Warhead (Vũ khí) - SAME FOR BOTH
        FOREACH (w_name IN row.warheads |
            MERGE (w:Warhead {name: w_name})
            MERGE (m)-[:CONTAINS_WARHEAD]->(w)
        )
        
        // 5. Tạo MoA (Cơ chế) - SAME FOR BOTH
        MERGE (moa:MoA {name: row.moa})
        MERGE (m)-[:ACTS_VIA]->(moa)
        
        // 6. Tạo Target (Mục tiêu - Đã suy luận) - SAME FOR BOTH
        MERGE (t:Target {name: row.target})
        MERGE (m)-[:TESTED_AGAINST]->(t)
        
        // 7. Tạo POTENT_AGAINST ONLY for Active Experimental Molecules
        // Virtual molecules NEVER get this relationship
        FOREACH (_ IN CASE 
            WHEN NOT row.is_virtual AND row.activity = 1 
            THEN [1] 
            ELSE [] 
        END |
            MERGE (m)-[:POTENT_AGAINST]->(t)
        )
        """
        with self.driver.session() as session:
            session.run(query, batch=batch)

# =============================================================================
# 5. MAIN PROGRAM - MODIFIED TO HANDLE BOTH DATA TYPES
# =============================================================================
def process_experimental_molecules(df, builder):
    """Xử lý phân tử experimental (từ data_end.csv)"""
    print(f"📊 Đang xử lý {len(df)} phân tử EXPERIMENTAL...")
    
    batch_size = 1000
    batch_data = []
    
    for idx, row in df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if not mol: continue
        
        # --- CHEMISTRY LOGIC (SAME) ---
        scaffold = get_scaffold(mol)
        warheads, moa = identify_warhead_and_moa(mol)
        functional_prompts = get_functional_prompts(mol)
        target = predict_target(mol)
        
        # --- EXPERIMENTAL-SPECIFIC PROPERTIES ---
        label = str(row['Final_Label']).lower()
        activity = 1 if label == 'active' else 0
        ic50 = float(row.get('IC50 value(nM)', 0))
        
        # Đóng gói
        item = {
            'smiles': smiles,
            'is_virtual': False,  # <- Experimental molecule
            'source': 'Experimental',
            'activity': activity,
            'ic50': ic50,
            'docking_affinity': None,  # Not applicable
            'ligand_id': None,  # Not applicable
            'scaffold': scaffold,
            'warheads': warheads,
            'moa': moa,
            'functional_prompts': functional_prompts,
            'target': target
        }
        
        batch_data.append(item)
        
        if len(batch_data) >= batch_size:
            builder.import_batch(batch_data)
            print(f"   ✅ Nạp experimental: {idx + 1}/{len(df)}")
            batch_data = []
            
    # Nạp nốt phần còn dư
    if batch_data:
        builder.import_batch(batch_data)
    
    print(f"✅ Hoàn tất {len(df)} phân tử experimental!")

def process_denovo_molecules(df, builder):
    """Xử lý phân tử de novo (từ DeNovo_Molecule.csv)"""
    print(f"🧪 Đang xử lý {len(df)} phân tử DE NOVO (Virtual)...")
    
    batch_size = 1000
    batch_data = []
    
    for idx, row in df.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)
        
        if not mol: continue
        
        # --- CHEMISTRY LOGIC (SAME AS EXPERIMENTAL) ---
        scaffold = get_scaffold(mol)
        warheads, moa = identify_warhead_and_moa(mol)
        functional_prompts = get_functional_prompts(mol)
        target = predict_target(mol)
        
        # --- VIRTUAL-SPECIFIC PROPERTIES ---
        # NO activity label, NO ic50
        ligand_id = int(row['ligand_id'])
        docking_affinity = float(row['affinity'])
        
        # Đóng gói
        item = {
            'smiles': smiles,
            'is_virtual': True,  # <- Virtual molecule
            'source': 'DiffSBDD',
            'activity': 0,  # Placeholder (not used)
            'ic50': None,  # Not applicable
            'docking_affinity': docking_affinity,
            'ligand_id': ligand_id,
            'scaffold': scaffold,
            'warheads': warheads,
            'moa': moa,
            'functional_prompts': functional_prompts,
            'target': target
        }
        
        batch_data.append(item)
        
        if len(batch_data) >= batch_size:
            builder.import_batch(batch_data)
            print(f"   ✅ Nạp de novo: {idx + 1}/{len(df)}")
            batch_data = []
            
    # Nạp nốt phần còn dư
    if batch_data:
        builder.import_batch(batch_data)
    
    print(f"✅ Hoàn tất {len(df)} phân tử de novo!")

def main():
    # Khởi tạo DB
    builder = KANO_KG_Builder(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
    builder.nuke_and_prepare_db()
    
    # --- PROCESS EXPERIMENTAL MOLECULES ---
    print("\n" + "="*80)
    print("BƯỚC 1: NẠP DỮ LIỆU EXPERIMENTAL")
    print("="*80)
    try:
        df_experimental = pd.read_csv(EXPERIMENTAL_CSV_PATH)
        process_experimental_molecules(df_experimental, builder)
    except FileNotFoundError:
        print(f"⚠️  Không tìm thấy {EXPERIMENTAL_CSV_PATH}, bỏ qua experimental data.")
    
    # --- PROCESS DE NOVO MOLECULES ---
    print("\n" + "="*80)
    print("BƯỚC 2: NẠP DỮ LIỆU DE NOVO (VIRTUAL)")
    print("="*80)
    try:
        df_denovo = pd.read_csv(DENOVO_CSV_PATH)
        process_denovo_molecules(df_denovo, builder)
    except FileNotFoundError:
        print(f"⚠️  Không tìm thấy {DENOVO_CSV_PATH}, bỏ qua de novo data.")
    
    builder.close()
    
    print("\n" + "="*80)
    print("🏆 HOÀN TẤT! Knowledge Graph đã được xây dựng với:")
    print("   ✅ Experimental molecules (với IC50, activity)")
    print("   ✅ De novo virtual molecules (với docking affinity)")
    print("   ✅ POTENT_AGAINST chỉ áp dụng cho experimental actives")
    print("="*80)

if __name__ == "__main__":
    main()                 