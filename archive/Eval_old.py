import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, roc_curve

# ================= CẤU HÌNH =================
DATA_PATH = 'Data/final_clean_dataset.csv'  # File dữ liệu sạch
RANDOM_SEED = 42

# ================= HÀM TẠO FINGERPRINT (Giữ nguyên) =================
def generate_ecfp4(smiles, n_bits=1024):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
            return np.array(fp)
    except:
        return np.zeros(n_bits)
    return np.zeros(n_bits)

def generate_maccs(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            fp = MACCSkeys.GenMACCSKeys(mol)
            return np.array(fp)
    except:
        return np.zeros(167)
    return np.zeros(167)

# ================= HÀM BENCHMARK MLP =================
def run_mlp_benchmark(X, y, name="Method"):
    print(f"\n--- Đang training MLP cho: {name} ---")
    
    # Chia tập train/test (80/20)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=RANDOM_SEED, stratify=y
    )
    
    # Cấu hình MLP (Mạng Neural)
    clf = MLPClassifier(
        hidden_layer_sizes=(512, 128, 64),  # 3 lớp ẩn: to -> nhỏ dần
        activation='relu',                  # Hàm kích hoạt chuẩn
        solver='adam',                      # Bộ tối ưu hóa
        alpha=0.0001,                       # L2 Regularization (tránh overfitting)
        batch_size=128,                     # Kích thước batch
        learning_rate_init=0.001,           # Tốc độ học
        max_iter=200,                       # Số vòng lặp tối đa
        early_stopping=True,                # Dừng sớm nếu không học thêm được (QUAN TRỌNG)
        validation_fraction=0.1,            # Dành 10% tập train để validate
        n_iter_no_change=10,                # Nếu 10 vòng ko tiến bộ thì dừng
        random_state=RANDOM_SEED,
        verbose=True                        # In quá trình training
    )
    
    clf.fit(X_train, y_train)
    
    # Dự đoán
    y_pred = clf.predict(X_test)
    y_prob = clf.predict_proba(X_test)[:, 1]
    
    # Tính Metrics
    acc = accuracy_score(y_test, y_pred)
    auc_score = roc_auc_score(y_test, y_prob)
    f1 = f1_score(y_test, y_pred)
    
    print(f"✅ Kết quả {name}:")
    print(f"   - Accuracy: {acc:.4f}")
    print(f"   - AUC-ROC:  {auc_score:.4f}")
    print(f"   - F1-Score: {f1:.4f}")
    
    return y_test, y_prob, auc_score

# ================= MAIN =================
def main():
    # 1. Load dữ liệu
    df = pd.read_csv(DATA_PATH)
    
    # Xử lý nhãn
    if df['Final_Label'].dtype == object:
        df['Label_Num'] = df['Final_Label'].apply(lambda x: 1 if str(x).lower() == 'active' else 0)
    else:
        df['Label_Num'] = df['Final_Label']
    y = df['Label_Num'].values

    # 2. Tạo Feature
    print("Đang tạo Features...")
    X_ecfp4 = np.array([generate_ecfp4(s) for s in df['SMILES']])
    X_maccs = np.array([generate_maccs(s) for s in df['SMILES']])
    
    # 3. Chạy Benchmark với MLP
    y_test_e, y_prob_e, auc_e = run_mlp_benchmark(X_ecfp4, y, "MLP + ECFP4")
    y_test_m, y_prob_m, auc_m = run_mlp_benchmark(X_maccs, y, "MLP + MACCS")
    
    # 4. Vẽ biểu đồ
    plt.figure(figsize=(8, 6))
    fpr_e, tpr_e, _ = roc_curve(y_test_e, y_prob_e)
    plt.plot(fpr_e, tpr_e, label=f'MLP + ECFP4 (AUC = {auc_e:.3f})', lw=2)
    
    fpr_m, tpr_m, _ = roc_curve(y_test_m, y_prob_m)
    plt.plot(fpr_m, tpr_m, label=f'MLP + MACCS (AUC = {auc_m:.3f})', lw=2)
    
    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('MLP Benchmark: ECFP4 vs MACCS')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.show()

if __name__ == "__main__":
    main()