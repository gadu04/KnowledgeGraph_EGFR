import pandas as pd
import numpy as np

def force_label_uncertains(input_file, output_file, threshold_nm=1000.0):
    print(f"📖 Đang đọc file: {input_file}...")
    df = pd.read_csv(input_file)
    
    # Đếm số lượng Uncertain còn sót lại
    uncertain_mask = df['Final_Label'] == 'Uncertain'
    n_uncertain = uncertain_mask.sum()
    print(f"⚠️ Phát hiện {n_uncertain} chất vẫn còn nhãn 'Uncertain'.")
    
    if n_uncertain == 0:
        print("✅ Không còn chất Uncertain nào. Dữ liệu đã sạch!")
        return

    print(f"🛠 Đang thực hiện gán nhãn cưỡng bức dựa trên IC50 (Ngưỡng = {threshold_nm} nM)...")

    # Hàm xử lý logic
    def resolve_uncertain(row):
        # Nếu đã có nhãn xịn (Active/Inactive) thì giữ nguyên
        if row['Final_Label'] != 'Uncertain':
            return row['Final_Label']
        
        # Nếu là Uncertain, nhìn vào IC50
        ic50 = float(row['IC50 value(nM)'])
        
        if ic50 <= threshold_nm:
            return 'Active'  # IC50 thấp -> Hoạt tính mạnh
        else:
            return 'Inactive' # IC50 cao -> Hoạt tính yếu

    # Áp dụng
    df['Final_Label'] = df.apply(resolve_uncertain, axis=1)
    
    # Kiểm tra lại
    print("\n📊 Phân bố nhãn cuối cùng (Final Distribution):")
    print(df['Final_Label'].value_counts())
    
    # Lưu file
    df.to_csv(output_file, index=False)
    print(f"\n✅ Đã lưu bộ dữ liệu ĐẦY ĐỦ (Không xóa dòng nào) vào: {output_file}")
    print("👉 Bạn hãy dùng file này để Build Knowledge Graph và Train Model.")

# --- CHẠY ---
# Input là file kết quả của bước Similarity trước đó
force_label_uncertains(
    input_file='data_final.csv', 
    output_file='data_end.csv',
    threshold_nm=30.0 # Bạn có thể đổi thành 500 nếu muốn khắt khe hơn
)