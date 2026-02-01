import pandas as pd

# 1. Đọc dữ liệu từ file CSV
df = pd.read_csv('Data/docking_results.csv')

# 2. Sắp xếp theo Affinity tăng dần 
# Trong docking, giá trị năng lượng (kcal/mol) càng âm thì ái lực liên kết càng mạnh.
top_30_egfr = df.sort_values(by='affinity', ascending=True).head(30)

# 3. Lưu kết quả ra file CSV mới
top_30_egfr.to_csv('egfr.csv', index=False)

# Hiển thị 5 dòng đầu tiên của kết quả
print(top_30_egfr.head())
print(f"\nĐã trích xuất thành công {len(top_30_egfr)} chất vào file 'egfr.csv'")