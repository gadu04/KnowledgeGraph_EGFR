# Drug Knowledge Graph for EGFR Inhibitor Prediction

## 🚀 Quick Start

### Prerequisites
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) hoặc [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Neo4j Database (sử dụng Docker hoặc local)

### Installation

**Option 1: Sử dụng Conda (Khuyến nghị)**

```bash
# Clone repository
git clone https://github.com/your-username/Drug-KG.git
cd Drug-KG

# Tạo môi trường từ file environment.yml
conda env create -f environment.yml

# Kích hoạt môi trường
conda activate egfr_ml
```

**Option 2: Sử dụng pip**

```bash
# Tạo môi trường mới
conda create -n egfr_ml python=3.10
conda activate egfr_ml

# Cài RDKit qua conda (quan trọng!)
conda install -c conda-forge rdkit

# Cài các package còn lại
pip install -r requirements.txt
```

### Setup Environment

1. **Tạo file `.env`:**
```bash
cp .env.example .env
```

2. **Chỉnh sửa `.env`:**
```env
NEO4J_URI=bolt://localhost:7687
NEO4J_USER=neo4j
NEO4J_PASSWORD=your_password_here
```

3. **Khởi động Neo4j:**
```bash
docker-compose up -d
```

### Build Knowledge Graph

```bash
python BuildKG.py
```

### Run Experiments

```bash
# Mở Jupyter Notebook
jupyter notebook

# Hoặc chạy stability test trực tiếp
jupyter nbconvert --to notebook --execute stability_test_hgt.ipynb
```

## 📂 Project Structure

```
Drug-KG/
├── BuildKG.py              # Xây dựng Knowledge Graph
├── bench.ipynb            # Benchmark so sánh ECFP4 vs KG
├── stability_test_hgt.ipynb      # Test ổn định HGT model
├── stability_test_graphsage.ipynb # Test ổn định GraphSAGE
├── Data/
│   ├── data_end.csv       # Dữ liệu experimental
│   └── DeNovo_Molecule.csv # Dữ liệu de novo
├── environment.yml         # Conda environment
├── requirements.txt        # Pip requirements
└── .env                   # Cấu hình (không commit)
```

## 🔐 Security Note

⚠️ **KHÔNG BAO GIỜ commit file `.env` lên Git!**

File [`.gitignore`](.gitignore) đã được cấu hình để bỏ qua file này.

## 📊 Results

Kết quả thí nghiệm được lưu trong:
- `multi_seed_results.csv` - HGT results
- `multi_seed_results_graphsage.csv` - GraphSAGE results
- `multi_seed_results_corrected.csv` - Corrected results

## 🐛 Troubleshooting

**Lỗi RDKit import:**
```bash
conda install -c conda-forge rdkit
```

**Lỗi Neo4j connection:**
- Kiểm tra Docker: `docker ps`
- Kiểm tra `.env` file
- Kiểm tra port 7687 không bị chiếm

**Lỗi PyTorch Geometric:**
```bash
pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
pip install torch-geometric
```

## 📝 Citation

Nếu sử dụng code này, vui lòng cite:

```bibtex
@software{drug_kg_2024,
  author = {Your Name},
  title = {Drug Knowledge Graph for EGFR Inhibitor Prediction},
  year = {2024},
  url = {https://github.com/gadu04/KnowledgeGraph_EGFR}
}
```

## 📄 License

MIT License - xem file LICENSE để biết chi tiết.