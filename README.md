# Drug Knowledge Graph for EGFR Inhibitor Prediction

## 🎯 Giới thiệu

Dự án xây dựng **Knowledge Graph** cho dự đoán hoạt tính của chất ức chế EGFR, kết hợp:
- **KANO Architecture**: Knowledge-Aware Neural Operator
- **HGT Model**: Heterogeneous Graph Transformer
- **GraphSAGE**: Graph Sample and Aggregate
- **Neo4j**: Graph Database

---

## 📂 Cấu trúc thư mục (REFACTORED)

```
KnowledgeGraph_EGFR/
├── src/                           # Source code
│   ├── config/                    # Cấu hình toàn cục
│   ├── kg/                        # Knowledge Graph builder
│   ├── models/                    # ML Models (HGT, GraphSAGE)
│   ├── preprocessing/             # Data preprocessing
│   ├── evaluation/                # Model evaluation
│   └── utils/                     # Utilities (chemistry.py)
├── scripts/                       # Executable scripts
│   └── build_kg.py                # Build Knowledge Graph
├── notebooks/                     # Jupyter notebooks
│   ├── experiments/               # Stability tests
│   └── exploratory/               # Benchmark analysis
├── data/
│   ├── processed/                 # data_end.csv, DeNovo_Molecule.csv
│   └── results/                   # multi_seed_results*.csv
├── tests/                         # Unit tests
├── archive/                       # Backup code cũ
└── neo4j_data/                    # Neo4j database (gitignored)
```

---

## 🚀 Cài đặt

### Prerequisites
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) hoặc [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Docker](https://www.docker.com/) (cho Neo4j)

### Bước 1: Clone repository

```bash
git clone https://github.com/gadu04/KnowledgeGraph_EGFR.git
cd KnowledgeGraph_EGFR
```

### Bước 2: Tạo môi trường Conda

```bash
conda env create -f environment.yml
conda activate egfr_ml
```

### Bước 3: Cấu hình môi trường

```bash
cp .env.example .env
# Chỉnh sửa .env với mật khẩu Neo4j của bạn
```

### Bước 4: Khởi động Neo4j

```bash
docker-compose up -d
```

---

## 📊 Sử dụng

### 1. Build Knowledge Graph

```bash
python scripts/build_kg.py
```

### 2. Chạy thí nghiệm

**HGT Model:**
```bash
jupyter notebook notebooks/experiments/stability_test_hgt.ipynb
```

**GraphSAGE Model:**
```bash
jupyter notebook notebooks/experiments/stability_test_graphsage.ipynb
```

**Benchmark:**
```bash
jupyter notebook notebooks/exploratory/bench.ipynb
```

---

## 🧪 Kiến trúc Knowledge Graph

### Node Types
- **Molecule**: Phân tử (experimental/virtual)
- **Scaffold**: Murcko scaffold
- **Target**: EGFR_WT, EGFR_T790M, EGFR_Generic
- **Warhead**: Acrylamide, Propynamide, etc.
- **MoA**: Covalent/Reversible Inhibitor
- **FunctionalGroup**: Quinazoline_Core, Aniline_Group, etc.

### Relationships
- `(Molecule)-[:HAS_SCAFFOLD]->(Scaffold)`
- `(Molecule)-[:TESTED_AGAINST]->(Target)`
- `(Molecule)-[:POTENT_AGAINST]->(Target)` (active only)
- `(Molecule)-[:CONTAINS_WARHEAD]->(Warhead)`
- `(Molecule)-[:ACTS_VIA]->(MoA)`
- `(Molecule)-[:HAS_FUNCTIONAL_GROUP]->(FunctionalGroup)`

---

## 📈 Kết quả

Kết quả thí nghiệm trong [`data/results/`](data/results):
- [`multi_seed_results.csv`](data/results/multi_seed_results.csv) - HGT (10 seeds)
- [`multi_seed_results_graphsage.csv`](data/results/multi_seed_results_graphsage.csv) - GraphSAGE
- [`multi_seed_results_corrected.csv`](data/results/multi_seed_results_corrected.csv) - Corrected

**Metrics:** Accuracy, Precision, Recall, F1-score, ROC-AUC

---

## 🔐 Bảo mật

⚠️ **KHÔNG BAO GIỜ commit file `.env` lên Git!**

File [`.gitignore`](.gitignore) đã được cấu hình để bỏ qua file này.

---

## 🐛 Troubleshooting

**Lỗi RDKit import:**
```bash
conda install -c conda-forge rdkit
```

**Lỗi Neo4j connection:**
```bash
docker ps
docker logs <container_id>
docker-compose restart
```

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