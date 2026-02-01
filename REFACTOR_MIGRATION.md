# 🚀 REFACTOR MIGRATION GUIDE

## ✅ ĐÃ HOÀN THÀNH

### 1. Cấu trúc thư mục mới

```
KnowledgeGraph_EGFR/
├── src/                    ✅ Created
│   ├── config/            ✅ Settings & configuration
│   ├── kg/                ✅ Knowledge Graph builder
│   ├── models/            ✅ ML models (ready for HGT/GraphSAGE)
│   ├── preprocessing/     ✅ Data preprocessing
│   ├── evaluation/        ✅ Model evaluation
│   └── utils/             ✅ Chemistry utilities
├── scripts/               ✅ Executable scripts
├── notebooks/
│   ├── experiments/       ✅ Stability tests moved
│   └── exploratory/       ✅ Bench notebook moved
├── data/
│   ├── processed/         ✅ data_end.csv, DeNovo_Molecule.csv
│   └── results/           ✅ multi_seed_results*.csv
├── tests/                 ✅ Ready for unit tests
└── archive/               ✅ Backup files moved
```

### 2. Files đã tạo

**Config:**
- ✅ `src/config/settings.py` - Cấu hình toàn cục
- ✅ `src/config/__init__.py` - Export config variables

**Knowledge Graph:**
- ✅ `src/kg/builder.py` - KnowledgeGraphBuilder class
- ✅ `src/kg/__init__.py` - Export builder

**Utilities:**
- ✅ `src/utils/chemistry.py` - RDKit chemistry functions
- ✅ `src/utils/__init__.py` - Export chemistry utils

**Scripts:**
- ✅ `scripts/build_kg.py` - Entry point để build KG

**Config:**
- ✅ `.env.example` - Template cho environment variables

**Documentation:**
- ✅ `README.md` - Updated với cấu trúc mới

### 3. Files đã di chuyển

**Data:**
- ✅ `Data/data_end.csv` → `data/processed/data_end.csv`
- ✅ `Data/DeNovo_Molecule.csv` → `data/processed/DeNovo_Molecule.csv`
- ✅ `multi_seed_results*.csv` → `data/results/`
- ✅ `output/*.csv` → `data/results/`

**Notebooks:**
- ✅ `bench.ipynb` → `notebooks/exploratory/`
- ✅ `stability_test_hgt.ipynb` → `notebooks/experiments/`
- ✅ `stability_test_graphsage.ipynb` → `notebooks/experiments/`

**Backup:**
- ✅ `backup/*` → `archive/`

---

## 🔄 BƯỚC TIẾP THEO

### 1. Xóa files deprecated (OLD code)

```powershell
# Di chuyển old files vào archive
Move-Item "BuildKG.py" "archive/BuildKG_old.py" -Force
Move-Item "data.py" "archive/data_old.py" -Force
Move-Item "Eval.py" "archive/Eval_old.py" -Force

# Xóa thư mục trống
Remove-Item "backup" -Force -ErrorAction SilentlyContinue
Remove-Item "output" -Force -ErrorAction SilentlyContinue
Remove-Item "Data" -Force -ErrorAction SilentlyContinue
```

### 2. Test script mới

```bash
# Test build KG
python scripts/build_kg.py
```

### 3. Update notebooks (nếu cần)

Notebooks cần update import paths:

**OLD:**
```python
from BuildKG import *
```

**NEW:**
```python
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.kg.builder import KnowledgeGraphBuilder
from src.config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD
```

### 4. Commit refactored code

```bash
# Stage new files
git add src/ scripts/ notebooks/ data/ .env.example README.md

# Remove old files from git
git rm BuildKG.py data.py Eval.py

# Commit
git commit -m "refactor: Reorganize code with professional structure

- src/config: Centralized configuration
- src/kg: Knowledge Graph builder (KANO architecture)
- src/utils: Reusable chemistry utilities
- scripts: Executable entry points
- notebooks: Organized experiments & exploratory
- data: Structured data files (processed & results)
- Clean separation of concerns"

# Push
git push origin main
```

---

## 📊 Lợi ích của refactor

### 1. **Code Organization**
- ✅ Separation of concerns: Mỗi module có trách nhiệm riêng
- ✅ Reusable utilities: Chemistry functions có thể dùng ở nhiều nơi
- ✅ Easy imports: `from src.kg.builder import KnowledgeGraphBuilder`

### 2. **Maintainability**
- ✅ Dễ tìm code: Biết chính xác file nào ở đâu
- ✅ Dễ test: Unit tests cho từng module
- ✅ Dễ extend: Thêm models/evaluation mới vào đúng thư mục

### 3. **Professional**
- ✅ Follow Python best practices
- ✅ Cấu trúc giống các open-source projects
- ✅ Dễ onboard developers mới

### 4. **Scalability**
- ✅ Dễ thêm features mới
- ✅ Dễ refactor từng phần
- ✅ Dễ tạo CI/CD pipeline

---

## 🧪 Verification Checklist

- [x] Cấu trúc thư mục đã tạo
- [x] Files config đã tạo
- [x] KG Builder đã refactor
- [x] Chemistry utils đã tách riêng
- [x] Entry point script đã tạo
- [x] Data files đã di chuyển
- [x] Notebooks đã di chuyển
- [x] Backup files đã lưu
- [x] README đã update
- [x] .env.example đã tạo
- [x] Test imports thành công
- [ ] Notebooks update imports (TODO)
- [ ] Old files đã xóa (TODO)
- [ ] Git commit (TODO)
- [ ] Git push (TODO)

---

## 📝 Notes

### Import trong Notebooks

Thêm vào đầu notebook:

```python
import sys
from pathlib import Path

# Add project root to path
project_root = Path.cwd().parent if 'notebooks' in str(Path.cwd()) else Path.cwd()
sys.path.insert(0, str(project_root))

# Now can import
from src.config import NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD
from src.kg.builder import KnowledgeGraphBuilder
from src.utils.chemistry import get_ecfp4, predict_target
```

### Next Refactoring Tasks

1. **Models Module** (`src/models/`)
   - Extract HGT model từ notebooks
   - Extract GraphSAGE model từ notebooks
   - Tạo base model class

2. **Evaluation Module** (`src/evaluation/`)
   - Extract evaluation metrics
   - Tạo evaluator class

3. **Preprocessing Module** (`src/preprocessing/`)
   - Extract data preprocessing logic
   - Tạo data loader class

4. **Tests** (`tests/`)
   - Unit tests cho chemistry utils
   - Unit tests cho KG builder
   - Integration tests

---

## 🎯 Summary

**Trước refactor:**
```
BuildKG.py (800+ lines)
data.py (messy)
Eval.py (mixed concerns)
```

**Sau refactor:**
```
src/
  config/settings.py (50 lines - config only)
  kg/builder.py (150 lines - KG logic only)
  utils/chemistry.py (180 lines - chemistry only)
scripts/build_kg.py (40 lines - entry point only)
```

✅ **Cleaner, more maintainable, professional!**
