# KG-MF: Drug Knowledge Graph for EGFR Activity Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/) [![Neo4j](https://img.shields.io/badge/Neo4j-Graph%20DB-brightgreen)](https://neo4j.com/)

Hybrid scaffold-hopping pipeline for EGFR activity prediction. This repository combines:

- MolDeBERTa-FT — SMILES-fine-tuned MolDeBERTa with focal loss and SMILES enumeration augmentation.
- KG-MF — Rank-average fusion of MolDeBERTa-FT scores with a Neo4j knowledge-graph branch using feature-based KNN projection.
- HGT — Heterogeneous Graph Transformer experiments for KG-aware GNN baselines.
- GraphSAGE — GraphSAGE baselines and stability experiments.

Table of Contents
-----------------

1. [Overview](#overview)
2. [Repository Structure](#repository-structure)
3. [Quickstart](#quickstart)
4. [Usage](#usage)
5. [Knowledge Graph Schema](#knowledge-graph-schema)
6. [Results & Metrics](#results--metrics)
7. [Troubleshooting](#troubleshooting)
8. [Citation](#citation)
9. [License](#license)

Overview
--------

The project follows three design principles:

1. Scaffold-based splitting to avoid leakage.
2. A feature-based KG branch (not an end-to-end GNN) to preserve interpretability.
3. Validation-only threshold tuning with untouched test evaluation.

Repository Structure
--------------------

Simplified view (top-level):

```
.
├── src/                 # Source code: config, KG builder, models, utils
├── scripts/             # Scripts (e.g., build_kg.py)
├── notebooks/           # Notebooks for experiments and benchmarks
├── data/                # Processed data and experiment results
└── neo4j_data/          # Local Neo4j DB files (gitignored)
```

Quickstart
----------

Prerequisites

- Conda or Miniconda
- Docker (for Neo4j)

Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate egfr_ml
```

Copy and edit environment variables (set Neo4j password):

```bash
cp .env.example .env
# Edit .env to set NEO4J_PASSWORD and other vars
```

Start Neo4j (Docker):

```bash
docker-compose up -d
```

Default endpoints used by the project:

- Bolt: `bolt://localhost:7688`
- Browser: `http://localhost:7475`

Usage
-----

Build the knowledge graph from processed input data:

```bash
python scripts/build_kg.py
```

Run experiments (examples):

- HGT stability tests: open `notebooks/experiments/stability_test_hgt.ipynb`
- GraphSAGE stability tests: open `notebooks/experiments/stability_test_graphsage.ipynb`
- Benchmarks: open `notebooks/exploratory/bench.ipynb`

Knowledge Graph Schema
----------------------

Primary node types:

- Molecule (experimental / virtual)
- Scaffold (Murcko scaffold)
- Target (e.g., EGFR_WT, EGFR_T790M)
- Warhead (e.g., Acrylamide)
- MoA (mechanism of action: covalent / reversible)
- FunctionalGroup (e.g., Quinazoline_Core)

Core relationships:

- (Molecule)-[:HAS_SCAFFOLD]->(Scaffold)
- (Molecule)-[:TESTED_AGAINST]->(Target)
- (Molecule)-[:POTENT_AGAINST]->(Target)  // active labels
- (Molecule)-[:CONTAINS_WARHEAD]->(Warhead)
- (Molecule)-[:ACTS_VIA]->(MoA)
- (Molecule)-[:HAS_FUNCTIONAL_GROUP]->(FunctionalGroup)

Results & Metrics
-----------------

Experiment outputs are stored under `data/results/`:

- `multi_seed_results.csv` — HGT (10 seeds)
- `multi_seed_results_graphsage.csv` — GraphSAGE
- `multi_seed_results_corrected.csv` — Corrected/merged results

Reported metrics include: Accuracy, Precision, Recall, F1-score, ROC-AUC.

Troubleshooting
---------------

- RDKit import error:

```bash
conda install -c conda-forge rdkit
```

- Neo4j connection issues:

```bash
docker ps
docker logs <container_id>
docker-compose restart
```

- PyTorch Geometric installation (CPU example):

```bash
pip install torch-scatter torch-sparse torch-cluster -f https://data.pyg.org/whl/torch-2.0.0+cpu.html
pip install torch-geometric
```

Citation
--------

If you use this code, please cite:

```bibtex
@software{drug_kg_2024,
  author = {gadu04},
  title = {KG-MF: Drug Knowledge Graph for EGFR Activity Prediction},
  year = {2024},
  url = {https://github.com/gadu04/KnowledgeGraph_EGFR}
}
```

License
-------

This project is released under the MIT License. See the `LICENSE` file for details.