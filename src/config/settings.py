"""Cấu hình toàn cục cho project"""
import os
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Project root
PROJECT_ROOT = Path(__file__).parent.parent.parent

# Neo4j Configuration
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD")

# Data Paths
DATA_DIR = PROJECT_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
RESULTS_DIR = DATA_DIR / "results"

# Input files
EXPERIMENTAL_CSV = PROCESSED_DATA_DIR / "data_end.csv"
DENOVO_CSV = PROCESSED_DATA_DIR / "DeNovo_Molecule.csv"

# Model Configuration
RANDOM_SEED = 42
TEST_SIZE = 0.2
EMBED_DIM = 128
BATCH_SIZE = 1000

# HGT Configuration
HGT_HIDDEN_CHANNELS = 128
HGT_NUM_LAYERS = 2
HGT_NUM_HEADS = 4
HGT_EPOCHS = 50

# GraphSAGE Configuration
SAGE_HIDDEN_CHANNELS = 128
SAGE_NUM_LAYERS = 2
SAGE_EPOCHS = 100

# Random Forest Configuration
RF_N_ESTIMATORS = 100
RF_RANDOM_STATE = RANDOM_SEED
