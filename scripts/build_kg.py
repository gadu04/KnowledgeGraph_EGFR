"""Build Knowledge Graph from CSV data"""
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
from src.kg.builder import KnowledgeGraphBuilder
from src.config import EXPERIMENTAL_CSV, DENOVO_CSV


def main():
    """Main function to build knowledge graph"""
    print("\n" + "="*80)
    print("🏗️  BUILDING EGFR KNOWLEDGE GRAPH")
    print("="*80)
    
    # Initialize builder
    with KnowledgeGraphBuilder() as builder:
        builder.nuke_and_prepare_db()
        
        # Process experimental molecules
        print("\n" + "="*80)
        print("STEP 1: IMPORT EXPERIMENTAL DATA")
        print("="*80)
        try:
            df_experimental = pd.read_csv(EXPERIMENTAL_CSV)
            builder.process_experimental_molecules(df_experimental)
        except FileNotFoundError:
            print(f"⚠️  File not found: {EXPERIMENTAL_CSV}")
        
        # Process de novo molecules
        print("\n" + "="*80)
        print("STEP 2: IMPORT DE NOVO DATA (VIRTUAL)")
        print("="*80)
        try:
            df_denovo = pd.read_csv(DENOVO_CSV)
            builder.process_denovo_molecules(df_denovo)
        except FileNotFoundError:
            print(f"⚠️  File not found: {DENOVO_CSV}")
    
    print("\n" + "="*80)
    print("✅ KNOWLEDGE GRAPH BUILD COMPLETED!")
    print("="*80)


if __name__ == "__main__":
    main()
