# PANGOLIN: Pan-cancer Analysis of Gene regulatory landscape Of LIoness Networks

## Overview

PANGOLIN is a comprehensive Snakemake pipeline for pan-cancer analysis of gene regulatory networks using TCGA data, with a specific focus on PD-1 pathway analysis. This pipeline reproduces all analyses and figures from our manuscript investigating gene regulatory landscapes across cancer types.

## 🚀 Quick Start

### Prerequisites

- **Snakemake** ≥ 7.23.1
- **Conda/Mamba** for environment management
- **R** ≥ 4.2.1 with Bioconductor packages
- **Python** ≥ 3.8
- **Sufficient disk space**: ~700GB for full workflow, ~50GB for precomputed analysis

### Installation

```bash
# Clone the repository
git clone https://github.com/[username]/PANGOLIN.git
cd PANGOLIN
```

### Configuration

1. **Choose analysis type** in `config.yaml`:
   ```yaml
   # For complete reproduction (takes several days/weeks due to the single-sample network reconstruction for over 9 000 samples
   analysis_type: "full_workflow"
   
   # For rapid figure reproduction using precomputed data (RECOMMENDED)
   analysis_type: "precomputed"
   ```

2. **Set up Zenodo credentials** (for precomputed workflow):
   ```yaml
   zenodo_record_id: "YOUR_ZENODO_RECORD_ID"
   ```

### Execution

```bash
# Load required modules (on HPC systems)
module load snakemake/7.23.1-foss-2022a
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1

# Dry run to check workflow
snakemake --cores 1 -np

# Execute pipeline
snakemake --use-conda --conda-frontend conda --cores 1
```

## 📊 Analysis Workflow

### Data Processing Pipeline

1. **📥 Data Acquisition**
   - Downloads TCGA expression data for 33 cancer types
   - Retrieves clinical information
   - Processes batch information and produces batch figures (this is rather time consuming step)

2. **🔧 Expression Data Normalization**
   - Combines multi-cancer expression matrices
   - Applies qsmooth normalization using PySNAIL
   - Performs batch effect detection and correction

3. **🕸️ Network Inference** (Full workflow only)
   - Constructs gene regulatory network using PANDA
   - Generates patient-specific networks with LIONESS
   - Calculates network-based features (in-degree)

4. **📈 Dimensionality Reduction & Pathway Analysis**
   - t-SNE analysis of expression and network features
   - PORCUPINE pathway heterogeneity analysis 
   - Reactome pathway enrichment

5. **🧬 PD-1 Pathway Analysis**
   - Extracts PD-1 pathway components and scores
   - Correlates with immune infiltration (CIBERSORTx)
   - Clinical association analysis (survival, molecular features)

6. **🎯 Clustering Analysis**
   - Consensus clustering using Cola
   - Cancer type-specific cluster characterization
   - Survival analysis of identified clusters
   - Comparison of clusters for PRAD

## 📁 Output Structure

```
results/
├── all_cancers/                    # Pan-cancer results
│   ├── gdc_data/                  # Raw TCGA data
│   ├── batch_analysis/            # Batch effect analysis
│   ├── porcupine_results/         # Pathway analysis
│   ├── clinical_associations_PD1/ # PD-1 clinical associations
│   └── combined_gdc_data/         # Normalized expression
├── individual_cancers/            # Cancer-specific results
│   └── [CANCER]/
│       ├── pd1_data/             # PD-1 scores and mappings
│       ├── cox/                  # Survival analysis
│       ├── clustering/           # Consensus clustering
│       └── clinical_associations/ # Clinical correlations
└── figures/                       # Publication-ready figures
    ├── MBatch_DSC.pdf            # Batch effect summary
    ├── TSNE_all_cancers*.pdf     # t-SNE visualizations
    ├── PC_immune_correlations.png # Immune correlations
    └── summary_table_PD1.html    # Results summary
```

## 🎯 Key Features

### Analysis Types

- **Full Workflow**: Complete analysis from raw data (~3-7 days runtime)
  - Downloads and processes all TCGA data
  - Constructs patient-specific regulatory networks
  - Performs all downstream analyses

- **Precomputed Workflow**: Rapid reproduction using intermediate files (~2-4 hours)
  - Downloads precomputed network features from Zenodo
  - Peforms most of the downstream analysis (exluding network generation)
  - Focuses on statistical analysis and figure generation


### Cancer Types Analyzed

33 TCGA cancer types: ACC, BLCA, BRCA, CESC, CHOL, COAD, DLBC, ESCA, GBM, HNSC, KICH, KIRC, KIRP, LAML, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS, UVM

### Software Versions

- Snakemake: 7.23.1
- R: 4.2.1
- Python: 3.8+
- Key R packages: survival, dplyr, data.table, ComplexHeatmap
- Key Python packages: pandas, numpy, scipy


## 📊 Generated Figures

The pipeline generates all manuscript figures:

1. **Figure 1**: Pan-cancer t-SNE analysis
2. **Figure 2**: Batch effect analysis and correction
3. **Figure 3**: PD-1 pathway activity across cancers
4. **Figure 4**: Immune infiltration correlations
5. **Figure 5**: Clinical association heatmaps
6. **Supplementary Figures**: Additional analyses and validations

## 📚 Citation

If you use PANGOLIN in your research, please cite:

```bibtex
@article{.....,
  title={PPan-cancer analysis of patient-specific gene regulatory landscapes identifies recurrent PD-1 pathway dysregulation},
  author={Belova et al.},
  journal={Journal Name},
  year={2025},
  doi={10.xxxx/xxxxx}
}
```

