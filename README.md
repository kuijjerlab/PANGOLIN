# PANGOLIN: Pan-cancer Analysis of Gene regulatory landscape Of LIoness Networks

## Overview

PANGOLIN is a comprehensive Snakemake pipeline for pan-cancer analysis of gene regulatory networks using TCGA data, with a specific focus on PD-1 pathway analysis. This pipeline reproduces all analyses and figures from our manuscript investigating gene regulatory landscapes across 33 cancer types.

## 🚀 Quick Start

### Prerequisites

- **Snakemake** ≥ 7.23.1
- **Conda** for environment management
- **R** ≥ 4.2.1 with Bioconductor packages
- **Python** ≥ 3.8
- **Sufficient disk space**: ~1TB for full workflow, ~100GB for precomputed analysis

### Installation

```bash
# Clone the repository
git clone https://github.com/tatianub/PANGOLIN.git
cd PANGOLIN
```

### Configuration

1. **Choose analysis type** in `config.yaml`:
   ```yaml
   # For complete reproduction (takes several days/weeks due to the single-sample network reconstruction for over 9 000 samples)
   analysis_type: "full_workflow"
   
   # For a more rapid analysis and figure reproduction using precomputed data (RECOMMENDED)
   analysis_type: "precomputed"
   ```

2. **Set up Zenodo credentials** (for precomputed workflow):
   ```yaml
   # This Zenodo ID contains precomputed network for PD1 pathway analysis
   zenodo_record_id: "10.5281/zenodo.17232919"
   ```

### Execution

```bash

# Dry run to check workflow
snakemake --use-conda --conda-frontend conda --cores 1 -np

# Execute pipeline
snakemake --use-conda --conda-frontend conda --cores 1
```

## 📊 Analysis Workflow

### Data Processing Pipeline

1. **📥 Data Acquisition**
   - Downloads TCGA expression data for 33 cancer types using TCGAbiolinks package
   - Retrieves clinical information
   - Processes batch information and produces batch figure 

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
├── data_all/                              # Pan-cancer results
│   ├── gdc_data/                            # Raw TCGA data
│   ├── batch_analysis/                      # Batch effect analysis
│   ├── batch_corrected_expression/          # Batch-corrected expression data
│   ├── clinical_associations_PD1/           # PD-1 clinical associations
│   ├── cola_consensus_clustering/           # all cancers consensus clustering
│   ├── combined_gdc_data/                   # Combined normalized expression
│   ├── cox_results_all/                     # survival analysis
│   ├── porcupine/                           # PORCUPINE pathway analysis
│   ├── pysnail_normalized_individual_cancer_expression/ # Normalized data per cancer
│   ├── tsne_results/                        # t-SNE dimensionality reduction
│   └── logs/                                # Processing logs
├── data_individual_cancers/                 # Cancer-specific results
│   └── [CANCER]/                            # Individual cancer directories
│       ├── pd1_data/                          # PD-1 scores and mappings
│       ├── cox/                               # Survival analysis
│       ├── consensus_clustering/              # Cancer-specific clustering
│       ├── final_clusters/                    # Final cluster assignments
│       ├── clinical_associations/             # Clinical correlations
│       ├── indegrees_norm/                    # Network in-degree features
│       ├── porcupine/                         # Pathway analysis results
│       └── clinical/                          # Clinical data processing
├── panda_input/                           # PANDA network input files
└── figs/                                  # Publication-ready figures
    ├── MBatch_DSC.pdf                       # Batch effect summary
    ├── TSNE_*.pdf                           # t-SNE visualizations
    ├── PC_immune_correlations_cibersort.png # Immune correlations
    ├── cox_results_final_clusters_*.pdf     # Survival analysis plots
    ├── PRAD_clusters_*.pdf                  # PRAD-specific analyses
    ├── pathways_intersection_pcp.pdf        # Pathway intersection analysis
    ├── sankey_plot_indegree_expression.pdf  # Sankey diagrams
    └── summary_table_PD1.html               # Results summary table
```

## 🎯 Key Features

### Analysis Types

- **Full Workflow**: Complete analysis from raw data (a very long runtime)
  - Downloads and processes all TCGA data
  - Constructs patient-specific regulatory networks
  - Performs all downstream analyses

- **Precomputed Workflow**: Rapid reproduction using intermediate files (~4 hours)
  - Downloads precomputed expression data from Zenodo
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

The pipeline generates all manuscript figures

## 📚 Citation

If you use PANGOLIN in your research, please cite:

```bibtex
@article{.....,
  title={Pan-cancer analysis of patient-specific gene regulatory landscapes identifies recurrent PD-1 pathway dysregulation},
  author={Belova et al.},
  journal={Journal Name},
  year={2025},
  doi={10.xxxx/xxxxx}
}
```

