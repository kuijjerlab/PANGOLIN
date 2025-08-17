#!/bin/bash

# List of cancer type abbreviations
# CANCER_TYPES=("BLCA" "CESC" "BRCA" "COAD" "ESCA" "GBM" "HNSC" "KICH" "KIRC" "KIRP" "LAML" "LGG" "LIHC" "LUAD" "LUSC" "MESO" "OV" "PAAD" "PCPG" "PRAD" "READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS" "UVM")
CANCER_TYPES=("ACC" "CHOL" "DLBC")
# Loop through each cancer type and move the corresponding file
for CANCER in "${CANCER_TYPES[@]}"; do
    mv "DOBERMAN/indegrees_norm_edges/indegree_norm_${CANCER}.RData" "PANGOLIN/data_individual_cancers/${CANCER}/indegrees_norm/"
done

echo "Files moved successfully."
