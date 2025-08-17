#!/bin/bash

# List of cancer type abbreviations
CANCER_TYPES=("ACC" "CHOL" "DLBC" "BLCA" "CESC" "BRCA" "COAD" "ESCA" "GBM" "HNSC" "KICH" "KIRC" "KIRP" "LAML" "LGG" "LIHC" "LUAD" "LUSC" "MESO" "OV" "PAAD" "PCPG" "PRAD" "READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS" "UVM")

# Loop through each cancer type and move the corresponding file
for CANCER in "${CANCER_TYPES[@]}"; do
    mv "DOBERMAN/pan_cancer_results/PORCUPINE_results_modified/pcp_results_with_variance_${CANCER}.txt" "PANGOLIN/data_individual_cancers/${CANCER}/porcupine/"
done

echo "Files moved successfully."
