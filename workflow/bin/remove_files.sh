#!/bin/bash

# List of cancer type abbreviations
CANCER_TYPES=("ACC" "CHOL" "DLBC" "BLCA" "CESC" "BRCA" "COAD" "ESCA" "GBM" "HNSC" "KICH" "KIRC" "KIRP" "LAML" "LGG" "LIHC" "LUAD" "LUSC" "MESO" "OV" "PAAD" "PCPG" "PRAD" "READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM" "UCEC" "UCS" "UVM")

# Loop through each cancer type and move the corresponding file
for CANCER in "${CANCER_TYPES[@]}"; do
    rm "PANGOLIN/data_individual_cancers/${CANCER}/pd1_data/pd1_net_norm_${CANCER}.RData"
    mv "PANGOLIN/data_individual_cancers/${CANCER}/pd1_data/pd1_net_norm_named_${CANCER}.RData" "PANGOLIN/data_individual_cancers/${CANCER}/pd1_data/pd1_net_norm_${CANCER}.RData"
done

echo "Files moved successfully."
