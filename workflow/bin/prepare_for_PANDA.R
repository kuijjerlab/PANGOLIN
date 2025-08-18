require(data.table)
require("tidyr")
library("plyr")
panda_dir_primary <- "/storage/kuijjerarea/tatiana/DOBERMAN/pan_cancer_results/panda_input_primary"
project_dir <- "/storage/kuijjerarea/tatiana/DOBERMAN/"
data_dir <- file.path(project_dir, "data/")
panda_dir <- file.path(project_dir, "new_panda_input/")
ppi_file <- file.path(panda_dir, "/newest_priors/ppi_prior.tsv")
motif_prior <- file.path(panda_dir, "/newest_priors/tf_prior_names_fixed.tsv")

group_file <- file.path(data_dir, "/hg38_sample_groups.tsv")
batch_file <- file.path(data_dir, "/batch_information.txt")
features_file <- file.path(data_dir, "hg38_features.RData")
log2_exp_file <- file.path(data_dir, "log2_snail_norm_hg38.RData")
samples_clean_file <- file.path(panda_dir, "samples_cancers_clean.tsv")
combat_dir <- file.path(project_dir, "/combat_data/")
source(file.path(project_dir, "/scripts/analyze_batch_fn.R"))
# clean samples

samples_clean <- fread(samples_clean_file)
# load features file
load(features_file, features <- new.env())
features <- features[['features']]

# load batch info
batch_info <- load_batch(batch_file)
batch_info <- batch_info[batch_info$Tissues %in% c("cancer")] #10293
dim(batch_info)

# load expression file
load(log2_exp_file, exp <- new.env())
ls(exp)
log2exp <- exp[['log2exp']]
table(colnames(log2exp) %in% batch_info$Samples)

# read in info on cancers
groups <- fread(group_file, head = F)
head(groups)
# cancers that were combat corrected or some samples were removed as outliers
cancers_to_replace <- c("TCGA-COAD", "TCGA-DLBC", "TCGA-LUAD",
                        "TCGA-PRAD", "TCGA-READ", "TCGA-UCEC")

data_to_replace <- groups[groups$V2 %in% cancers_to_replace, ]
log2exp <- log2exp[, !colnames(log2exp) %in% data_to_replace$V1]

# include only cancer, remove the healthy tissues
log2exp <- log2exp[, colnames(log2exp) %in% batch_info$Samples]

# read in corrected data
exp_coad <- readRDS(file.path(combat_dir, "COAD_combat_exp.RData"))
exp_dlbc <- readRDS(file.path(combat_dir, "DLBC_combat_exp.RData"))
exp_luad <- readRDS(file.path(combat_dir, "LUAD_combat_exp.RData"))
exp_prad <- readRDS(file.path(combat_dir, "PRAD_nobatch_exp.RData"))
exp_read <- readRDS(file.path(combat_dir, "READ_combat_exp.RData"))
exp_ucec <- readRDS(file.path(combat_dir, "UCEC_combat_exp.RData"))

# bind with corrected data
log2exp_new <- cbind(log2exp, exp_coad, exp_dlbc, exp_luad,
                        exp_prad, exp_read, exp_ucec)

# extract only protein coding genes
log2exp_clean <- log2exp_new[features$gene_type %in% c("protein_coding"), ]
features_clean <- features[features$gene_type %in% c("protein_coding"), ]
head(features_clean)
length(features_clean$gene_id)==length(unique(features_clean$gene_id))

# there are some duplicate gene names
duplicates <- features_clean[duplicated(features_clean$gene_name) |
                        duplicated(features_clean$gene_name, fromLast = TRUE), ]

length(unique(duplicates$gene_name)) # 24 unique
#par_y duplicates to remove
par_y <- duplicates[grep("_PAR_Y", duplicates$gene_id), ]
# when par_y removed, what are remaining duplicates
duplicates <- duplicates[-grep("_PAR_Y", duplicates$gene_id), ]

duplicates <- duplicates[duplicated(duplicates$gene_name) |
                        duplicated(duplicates$gene_name, fromLast = TRUE), ]
length(unique(duplicates$gene_name)) # only 6 remains
# at the end we remove PAR_Y and the other duplicates (only 6)
genes_to_exclude <- c(par_y$gene_id, duplicates$gene_id)
genes_to_exclude

log2exp_clean <- log2exp_clean[!features_clean$gene_id %in% genes_to_exclude, ]

rownames(log2exp_clean) <-
        features_clean$gene_name[!features_clean$gene_id %in% genes_to_exclude]

dim(log2exp_clean) #19932

log2exp_clean <- log2exp_clean[, colnames(log2exp_clean) %in% samples_clean$sample_id]
log2exp_clean <- log2exp_clean[rowSums(log2exp_clean) > 0, ] #19556 genes, 10283 samples
non_zero <- apply(log2exp_clean, 1, function(i) sum(i > 0))
cancers <- data.table(groups[match(colnames(log2exp_clean), groups$V1),])
table(cancers$V2)
# the least number of samples is in TCGA-CHOL (35 samples)
table(non_zero <= 20) #remove genes present in less than 20 samples
log2exp_clean <- log2exp_clean[non_zero >=20, ]

# PPI 
ppi <- fread(ppi_file)
ppis <- unique(c(ppi$V1, ppi$V2))
head(ppi)
# motif prior
motif <- fread(motif_prior)
colnames(motif) <- c("tf", "gene", "prior")
motif_genes <- unique(motif$gene)
motif_tfs <- unique(motif$tf)
length(motif_genes) #61 050


table(motif_genes %in% rownames(log2exp_clean)) 
table(motif_tfs %in% rownames(log2exp_clean)) 


motif_clean <- motif[motif$gene %in% rownames(log2exp_clean), ]
motif_clean <- motif_clean[motif_clean$tf %in% rownames(log2exp_clean), ]

write.table(motif_clean, file.path(panda_dir_primary, "motif_tcga_primary.tsv"),
            col.names = F, sep = "\t", row.names = F, quote = F)

table(ppis %in% rownames(log2exp_clean)) # 631
ppi_clean <- ppi[ppi$V1 %in% rownames(log2exp_clean), ]
ppi_clean <- ppi_clean[ppi_clean$V2 %in% rownames(log2exp_clean), ]
write.table(ppi_clean, file.path(panda_dir_primary, "ppi_tcga_primary.tsv"),
            col.names = F, sep = "\t", row.names = F, quote = F)

samples <- data.table(colnames(log2exp_clean))
write.table(samples, file.path(panda_dir_primary, "samples_primary.tsv"),
            col.names = F, sep = "\t", row.names = T, quote = F)


genes_to_include <- c(motif_clean$tf, motif_clean$gene,
                        ppi_clean$V1, ppi_clean$V2)
genes_to_include <- unique(genes_to_include)

log2exp_clean <- log2exp_clean[rownames(log2exp_clean) %in% genes_to_include, ]
write.table(log2exp_clean, file.path(panda_dir_primary, "exp_tcga_primary.tsv"),
            col.names = F, sep = "\t", row.names = T, quote = F)

cancers <- data.table(groups[match(colnames(log2exp_clean), groups$V1), ])
colnames(cancers) <- c("sample_id", "cancer")
write.table(cancers, file.path(panda_dir_primary, "samples_cancers_primary.tsv"),
            col.names = T, sep = "\t", row.names = F, quote = F)

