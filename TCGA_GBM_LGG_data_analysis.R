# Load packages
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("SummarizedExperiment")
library("gProfileR")
library("genefilter")
library('sesame')
library('DT')
library('GDCquery_clinic')



# A total of 2.27 GB
query_exp_lgg <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp_lgg)
exp_lgg <- GDCprepare(
  query = query_exp_lgg
)

query_exp_gbm <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp_gbm)
exp_gbm <- GDCprepare(
  query = query_exp_gbm
)

# The following clinical data is not available in GBM
missing_cols <- setdiff(colnames(colData(exp_lgg)),colnames(colData(exp_gbm)))
for(i in missing_cols){
  exp_lgg[[i]] <- NULL
}



library(SummarizedExperiment)

# Load object from TCGAWorkflowData package
# This object will be created in subsequent sections for enhanced clarity and understanding.


# get genes information
genes.info <- rowRanges(exp_gbm)


# get sample information
sample.info <- colData(exp_gbm)
datatable(
  data = as.data.frame(sample.info), 
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# get  clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(
  project = "TCGA-GBM", 
  type = "Clinical"
)

# get clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(
  project = "TCGA-LGG", 
  type = "Clinical"
)


gbm_res = getResults(query_exp_gbm) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(gbm_res) # columns present in the table

lgg_res = getResults(query_exp_lgg) # make results as table
# head(lihc_res) # data of the first 6 patients.
colnames(lgg_res) # columns present in the table



head(gbm_res$sample_type) # first 6 types of tissue.
head(lgg_res$sample_type) # first 6 types of tissue.

summary(factor(gbm_res$sample_type)) # summary of distinct tissues types present in this study
summary(factor(lgg_res$sample_type)) # summary of distinct tissues types present in this study


query_exp_gbm = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

query_exp_lgg = GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query = query_exp_gbm)
GDCdownload(query = query_exp_lgg)


tcga_data_gbm = GDCprepare(query_exp_gbm)
tcga_data_LGG = GDCprepare(query_exp_lgg)

dim(tcga_data_gbm)
dim(tcga_data_LGG)


# Convert to dataframe
tcga_data_gbm_df <- as.data.frame(assay(tcga_data_gbm))
# Add rownames as a column
tcga_data_gbm_df$Gene <- rownames(tcga_data_gbm_df)
write.csv(tcga_data_gbm_df, "tcga_data_gbm.csv")                 # download exp df
write.csv(query_exp_gbm[[1]][[1]], "tcga_gbm_sampleinfo.csv")    # download sampleinfo df
write.csv(gbm_clin, "gbm_clinical.csv")


# Convert to dataframe
tcga_data_lgg_df <- as.data.frame(assay(tcga_data_LGG))
# Add rownames as a column
tcga_data_lgg_df$Gene <- rownames(tcga_data_LGG_df)
write.csv(tcga_data_lgg_df, "tcga_data_lgg.csv")                # download exp df
write.csv(query_exp_lgg[[1]][[1]], "tcga_lgg_sampleinfo.csv")  # download sampleinfo df
write.csv(lgg_clin, "lgg_clinical.csv")






###### Diffential expression of TMBIM6 in GBM compared to normal
# load library
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(ggpubr)
library(edgeR)

# Step 1: preparing count data--------

# read in counts data
dat_gbm <- read.csv('tcga_data_gbm.csv', row.names = "gene")
head(dat)
# read in sample info
metadat_gbm <- read.csv('col_data_GBM.csv', row.names = "sample")
head(metadat)



# making sure the row names in colData matches to column names in counts_data
all(colnames(dat_gbm) %in% rownames(metadat_gbm))
# are they in the same order?
all(colnames(dat_gbm) == rownames(metadat_gbm))

# Standardize format in dat
colnames(dat_gbm) <- gsub("\\.", "-", colnames(dat_gbm))

# Check if they are in the same order
all(colnames(dat_gbm) == rownames(metadat_gbm))


# Create DGEList object
dge <- DGEList(counts = dat_gbm)

# Calculate normalization factors using TMM
norm_factors <- calcNormFactors(dge)

# Ensure that the expression data is numeric
expression_data_gbm <- as.matrix(dat)

# Extract the normalization factors from the list
norm_factors_gbm <- as.numeric(norm_factors$samples$norm.factors)

# Normalize the data
normalized_data_gbm <- cpm(dge) / norm_factors


# Convert the matrix to a data frame
normalized_data_df_gbm <- as.data.frame(normalized_data_gbm)

# Add a 'gene' column based on row names
normalized_data_df_gbm$gene <- rownames(normalized_data_df_gbm)

# reshaping data - from wide to long--------
dat.long <- normalized_data_df_gbm %>%
  gather(key = 'samples', value = 'TMM', -gene)


gene_names <- rownames(normalized_data_gbm)
# Transfer row names to the first column
normalized_data_gbm <- cbind(gene = gene_names, normalized_data_gbm)
# Convert the matrix to a data frame
normalized_data_df_gbm <- as.data.frame(normalized_data_gbm)

# Reshape the data using gather
dat.long_gbm <- tidyr::gather(normalized_data_df_gbm, key = 'samples', value = 'TMM', -gene)

colnames(dat.long-gbm)
colnames(metadat-gbm)
str(dat.long_gbm$samples)
str(metadat-gbm$samples)
common_samples_gbm <- intersect(dat.long$samples, metadat$samples)
common_samples-gbm

# combine the dat.long and metadat.gbm
combined_data_gbm <- merge(dat.long, metadat, by = 'samples', all.x = TRUE)


## TMBIM6 Expression in GBM and Lgg compared to normal
# Specify the gene you want to plot
selected_gene_gbm <- "ENSG00000139644"

# Filter data for the selected gene
selected_gene_data_gbm <- combined_data_gbm %>% 
  filter(gene == selected_gene_gbm)

# Convert 'TMM' to numeric
selected_gene_data_gbm$TMM <- as.numeric(selected_gene_data_gbm$TMM)

# Boxplot for the specific gene using ggplot2
p_gene-gbm <- ggplot(selected_gene_data_gbm, aes(x = group, y = TMM, fill = group)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), 
               color = "black",      # Border color of the boxplot
               alpha = 0.9,          # Transparency of the boxplot
               outlier.shape = NA,   # Remove outliers
               outlier.colour = "green") +  # Outlier color
  geom_jitter(aes(group = group), width = 0.3, alpha = 0.5, color = "orange") + # Add jittered points
  labs(y = "TMBIM6 expression value", x = "TCGA-GBM") +
  scale_fill_manual(values = c("GBM" = "red", "Normal" = "blue")) 

# Add significance level using ggpubr
p_gene_gbm <- p_gene_gbm + stat_compare_means(
  comparisons = list(c("GBM", "Normal")),
  method = "wilcox.test",  # Use Mann-Whitney U test
  label = "p.format",
  position = position_dodge(width = 0.9),
  show.legend = TRUE,            # Show the legend for significance levels   
)

# Print boxplot for the specific gene with significance levels
print(p_gene_gbm)






###### Diffential expression of TMBIM6 in LGG compared to normal
# read in counts data
dat_lgg <- read.csv('tcga_data_lgg.csv', row.names = "gene")
head(dat_lgg)
# read in sample info
metadat_lgg <- read.csv('col_data_LGG.csv', row.names = "sample")
head(metadat_lgg)



# making sure the row names in colData matches to column names in counts_data
all(colnames(dat_lgg) %in% rownames(metadat_lgg))
# are they in the same order?
all(colnames(dat_lgg) == rownames(metadat_lgg))

# Standardize format in dat
colnames(dat_lgg) <- gsub("\\.", "-", colnames(dat_lgg))

# Check if they are in the same order
all(colnames(dat_lgg) == rownames(metadat_lgg))

library(edgeR)
# Create DGEList object
dge <- DGEList(counts = dat_lgg)

# Calculate normalization factors using TMM
norm_factors <- calcNormFactors(dge)

# Ensure that the expression data is numeric
expression_data_lgg <- as.matrix(dat)

# Extract the normalization factors from the list
norm_factors_lgg <- as.numeric(norm_factors$samples$norm.factors)

# Normalize the data
normalized_data_lgg <- cpm(dge) / norm_factors



# Convert the matrix to a data frame
normalized_data_df_lgg <- as.data.frame(normalized_data)

# Add a 'gene' column based on row names
normalized_data_df_lgg$gene <- rownames(normalized_data_df)

# reshaping data - from wide to long--------
dat.long_lgg <- normalized_data_df_lgg %>%
  gather(key = 'samples', value = 'TMM', -gene)


gene_names_lgg <- rownames(normalized_data)
# Transfer row names to the first column
normalized_data_lgg <- cbind(gene = gene_names, normalized_data)
# Convert the matrix to a data frame
normalized_data_df_lgg <- as.data.frame(normalized_data)


# Reshape the data using gather
dat.long_lgg <- tidyr::gather(normalized_data_df, key = 'samples', value = 'TMM', -gene)


colnames(dat.long_lgg)
colnames(metadat_lgg)
str(dat.long_lgg$samples)
str(metadat_lgg$samples)
common_samples_lgg <- intersect(dat.long$samples, metadat$samples)
common_samples_lgg

combined_data_lgg <- merge(dat.long_lgg, metadat_lgg, by = 'samples', all.x = TRUE)


##Only for TMBIM6
# Specify the gene you want to plot
selected_gene_lgg <- "ENSG00000139644.13"

# Filter data for the selected gene
selected_gene_data_lgg <- combined_data %>% 
  filter(gene == selected_gene)

# Convert 'TMM' to numeric
selected_gene_data_lgg$TMM <- as.numeric(selected_gene_data$TMM)

# Boxplot for the specific gene using ggplot2
p_gene_lgg <- ggplot(selected_gene_data_lgg, aes(x = group, y = TMM, fill = group)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), 
               color = "black",      # Border color of the boxplot
               alpha = 0.9,          # Transparency of the boxplot
               outlier.shape = NA,   # Remove outliers
               outlier.colour = "green") +
  geom_jitter(aes(group = group), width = 0.3, alpha = 0.5, color = "orange") + # Add jittered points# Outlier color
  labs(y = "TMBIM6 expression value", x = "TCGA-LGG") +
  scale_fill_manual(values = c("LGG" = "red", "Normal" = "blue")) 

# Add significance level using ggpubr
p_gene_lgg <- p_gene + stat_compare_means(
  comparisons = list(c("LGG", "Normal")),
  method = "wilcox.test",  # Use Mann-Whitney U test
  label = "p.format",
  position = position_dodge(width = 0.9),
  show.legend = TRUE,            # Show the legend for significance levels   
)

# Print boxplot for the specific gene with significance levels
print(p_gene_lgg)
