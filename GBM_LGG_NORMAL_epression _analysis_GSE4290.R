# load library
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(ggpubr)
getwd()


# Step 1: preparing count data--------

# read in counts data
dat <- read.csv('counts_data_GBM_LGG_normal.csv', row.names = "gene")
head(dat)
# read in sample info
metadat <- read.csv('col_data.csv', row.names = "sample")
head(metadat)



# making sure the row names in colData matches to column names in counts_data
all(colnames(dat) %in% rownames(metadat))
# are they in the same order?
all(colnames(dat) == rownames(metadat))

# Standardize format in dat
colnames(dat) <- gsub("\\.", "-", colnames(dat))

# Check if they are in the same order
all(colnames(dat) == rownames(metadat))



library(edgeR)

# Add a 'gene' column based on row names
dat$gene <- rownames(dat)

# Reshaping data - from wide to long
dat.long <- dat %>%
  gather(key = 'samples', value = 'log2(CPM+1)', -gene)

#make a column of gene in normalized_data
gene_names <- rownames(dat)
# Transfer row names to the first column
dat <- cbind(gene = gene_names, dat)
# Convert the matrix to a data frame
df <- as.data.frame(dat)

# Reshape the data using gather
dat.long <- tidyr::gather(dat, key = 'samples', value = 'log2(CPM+1)', -gene)

combined_data <- merge(dat.long, metadat, by = 'samples', all.x = TRUE)

# Remove the column you want to exclude
combined_data <- combined_data %>%
  select(-3)


##Only for TMBIM6
selected_gene <- "TMBIM6"

# Filter data for the selected gene
selected_gene_data <- combined_data %>% 
  filter(gene == selected_gene)


library(ggplot2)


# Boxplot for the specific gene using ggplot2
p_gene <- ggplot(selected_gene_data, aes(x = group, y = `log2(CPM+1)`, fill = group)) +
  geom_boxplot(width = 0.7, position = position_dodge(width = 0.8), 
               color = "black",      # Border color of the boxplot
               alpha = 0.9,          # Transparency of the boxplot
               outlier.shape = NA) + # Remove outliers
  labs(y = "TMBIM6 expression value", x = "GSE4290") +
  geom_jitter(aes(group = group), width = 0.3, alpha = 0.5, color = "orange") + # Add jittered points# Outlier color
  scale_fill_manual(values = c("GBM" = "red", "LGG" = "darkred", "normal" = "blue"))  # Specify uppercase for

# Add significance level using ggpubr
p_gene <- p_gene + stat_compare_means(
  comparisons = list(
    c("GBM", "normal"),  # Comparison between GBM and normal
    c("LGG", "normal")   # Comparison between LGG and normal
  ),
  method = "wilcox.test",  # Use Mann-Whitney U test
  label = "p.format",
  position = position_dodge(width = 0.8),
  show.legend = TRUE,    # Show the legend for significance levels   
)  # Adjust the position of labels

# Print boxplot for the specific gene with significance levels
print(p_gene)




getwd()






