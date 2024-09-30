# Load necessary libraries

library(Biobase)
library(GEOquery)
library(limma)

# Load the GEO data from your local file

gse <- getGEO(filename="file.txt")

# Extract expression data and phenotype data

exprs_data <- exprs(gse)
pdata <- pData(phenoData(gse))

# Check unique values in 'status:ch1'

unique(pdata$status:ch1)

# Create a 'condition' column based on 'status:ch1'

pdata$condition <- factor(pdata$status:ch1, levels = c("CTL", "MCI", "AD"))

# Design matrix for limma (comparing AD and MCI vs Control)

design <- model.matrix(~ 0 + pdata$condition)
colnames(design) <- levels(pdata$condition)

# Fit the linear model

fit <- lmFit(exprs_data, design)

# Define contrasts for AD vs Control and MCI vs Control

contrast_matrix <- makeContrasts(
  ADvsControl = AD - CTL,
  MCIvsControl = MCI - CTL,
  levels = design
)

# Fit the contrasts

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract differentially expressed genes (DEGs) for AD vs Control

top_table_AD <- topTable(fit2, coef = "ADvsControl", adjust = "fdr", sort.by = "P", number = Inf)

# Extract differentially expressed genes (DEGs) for MCI vs Control

top_table_MCI <- topTable(fit2, coef = "MCIvsControl", adjust = "fdr", sort.by = "P", number = Inf)

head(top_table_AD)

head(top_table_MCI)
