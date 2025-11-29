#------------------------------------
# xCell2 Cell Type Enrichment Analysis
#------------------------------------

# Load necessary R packages
library(xCell2)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(dplyr)

dir.create("./tissue1/Resultados/XCELL2")
setwd("./tissue1/Resultados/XCELL2")

#------------------------------------------
# 1. Data Preparation
#------------------------------------------

# Connect to the Ensembl database for human genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch gene information, including ID, positions, and symbol
gene_info <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position", "hgnc_symbol"),
                   mart = mart)

# Calculate gene lengths and name the list with their Ensembl IDs
gene_lengths <- abs(gene_info$end_position - gene_info$start_position)
names(gene_lengths) <- gene_info$ensembl_gene_id

# Load the raw count matrix
counts <- read.csv("Matriz_conteo_align.csv", header = TRUE, row.names = "GeneID")
counts
# Remove the "Length" column
counts$Length <- NULL

# Remove the version number from Ensembl IDs in the row names
counts$EnsemblID <- sub("\\..*", "", rownames(counts))

# Map Ensembl IDs to gene symbols
matched_symbols <- gene_info$hgnc_symbol[match(counts$EnsemblID, gene_info$ensembl_gene_id)]
counts$Symbol <- matched_symbols

# Filter to get only genes with unique and valid symbols
counts_unique <- counts[!is.na(counts$Symbol) & counts$Symbol != "", ]
counts_unique <- counts_unique[!duplicated(counts_unique$Symbol), ]

# Set row names to the gene symbols
rownames(counts_unique) <- counts_unique$Symbol

# Remove temporary columns
counts_unique$EnsemblID <- NULL
counts_unique$Symbol <- NULL

# Create a gene length vector aligned with the count matrix
ensembl_ids_in_counts <- gene_info$ensembl_gene_id[match(rownames(counts_unique), gene_info$hgnc_symbol)]
gene_lengths_filtered <- gene_lengths[ensembl_ids_in_counts]

# Remove genes without a corresponding length
na_genes <- is.na(gene_lengths_filtered)
if (any(na_genes)) {
  counts_unique <- counts_unique[!na_genes, ]
  gene_lengths_filtered <- gene_lengths_filtered[!na_genes]
}

#------------------------------------------
# 2. TPM Normalization
#------------------------------------------

# Function for TPM normalization
counts_to_tpm <- function(counts, gene_lengths) {
  rpk <- counts / gene_lengths
  per_million_scaling_factor <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, per_million_scaling_factor, FUN = "/")
  return(tpm)
}

# Perform TPM normalization
tpm_matrix <- counts_to_tpm(counts_unique, gene_lengths_filtered)

#------------------------------------------
# 3. xCell2 Analysis
#------------------------------------------
#---------------------------------------------------|

#  BlueprintEncode.xCell2Ref        | 43  CellTypes |

#  ImmGenData.xCell2Ref             | 19  CellTypes |

#  ImmuneCompendium.xCell2Ref       | 40  CellTypes |

#  TMECompendium.xCell2Ref          | 25  CellTypes |

#  PanCancer.xCell2Ref              | 29  CellTypes |

#  TabulaSapiensBlood.xCell2Ref     | 21  CellTypes |

#  LM22.xCell2Ref                   | 22  CellTypes |

#---------------------------------------------------|

#  TabulaMurisBlood.xCell2Ref       | 6   CellTypes |

#  MouseRNAseqData.xCell2Ref        | 18  CellTypes |

#  DICE_demo.xCell2Ref              | x   CellTypes |

#---------------------------------------------------|


# Load the xCell2 reference gene set
data("ImmuneCompendium.xCell2Ref", package = "xCell2")
data("BlueprintEncode.xCell2Ref", package = "xCell2")
data("ImmGenData.xCell2Ref", package = "xCell2")
data("TMECompendium.xCell2Ref", package = "xCell2")
data("PanCancer.xCell2Ref", package = "xCell2")
data("LM22.xCell2Ref", package = "xCell2")
data("TabulaSapiensBlood.xCell2Ref", package = "xCell2")
# Check gene overlap
genes_ref <- getGenesUsed(TMECompendium.xCell2Ref)
genes_mix <- rownames(tpm_matrix)
overlap_percentage <- round(length(intersect(genes_mix, genes_ref)) / length(genes_ref) * 100, 2)
print(paste0("Overlap between your genes and the reference is: ", overlap_percentage, "%"))

# Perform cell type enrichment analysis
xcell2_results <- xCell2Analysis(
  mix = tpm_matrix,
  xcell2object = TMECompendium.xCell2Ref,
  minSharedGenes = 0.8
)

# Save the results to a CSV file
write.csv(xcell2_results, file = "xCell2_enrichment_scores_TMECompendium.csv", row.names = TRUE)

# Display the structure and first few rows of the results
str(xcell2_results)
head(xcell2_results)




#------------------------------------------
# 4. Visualization
#------------------------------------------
# Convert results to a data frame and add cell type as a column
xcell2_df <- as.data.frame(xcell2_results)
xcell2_df$CellType <- rownames(xcell2_df)

# Convert to long format for ggplot2
xcell2_long <- pivot_longer(xcell2_df,
                            cols = -CellType,
                            names_to = "Sample",
                            values_to = "Score")

# Cargar la información del experimento para asignar grupos
experiment_info <- read.csv("../../experimento.csv", header = TRUE, row.names = 1)

# Asegurar que los nombres de muestra coincidan
if (!all(xcell2_long$Sample %in% rownames(experiment_info))) {
  stop("Algunas muestras en xCell2 no están presentes en el archivo experimento.csv")
}

# Asignar grupo desde el archivo de experimento
xcell2_long$Group <- experiment_info[xcell2_long$Sample, "Condition"]

# Calculate the mean score for each cell type per group
grouped_scores <- xcell2_long %>%
  group_by(Group, CellType) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE))
grouped_scores
# Generate the final stacked bar chart with only two bars
ggplot(grouped_scores, aes(x = Group, y = MeanScore, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Stacked Bar Chart of Mean Cell Type Scores by Group",
       x = "Group",
       y = "Mean Enrichment Score",
       fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))



#--------------------------------------------
# 5. STATISTICAL VALIDATION
#--------------------------------------------

# 5. Statistical Validation

# Load necessary libraries for data manipulation and statistical tests
library(dplyr)

# Get the list of unique cell types
cell_types <- unique(xcell2_long$CellType)

# Initialize an empty data frame to store results
p_values_df <- data.frame(
  CellType = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each cell type and perform a t-test
for (ct in cell_types) {
  # Filter data for the current cell type
  subset_data <- xcell2_long %>%
    filter(CellType == ct)
  
  # Check if both groups have data for the cell type
  if (length(unique(subset_data$Group)) == 2) {
    # Perform a two-sample t-test
    ttest_result <- t.test(Score ~ Group, data = subset_data)
    
    # Store the result
    p_values_df <- rbind(p_values_df, data.frame(
      CellType = ct,
      p_value = ttest_result$p.value
    ))
  } else {
    # If one group is missing data, record an NA p-value
    p_values_df <- rbind(p_values_df, data.frame(
      CellType = ct,
      p_value = NA
    ))
  }
}

# Adjust p-values for multiple comparisons using the Benjamini-Hochberg method
p_values_df$q_value <- p.adjust(p_values_df$p_value, method = "BH")

# Order the results by adjusted p-value (q-value)
p_values_df <- p_values_df %>%
  arrange(q_value)

# Display the final results, showing which cell types are significantly different
# A q-value < 0.05 indicates a significant difference between groups
print(p_values_df)

# To visualize the results, you can filter for significant cell types
significant_cell_types <- p_values_df %>%
  filter(q_value < 0.05) %>%
  pull(CellType)

print("Significantly different cell types (q-value < 0.05):")
print(significant_cell_types)



#------------------------------------------
# 6. Box Plot for a Specific Cell Type with Significance
#------------------------------------------

# Set the cell type you want to visualize
# Change "neutrophil" to any cell type from the LM22.xCell2Ref data set
target_cell_type <- "neutrophil"

# Filter the data for the target cell type
plot_data <- xcell2_long %>%
  filter(CellType == target_cell_type)

# Get the q-value for the target cell type
target_q_value <- p_values_df %>%
  filter(CellType == target_cell_type) %>%
  pull(q_value)

# Format the p-value for the plot title
formatted_p_value <- ifelse(target_q_value < 0.001, "p < 0.001", paste0("p = ", round(target_q_value, 4)))

# Create the box plot
ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  theme_minimal() +
  labs(title = paste("Enrichment Score for", target_cell_type),
       subtitle = paste("Significance:", formatted_p_value),
       x = "Group",
       y = "Enrichment Score",
       fill = "Group") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Prostate Cancer" = "coral1", "Healthy" = "skyblue"))

#############################################################################################







#################################################################
##
##  xCell2 Cell Type Enrichment Analysis Workflow for Tissue 1
##
#################################################################


# ===============================================================
# SECTION 0: SETUP & LIBRARIES
# ===============================================================

# --- Load all necessary R packages ---
library(xCell2)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(dplyr)

# --- Set Working Directory ---
# Creates the directory if it doesn't exist and sets it as the working directory.
dir.create("./tissue1/Resultados/XCELL2", showWarnings = FALSE, recursive = TRUE)
setwd("./tissue1/Resultados/XCELL2")


# ===============================================================
# SECTION 1: DATA PREPARATION & GENE ANNOTATION
# ===============================================================

# --- Connect to Ensembl to fetch gene information ---
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position", "hgnc_symbol"),
  mart = mart
)

# --- Calculate gene lengths ---
gene_lengths <- abs(gene_info$end_position - gene_info$start_position)
names(gene_lengths) <- gene_info$ensembl_gene_id

# --- Load and process the raw count matrix ---
# Ensure the CSV file for Tissue 1 is in the correct path.
counts <- read.csv("Matriz_conteo_align.csv", header = TRUE, row.names = "GeneID")
counts$Length <- NULL # Remove the 'Length' column if it exists.

# --- Map Ensembl IDs to Gene Symbols ---
# Remove version numbers from Ensembl IDs (e.g., ENSG00000000003.14 -> ENSG00000000003).
counts$EnsemblID <- sub("\\..*", "", rownames(counts))
counts$Symbol <- gene_info$hgnc_symbol[match(counts$EnsemblID, gene_info$ensembl_gene_id)]

# --- Filter for unique and valid gene symbols ---
counts_unique <- counts[!is.na(counts$Symbol) & counts$Symbol != "", ]
counts_unique <- counts_unique[!duplicated(counts_unique$Symbol), ]
rownames(counts_unique) <- counts_unique$Symbol
counts_unique$EnsemblID <- NULL
counts_unique$Symbol <- NULL

# --- Align gene lengths with the final count matrix ---
ensembl_ids_in_counts <- gene_info$ensembl_gene_id[match(rownames(counts_unique), gene_info$hgnc_symbol)]
gene_lengths_filtered <- gene_lengths[ensembl_ids_in_counts]

# --- Remove genes that do not have a corresponding length ---
na_genes <- is.na(gene_lengths_filtered)
if (any(na_genes)) {
  counts_unique <- counts_unique[!na_genes, ]
  gene_lengths_filtered <- gene_lengths_filtered[!na_genes]
}


# ===============================================================
# SECTION 2: TPM NORMALIZATION
# ===============================================================

# --- Define the TPM normalization function ---
counts_to_tpm <- function(counts, gene_lengths) {
  # Divide by gene length to get Reads Per Kilobase (RPK)
  rpk <- counts / gene_lengths
  # Calculate the "per million" scaling factor for each sample
  per_million_scaling_factor <- colSums(rpk) / 1e6
  # Divide RPK by the scaling factor to get TPM
  tpm <- sweep(rpk, 2, per_million_scaling_factor, FUN = "/")
  return(tpm)
}

# --- Perform TPM normalization on the count data ---
tpm_matrix <- counts_to_tpm(counts_unique, gene_lengths_filtered)


# ===============================================================
# SECTION 3: XCELL2 ANALYSIS
# ===============================================================

# --- Available xCell2 Reference Databases ---
#---------------------------------------------------|
#  BlueprintEncode.xCell2Ref       | 43  CellTypes |
#  ImmGenData.xCell2Ref            | 19  CellTypes |
#  ImmuneCompendium.xCell2Ref      | 40  CellTypes |
#  TMECompendium.xCell2Ref         | 25  CellTypes |
#  PanCancer.xCell2Ref             | 25  CellTypes |
#  TabulaSapiensBlood.xCell2Ref    | 21  CellTypes |
#  LM22.xCell2Ref                  | 22  CellTypes |
#---------------------------------------------------|
#  TabulaMurisBlood.xCell2Ref      | 6   CellTypes |
#  MouseRNAseqData.xCell2Ref       | 18  CellTypes |
#---------------------------------------------------|

# --- Load the desired reference gene sets ---
data("BlueprintEncode.xCell2Ref", package = "xCell2")
data("ImmGenData.xCell2Ref", package = "xCell2")
data("ImmuneCompendium.xCell2Ref", package = "xCell2")
data("TMECompendium.xCell2Ref", package = "xCell2")
data("PanCancer.xCell2Ref", package = "xCell2")
data("TabulaSapiensBlood.xCell2Ref", package = "xCell2")
data("LM22.xCell2Ref", package = "xCell2")
# ... add other 'data()' calls here if needed ...

# --- Select your reference object for the analysis ---
# Choose one of the loaded references, for example, TMECompendium.xCell2Ref.
reference_object <- LM22.xCell2Ref

# --- Check gene overlap with the selected reference ---
genes_ref <- getGenesUsed(reference_object)
overlap_percentage <- round(length(intersect(rownames(tpm_matrix), genes_ref)) / length(genes_ref) * 100, 2)
print(paste0("Gene overlap with reference is: ", overlap_percentage, "%"))

# --- Perform the cell type enrichment analysis ---
xcell2_results <- xCell2Analysis(
  mix = tpm_matrix,
  xcell2object = reference_object,
  minSharedGenes = 0.8 # Minimum gene overlap threshold.
)

# --- Save and display the raw enrichment scores ---
write.csv(xcell2_results, file = "xCell2_enrichment_scores_LM22.csv", row.names = TRUE)
print("xCell2 analysis complete. First few results:")
print(head(xcell2_results))


# ===============================================================
# SECTION 4: GROUP ASSIGNMENT & DATA RESHAPING
# ===============================================================

# --- Convert the results matrix to a long-format data frame for plotting ---
xcell2_df <- as.data.frame(xcell2_results)
xcell2_df$CellType <- rownames(xcell2_df)
xcell2_long <- pivot_longer(
  xcell2_df,
  cols = -CellType,
  names_to = "Sample",
  values_to = "Score"
)

# --- Load experiment metadata and assign groups to samples ---
# This uses the specific metadata file for Tissue 1.
experiment_info <- read.csv("../../experimento.csv", header = TRUE, row.names = 1)

# --- Verify that all sample names match between the results and the metadata ---
# Note: This version assumes sample names already match perfectly.
if (!all(xcell2_long$Sample %in% rownames(experiment_info))) {
  stop("Critical Error: Sample names do not match between results and metadata file.")
}

# --- Assign the 'Group' column based on the metadata ---
xcell2_long$Group <- experiment_info[xcell2_long$Sample, "Condition"]


# ===============================================================
# SECTION 5: STATISTICAL VALIDATION
# ===============================================================

# --- Initialize a data frame to store statistical results ---
p_values_df <- data.frame(
  CellType = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# --- Loop through each cell type and perform a two-sample t-test ---
for (ct in unique(xcell2_long$CellType)) {
  subset_data <- xcell2_long %>% filter(CellType == ct)
  
  # Check if both groups ('Case' and 'Control') have data for the test.
  if (length(unique(subset_data$Group)) == 2) {
    ttest_result <- t.test(Score ~ Group, data = subset_data)
    p_values_df <- rbind(p_values_df, data.frame(
      CellType = ct,
      p_value = ttest_result$p.value
    ))
  } else {
    # If only one group is present, record NA as the p-value.
    p_values_df <- rbind(p_values_df, data.frame(
      CellType = ct,
      p_value = NA
    ))
  }
}

# --- Adjust p-values for multiple comparisons (Benjamini-Hochberg method) ---
# This is essential to control for false positives when running many tests.
p_values_df$q_value <- p.adjust(p_values_df$p_value, method = "BH")

# --- Order results by the adjusted p-value (q-value) ---
p_values_df <- p_values_df %>% arrange(q_value)

# --- Display the final statistical results ---
print("Statistical Validation Results (t-test):")
print(p_values_df)

# --- Filter and display the significantly different cell types ---
significant_cell_types <- p_values_df %>%
  filter(q_value < 0.05) %>%
  pull(CellType)

print("Significantly different cell types (q-value < 0.05):")
if (length(significant_cell_types) > 0) {
  print(significant_cell_types)
} else {
  print("No significant differences found.")
}


# ===============================================================
# SECTION 6: VISUALIZATION
# ===============================================================

# --- Stacked Bar Chart of Mean Scores ---
grouped_scores <- xcell2_long %>%
  group_by(Group, CellType) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = 'drop')

ggplot(grouped_scores, aes(x = Group, y = MeanScore, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Mean Cell Type Scores by Group",
    x = "Group",
    y = "Mean Enrichment Score",
    fill = "Cell Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# --- Box Plot for a Specific Cell Type ---
# You can change this to any cell type name from your results.
# For example, use one from the 'significant_cell_types' list.
target_cell_type <- "neutrophil"

# Filter data for the target cell type
plot_data <- xcell2_long %>%
  filter(CellType == target_cell_type)

# Get the adjusted p-value for the plot title
target_q_value <- p_values_df %>%
  filter(CellType == target_cell_type) %>%
  pull(q_value)

# Format the q-value for display.
formatted_q_value <- ifelse(is.na(target_q_value), "NA",
                            ifelse(target_q_value < 0.001, 
                                   "q < 0.001", 
                                   paste0("q = ", round(target_q_value, 4))
                            )
)

# Create the box plot
ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  theme_minimal() +
  labs(
    title = paste("Enrichment Score for", target_cell_type),
    subtitle = paste("Significance (BH-adjusted):", formatted_q_value),
    x = "Group",
    y = "Enrichment Score"
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Case" = "coral1", "Control" = "skyblue"))