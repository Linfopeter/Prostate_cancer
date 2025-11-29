#################################################################
##
##  xCell2 Cell Type Enrichment Analysis Workflow
##
#################################################################


# ===============================================================
# SECTION 0: SETUP
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
# dir.create("./tissue1/Resultados/XCELL2", showWarnings = FALSE, recursive = TRUE)
# setwd("./tissue1/Resultados/XCELL2")


# ===============================================================
# SECTION 1: DATA PREPARATION
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
# Ensure the CSV file is in your working directory.
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
#  DICE_demo.xCell2Ref             | x   CellTypes |
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
# Make sure the path to your experiment file is correct.
experiment_info <- read.csv("experimento.csv", header = TRUE, row.names = 1)

# --- Clean sample names from both sources to ensure they match ---
# This step is crucial and handles potential formatting issues like extra spaces.
samples_for_lookup <- trimws(xcell2_long$Sample)
rownames(experiment_info) <- trimws(rownames(experiment_info))

# --- Verify that all sample names match between the results and the metadata ---
if (!all(samples_for_lookup %in% rownames(experiment_info))) {
  stop("Critical Error: Sample names do not match between results and metadata file after cleaning.")
}

# --- Assign the 'Group' column based on the metadata ---
xcell2_long$Group <- experiment_info[samples_for_lookup, "Condition"]


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
formatted_q_value <- ifelse(target_q_value < 0.001, 
                            "q < 0.001", 
                            paste0("q = ", round(target_q_value, 4)))

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
