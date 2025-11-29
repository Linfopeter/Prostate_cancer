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

dir.create("./pbmcs/Resultados/XCELL2")
setwd("./pbmcs/Resultados/XCELL2")

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

#  TabulaSapiensBlood.xCell2Ref     | 18  CellTypes |

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
genes_ref <- getGenesUsed(LM22.xCell2Ref)
genes_mix <- rownames(tpm_matrix)
overlap_percentage <- round(length(intersect(genes_mix, genes_ref)) / length(genes_ref) * 100, 2)
print(paste0("Overlap between your genes and the reference is: ", overlap_percentage, "%"))

# Perform cell type enrichment analysis
xcell2_results <- xCell2Analysis(
  mix = tpm_matrix,
  xcell2object = LM22.xCell2Ref,
  minSharedGenes = 0.9
)

# Save the results to a CSV file
write.csv(xcell2_results, file = "xCell2_enrichment_scores_LM22.csv", row.names = TRUE)

# Display the structure and first few rows of the results
str(xcell2_results)
head(xcell2_results)




#------------------------------------------
# 4. Visualization
#------------------------------------------
# Convert results to a data frame
xcell2_df <- as.data.frame(xcell2_results)

# Ensure no existing 'CellType' column to avoid duplication
if ("CellType" %in% colnames(xcell2_df)) {
  xcell2_df$CellType <- NULL  # Remove if exists
}

library(tidyverse)

# Convert row names (cell types) into a proper column, then pivot to long format
xcell2_long <- xcell2_df %>%
  tibble::rownames_to_column(var = "CellType") %>%   # Now safe: adds CellType from row names
  pivot_longer(
    cols = -CellType,
    names_to = "Sample",
    values_to = "Score"
  )

# Assign group: Healthy (non-P) vs Prostate Cancer (starts with P)
xcell2_long <- xcell2_long %>%
  mutate(Group = ifelse(grepl("^P", Sample), "Prostate Cancer", "Healthy"))

# Calculate mean score per cell type per group
grouped_scores <- xcell2_long %>%
  group_by(Group, CellType) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = "drop")

# Optional: Inspect the results
print(grouped_scores)
# View all unique cell types:
unique(grouped_scores$CellType)

# Define the custom color vector in the same order as your 22 cell types
custom_colors <- c(
  "alternatively activated macrophage" = "#50FA7B",
  "CD8-positive, alpha-beta T cell" = "#33A02C",
  "Dendritic cells activated" = "#A6CEE3",
  "Dendritic cells resting" = "#1F78B4",
  "eosinophil" = "#E31A1C",
  "gamma-delta T cell" = "#FB9A99",
  "inflammatory macrophage" = "#B8860B",
  "macrophage" = "#A52A2A",
  "Mast cells activated" = "#DDA0DD",
  "Mast cells resting" = "#6A3D9A",
  "memory B cell" = "#FF9900",
  "monocyte" = "#B15928",
  "naive B cell" = "#8DD3C7",
  "naive thymus-derived CD4-positive, alpha-beta T cell" = "#7F7F7F",
  "neutrophil" = "#0066FF",
  "NK cells activated" = "#0088A8",
  "NK cells resting" = "#6495ED", ##6495ED
  "plasma cell" = "#FFD700",
  "regulatory T cell" = "#4DAF4A",
  "T cells CD4 memory activated" = "#F08080",
  "T cells CD4 memory resting" = "#9370DB",
  "T follicular helper cell" = "#000080"
)
#Color palette
#c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#8DD3C7", "#FFFF33", "#333333")
#c("#FF3333", "#FF9900", "#CCFF00", "#33CC00", "#00CCCC", "#0066FF", "#6600FF", "#CC00FF", "#FF0066", "#FF6633")
#c("#8B4000", "#D2691E", "#B8860B", "#556B2F", "#2E8B57", "#006400", "#808000", "#A0522D", "#CD853F", "#F4A460")
#c("#1E90FF", "#4169E1", "#0000CD", "#000080", "#4682B4", "#6495ED", "#7B68EE", "#9370DB", "#8A2BE2", "#9932CC")
#c("#FFB6C1", "#FFA07A", "#FFE4B5", "#7FFFD4", "#E0FFFF", "#D8BFD8", "#DDA0DD", "#F0E68C", "#B0E0E6", "#98FB98")
#c("#FF5555", "#50FA7B", "#8BE9FD", "#BD93F9", "#FFB86C", "#F1FA8C", "#FF79C6", "#50FA7B", "#FF55BB", "#A89984")
#c("#D62728", "#1F77B4", "#2CA02C", "#FF7F0E", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
#c("#FF4500", "#FF6347", "#FF7F50", "#FFA07A", "#F08080", "#CD5C5C", "#DC143C", "#B22222", "#8B0000", "#A52A2A")
#c("#003366", "#005A8C", "#0088A8", "#00B2A9", "#00D3A6", "#00E096", "#38EC8D", "#7AF588", "#C1FA82", "#FFFF7C")
#c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FF9900", "#9933FF", "#33CCFF", "#FF6666")
#
X11()
# Generate the final stacked bar chart with custom colors
ggplot(grouped_scores, aes(x = Group, y = MeanScore, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "", #"Stacked Bar Chart of Mean Cell Type Scores by Group",
    x = "Group",
    y = "Mean Enrichment Score",
    fill = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 5, hjust = 1, size = 10),
    legend.position = "right"
  ) +
  scale_fill_manual(values = custom_colors)  # Use custom colors

# ------------------------------------------------------------------------------
# Guardar el gráfico en un archivo PNG de alta resolución (no formateado)
# ------------------------------------------------------------------------------

# Si no asignaste el ggplot a una variable (ej. final_plot), usa last_plot()
output_filename <- file.path("XCELL2_pbmcs.png")

ggsave(
  filename = output_filename,
  plot = last_plot(),    # o usa 'final_plot' si lo guardaste en una variable
  device = "png",
  width = 12, height = 10, units = "in",
  dpi = 300,
  bg = "white"
)

message("Gráfico guardado en: ", output_filename)


# ------------------------------------------------------------------------------
# Ajustes para texto más grande y leyenda en una sola columna (formateado)
# ------------------------------------------------------------------------------

final_plot <- ggplot(grouped_scores, aes(x = Group, y = MeanScore, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 18) +  # Aumenta el tamaño base de texto
  labs(
    x = "Group",
    y = "Mean Enrichment Score",
    fill = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 5, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key.size = unit(1.1, "lines"),
    legend.position = "right",
    legend.box = "vertical",      # mantiene la leyenda en una sola columna
    legend.spacing.y = unit(0.3, "cm")
  ) +
  guides(
    fill = guide_legend(ncol = 1) # fuerza una sola columna
  ) +
  scale_fill_manual(values = custom_colors)

# ------------------------------------------------------------------------------
# Guardar gráfico en PNG con alta resolución
# ------------------------------------------------------------------------------
output_filename <- file.path("XCELL2_CORRECTED.png")

ggsave(
  filename = output_filename,
  plot = final_plot,
  device = "png",
  width = 14, height = 10, units = "in",
  dpi = 400, bg = "white"
)

message("Gráfico guardado en: ", output_filename)

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
x11()
# Create the box plot
ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  theme_minimal() +
  labs(title = paste("Enrichment Score for", target_cell_type),
       #subtitle = paste("Significance:", formatted_p_value),
       x = "Group",
       y = "Enrichment Score",
       fill = "Group") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("Prostate Cancer" = "coral1", "Healthy" = "skyblue"))



# ------------------------------------------------------------------------------
# Boxplot con texto más grande y colores personalizados
# ------------------------------------------------------------------------------

box_plot <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.size = 2, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 1.8) +
  scale_fill_manual(values = c("Prostate Cancer" = "coral1", "Healthy" = "skyblue")) +
  theme_minimal(base_size = 18) +
  labs(
    title = paste("Enrichment Score for", target_cell_type),
    # subtitle = paste("Significance:", formatted_p_value),
    x = "Group",
    y = "Enrichment Score",
    fill = "Group"
  ) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none",   # sin leyenda, ya que los grupos están en el eje X
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# Guardar el gráfico en un archivo PNG de alta resolución
# ------------------------------------------------------------------------------

# Asegurarse de que la carpeta de salida existe
out_dir <- "figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

output_filename <- file.path(out_dir, paste0("Boxplot_Enrichment_", target_cell_type, ".png"))

ggsave(
  filename = output_filename,
  plot = box_plot,
  device = "png",
  width = 10, height = 8, units = "in",
  dpi = 400, bg = "white"
)

message("✅ Gráfico de boxplot guardado en: ", output_filename)
