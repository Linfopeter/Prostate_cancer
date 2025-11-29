#-----------------------------------|
# IMMUNE DECONVOLUTION "GRANULATOR" |
#-----------------------------------|
# 1. Install and load the granulator package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("granulator")

library(granulator)
library(biomaRt)
setwd("../../")
#DIRECTORI ASIGNEMENT
#dir.create("./pbmcs/Resultados/granulator")
setwd("./pbmcs/Resultados/granulator")


# 2. Load the raw count data (ADD MANUALLY BEFORE)
counts <- read.csv("Matriz_conteo_align.csv", header = TRUE, row.names = "GeneID")
counts
# Remove the "Length" column
counts$Length <- NULL

# 3. Filtrar genes con baja expresión
# Mantener genes con conteo >= 10 en al menos el 20% de las muestras
min_samples <- ceiling(0.2 * ncol(counts))
counts_filtered <- counts[rowSums(counts >= 10) >= min_samples, ]


# Connect to the Ensembl database for human genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch gene information, including ID, positions, and symbol
gene_info <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position", "hgnc_symbol"),
                   mart = mart)

# Calculate gene lengths and name the list with their Ensembl IDs
gene_lengths <- abs(gene_info$end_position - gene_info$start_position)
names(gene_lengths) <- gene_info$ensembl_gene_id


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


# 4. Convert raw counts to TPM (Transcripts Per Million)
# First, convert the filtered counts data frame to a matrix
counts_matrix <- as.matrix(counts_unique)

# Now use the `get_TPM` function with the matrix and gene lengths
tpm_matrix <- get_TPM(counts = counts_matrix, effLen = gene_lengths_filtered)

# 5. Load and align the reference signature matrix
# It's important to ensure your genes are also in the signature matrix.
# We will use the `granulator` demo data as an example.
load_ABIS()
signature_matrix <- sigMatrix_ABIS_S0

# Filter both matrices to keep only genes they have in common
common_genes <- intersect(rownames(tpm_matrix), rownames(signature_matrix))
tpm_matrix_filtered <- tpm_matrix[common_genes, ]
signature_matrix_filtered <- signature_matrix[common_genes, ]

# 6. Perform deconvolution
# `granulator`'s `deconvolute` function requires the `sigMatrix` to be a named list.
deconvolution_results <- deconvolute(
  m = tpm_matrix_filtered,
  sigMatrix = list("S0" = signature_matrix_filtered),
  methods = get_decon_methods(), # Runs all available methods
)

# 7. Access and analyze the results
# The output is a list containing estimated coefficients and proportions.
# The proportions are what you typically use for downstream analysis.
head(deconvolution_results$proportions)

# 8. Visualize the deconvolution results
# `granulator` has built-in plotting functions to visualize the output.
# We will use "nnls" as an example method and "S0" as the signature.
plot_proportions(
  deconvoluted = deconvolution_results,
  #method = "ols",
  #method = "nnls", #Pick the desired method by uncommenting
  #method = "qprog",
  #method = "qprogwc", # <-
  method = "dtangle", # <-
  #method = "rls",
  #method = "svr",
  signature = "S0"
)

# You can also use `plot_deconvolute` to compare different methods.
plot_deconvolute(deconvoluted = deconvolution_results)

# 9. Save the results
# It's a good practice to save the proportions as a CSV file.
write.csv(deconvolution_results$proportions, "deconvolution_proportions.csv", row.names = FALSE)

#




#CHECK FOR:
deconvolution_results
# 1. Prepare the data for statistical analysis

# Correctly extract the desired data frame from the list
# The 'deconvolution_results$proportions' list contains multiple data frames.
# You need to specify which one you want to use for analysis.
# We'll use the 'dtangle_S0' method as an example.
proportions_df <- deconvolution_results$proportions$dtangle_S0 #CHANGE DESIRED ALGORITHM

# Load the required packages
library(tidyr)
library(dplyr)
library(ggplot2)

# Pivot the data to a long format suitable for plotting and testing
# The 'cell_type' column needs to be created from the row names
proportions_df <- proportions_df %>%
  tibble::rownames_to_column(var = "SampleID")

proportions_long <- proportions_df %>%
  pivot_longer(
    cols = -SampleID, # This is the corrected line. We now exclude the SampleID column
    names_to = "cell_type",
    values_to = "Proportion"
  )

# Add a new 'Group' column based on the SampleID prefix
proportions_long <- proportions_long %>%
  mutate(
    Group = case_when(
      startsWith(SampleID, "C") | startsWith(SampleID, "S") ~ "Healthy",
      startsWith(SampleID, "P") ~ "Cancer",
      TRUE ~ "Unknown" # For any samples that don't match
    )
  )


proportions_long$Group <- factor(proportions_long$Group, levels = c("Healthy", "Cancer", "Unknown"))


# 2. Perform statistical tests (Wilcoxon rank-sum test)
cell_types <- unique(proportions_long$cell_type)
p_values <- list()

for (cell in cell_types) {
  subset_data <- proportions_long %>%
    filter(cell_type == cell)
  
  if (n_distinct(subset_data$Group) == 2 && all(table(subset_data$Group) >= 3)) {
    test_result <- wilcox.test(Proportion ~ Group, data = subset_data)
    p_values[[cell]] <- test_result$p.value
  } else {
    p_values[[cell]] <- NA
  }
}

p_values_df <- data.frame(
  cell_type = names(p_values),
  p_value = unlist(p_values)
)

# 3. Adjust p-values for multiple comparisons (Benjamini-Hochberg method)
p_values_df <- p_values_df %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
  arrange(adjusted_p_value)

print("Statistical results:")
print(p_values_df)

# 4. Visualize the most significant differences with a box plot
top_cell_types <- p_values_df %>%
  arrange(adjusted_p_value) %>%
  top_n(5, -adjusted_p_value) %>%
  pull(cell_type)

proportions_to_plot_box <- proportions_long %>%
  filter(cell_type %in% top_cell_types)

# Verificar si hay tipos celulares significativos antes de graficar
if (length(top_cell_types) > 0) {
  proportions_to_plot_box <- proportions_long %>%
    filter(cell_type %in% top_cell_types)
  
  ggplot(proportions_to_plot_box, aes(x = Group, y = Proportion, fill = Group)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(
      title = "Significant Immune Cell Proportions: Healthy vs. Cancer",
      y = "Estimated Cell Proportion",
      x = "Group"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
} else {
  message("No significant cell types found for plotting.")
}

# 5. Visualize all proportions with a stacked bar plot
healthy_cancer_proportions_bar <- proportions_long %>%
  filter(Group %in% c("Healthy", "Cancer"))

ggplot(healthy_cancer_proportions_bar, aes(x = SampleID, y = Proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  labs(
    title = "Overall Immune Cell Proportions: Healthy vs. Cancer",
    y = "Estimated Cell Proportion",
    x = "Sample ID",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#
# TWO GROUPS BAR CHART
#

# Calculate the mean proportion for each cell type in each group
pooled_proportions <- proportions_long %>%
  group_by(Group, cell_type) %>%
  summarise(
    MeanProportion = mean(Proportion),
    .groups = 'drop'
  ) %>%
  filter(Group != "Unknown") # Exclude the 'Unknown' group if it exists

# Create the stacked bar plot with only two bars
ggplot(pooled_proportions, aes(x = Group, y = MeanProportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Pooled Immune Cell Proportions: Healthy vs. Cancer",
    y = "Mean Estimated Cell Proportion",
    x = "Group",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Center the group labels
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )









##



#-----------------------------------------------------|
# IMMUNE DECONVOLUTION "GRANULATOR"  > immunoStates < |
#-----------------------------------------------------|
# 1. Install and load the granulator package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("granulator")

library(granulator)
library(biomaRt)

#DIRECTORI ASIGNEMENT
dir.create("./pbmcs/Resultados/granulator")
setwd("../")
setwd("./pbmcs/Resultados/granulator")


# 2. Load the raw count data (ADD MANUALLY BEFORE)
counts <- read.csv("Matriz_conteo_align.csv", header = TRUE, row.names = "GeneID")
counts
# Remove the "Length" column
counts$Length <- NULL

# Filtrar genes con baja expresión: mantener genes con conteo >= 10 en al menos el 20% de las muestras
min_samples <- ceiling(0.2 * ncol(counts))
counts <- counts[rowSums(counts >= 10) >= min_samples, ]

# Connect to the Ensembl database for human genes
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Fetch gene information, including ID, positions, and symbol
gene_info <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position", "hgnc_symbol"),
                   mart = mart)

# Calculate gene lengths and name the list with their Ensembl IDs
gene_lengths <- abs(gene_info$end_position - gene_info$start_position)
names(gene_lengths) <- gene_info$ensembl_gene_id


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


# 4. Convert raw counts to TPM (Transcripts Per Million)
# First, convert the filtered counts data frame to a matrix
counts_matrix <- as.matrix(counts_unique)

# Now use the `get_TPM` function with the matrix and gene lengths
tpm_matrix <- get_TPM(counts = counts_matrix, effLen = gene_lengths_filtered)

# Leer el archivo sin usar la primera fila como encabezado
raw_matrix <- read.csv("immunoStates.csv", header = FALSE, stringsAsFactors = FALSE)

# Usar la tercera fila como nombres de columna
colnames(raw_matrix) <- as.character(unlist(raw_matrix[3, ]))

# Eliminar las tres primeras filas (metadatos o encabezados duplicados)
signature_matrix <- raw_matrix[-c(1, 2, 3), ]

# Extraer los nombres de los genes antes de convertir a numérico
gene_names <- signature_matrix[, 1]

# Eliminar la columna de genes para convertir el resto a numérico
signature_matrix <- signature_matrix[, -1]

# Convertir todos los valores a numéricos
signature_matrix <- apply(signature_matrix, 2, as.numeric)

# Restaurar los nombres de fila con los nombres de los genes
rownames(signature_matrix) <- gene_names

# Convertir a data frame
signature_matrix <- as.data.frame(signature_matrix)

# Mostrar la matriz final
signature_matrix

# Filter both matrices to keep only genes they have in common
common_genes <- intersect(rownames(tpm_matrix), rownames(signature_matrix))
tpm_matrix_filtered <- tpm_matrix[common_genes, ]
signature_matrix_filtered <- signature_matrix[common_genes, ]
# Convertir la firma filtrada a matriz
signature_matrix_filtered <- as.matrix(signature_matrix_filtered)
# 6. Perform deconvolution
# `granulator`'s `deconvolute` function requires the `sigMatrix` to be a named list.
deconvolution_results <- deconvolute(
  m = tpm_matrix_filtered,
  sigMatrix = list("ImmunoStates" = signature_matrix_filtered),
  methods = get_decon_methods(), # Runs all available methods
)

# 7. Access and analyze the results
# The output is a list containing estimated coefficients and proportions.
# The proportions are what you typically use for downstream analysis.
head(deconvolution_results$proportions)

# 8. Visualize the deconvolution results
# `granulator` has built-in plotting functions to visualize the output.
# We will use "nnls" as an example method and "S0" as the signature.
plot_proportions(
  deconvoluted = deconvolution_results,
  #method = "ols",
  method = "nnls", #Pick the desired method by uncommenting
  #method = "qprog",
  #method = "qprogwc", # <-
  #method = "dtangle", # <-
  #method = "rls",
  #method = "svr",
  signature = "ImmunoStates"
)

# You can also use `plot_deconvolute` to compare different methods.
plot_deconvolute(deconvoluted = deconvolution_results)

# 9. Save the results
# It's a good practice to save the proportions as a CSV file.
write.csv(deconvolution_results$proportions, "deconvolution_proportions.csv", row.names = FALSE)

#
# 7. Access and analyze the results
# Your code uses dtangle, which is a solid choice.
proportions_df <- deconvolution_results$proportions$dtangle_ImmunoStates

# Load the required packages
library(tidyr)
library(dplyr)
library(ggplot2)

# Pivot the data to a long format suitable for plotting and testing
proportions_long <- proportions_df %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "cell_type",
    values_to = "Proportion"
  )

# Add a new 'Group' column based on the SampleID prefix
proportions_long <- proportions_long %>%
  mutate(
    Group = case_when(
      startsWith(SampleID, "C") | startsWith(SampleID, "S") ~ "Healthy",
      startsWith(SampleID, "P") ~ "Cancer",
      TRUE ~ "Unknown"
    )
  )

proportions_long$Group <- factor(proportions_long$Group, levels = c("Healthy", "Cancer", "Unknown"))

# Corrected Statistical Test section:
# Check and remove groups with insufficient samples before the test
proportions_filtered <- proportions_long %>%
  group_by(cell_type, Group) %>%
  mutate(n = n()) %>%
  filter(n >= 3) %>%
  ungroup()

# Get unique cell types from the filtered data
cell_types <- unique(proportions_filtered$cell_type)
p_values <- list()

for (cell in cell_types) {
  subset_data <- proportions_filtered %>%
    filter(cell_type == cell)
  
  # Check again if there are exactly two groups for the test
  if (n_distinct(subset_data$Group) == 2) {
    test_result <- wilcox.test(Proportion ~ Group, data = subset_data)
    p_values[[cell]] <- test_result$p.value
  } else {
    p_values[[cell]] <- NA
  }
}

p_values_df <- data.frame(
  cell_type = names(p_values),
  p_value = unlist(p_values)
)

# 3. Adjust p-values for multiple comparisons
p_values_df <- p_values_df %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
  arrange(adjusted_p_value)

print("Statistical results:")
print(p_values_df)

# 4. Visualize the most significant differences with a box plot
top_cell_types <- p_values_df %>%
  arrange(adjusted_p_value) %>%
  top_n(5, -adjusted_p_value) %>%
  pull(cell_type)

proportions_to_plot_box <- proportions_long %>%
  filter(cell_type %in% top_cell_types)

# Verificar si hay tipos celulares significativos antes de graficar
if (length(top_cell_types) > 0) {
  proportions_to_plot_box <- proportions_long %>%
    filter(cell_type %in% top_cell_types)
  
  ggplot(proportions_to_plot_box, aes(x = Group, y = Proportion, fill = Group)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5) +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(
      title = "Significant Immune Cell Proportions: Healthy vs. Cancer",
      y = "Estimated Cell Proportion",
      x = "Group"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
} else {
  message("No significant cell types found for plotting.")
}

# 5. Visualize all proportions with a stacked bar plot
healthy_cancer_proportions_bar <- proportions_long %>%
  filter(Group %in% c("Healthy", "Cancer"))

ggplot(healthy_cancer_proportions_bar, aes(x = SampleID, y = Proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  labs(
    title = "Overall Immune Cell Proportions: Healthy vs. Cancer",
    y = "Estimated Cell Proportion",
    x = "Sample ID",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#
# TWO GROUPS BAR CHART
#

# Calculate the mean proportion for each cell type in each group
pooled_proportions <- proportions_long %>%
  group_by(Group, cell_type) %>%
  summarise(
    MeanProportion = mean(Proportion),
    .groups = 'drop'
  ) %>%
  filter(Group != "Unknown") # Exclude the 'Unknown' group if it exists

# Create the stacked bar plot with only two bars
ggplot(pooled_proportions, aes(x = Group, y = MeanProportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Pooled Immune Cell Proportions: Healthy vs. Cancer",
    y = "Mean Estimated Cell Proportion",
    x = "Group",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), # Center the group labels
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )



