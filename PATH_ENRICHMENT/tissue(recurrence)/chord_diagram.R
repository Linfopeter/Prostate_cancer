# ==============================================================================
# SCRIPT INTEGRADO:
# 1. Análisis de Enriquecimiento de Vías con pathfindR
# 2. Visualización de Red Gen-Vía con ggraph
# ==============================================================================

# ------------------------------------------------------------------------------
# PASO 1: Cargar todas las librerías necesarias
# ------------------------------------------------------------------------------
# Se han combinado las librerías de ambos scripts. No hay conflictos de
# compatibilidad conocidos entre estos paquetes.
message("Cargando librerías...")
library(tidyverse)
library(pathfindR)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(ggrepel)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# PASO 2: Configuración de Rutas y Carga de Datos
# ------------------------------------------------------------------------------
getwd()
# Cargar los resultados del análisis de expresión diferencial
path_to_results <- file.path('Res_DS2_align_annot.csv') # Estandarizado el nombre del archivo
message("Cargando el archivo de resultados anotado: ", path_to_results)

if (!file.exists(path_to_results)) {
  stop("Error: El archivo 'ResDS2align_annot.csv' no se encontró en la ruta esperada.")
}
ds_volcano <- read.csv(path_to_results, header = TRUE, row.names = 1)

# ------------------------------------------------------------------------------
# PASO 3: Preparación de datos y Análisis de Enriquecimiento con pathfindR
# ------------------------------------------------------------------------------
message("Preparando datos para pathfindR...")
# Renombrar y seleccionar columnas (Genesymbol, logFC, FDRp) y eliminar NAs
datos_pathfindr <- ds_volcano %>%
  select(external_gene_name, log2FoldChange, padj) %>%
  na.omit() %>%
  rename(Genesymbol = external_gene_name, logFC = log2FoldChange, FDR_p = padj)

# Verificar el formato de entrada
input_testing(datos_pathfindr, p_val_threshold = 0.05)

message("\n--- Iniciando análisis de enriquecimiento de vías KEGG con pathfindR ---")
# Correr pathfindR para el conjunto de genes KEGG
# plot_enrichment_chart = FALSE porque generaremos un gráfico de red más avanzado
pathfindr_results <- run_pathfindR(
  datos_pathfindr,
  gene_sets = "KEGG",
  min_gset_size = 10,
  max_gset_size = 300,
  pin_name_path = "KEGG", # Usar la red de interacciones de KEGG
  p_val_threshold = 0.05,
  enrichment_threshold = 0.05,
  output_dir = getwd(),
  plot_enrichment_chart = FALSE
)

# PASO 4: Adaptar los Resultados de pathfindR para la Visualización (CORREGIDO)
message("Adaptando los resultados de pathfindR para el formato del gráfico de red...")

kegg_results_for_graphing <- pathfindr_results %>%
  # Une las columnas de genes Up y Down, manejando correctamente los casos vacíos
  tidyr::unite("geneID", c(Up_regulated, Down_regulated), sep = ",", na.rm = TRUE, remove = TRUE) %>%
  
  # Reemplaza las comas restantes por slashes ("/")
  mutate(geneID = str_replace_all(geneID, ",", "/")) %>%
  
  # Calcula la columna 'Count' (número de genes en la vía)
  mutate(Count = sapply(strsplit(geneID, "/"), length)) %>%
  
  # Renombrar columnas para que coincidan con el script de graficación
  # ❗ CORRECCIÓN: Se usa 'lowest_p' en lugar de 'adj_p_value'
  rename(
    Description = Term_Description,
    p.adjust = lowest_p 
  ) %>%
  
  # Seleccionar las columnas necesarias para mantener la compatibilidad
  select(ID, Description, p.adjust, Count, geneID)

# Imprime las primeras filas para verificar el resultado
print(head(kegg_results_for_graphing))

# ------------------------------------------------------------------------------
# PASO 5: Generar el Gráfico de Red con ggraph (usando los datos de pathfindR)
# ------------------------------------------------------------------------------
message("Construyendo el gráfico de red para las 5 vías más significativas...")

# 1) Seleccionar las 5 vías más significativas
top_kegg_df <- kegg_results_for_graphing %>%
  arrange(p.adjust) %>%
  head(6) %>%

  # --- AÑADE ESTA LÍNEA ---
  # Convierte la columna 'Description' en un factor para preservar el orden por significancia
  mutate(Description = factor(Description, levels = unique(Description)))

# 2) Crear los arcos (edges) entre vías y genes (VERSIÓN CON FILTRO ROBUSTO)
edge_data <- top_kegg_df %>%
  select(ID, geneID, Description) %>%
  separate_rows(geneID, sep = "/") %>%
  
  # Limpiamos la columna geneID directamente
  mutate(geneID = trimws(geneID)) %>% 
  
  # --- CORRECCIÓN FINAL ---
  # Filtro robusto que elimina tanto NA como cadenas de texto vacías
  filter(!is.na(geneID) & geneID != "") %>%
  
  # Renombramos las columnas finales correctamente
  rename(from = ID, to = geneID, category = Description)

# 3) Crear los nodos (nodes) para las vías
pathway_nodes <- top_kegg_df %>%
  select(ID, Description, Count) %>%
  rename(name = ID, label = Description, size = Count) %>%
  mutate(type = "pathway", category = label)

# 4) Crear los nodos para los genes (VERSIÓN FINAL Y DEFINITIVA)
# Se usa 'datos_pathfindr' que tiene la información original de logFC por gen
gene_nodes <- data.frame(name = unique(edge_data$to)) %>%
  mutate(type = "gene", size = 3) %>%
  left_join(
    datos_pathfindr %>% select(Genesymbol, logFC),
    by = c("name" = "Genesymbol")
  ) %>%
  # --- INICIO DE LA CORRECCIÓN ---
  # Renombramos solo logFC y creamos la columna 'label' a partir de 'name',
  # ¡sin eliminar la columna 'name' original!
  rename(foldChange = logFC) %>%
  mutate(label = name)
# --- FIN DE LA CORRECCIÓN ---

# 5) Unir todos los nodos y crear el objeto grafo
node_data <- bind_rows(pathway_nodes, gene_nodes)
graph_data <- tbl_graph(nodes = node_data, edges = edge_data, directed = FALSE)

# 6) Definir el layout circular y la paleta de colores
layout <- create_layout(graph_data, layout = "linear", circular = TRUE)
custom_colors <- c("#0066FF", "#FF9900", "#FF00FF", "#4f9f3b", "#e93e3f","#6A3D9A")
names(custom_colors) <- unique(pathway_nodes$label)

# 7) Generar el gráfico final
message("Renderizando el gráfico final...")
final_plot <- ggraph(layout) +
  geom_edge_arc(aes(color = category), strength = 0.5, alpha = 0.8, width = 0.8) +
  # --- AÑADE ESTA LÍNEA NUEVA ---
  # Esta capa dibuja el borde negro (ligeramente más grande, sin relleno)
  geom_node_point(aes(size = size, filter = type == "gene"), 
                  shape = 1, color = "black", stroke = 0.7) + # shape=1 es un círculo vacío con borde
  # --- FIN DE LA LÍNEA NUEVA ---
  geom_node_point(aes(size = size, color = foldChange, filter = type == "gene")) +
  geom_node_point(aes(size = size, fill = category, filter = type == "pathway"),
                  shape = 21, color = "black", stroke = 1) +
  geom_node_text(aes(label = label, filter = type == "gene"),
                 repel = TRUE, size = 4, max.overlaps = 100) +
  scale_color_gradient2(
    low = "#008000", mid = "white", high = "#FF0000", midpoint = 0,
    na.value = "grey70", name = "Log2 FC"
  ) +
  scale_fill_manual(values = custom_colors, name = "category") +
  scale_edge_color_manual(values = custom_colors, name = "category") +
  scale_size_continuous(range = c(2, 10), name = "size", breaks = c(14, 20, 26)) +
  theme_void() +
  theme(
    legend.position = "right",
    # --- AÑADE ESTA LÍNEA ---
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Margen de 1 cm en todos los lados
  )


# 8) Mostrar y guardar el gráfico
# Para mostrar en una ventana emergente (si usas RStudio, se muestra en el panel Plots)
x11()
print(final_plot)

# Guardar el gráfico en un archivo
output_filename <- file.path("KEGG_Network_from_pathfindR_pbmcs.png")
ggsave(
  filename = output_filename,
  plot = final_plot,
  width = 12, height = 10, units = "in", dpi = 300, bg = "white"
)

message(paste0("✅ ¡Proceso completado! El gráfico de red se ha guardado en: ", output_filename))

















# --- CÓDIGO DE DIAGNÓSTICO ---
# (Ejecutar solo si el error de 'tbl_graph' vuelve a aparecer)

message("--- Resultados del Diagnóstico ---")

message("\n1. Verificando valores inválidos en la columna 'to' de edge_data:")
print(paste("Número de NAs:", sum(is.na(edge_data$to))))
print(paste("Número de textos vacíos:", sum(edge_data$to == "")))

message("\n2. Verificando que todos los nodos del grafo existan en la tabla de nodos:")
# Esta línea debe devolver TRUE. Si devuelve FALSE, aquí está el problema.
nodos_en_edges <- unique(c(edge_data$from, edge_data$to))
nodos_en_nodes <- unique(node_data$name)
print(paste("Todos los nodos de 'edges' están en 'nodes':", all(nodos_en_edges %in% nodos_en_nodes)))

message("\n3. Mostrando algunas filas de los datos justo antes del error:")
message("--- head(edge_data) ---")
print(head(edge_data))
message("\n--- head(node_data) ---")
print(head(node_data))

message("--- Fin del Diagnóstico ---")











"#B2DF8A"
"#33A02C"
"#A6CEE3"
"#1F78B4"
"#E31A1C"
"#FB9A99"
"#FF7F00"
"#FDBF6F"
"#CAB2D6"
"#6A3D9A"
"#FFFF99"
"#B15928"
"#8DD3C7"
"#BEBADA"
"#FB8072"
"#80B1D3"
"#FDB462"
"#FFD700"
"#ADFF2F" 
"#6495ED"
"#FFA07A"
"#DDA0DD"