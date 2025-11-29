# ==============================================================================
#      PathfindR Enrichment Analysis
#      Bubble charts for every gene_set
# ==============================================================================

# Cargar las librerías necesarias
library(tidyverse)
library(pathfindR)
library(ggplot2)

## --- Configuración de Rutas y Carga de Datos ---
# Se asume que getwd() es el directorio "RNAseq"
dir_proyecto <- getwd()
dir_pbmcs <- file.path(dir_proyecto, "pbmcs")
dir_resultados <- file.path(dir_pbmcs, "Resultados")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")
dir_output_enrich <- file.path(dir_resultados, "Enrich")

# Directorio base para resultados de enriquecimiento
dir.create(dir_output_enrich, recursive = TRUE)

path_to_results <- file.path(dir_deseq_out, 'Res_DS2_align_annot.csv')

message("Cargando el archivo de resultados anotado...")
if (!file.exists(path_to_results)) {
  stop("Error: El archivo 'Res_DS2_align_annot.csv' no se encontró en la ruta esperada.")
}
ds_volcano <- read.csv(path_to_results, header = TRUE, row.names = 1)

## --- Preparación de datos para PathfindR ---
# Renombrar y seleccionar las columnas requeridas (Gene_symbol, logFC, FDR_p)
# y eliminar las filas con valores NA.
datos_pathfindr <- ds_volcano %>%
  select(external_gene_name, log2FoldChange, padj) %>%
  na.omit() %>%
  rename(Gene_symbol = external_gene_name, logFC = log2FoldChange, FDR_p = padj)

message("Verificando el formato de los datos de entrada...")
input_testing(datos_pathfindr, p_val_threshold = 0.05)
message("El formato de los datos de entrada es correcto.")

## --- Función para realizar el análisis de enriquecimiento y generar el gráfico de burbujas ---
perform_enrichment <- function(data, gene_set_name) {
  message(paste0("\n--- Iniciando análisis para vías ", gene_set_name, " ---"))
  
  # Definir el directorio específico para los resultados
  output_dir <- file.path(dir_output_enrich, gene_set_name)
  dir.create(output_dir, recursive = TRUE)
  
  # Correr pathfindR con el conjunto de genes especificado
  enrichment_results <- run_pathfindR(
    data,
    gene_sets = gene_set_name,
    min_gset_size = 10,
    max_gset_size = 300,
    pin_name_path = "KEGG",
    p_val_threshold = 0.05,
    enrichment_threshold = 0.05,
    output_dir = output_dir,
    plot_enrichment_chart = FALSE # Se generará el gráfico manualmente
  )
  
  if (!is.null(enrichment_results) && nrow(enrichment_results) > 0) {
    # Generar el gráfico de burbujas para las top 20 vías
    message(paste0("\nGenerando y mostrando el gráfico de burbujas de ", gene_set_name, "."))
    bubble_chart <- enrichment_chart(
      enrichment_results,
      top_terms = 20,
      plot_by_cluster = FALSE,
      even_breaks = TRUE
    )
    print(bubble_chart)
    
    # Guardar el gráfico de burbujas como un archivo .png en el directorio de salida
    chart_filename <- file.path(output_dir, paste0("Bubble_Chart_", gene_set_name, ".png"))
    ggsave(filename = chart_filename, plot = bubble_chart, width = 10, height = 8, dpi = 300)
    message(paste0("Gráfico de burbujas guardado en: ", chart_filename))
    
  } else {
    message(paste0("No se encontraron vías enriquecidas para ", gene_set_name, "."))
  }
}

## --- Ejecutar los análisis para KEGG, GO-All y Biogrid ---
perform_enrichment(datos_pathfindr, "KEGG")
perform_enrichment(datos_pathfindr, "GO-All")
perform_enrichment(datos_pathfindr, "Reactome")
perform_enrichment(datos_pathfindr, "BioCarta")
perform_enrichment(datos_pathfindr, "GO-CC")
perform_enrichment(datos_pathfindr, "GO-BP")
perform_enrichment(datos_pathfindr, "GO-MF")
perform_enrichment(datos_pathfindr, "cell_markers")

message("\nAnálisis de enriquecimiento finalizado. Revisa el panel de 'Plots' y las carpetas 'Enrich/KEGG', 'Enrich/GO-All' y 'Enrich/Biogrid' para los resultados.")





# X11() GRAPH PERSONALIZATION


perform_enrichment <- function(data, gene_set_name) {
  message(paste0("\n--- Iniciando análisis para vías ", gene_set_name, " ---"))
  
  output_dir <- file.path(dir_output_enrich, gene_set_name)
  dir.create(output_dir, recursive = TRUE)
  
  enrichment_results <- run_pathfindR(
    data,
    gene_sets = gene_set_name,
    min_gset_size = 10,
    max_gset_size = 300,
    pin_name_path = "KEGG",
    p_val_threshold = 0.05,
    enrichment_threshold = 0.05,
    output_dir = output_dir,
    plot_enrichment_chart = FALSE
  )
  
  if (!is.null(enrichment_results) && nrow(enrichment_results) > 0) {
    message(paste0("\nMostrando gráfico de burbujas para ", gene_set_name, "."))
    
    bubble_chart <- enrichment_chart(
      enrichment_results,
      top_terms = 20,
      plot_by_cluster = FALSE,
      even_breaks = TRUE
    )
    
    # Abrir ventana interactiva
    if (.Platform$OS.type == "windows") {
      x11()
    } else if (Sys.info()["sysname"] == "Darwin") {
      quartz()
    } else {
      x11()
    }
    
    print(bubble_chart)
    
    message("Personaliza y guarda el gráfico manualmente desde la ventana abierta.")
    message("Presiona [Enter] en la consola cuando termines para continuar con el siguiente análisis.")
    readline()
    
  } else {
    message(paste0("No se encontraron vías enriquecidas para ", gene_set_name, "."))
  }
}


## --- Ejecutar los análisis para KEGG, GO-All y Biogrid ---
perform_enrichment(datos_pathfindr, "KEGG")
perform_enrichment(datos_pathfindr, "GO-All")
perform_enrichment(datos_pathfindr, "Reactome")
perform_enrichment(datos_pathfindr, "BioCarta")
perform_enrichment(datos_pathfindr, "GO-CC")
perform_enrichment(datos_pathfindr, "GO-BP")
perform_enrichment(datos_pathfindr, "GO-MF")
perform_enrichment(datos_pathfindr, "cell_markers")



