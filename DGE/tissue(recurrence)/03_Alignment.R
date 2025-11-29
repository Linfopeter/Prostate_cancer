# ==============================================================================
# Script Unificado: Alineamiento y Conteo de Datos Single-End
# Adaptado para el análisis en el directorio 'tissue2'
# ==============================================================================

# Cargar la librería necesaria
library(Rsubread)
rm(list = ls()) # Limpia el entorno

# 1. Configuración de Rutas y Creación de Directorios
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definir las rutas base
dir_proyecto <- "."
dir_analisis <- file.path(dir_proyecto, "tissue2")
dir_refdata <- file.path(dir_proyecto, "RefData")
dir_refgen <- file.path(dir_proyecto, "RefGen48")

# Crear los directorios de resultados si no existen
dir_resultados <- file.path(dir_analisis, "Resultados")
dir_align_out <- file.path(dir_resultados, "align")
dir_counts_out <- file.path(dir_resultados, "counts", "align")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")

# Asegura que todos los directorios se creen correctamente de forma recursiva 
dir.create(dir_align_out, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_counts_out, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_deseq_out, showWarnings = FALSE, recursive = TRUE)

# Directorio donde se encuentran tus datos raw
dir_data <- file.path(dir_analisis, "data")

# 2. Alineamiento con Rsubread::align() para datos Single-End
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Iniciando el proceso de alineamiento para 'tissue2'...")

# Listar archivos fastq.gz
archivos_fastq <- list.files(path = dir_data, pattern = ".fq.gz$", full.names = FALSE)

if (length(archivos_fastq) == 0) {
  stop("No se encontraron archivos .fq.gz en tissue2/data")
}

message("Archivos a procesar:")
print(archivos_fastq)

# Bucle para alinear cada muestra
for (i in seq_along(archivos_fastq)) {
  nombre_muestra <- sub(".fq.gz$", "", archivos_fastq[i])
  file1 <- file.path(dir_data, archivos_fastq[i])
  bam_output_file <- file.path(dir_align_out, paste0(nombre_muestra, ".BAM"))
  
  message(paste("Alineando muestra:", nombre_muestra))
  
  align(
    index = file.path(dir_refgen, "Hg38v48"),
    readfile1 = file1,
    type = "rna",
    output_format = "BAM",
    output_file = bam_output_file,
    phredOffset = 33,
    nthreads = 24,
    useAnnotation = TRUE,
    annot.ext = file.path(dir_refdata, "gencode.v48.primary_assembly.annotation.gtf.gz"),
    isGTF = TRUE
  )
}

# 3. Conteo de lecturas con Rsubread::featureCounts()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Iniciando el conteo de lecturas...")

# Listar todos los archivos BAM generados
archivos_bam <- list.files(path = dir_align_out, pattern = ".BAM$", full.names = TRUE)

if (length(archivos_bam) == 0) {
  stop("No se encontraron archivos .BAM en tissue2/Resultados/align")
}

message("Archivos BAM a contar:")
print(archivos_bam)

# Conteo de todas las muestras en una sola llamada (Correcto para análisis downstream)
cuentas <- featureCounts(
  files = archivos_bam,
  annot.ext = file.path(dir_refdata, "gencode.v48.primary_assembly.annotation.gtf.gz"),
  isGTFAnnotationFile = TRUE,
  nthreads = 24
)

# 4. Generación y Exportación de la Matriz de Conteo Final (.csv)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Generando la matriz de conteo final...")

# Combinar la anotación y los conteos en un único data frame
tabla_conteos <- data.frame(
  GeneID = cuentas$annotation$GeneID,
  Length = cuentas$annotation$Length,
  cuentas$counts,
  stringsAsFactors = FALSE
)

# Exportar la matriz de conteo a un archivo CSV
file_out <- file.path(dir_counts_out, "Matriz_conteo_align_single_end.csv")
write.csv(tabla_conteos, file = file_out, quote = FALSE, row.names = FALSE)

message(paste("Proceso completado. La matriz de conteo se ha guardado en:", file_out))

# Guardar el objeto de conteo de R para análisis posteriores
save(cuentas, file = file.path(dir_deseq_out, "cuentas_raw_single_end.RData"))