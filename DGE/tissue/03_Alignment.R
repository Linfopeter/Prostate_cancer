# ==============================================================================
# Script Unificado: Alineamiento y Conteo de Datos Paired-End
# Adaptado para el análisis en el directorio 'tissue1'
# ==============================================================================

# Cargar la librería necesaria
library(Rsubread)

# Limpiar las variables del sistema
rm(list = ls())

# ===========================================================================
# 1. Configuración de Rutas y Creación de Directorios
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Definir las rutas base
dir_proyecto <- "."
dir_analisis <- file.path(dir_proyecto, "tissue1") # Directorio específico para este análisis
dir_refdata <- file.path(dir_proyecto, "RefData")
dir_refgen <- file.path(dir_proyecto, "RefGen48") # Versión del genoma v48

# Crear los directorios de resultados si no existen
dir_resultados <- file.path(dir_analisis, "Resultados")
dir_align_out <- file.path(dir_resultados, "align")
dir_counts_out <- file.path(dir_resultados, "counts", "align")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")

# Asegura que todos los directorios se creen correctamente
# Usar recursive = TRUE para crear directorios anidados.
dir.create(dir_resultados, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_align_out, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_counts_out, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_deseq_out, showWarnings = FALSE, recursive = TRUE)

# Directorio donde se encuentran tus datos raw
dir_data <- file.path(dir_analisis, "data")

# ===========================================================================
# 2. Alineamiento con Rsubread::align() para datos Paired-End
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Iniciando el proceso de alineamiento para 'tissue1'...")

# Listar y emparejar archivos fastq.gz
archivos_1 <- list.files(path = dir_data, pattern = "_1.fq.gz$")
archivos_2 <- list.files(path = dir_data, pattern = "_2.fq.gz$")

# Ordenar para asegurar que los pares coincidan
archivos_1 <- sort(archivos_1)
archivos_2 <- sort(archivos_2)

# Verificar que el número de archivos coincida
if (length(archivos_1) != length(archivos_2)) {
  stop("El número de archivos _1 y _2 no coincide. Revisa tus datos.")
}

message("Archivos _1 a procesar:")
print(archivos_1)
message("Archivos _2 a procesar:")
print(archivos_2)

# Bucle para alinear cada par de muestras
for (i in seq_along(archivos_1)) {
  nombre_muestra <- sub("_1.fq.gz$", "", archivos_1[i])
  file1 <- file.path(dir_data, archivos_1[i])
  file2 <- file.path(dir_data, archivos_2[i])
  bam_output_file <- file.path(dir_align_out, paste0(nombre_muestra, ".BAM"))
  
  message(paste("Alineando muestra:", nombre_muestra))
  
  align(
    index = file.path(dir_refgen, "Hg38v48"),
    readfile1 = file1,
    readfile2 = file2,
    type = "rna",
    input_format = "gzFASTQ",
    output_format = "BAM",
    output_file = bam_output_file,
    phredOffset = 33,
    nsubreads = 10,
    TH1 = 3,
    TH2 = 1,
    maxMismatches = 3,
    unique = FALSE,
    nBestLocations = 1,
    indels = 5,
    complexIndels = FALSE,
    nTrim5 = 10,
    nTrim3 = 0,
    minFragLength = 50,
    maxFragLength = 600,
    PE_orientation = "fr",
    nthreads = 8,
    keepReadOrder = FALSE,
    sortReadsByCoordinates = TRUE,
    color2base = FALSE,
    DP_GapOpenPenalty = -1,
    DP_GapExtPenalty = 0,
    DP_MismatchPenalty = 0,
    DP_MatchScore = 2,
    detectSV = FALSE,
    useAnnotation = TRUE,
    annot.ext = file.path(dir_refdata, "gencode.v48.primary_assembly.annotation.gtf.gz"),
    isGTF = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    chrAliases = NULL
  )
}

# ==============================================================================
# 3. Conteo de lecturas con Rsubread::featureCounts()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Iniciando el conteo de lecturas...")

# Listar todos los archivos BAM generados
archivos_bam <- list.files(path = dir_align_out, pattern = ".BAM$", full.names = TRUE)
message("Archivos BAM a contar:")
print(archivos_bam)

# Conteo de todas las muestras en una sola llamada
cuentas <- featureCounts(
  files = archivos_bam,
  annot.ext = file.path(dir_refdata, "gencode.v48.primary_assembly.annotation.gtf.gz"),
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  nthreads = 22,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id"
)

# ==============================================================================
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
file_out <- file.path(dir_counts_out, "Matriz_conteo_align.csv")
write.csv(tabla_conteos, file = file_out, quote = FALSE, row.names = FALSE)

message(paste("Proceso completado. La matriz de conteo se ha guardado en:", file_out))

# Se recomienda guardar el objeto de conteo de R para análisis posteriores
save(cuentas, file = file.path(dir_deseq_out, "cuentas_raw.RData"))
