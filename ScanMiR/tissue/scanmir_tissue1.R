## =======================
## scanMiR — Significance-only (sites & pairs) + FORCED PARALLEL (PSOCK) — PNG/PDF
## =======================
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

bio_pkgs <- c("scanMiR","scanMiRData","Biostrings","GenomicRanges",
              "data.table","BiocParallel","ggplot2","pheatmap","matrixStats")
for (p in bio_pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
if (!requireNamespace("metap", quietly=TRUE)) install.packages("metap")  # Fisher

suppressPackageStartupMessages({
  library(scanMiR); library(scanMiRData); library(Biostrings); library(GenomicRanges)
  library(data.table); library(BiocParallel); library(ggplot2); library(pheatmap)
  library(matrixStats); library(metap)
})

## =======================
## Parámetros
## =======================
# setwd("./tissue1")                     # <-- ajusta si procede
FASTA_IN   <- "lncrnas_tissue1.fasta"      # <-- tu FASTA
SPECIES    <- "hsa"
N_WORKERS  <- 22                         # <-- nº de núcleos
ALPHA      <- 0.05                       # umbral FDR para sitios y pares
CHUNK_SIZE <- 200                        # nº de miRNAs por chunk (ajústalo si quieres)
OUTDIR     <- "scanmir_out"; dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
vizdir     <- file.path(OUTDIR, "plots"); dir.create(vizdir, showWarnings=FALSE)

## Temp local (ayuda a PSOCK en Windows)
if (.Platform$OS.type == "windows") {
  if (!dir.exists("C:/tmp")) dir.create("C:/tmp", showWarnings=FALSE)
  Sys.setenv(TMPDIR="C:/tmp")
}

## =======================
## Paralelismo PSOCK explícito
## =======================
BPP <- SnowParam(workers = N_WORKERS, type = "SOCK",
                 progressbar = TRUE, RNGseed = 1, stop.on.error = FALSE,
                 timeout = 1e7)
register(BPP, default = TRUE)
bpstart(BPP)
message(sprintf("[scanMiR] PSOCK online con %d workers.", bpworkers(bpparam())))

## (opcional) sanity check
.try <- try({
  pids <- bplapply(1:max(4, min(8, N_WORKERS)), function(i) Sys.getpid())
  message(sprintf("[check] workers únicos detectados: %d", length(unique(unlist(pids)))))
}, silent=TRUE)

## =======================
## 1) Cargar FASTA (U->T) y asegurar nombres únicos
## =======================
tx <- readDNAStringSet(FASTA_IN)
tx_chr <- as.character(tx)
if (any(grepl("U", tx_chr, fixed=TRUE))) {
  tx_chr <- chartr("U","T", tx_chr)
  tx <- DNAStringSet(tx_chr)
}
names(tx) <- make.unique(names(tx))
stopifnot(length(tx) > 0)

## =======================
## 2) Modelos Kd
## =======================
kmods <- getKdModels(SPECIES)
stopifnot(length(kmods) > 0)
mi_names <- names(kmods)

## =======================
## 3) ESCANEO FORZADO EN PARALELO POR CHUNKS (sin variables globales)
## =======================
## Armar lista de chunks con los KdModels ya recortados
split_idx   <- split(seq_along(mi_names), ceiling(seq_along(mi_names) / CHUNK_SIZE))
seed_chunks <- lapply(split_idx, function(idxs) kmods[ mi_names[idxs] ])

message(sprintf("[scan] %d miRNAs en %d chunks de ~%d c/u",
                length(mi_names), length(seed_chunks), CHUNK_SIZE))

## Función que corre en cada worker: recibe todo como argumentos
chunk_fun <- function(seeds_sub, tx_obj) {
  suppressPackageStartupMessages(library(scanMiR))  # asegurar namespace en el worker
  findSeedMatches(seqs = tx_obj, seeds = seeds_sub, useTmpFiles = TRUE)
}

## Ejecuta chunks en paralelo
hits_list <- bplapply(seed_chunks, chunk_fun, tx_obj = tx)

## Combinar resultados (GRanges)
## --- Combinar resultados (robusto) ---
stopifnot(length(hits_list) > 0)

# Diagnóstico rápido
lens <- vapply(hits_list, function(x) if (is(x, "GRanges")) length(x) else 0L, integer(1))
message(sprintf("[combine] chunks con ≥1 hit: %d / %d; total hits=%d",
                sum(lens > 0), length(lens), sum(lens)))

# Forzar a GRanges (evita GRangesList / clases raras)
hits <- suppressWarnings(do.call(c, lapply(hits_list, function(x) {
  if (is(x, "GRanges")) x else as(x, "GRanges")
})))
if (!is(hits, "GRanges")) {
  hits <- unlist(GRangesList(hits_list), use.names = FALSE)
}

if (length(hits) == 0L) {
  warning("No se encontraron matches en ningún chunk. Revisa 'warnings()' y verifica FASTA/kmods.")
  print(warnings())
  quit(save="no", status=0)  # o salta a fin si prefieres
}


## Agregación
agg <- aggregateMatches(hits, keepSiteInfo = TRUE)

## Guardar crudos
hits_dt_raw <- as.data.table(as.data.frame(hits))
fwrite(hits_dt_raw, file.path(OUTDIR, sprintf("hits_%s.tsv.gz", SPECIES)), sep="\t")
fwrite(as.data.table(agg), file.path(OUTDIR, sprintf("agg_%s.tsv.gz", SPECIES)), sep="\t")

## =======================
## 4) Estadística de significancia (sitios y pares)
## =======================
hits_dt <- copy(hits_dt_raw)
if (!"type" %in% names(hits_dt)) hits_dt[, type := NA_character_]
if ("log_kd" %in% names(hits_dt)) setnames(hits_dt, "log_kd", "logKd")
stopifnot("logKd" %in% names(hits_dt))

hits_dt[, canonical := grepl("8mer|7mer", type, ignore.case=TRUE)]
hits_can <- hits_dt[canonical == TRUE & is.finite(logKd)]
if (nrow(hits_can) == 0L) stop("No hay sitios canónicos con logKd finito.")

hits_can[, T := -logKd]                      # mayor = mejor afinidad
Fhat <- ecdf(hits_can$T)
hits_can[, p_site := pmax(1 - Fhat(T) + 1e-12, 1e-12)]

has_qvalue <- requireNamespace("qvalue", quietly=TRUE)
if (has_qvalue) hits_can[, q_site := qvalue::qvalue(p_site)$qvalues] else hits_can[, q_site := p.adjust(p_site, method = "BH")]

hits_can[, transcript := as.character(seqnames)]
hits_sig <- hits_can[q_site < ALPHA]
fwrite(hits_sig, file.path(OUTDIR, sprintf("sites_SIGNIFICANT_FDR%.02f.tsv.gz", ALPHA)), sep="\t")

## A nivel par (Fisher) + FDR
pair_p <- hits_can[, {
  ps <- p_site
  X  <- -2 * sum(log(ps))
  k  <- length(ps)
  p_comb <- pchisq(X, df = 2*k, lower.tail = FALSE)
  .(p_pair = p_comb)
}, by = .(miRNA, transcript)]

pair_p[, q_pair := p.adjust(p_pair, method="BH")]
pairs_sig <- pair_p[q_pair < ALPHA]
fwrite(pairs_sig, file.path(OUTDIR, sprintf("pairs_SIGNIFICANT_FDR%.02f.tsv.gz", ALPHA)), sep="\t")

## =======================
## 5) Figuras sólo con significativos
## =======================
## 5.1 KdModels
sig_miRNAs <- unique(pairs_sig$miRNA)
if (length(sig_miRNAs) > 0) {
  if (length(sig_miRNAs) > 12) {
    freq <- pairs_sig[, .N, by=miRNA][order(-N)]
    sig_miRNAs <- freq$miRNA[1:12]
  }
  for (m in sig_miRNAs) {
    if (!is.null(kmods[[m]])) {
      p <- plotKdModel(kmods[[m]], what="both", n=10)
      if ("ggplot" %in% class(p)) p <- p + ggtitle(sprintf("%s (significant targets present; FDR<%.02f)", m, ALPHA))
      safe_m <- gsub("[^A-Za-z0-9_.-]", "_", m)
      ggsave(file.path(vizdir, sprintf("01_KdModel_%s.png", safe_m)), plot=p, width=8, height=5.5, dpi=300)
      ggsave(file.path(vizdir, sprintf("01_KdModel_%s.pdf", safe_m)), plot=p, width=8, height=5.5, device=cairo_pdf)
    }
  }
} else message("No hay miRNAs con pares significativos para KdModels.")

## 5.2 Alineamientos — sitios significativos (hasta 20)
if (nrow(hits_sig) > 0) {
  setorder(hits_sig, -T)
  hits_ix <- as.data.table(as.data.frame(hits))[
    , .(idx = .I, miRNA, seqnames = as.character(seqnames), start, end)
  ]
  hits_sig[, site_key := paste(miRNA, transcript, start, end, sep=":")]
  sig_sites <- hits_sig[!duplicated(site_key)][1:min(20, .N)]
  
  has_svglite <- requireNamespace("svglite", quietly = TRUE)
  for (i in seq_len(nrow(sig_sites))) {
    row <- sig_sites[i]
    idx <- hits_ix[
      miRNA == row$miRNA & seqnames == row$transcript &
        start == row$start & end == row$end, idx][1]
    if (length(idx) == 1 && !is.na(idx) && !is.null(kmods[[row$miRNA]])) {
      gr <- hits[idx]
      p  <- viewTargetAlignment(m = gr, miRNA = kmods[[row$miRNA]], seqs = tx,
                                min3pMatch = 3L, outputType = "ggplot")
      if ("ggplot" %in% class(p)) {
        tag <- sprintf("%s_on_%s_%s-%s_qsite=%.3g",
                       row$miRNA, row$transcript, row$start, row$end, row$q_site)
        tag <- gsub("[^A-Za-z0-9_.-]", "_", tag)
        ggsave(file.path(vizdir, sprintf("02_align_%02d_%s.png", i, tag)),
               plot=p, width=9, height=3.2, dpi=300)
        if (has_svglite) ggsave(file.path(vizdir, sprintf("02_align_%02d_%s.svg", i, tag)),
                                plot=p, width=9, height=3.2, device="svg")
        else ggsave(file.path(vizdir, sprintf("02_align_%02d_%s.pdf", i, tag)),
                    plot=p, width=9, height=3.2, device=cairo_pdf)
      }
    }
  }
} else message("No hay sitios significativos (q_site < ALPHA) para alinear.")

## 5.3 Heatmap — pares significativos
agg_dt <- as.data.table(agg)
score_col <- if ("repression" %in% names(agg_dt)) "repression" else "score"

if (nrow(pairs_sig) > 0) {
  setkey(pairs_sig, miRNA, transcript); setkey(agg_dt, miRNA, transcript)
  agg_sig <- agg_dt[pairs_sig, nomatch=0L]
  if (nrow(agg_sig) > 0) {
    mat <- dcast(agg_sig, miRNA ~ transcript, value.var = score_col, fill = 0)
    rn <- mat$miRNA; mat$miRNA <- NULL; m <- as.matrix(mat); rownames(m) <- rn
    
    nonzero_rows <- rowSums(m != 0) > 0; nonzero_cols <- colSums(m != 0) > 0
    m <- m[nonzero_rows, nonzero_cols, drop=FALSE]
    
    if (nrow(m) >= 2 && ncol(m) >= 2) {
      topR <- min(30L, nrow(m)); topC <- min(40L, ncol(m))
      vr <- matrixStats::rowVars(m); vc <- matrixStats::colVars(m)
      m  <- m[order(vr, decreasing=TRUE)[seq_len(topR)], order(vc, decreasing=TRUE)[seq_len(topC)], drop=FALSE]
      
      m_z <- t(scale(t(m))); m_z[is.na(m_z)] <- 0
      cap <- 2; m_z[m_z > cap] <- cap; m_z[m_z < -cap] <- -cap
      
      rownames(m_z) <- gsub("^hsa-", "", rownames(m_z))
      rownames(m_z) <- gsub("-(3p|5p)$", "", rownames(m_z))
      colnames(m_z) <- substr(colnames(m_z), 1, 24)
      
      pal <- colorRampPalette(c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
                                "#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026"))
      main_t <- sprintf("Significant miRNA interactions with predicted lncRNAs (row-wise z-scores; FDR<%.02f)", ALPHA)
      
      png(file.path(vizdir, "03_heatmap_repression_SIGNIFICANT.png"), width=2400, height=1600, res=300)
      pheatmap(m_z, color=pal(255), breaks=seq(-cap, cap, length.out=256),
               cluster_rows=TRUE, cluster_cols=TRUE,
               border_color="grey85", legend=TRUE,
               fontsize_row=8, fontsize_col=8, angle_col=45,
               main=main_t)
      dev.off()
      
      pdf(file.path(vizdir, "03_heatmap_repression_SIGNIFICANT.pdf"), width=12, height=8.2)
      pheatmap(m_z, color=pal(255), breaks=seq(-cap, cap, length.out=256),
               cluster_rows=TRUE, cluster_cols=TRUE,
               border_color="grey85", legend=TRUE,
               fontsize_row=8, fontsize_col=8, angle_col=45,
               main=main_t)
      dev.off()
    } else message("Tras FDR, la matriz no alcanza 2×2 para heatmap.")
  } else message("No hay valores agregados para los pares significativos.")
} else message("No hay pares significativos (q_pair < ALPHA).")

## =======================
## Cierre
## =======================
bpstop(BPP)
message(sprintf("✅ Listo. Sitios significativos: %d; Pares significativos: %d. Salida en: %s",
                nrow(hits_sig), nrow(pairs_sig), normalizePath(OUTDIR)))
