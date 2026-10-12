# =============================================================================
# Plasmodium falciparum persister parasites — guided scRNA-seq QC and exploration
# PURPOSE: inspect 10 Cell Ranger matrices, plot at every stage, audit organellar
#          features, optionally remove VERIFIED mitochondrial genes, optionally
#          filter cells, then perform per-library PCA and unclustered UMAP.
# IMPORTANT: No biological replicates, treatment groups or persister states are
#            inferred from folder names. No clustering / differential expression.
# RUN: open this file in RStudio, press Source; inspect the Plots pane on right.
# =============================================================================

# 0. USER SETTINGS ------------------------------------------------------------
data_root <- path.expand('~/OneDrive/scRNA_data')
sample_ids <- c('6_1','6_2','7_1','7_2','8_1','8_2','9_1','9_2','10_1','10_2')
output_root <- file.path(data_root, 'Pf_persister_guided_results')

# Possible annotations only: these are NOT proof a gene is mitochondrial.
mt_candidate_patterns <- c('^MT-', 'mitochond', '^cox[123]$', '^cob$',
                           'cytochrome b', '^PF3D7_MIT', '^MAL_MITO')
# Set to a text file with ONE verified feature identifier per line after review.
# Can contain ORIGINAL HDF5 names or the names shown in feature_name_mapping.csv.
verified_mt_file <- NULL
# Optional sample-level QC thresholds CSV; no cell threshold by default.
# Expected columns: sample_id,min_features,max_features,min_umis,max_umis,
#                   max_percent_mt (blank values mean no cutoff).
filter_thresholds_file <- NULL
run_exploratory <- TRUE  # exploratory PCA/UMAP, without biological interpretation
show_plots_in_rstudio <- TRUE  # print ggplots to the RStudio Plots pane
pause_after_plot <- FALSE  # set TRUE for optional interactive inspection
seed <- 42L

# 1. PREREQUISITES AND OUTPUTS ------------------------------------------------
# WHY: fail clearly before processing; do not reinstall packages on each run.
# EXPECTED: output directory and a manifest, without touching the input files.
required <- c('Seurat', 'SeuratObject', 'hdf5r', 'Matrix', 'ggplot2')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Missing packages: ', paste(missing, collapse=', '),
                          '. Install once with install.packages(c(',
                          paste(sprintf('"%s"', missing), collapse=', '), '))')
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
if (!dir.exists(data_root)) stop('Input directory not found: ', data_root)
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(output_root)) stop('Unable to create output directory: ', output_root)
set.seed(seed)

write_csv <- function(x, f) utils::write.csv(x, f, row.names=FALSE, na='')
step <- function(n, why, expected) {
  message('\n', strrep('=', 66), '\nSTEP ', n, ' | WHY: ', why,
          '\nEXPECTED: ', expected, '\n', strrep('=', 66))
}
plot_and_save <- function(p, path_no_ext, w=10, h=5) {
  # print(), not just ggsave(), is essential when using Source in RStudio.
  if (show_plots_in_rstudio && interactive()) {
    print(p)
    if (pause_after_plot) readline('Inspect the RStudio Plots pane; press Enter to continue...')
  }
  ggplot2::ggsave(paste0(path_no_ext, '.png'), p, width=w, height=h, dpi=300, bg='white')
  ggplot2::ggsave(paste0(path_no_ext, '.pdf'), p, width=w, height=h, bg='white')
  message('Plot saved: ', path_no_ext, ' [.png, .pdf]')
}
logline <- function(id, stage, status, note='') data.frame(sample_id=id, stage=stage,
                      status=status, note=as.character(note), stringsAsFactors=FALSE)

manifest_file <- file.path(output_root, 'sample_manifest_TO_ANNOTATE.csv')
if (!file.exists(manifest_file)) {
  manifest <- data.frame(sample_id=sample_ids,
    input_file=file.path(data_root, sample_ids, 'filtered_feature_bc_matrix.h5'),
    condition=NA_character_, biological_replicate=NA_character_,
    technical_replicate=NA_character_, timepoint=NA_character_,
    sequencing_batch=NA_character_, notes=NA_character_)
  write_csv(manifest, manifest_file)
}
# Never overwrite a manifest the researcher has already annotated.

verified_mt <- character()
if (!is.null(verified_mt_file)) {
  if (!file.exists(verified_mt_file)) stop('Verified MT list missing: ', verified_mt_file)
  verified_mt <- trimws(readLines(verified_mt_file, warn=FALSE))
  verified_mt <- unique(verified_mt[nzchar(verified_mt) & !grepl('^#', verified_mt)])
  if (!length(verified_mt)) stop('Verified mitochondrial feature list is empty')
}

thresholds <- NULL
if (!is.null(filter_thresholds_file)) {
  if (!file.exists(filter_thresholds_file)) stop('QC thresholds file missing')
  thresholds <- read.csv(filter_thresholds_file, stringsAsFactors=FALSE)
  cols <- c('sample_id','min_features','max_features','min_umis','max_umis','max_percent_mt')
  if (!all(cols %in% names(thresholds))) stop('Threshold CSV requires: ', paste(cols,collapse=', '))
  if (anyDuplicated(thresholds$sample_id) || nrow(thresholds)!=length(sample_ids) ||
      !setequal(thresholds$sample_id, sample_ids)) stop('Exactly one threshold row per sample is required')
  for (nm in setdiff(cols, 'sample_id')) {
    original <- thresholds[[nm]]
    thresholds[[nm]] <- suppressWarnings(as.numeric(original))
    bad <- !is.na(original) & nzchar(trimws(as.character(original))) & is.na(thresholds[[nm]])
    if (any(bad)) stop('Non-numeric threshold detected in: ', nm)
    if (any(!is.na(thresholds[[nm]]) & !is.finite(thresholds[[nm]]))) stop('Non-finite threshold: ', nm)
  }
  if (any(thresholds$min_features > thresholds$max_features, na.rm=TRUE) ||
      any(thresholds$min_umis > thresholds$max_umis, na.rm=TRUE)) stop('Min threshold exceeds max')
  if (any(!is.na(thresholds$max_percent_mt)) && !length(verified_mt))
    stop('Mitochondrial percentage filtering requires a verified mitochondrial feature list')
}

# 2. IMPORT, VALIDATE AND PLOT RAW COUNTS -------------------------------------
# WHY: earlier scripts failed with an empty MT-match table even though HDF5
#      imports worked. Keep raw barcode data and never impose min.features=200.
# EXPECTED: 10 raw-barcode CSVs, feature audits, QC scatterplots and diagnostics.
step(2, 'Validate each Cell Ranger gene-expression matrix',
     'Raw per-barcode metrics, feature names, and nonempty QC plots')
all_raw <- list(); summary_rows <- list(); status_rows <- list()

for (id in sample_ids) {
  message('\n--- SAMPLE ', id, ' ---')
  dest <- file.path(output_root, id)
  dir.create(dest, showWarnings=FALSE, recursive=TRUE)
  status <- list()
  tryCatch({
    file <- file.path(data_root, id, 'filtered_feature_bc_matrix.h5')
    if (!file.exists(file)) stop('Missing HDF5 input: ', file)
    m <- Seurat::Read10X_h5(file, use.names=TRUE, unique.features=TRUE)
    if (is.list(m)) {
      if (!'Gene Expression' %in% names(m)) stop('No Gene Expression modality; found: ', paste(names(m),collapse=', '))
      m <- m[['Gene Expression']]
    }
    if (!nrow(m) || !ncol(m)) stop('Empty gene-expression matrix')
    if (anyDuplicated(rownames(m)) || anyDuplicated(colnames(m))) stop('Duplicate features or cell barcodes')
    vals <- if (inherits(m,'sparseMatrix')) m@x else as.vector(m)
    if (anyNA(vals) || any(vals<0) || any(abs(vals-round(vals))>1e-8)) stop('Counts are missing, negative, or non-integers')
    counts <- as.numeric(Matrix::colSums(m))
    genes <- as.numeric(Matrix::colSums(m>0))
    raw <- data.frame(sample_id=id, barcode=colnames(m), nCount_raw=counts, nFeature_raw=genes)
    write_csv(raw, file.path(dest, 'raw_barcode_metrics.csv'))
    all_raw[[id]] <- raw
    status[[length(status)+1L]] <- logline(id,'01_import','OK',paste(nrow(m),'features;',ncol(m),'barcodes'))

    p <- ggplot2::ggplot(raw, ggplot2::aes(nCount_raw,nFeature_raw)) +
      ggplot2::geom_point(alpha=.24,size=.55) +
      ggplot2::scale_x_continuous(trans='log1p') + ggplot2::scale_y_continuous(trans='log1p') +
      ggplot2::labs(title=paste('Sample',id,'— raw cell complexity'), x='UMIs / barcode (log1p)',
                    y='Detected genes / barcode (log1p)') + ggplot2::theme_bw()
    plot_and_save(p,file.path(dest,'01_raw_UMI_vs_genes'),7,5)

    # 3. MITOCHONDRIAL FEATURE NAME INSPECTION --------------------------------
    # WHY: ^MT- matched zero of the 5,525 features in earlier diagnostics.
    # EXPECTED: full names, candidate annotations, and raw-to-Seurat mapping.
    step(paste('3, sample',id), 'Audit parasite mitochondrial gene annotations',
         'Feature-name mapping and candidate list; NO automatic deletion')
    raw_features <- rownames(m)
    hits <- vapply(raw_features, function(nm) any(vapply(mt_candidate_patterns,
             function(pattern) grepl(pattern,nm,ignore.case=TRUE,perl=TRUE),logical(1))),logical(1))
    write_csv(data.frame(raw_feature=raw_features, mitochondrial_candidate=hits),
              file.path(dest,'all_feature_names_audit.csv'))
    # Explicit zero-row frame avoids the earlier `arguments imply differing number of rows: 1, 0` error.
    write_csv(data.frame(raw_feature=raw_features[hits]),file.path(dest,'mitochondrial_CANDIDATES_review.csv'))
    message('Candidate mitochondrial features: ',sum(hits),' (review before removal)')

    # No automatic 200-gene cutoff: sample 8_2 may contain low-RNA cells.
    obj <- Seurat::CreateSeuratObject(m, project=paste0('Pf_',id),min.cells=0,min.features=0)
    obj$sample_id <- id
    mapping <- data.frame(raw_feature=raw_features,seurat_feature=rownames(obj))
    if (nrow(mapping)!=length(raw_features)) stop('Seurat feature mapping length mismatch')
    write_csv(mapping,file.path(dest,'feature_name_mapping.csv'))

    # 4. VERIFIED MITOCHONDRIAL REMOVAL (OPTIONAL) -----------------------------
    # WHY: organellar features can bias QC, but pattern hits do not establish identity.
    # EXPECTED: original percent.mt stored before removal; none removed by default.
    step(paste('4, sample',id),'Optionally exclude verified mitochondrial features',
         'Retained QC object; percent.mt only when a verified list is supplied')
    mt_seurat <- character()
    if (length(verified_mt)) {
      selected <- mapping$seurat_feature[mapping$raw_feature %in% verified_mt |
                                         mapping$seurat_feature %in% verified_mt]
      mt_seurat <- unique(selected)
      if (!length(mt_seurat)) stop('No verified mitochondrial IDs found in sample ',id)
      write_csv(data.frame(feature=mt_seurat),file.path(dest,'verified_mt_removed.csv'))
      obj$percent.mt <- Seurat::PercentageFeatureSet(obj,features=mt_seurat)
    }
    before <- obj[[]]
    write_csv(data.frame(barcode=rownames(before),before,check.names=FALSE),
              file.path(dest,'QC_metadata_before_MT.csv'))
    if (length(mt_seurat)) {
      keep <- setdiff(rownames(obj),mt_seurat)
      if (!length(keep)) stop('No genes left after mitochondrial exclusion')
      obj <- subset(obj,features=keep)
      cm <- SeuratObject::LayerData(obj,assay='RNA',layer='counts')
      obj$nCount_RNA <- as.numeric(Matrix::colSums(cm))
      obj$nFeature_RNA <- as.numeric(Matrix::colSums(cm>0))
      stopifnot(!any(mt_seurat %in% rownames(obj)))
      status[[length(status)+1L]] <- logline(id,'02_mitochondrial','REMOVED',length(mt_seurat))
    } else {
      status[[length(status)+1L]] <- logline(id,'02_mitochondrial','NOT_REMOVED','No verified features supplied')
    }
    after <- obj[[]]
    write_csv(data.frame(barcode=rownames(after),after,check.names=FALSE),
              file.path(dest,'QC_metadata_after_MT.csv'))
    saveRDS(obj,file.path(dest,'seurat_before_cell_filter.rds'))

    # PRE/POST QC plots use explicit ggplot metadata rather than relying on
    # potentially stale Seurat plotting layers.
    before$barcode <- rownames(before); after$barcode <- rownames(after)
    before$phase <- 'Before MT exclusion'; after$phase <- 'After MT step'
    combined <- rbind(before[,c('barcode','nFeature_RNA','nCount_RNA','phase')],
                      after[,c('barcode','nFeature_RNA','nCount_RNA','phase')])
    p <- ggplot2::ggplot(combined,ggplot2::aes(phase,nFeature_RNA)) +
         ggplot2::geom_violin(fill='grey85',scale='width') +
         ggplot2::geom_boxplot(width=.16,outlier.shape=NA) +
         ggplot2::theme_bw() + ggplot2::labs(title=paste(id,'— features before/after MT step'),
                                            x=NULL,y='Detected genes per barcode')
    plot_and_save(p,file.path(dest,'02_before_after_MT_features'),8,5)
    if ('percent.mt' %in% names(before)) {
      p <- ggplot2::ggplot(before,ggplot2::aes(x=percent.mt)) +
        ggplot2::geom_histogram(bins=35,fill='grey60',color='white') +
        ggplot2::theme_bw() + ggplot2::labs(title=paste(id,'— original mitochondrial counts (%)'),
                                     x='Percent of raw counts from verified MT genes',y='Barcodes')
      plot_and_save(p,file.path(dest,'03_original_percent_mt'),8,5)
    }

    # 5. OPTIONAL CELL QC FILTERING -------------------------------------------
    # WHY: avoid automatic loss of low-RNA persister candidates (notably 8_2).
    # EXPECTED: if no threshold CSV, retain all barcodes and label as unfiltered.
    step(paste('5, sample',id),'Apply only researcher-approved cell thresholds',
         'Cell retention metrics, with filtering explicitly marked if applied')
    cells_before <- ncol(obj)
    if (!is.null(thresholds)) {
      lim <- thresholds[match(id,thresholds$sample_id),]
      keep <- rep(TRUE,ncol(obj))
      if (!is.na(lim$min_features)) keep <- keep & obj$nFeature_RNA >= lim$min_features
      if (!is.na(lim$max_features)) keep <- keep & obj$nFeature_RNA <= lim$max_features
      if (!is.na(lim$min_umis)) keep <- keep & obj$nCount_RNA >= lim$min_umis
      if (!is.na(lim$max_umis)) keep <- keep & obj$nCount_RNA <= lim$max_umis
      if (!is.na(lim$max_percent_mt)) keep <- keep & obj$percent.mt <= lim$max_percent_mt
      keep[is.na(keep)] <- FALSE
      if (!any(keep)) stop('No cells pass the thresholds for ',id)
      obj <- subset(obj,cells=colnames(obj)[keep])
      status[[length(status)+1L]] <- logline(id,'03_cell_filter','FILTERED',paste(ncol(obj),'/',cells_before,'retained'))
    } else status[[length(status)+1L]] <- logline(id,'03_cell_filter','NOT_APPLIED','Exploratory QC only')
    saveRDS(obj,file.path(dest,'seurat_QC_ready.rds'))

    # 6. PER-SAMPLE NORMALIZATION, PCA, UNCLUSTERED UMAP ----------------------
    # WHY: exploratory within-library structure, NOT persister classification.
    # EXPECTED: PCA and UMAP on RStudio Plots pane; per-library RDS.
    step(paste('6, sample',id),'Explore transcriptomic heterogeneity within this library',
         'Normalized object, PCA, UMAP (if enough cells/features)')
    if (run_exploratory) {
      exploratory <- tryCatch({
        if (ncol(obj)<30L || nrow(obj)<50L) stop('Insufficient cells/features')
        obj <- Seurat::NormalizeData(obj,verbose=FALSE)
        obj <- Seurat::FindVariableFeatures(obj,nfeatures=min(2000L,nrow(obj)),verbose=FALSE)
        vars <- Seurat::VariableFeatures(obj)
        if (length(vars)<4L) stop('Insufficient variable genes')
        obj <- Seurat::ScaleData(obj,features=vars,verbose=FALSE)
        pcs <- min(20L,ncol(obj)-1L,length(vars)-1L)
        if (pcs<3L) stop('Insufficient PCA dimensions')
        obj <- Seurat::RunPCA(obj,npcs=pcs,features=vars,verbose=FALSE)
        p <- Seurat::DimPlot(obj,reduction='pca') + ggplot2::ggtitle(paste(id,'— exploratory PCA'))
        plot_and_save(p,file.path(dest,'04_PCA'),7,5)
        obj <- Seurat::RunUMAP(obj,dims=seq_len(min(10L,pcs)),seed.use=seed,verbose=FALSE)
        p <- Seurat::DimPlot(obj,reduction='umap') +
             ggplot2::ggtitle(paste(id,'— exploratory unclustered UMAP'))
        plot_and_save(p,file.path(dest,'05_UMAP'),7,5)
        saveRDS(obj,file.path(dest,'seurat_exploratory.rds'))
        logline(id,'04_exploration','OK','PCA and UMAP saved')
      }, error=function(e) logline(id,'04_exploration','ERROR',conditionMessage(e)))
      status[[length(status)+1L]] <- exploratory
    }
    summary_rows[[id]] <- data.frame(sample_id=id, original_features=nrow(m),
       original_barcodes=ncol(m), median_raw_genes=median(genes), median_raw_umis=median(counts),
       pct_barcodes_ge200_features=100*mean(genes>=200), mt_candidates=sum(hits),
       mt_verified_removed=length(mt_seurat), cells_after_optional_filter=ncol(obj),
       median_after_genes=median(obj$nFeature_RNA), median_after_umis=median(obj$nCount_RNA))
    status[[length(status)+1L]] <- logline(id,'pipeline','OK','Main QC stages completed')
  }, error=function(e) {
    message('Sample ',id,' failed: ',conditionMessage(e))
    status[[length(status)+1L]] <<- logline(id,'pipeline','ERROR',conditionMessage(e))
  })
  status_rows[[id]] <- do.call(rbind,status)
  invisible(gc(verbose=FALSE))
}

# 7. COMPARE ALL TEN SAMPLES --------------------------------------------------
# WHY: identify sample-level QC differences before choosing thresholds.
# EXPECTED: visible violin/scatter plots plus saved figures and CSV tables.
step(7,'Compare raw QC across libraries','Figures in RStudio Plots pane; combined sample summary')
status_table <- do.call(rbind,status_rows)
write_csv(status_table,file.path(output_root,'pipeline_status.csv'))
if (length(summary_rows)) {
  diag <- do.call(rbind,summary_rows); rownames(diag) <- NULL
  write_csv(diag,file.path(output_root,'all_samples_diagnostics.csv'))
  print(diag,row.names=FALSE)
}
if (length(all_raw)) {
  d <- do.call(rbind,all_raw); rownames(d) <- NULL
  d$sample_id <- factor(d$sample_id,levels=sample_ids)
  write_csv(d,file.path(output_root,'all_raw_barcode_metrics.csv'))
  p <- ggplot2::ggplot(d,ggplot2::aes(sample_id,nFeature_raw)) +
    ggplot2::geom_violin(fill='grey85',scale='width') +
    ggplot2::geom_boxplot(width=.12,outlier.shape=NA) + ggplot2::theme_bw() +
    ggplot2::labs(title='P. falciparum — raw genes per barcode',x='Library',y='Detected genes')
  plot_and_save(p,file.path(output_root,'all_samples_genes_violin'),11,5)
  p <- ggplot2::ggplot(d,ggplot2::aes(sample_id,nCount_raw)) +
    ggplot2::geom_violin(fill='grey85',scale='width') +
    ggplot2::geom_boxplot(width=.12,outlier.shape=NA) +
    ggplot2::scale_y_continuous(trans='log1p') + ggplot2::theme_bw() +
    ggplot2::labs(title='P. falciparum — raw UMIs per barcode',x='Library',y='UMIs (log1p)')
  plot_and_save(p,file.path(output_root,'all_samples_UMIs_violin'),11,5)
  p <- ggplot2::ggplot(d,ggplot2::aes(nCount_raw,nFeature_raw)) +
    ggplot2::geom_point(alpha=.18,size=.4) + ggplot2::facet_wrap(~sample_id,scales='free') +
    ggplot2::scale_x_continuous(trans='log1p') +
    ggplot2::scale_y_continuous(trans='log1p') + ggplot2::theme_bw() +
    ggplot2::labs(title='P. falciparum — transcript complexity',x='UMIs (log1p)',y='Genes (log1p)')
  plot_and_save(p,file.path(output_root,'all_samples_complexity'),12,8)
}
writeLines(capture.output(sessionInfo()),file.path(output_root,'session_info.txt'))
message('\nCOMPLETE: inspect ',file.path(output_root,'pipeline_status.csv'))
message('Plots appear in RStudio Plots pane when sourced interactively; use its arrows to review.')
message('No mitochondrial genes are removed without a verified list. No automatic QC cutoff.')
