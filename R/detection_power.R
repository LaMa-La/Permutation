#' Power DE call
#'
#' @param seurat_analysis Seurat object
#' @param method List of available options c("seurat_wilcox",
#' "linear_TDE",
#' "mixed_TDE",
#' "seurat_MAST",
#' "seurat_MAST_latent",
#' "MAST_longitudinal",
#' "MAST_longitudinal_RE",
#' "MAST_longitudinal_fixed")
#' @param donor.var Name of donor ID column in the metadata of the seurat object
#' @param time.var Name of time column in the metadata of the seurat object
#' @param comparison.group1 Define pair-wise comparisons
#' @param comparison.group2 Define pair-wise comparison control groups
#' @param seurat.test.use Passed to seurat FindMarkers test.use
#' @param min.pt Min percent filter
#' @param gene_filter Filter genes before model fitting
#' @param tde_method Can be "cell" or "pseudocell". See TDESeq for more information
#' @param cellWeights Cellweights from slingshot trajectory passed to TradeSeq FitGAM
#' @param nknots TradeSeq FitGAM nknots argument
#' @param save_DEG Logical. Save the output of the DE call as RDS
#' @param save_model Logical. Save the fitted models as RDS
#' @param outdir Define location where the output files are saved
#' @param verbose Logical
#' @param ncores Numbers of cores used
#'
#' @return List of DE output data frames
#' @export Power_DEcall
#'
#' @examples  \dontrun{
#' DEG.power <- Power_DEcall(seurat_analysis = seurat_analysis, method =
#' c("seurat_wilcox","linear_TDE", "mixed_TDE","seurat_MAST","seurat_MAST_latent",
#' "MAST_longitudinal","MAST_longitudinal_RE","MAST_longitudinal_fixed"),comparison.group1 = c("W3", "W4", "W5", "W6","W7"),
#' comparison.group2 = rep("W2", 5),save_DEG = T, ncores = 5, gene_filter = NULL, tde_method = "pseudocell", outdir = outdir)
#' }
Power_DEcall <- function(seurat_analysis,
                         method = c("seurat_wilcox","seurat_MAST", "seurat_MAST_latent"), #Does support Seurat, MAST, TDESeq and Tradeseq DE call
                         donor.var = "PID",
                         time.var = "time",
                         comparison.group1 = c("V2", "V3", "V4"),
                         comparison.group2 = rep("V1", 3),
                         seurat.test.use = list(),
                         min.pt = 0,
                         gene_filter = seurat_analysis@assays$RNA@var.features,
                         tde_method = c("linear"="cell", "mixed"="pseudocell"),
                         cellWeights = NULL,
                         nknots = 6,
                         save_DEG = TRUE,
                         save_model = FALSE,
                         load_model = FALSE,
                         model = list(),
                         outdir = NULL,
                         verbose = TRUE,
                         ncores = 1,
                         ebayes = T
){

  #Data prep
  ## For seurat
  if(any(grepl("seurat", method, ignore.case = T) | grepl("MAST", method, ignore.case = T)& !grepl("longitudinal", method, ignore.case = T))){
    #Define comparisons
    comparisons.df <- data.frame("group1" = comparison.group1,
                                 "group2" = comparison.group2)
    comparisons.df$comparison <- paste0(comparisons.df$group1, "_vs_", comparisons.df$group2)

    plan(multisession, workers = ncores)
  }

  ## For MAST
  if(any(!grepl("seurat", method, ignore.case = T) & grepl("MAST", method, ignore.case = T))){

    sca <- MAST::SceToSingleCellAssay(as.SingleCellExperiment(seurat_analysis))
    if(!is.null(gene_filter)){sca <- sca[gene_filter,]}
    sca <- MAST::filterLowExpressedGenes(sca, min.pt)
    options("mc.cores" = ncores)

  }

  ## For TDESeq
  if(any(grepl("TDE", method, ignore.case = T))){
    if(!is.null(gene_filter)){ tde_input <- suppressMessages(subset(seurat_analysis, features = gene_filter))}
    else {tde_input <- seurat_analysis}
    counts<- Seurat::GetAssayData(tde_input,'counts')
    norm.data<- Seurat::GetAssayData(tde_input,'data')
    meta.data<- tde_input@meta.data

    tde_power <- TDEseq::CreateTDEseqObject(counts = counts, data=norm.data, meta.data = meta.data)
    rm(tde_input)
  }

  ## TradeSeq
  if(any(grepl("tradeseq", method, ignore.case = T))){
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- ncores
  }

  DEG.power <- list()

  for (m in method){
    if(verbose){print(m)}
    DEG.power[[m]] <- data.frame()

    # Seurat based tests
    if(grepl("seurat", m, ignore.case = T)){
      #Determine the test to use based on the method string
      test <- ifelse(grepl("wilcox", m, ignore.case = T),  "wilcox",
                     ifelse(grepl("DESeq2", m, ignore.case = T),"DESeq2",
                            ifelse(grepl("MAST", m, ignore.case = T),"MAST", seurat.test.use[[m]])))

      # Include latent variable based on method string
      if (grepl("latent", m, ignore.case = TRUE)){
        latent.var <- donor.var} else {
          latent.var <- NULL }



      #Loop through comparisons
      for(c in comparisons.df$comparison){
        if(verbose){print(paste0("Performing DE call for ", c))}
        Idents(seurat_analysis) <- time.var
        tmp <- Seurat::FindMarkers(object = seurat_analysis,
                           ident.1 = comparisons.df[comparisons.df$comparison==c,]$group1,
                           ident.2 = comparisons.df[comparisons.df$comparison==c,]$group2,
                           only.pos = FALSE,
                           features = gene_filter,
                           verbose = verbose,
                           min.pct = min.pt,
                           logfc.threshold = 0,
                           test.use = test,
                           latent.vars =  latent.var)
        tmp <- tmp[!is.na(tmp$p_val),]
        if (nrow(tmp) > 1){
          tmp$comparison <- c
          tmp$gene <- row.names(tmp)
          tmp$regulation <- ifelse(tmp$avg_log2FC > 0, "up", "down")
          tmp$significance <- ifelse(tmp$p_val_adj < 0.05, "sign", "ns")
          tmp$comparison <- factor(tmp$comparison, levels=unique(tmp$comparison))
          tmp$comparison_regulation <- paste(tmp$comparison, tmp$regulation, sep="_")
          tmp$time <- str_split(tmp$comparison, "_", simplify = T)[,1]

          DEG.power[[m]] <- rbind(DEG.power[[m]], tmp)
        }}} #close seurat loop


    # longitudinal MAST based tests
    if(!grepl("seurat", m, ignore.case = T) & grepl("MAST", m, ignore.case = T)& grepl("longitudinal", m, ignore.case = T)){

      zlmTime_power <- list()
      if(load_model){zlmTime_power <- readRDS(model[[m]])}
      else{
      ifelse(!is.numeric(time.var), x <- paste0(time.var, "_num"), x <- time.var) # numerical input is required
      # Create formula based on method string
      formula <- as.formula(ifelse(grepl("fixed", m, ignore.case = T)|grepl("FE", m, ignore.case = T),  paste("~", x, "+", donor.var),
                                   ifelse(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T),paste("~", x, "+ (1|", donor.var, ")"),
                                          paste("~", x))))

      if(verbose){print(paste("Fitting model for", time.var, "with formula:", deparse(formula)))}

      # Fit the MAST model

      zlmTime_power[[time.var]] <- if(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T)){
        MAST::zlm(formula, sca, parallel = T, method = "glmer", ebayes = F)}
      else {
        MAST::zlm(formula, sca = sca, parallel = T, ebayes = ebayes)}
      if(save_model){saveRDS(zlmTime_power, file = paste0(outdir, "zlm_",m,"_power.rds"))}
      }




      if(verbose){print(paste("Doing DE call for", time.var))}
      ifelse(!is.numeric(time.var), x <- paste0(time.var, "_num"), x <- time.var)
      summaryCond_power <- MAST::summary(zlmTime_power[[time.var]], doLRT= x, parallel = T)


      summaryDt_power <- summaryCond_power$datatable
      tmp <- data.frame(merge(summaryDt_power[contrast==x & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                              summaryDt_power[contrast==x & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')) #logFC coefficients

      colnames(tmp) <- c("gene", "p_val", "avg_log2FC", "ci.hi", "ci.lo")
      tmp <- tmp[!is.na(tmp$p_val),]
      tmp$regulation <- ifelse(tmp$avg_log2FC > 0, "up", "down")
      rownames(tmp) <- tmp$gene
      tmp$p_val_adj <- p.adjust(tmp$p_val, 'fdr')
      tmp$comparison <- "longitudinal"


      DEG.power[[m]] <- rbind(DEG.power[[m]], tmp)
    }# closing MAST loop


    # MAST pairwise
    if(!grepl("seurat", m, ignore.case = T) & grepl("MAST", m, ignore.case = T)& !grepl("longitudinal", m, ignore.case = T)){

      zlmTime_power <- list()
      if(load_model){zlmTime_power <- readRDS(model[[m]])}
      else{
      for(control in unique(comparisons.df$group2)){
        ifelse(is.character(time.var), x <- time.var, x <-as.character(time.var))# categorical input is required
        colData(sca)[,x] <- factor(colData(sca)[,x])
        colData(sca)[,x] <- relevel(colData(sca)[,x],control)
        # Create formula based on method string
        formula <- as.formula(ifelse(grepl("fixed", m, ignore.case = T)|grepl("FE", m, ignore.case = T),  paste("~", x, "+", donor.var),
                                     ifelse(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T),paste("~", x, "+ (1|", donor.var, ")"),
                                            paste("~", x))))

        if(verbose){print(paste("Fitting model for reference level: ",control, "with formula:", deparse(formula)))}

        #Fit the MAST model
        zlmTime_power[[control]] <- if(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T)){
          MAST::zlm(formula, sca, parallel = T, method = "glmer", ebayes = F)}
        else{
          MAST::zlm(formula, sca = sca, parallel = T, ebayes = ebayes)}

      }
        if(save_model){saveRDS(zlmTime_power, file = paste0(outdir, "zlm_",m,"_power.rds"))}}

      for (c in unique(comparisons.df$comparison)){
        if(verbose){print(paste("Doing DE call for",c))}
        a <- paste0(time.var,comparisons.df[comparisons.df$comparison == c,]$group1)
        summaryCond_power <- MAST::summary(zlmTime_power[[comparisons.df[comparisons.df$comparison == c,]$group2]], doLRT= a , parallel = T)


        summaryDt_power <- summaryCond_power$datatable
        tmp <- data.frame(merge(summaryDt_power[contrast==a & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                summaryDt_power[contrast==a & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')) #logFC coefficients

        colnames(tmp) <- c("gene", "p_val", "avg_log2FC", "ci.hi", "ci.lo")
        tmp <- tmp[!is.na(tmp$p_val),]
        tmp$regulation <- ifelse(tmp$avg_log2FC > 0, "up", "down")
        rownames(tmp) <- tmp$gene
        tmp$p_val_adj <- p.adjust(tmp$p_val, 'fdr')
        tmp$comparison <- c


        DEG.power[[m]] <- rbind(DEG.power[[m]], tmp)
      }}# closing MAST pairwise loop


    # TDESeq tests
    if(grepl("TDE", m, ignore.case = T)){
      ifelse(grepl("linear", m, ignore.case = T), tde_method_use <- tde_method["linear"], tde_method_use <- tde_method["mixed"] )
      model <- ifelse(grepl("mixed", m, ignore.case = T), "lmm", "lm")

      tde_param <- list(sample.var = donor.var,
                        stage.var = time.var,
                        fit.model = model, # This defines the model to be linear or mixed
                        tde.thr = 0.05,
                        pct = min.pt,
                        mod = 'FastLMM')

      tde_power <- TDEseq::tdeseq(object = tde_power, tde.method=tde_method_use, tde.param=tde_param, num.core= ncores)
      tmp <- TDEseq::GetTDEseqAssayData(tde_power, slot='tde')
      tmp <- tmp[!is.na(tmp$pvalue),]
      tmp$regulation <- ifelse(tmp$patter %in% c("Recession","Trough"),"down",
                               ifelse(tmp$patter %in% c("None"), "None",
                                      "up"))
      tmp$comparison <- "longitudinal"
      DEG.power[[m]] <- rbind(DEG.power[[m]], tmp )
    }#closing TDESeq loop


    # Tradeseq
    if(grepl("tradeseq", m, ignore.case = T)){

      if (grepl("latent", m, ignore.case = TRUE)){
        PID <- seurat_analysis@meta.data[,donor.var]
        U <- model.matrix(~PID)} else { # (pseudo)time variable is added internally in the fitGAM function to the design matrix
          U <- NULL }

      threshold <- min.pt * ncol(seurat_analysis)
      genes_to_keep <- rowSums(seurat_analysis@assays$RNA@counts > 0) >= threshold
      seurat_use <- suppressMessages(subset(seurat_analysis, features = names(genes_to_keep[genes_to_keep])))

      if(!is.null(gene_filter)){seurat_use <- suppressMessages(subset(seurat_use, features = gene_filter))}

      set.seed(42)
      ifelse(!is.numeric(seurat_use@meta.data[,time.var]), x <- paste0(time.var, "_num"), x <- time.var) # numerical input is required
      ifelse(is.null(cellWeights),cellWeights <- data.frame("lineage" = rep(1, length(seurat_analysis@meta.data[,x])), row.names = names(seurat_analysis@meta.data[,x])), cellWeights <- cellWeights) # To define the lineages


      gam_power <- list()
      if(load_model){gam_power <- readRDS(model[[m]])}
      else{
      if(verbose){print(paste("Fitting model for", x))}
      gam_power[[time.var]] <- tradeSeq::fitGAM(counts = seurat_use@assays$RNA@counts,
                                                pseudotime = seurat_use[[x]],
                                                cellWeights = cellWeights,
                                                U = U,
                                                nknots = nknots,
                                                parallel = TRUE,
                                                BPPARAM = BPPARAM,
                                                verbose = verbose)

      if(save_model){saveRDS(gam_power, file = paste0(outdir,"GAM_", m, "_power.rds"))}
      }

      if(verbose){print(paste("Doing DE call for", time.var))}
      tmp <- tradeSeq::associationTest(gam_power[[time.var]])
      tmp <- tmp[!(is.na(tmp$pvalue)), ]
      tmp$padj <- p.adjust(tmp$pvalue, method = "fdr")
      tmp <- tmp[order(tmp$padj), ]
      tmp$gene <- rownames(tmp)
      tmp$regulation <- NA
      tmp$comparison <- "longitudinal"

      DEG.power[[m]] <- rbind(DEG.power[[m]], tmp)

    }# close tradeseq loop


    if(save_DEG){saveRDS(DEG.power[[m]], file = paste0(outdir, "DEG_",m,"_power.rds"))}
    invisible(gc())}#closing methods loop
  return(DEG.power)
}

#' Generate pseudobulk based on donor and time label
#'
#' @param seurat_analysis seurat object
#' @param gene_filter Filter genes by list
#' @param min.pt Min percent filter
#' @param time.var Name of the time variable column
#' @param assay Name of the Assay in the seurat object
#' @param design_matrix Design matrix for DESeq2
#' @param save_DDS Logical. Save output as rds
#' @param outdir Output location for saved files
#' @param verbose Logical
#'
#' @return DESeq2 opbject
#' @export generate_pseudobulk
#'
#' @examples \dontrun{dds <- generate_pseudobulk(seurat_analysis)}
generate_pseudobulk <- function(seurat_analysis,
                                gene_filter = seurat_analysis@assays$RNA@var.features,
                                min.pt = 0.1,
                                time.var = "time",
                                assay = "RNA",
                                design_matrix = ~ time + donor,
                                save_DDS = FALSE,
                                outdir = NULL,
                                verbose = TRUE){

  pseudobulk_counts_perm <- list()
  sample_table_perm <- list()
  dds_perm <- list()

  pseudobulk_counts_perm[[time.var]] <- AggregateExpression(seurat_analysis, assays = assay , slot = "counts" ,return.seurat = F, group.by = paste0("PID_", time.var))[[assay]]

  # Filter for NAs
  pseudobulk_counts_perm[[time.var]] <- pseudobulk_counts_perm[[time.var]][!is.na(rownames(pseudobulk_counts_perm[[time.var]])),]

  # Generate Sample table from column names
  tmp <- data.frame(row.names = seq(1, length(as.character(colnames(pseudobulk_counts_perm[[time.var]])))))
  # mandatory columns

  ## Ensure you generate a variable "ID" which matches the colnames in your pseudobulk counts matrix
  tmp$ID <- as.character(colnames(pseudobulk_counts_perm[[time.var]]))
  rownames(tmp) <- tmp$ID

  m <- stringr::str_split_fixed(as.character(colnames(pseudobulk_counts_perm[[time.var]])),pattern =  "_", n=2)
  # Ensure the presence of a column "condition" which will be used for the grouping of your samples
  tmp$condition <- factor(paste(m[,1], m[,2], sep="_"), levels = unique(paste(m[,1], m[,2], sep="_")))



  tmp$sample_name <- paste(m[,1], m[,2], sep="_")
  tmp$donor <- factor(m[,1], levels = sort(unique(m[,1])))


  # addtional columns
  tmp$time <- factor(m[,2], levels = sort(unique(m[,2])))


  # double-check your tmp includes all necessary information
  c("ID", "condition") %in% colnames(tmp)
  sample_table_perm[[time.var]] <- tmp

  dds_txi <- DESeqDataSetFromMatrix(countData = pseudobulk_counts_perm[[time.var]], colData = sample_table_perm[[time.var]] , design = as.formula(design_matrix))


  threshold <- min.pt * ncol(seurat_analysis)
  genes_to_keep <- rowSums(seurat_analysis@assays$RNA@counts > 0) >= threshold & rownames(seurat_analysis) %in% gene_filter
  dds <- dds_txi[genes_to_keep,]

  dds <- DESeq(dds)

  if(save_DDS){saveRDS(dds, file = paste0(outdir, "Pseudobulk_dds.rds"))}
  return(dds)
}

#' Power DE call pseudobulk
#' @param dds Output from generate_pseudobulk
#' @param method "Wald" or "LRT"
#' @param comparison.group1 Define pair-wise comparisons
#' @param comparison.group2 Define pair-wise comparison control groups
#' @param design_matrix Design matrix for the DE call
#' @param design_matrix_reduced reduced design matrix
#' @param donor.var Name of donor ID column in the metadata of the seurat object
#' @param time.var Name of time column in the metadata of the seurat object
#' @param save_DEG Logical. Save the output of the DE call as RDS
#' @param outdir Define location where the output files are saved
#' @param verbose Logical
#' @param ncores Numbers of cores used
#' @param alpha Significance Threshold
#' @param multiple_testing DE call parameter
#' @param pAdjustMethod Pvalue adjustment method
#' @param independentFiltering DE call parameter
#' @param shrinkage DE call parameter
#' @param shrinkType DE call parameter
#'
#'
#' @return Power DE call pseudobulk
#' @export Power_DEcall_pseudobulk
#'
#' @examples \dontrun{DEG.power <- Power_DEcall_pseudobulk(dds, method = c("pseudobulk_DESeq2_Wald","pseudobulk_DESeq2_LRT"),
#' comparison.group1 = c("W3", "W4", "W5", "W6","W7"),comparison.group2 = rep("W2", 5),save_DEG = T, outdir = outdir,
#' shrinkType = "apeglm")}

Power_DEcall_pseudobulk <- function(dds,
                                    method = c("pseudobulk_DESeq2_Wald","pseudobulk_DESeq2_LRT"),
                                    comparison.group1 = c("V2", "V3", "V4"),
                                    comparison.group2 = rep("V1", 3),
                                    design_matrix = ~ time + donor,
                                    design_matrix_reduced = ~ donor,
                                    time.var = "time",
                                    donor.var = "PID",
                                    save_DEG = TRUE,
                                    outdir = NULL,
                                    verbose = TRUE,
                                    ncores = 1,
                                    alpha = 0.05,
                                    pAdjustMethod = "BH",
                                    independentFiltering = TRUE,
                                    shrinkage = TRUE,
                                    shrinkType = "normal",
                                    multiple_testing="IHW"){

  if(ncores > 1){
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- ncores
    parallel <- TRUE}else{parallel <- FALSE}

  DEG.power <- list()

  for (m in method){
    if(verbose){print(m)}
    DEG.power[[m]] <- data.frame()

    if(grepl("Wald", m, ignore.case = T)){

      # Generate Comparison table
      comparison_table <- data.frame(comparison = comparison.group1,
                                     control = comparison.group2)



      design(dds) <- as.formula(design_matrix)

      norm_anno <- generate_norm_anno(dds_object = dds)

      # DE calculation
      DEresults <- data.frame()

      for (i in unique(comparison_table$control)){
        dds$time <- stats::relevel(dds$time, i)
        dds <- DESeq2::nbinomWaldTest(object = dds)

        for(x in 1:nrow(comparison_table)){
          if(verbose){print(paste("Doing DE call for", comparison_table$comparison[x] , "vs", comparison_table$control[x]))}
          if (multiple_testing=="IHW") {
            tmp <- DESeq2::results(dds,
                                   contrast = c(time.var,
                                                paste(comparison_table$comparison[x]),
                                                paste(comparison_table$control[x])),
                                   lfcThreshold = 0,
                                   alpha = alpha,
                                   filterFun = ihw,
                                   pAdjustMethod = pAdjustMethod,
                                   altHypothesis = "greaterAbs")
            # Independent Filtering
          }else {
            tmp <- DESeq2::results(dds,
                                   contrast = c(time.var,
                                                paste(comparison_table$comparison[x]),
                                                paste(comparison_table$control[x])),
                                   lfcThreshold = 0,
                                   alpha = alpha,
                                   independentFiltering = independentFiltering,
                                   altHypothesis = "greaterAbs",
                                   pAdjustMethod= pAdjustMethod)
          }
          if(shrinkage == TRUE){
            if(shrinkType %in% c("normal", "ashr")){

              tmp <- DESeq2::lfcShrink(dds,
                                       contrast = c(time.var,
                                                    paste(comparison_table$comparison[x]),
                                                    paste(comparison_table$control[x])),
                                       res= tmp,
                                       type = shrinkType)

            }else if(shrinkType == "apeglm"){

              tmp <- DESeq2::lfcShrink(dds,
                                       coef = paste0(time.var, "_",
                                                     comparison_table$comparison[x], "_vs_",
                                                     comparison_table$control[x]),
                                       res= tmp,
                                       type = shrinkType,
                                       returnList = F)
            }
          }
          tmp <- as.data.frame(tmp)

          tmp <- tmp[!is.na(tmp$pvalue),]
          tmp$comparison <- paste0(comparison_table$comparison[x],"_vs_",
                                   comparison_table$control[x])
          tmp$comparison <- factor(tmp$comparison, levels=unique(tmp$comparison))
          tmp$time <- stringr::str_split(tmp$comparison, "_", simplify = T)[,1]
          tmp$regulation <- ifelse(tmp$log2FoldChange > 0, "up", "down")
          tmp$significance <- ifelse(tmp$padj < 0.05, "sign", "ns")
          tmp$comparison_regulation <- paste(tmp$comparison, tmp$regulation, sep="_")
          tmp$gene <- rownames(tmp)

          DEresults <- rbind(DEresults, tmp)
        }}

      DEG.power[[m]] <- rbind(DEG.power[[m]], DEresults)
    } #close Wald loop

    if(grepl("LRT", m, ignore.case = T)){

      design(dds) <- as.formula(design_matrix)

      # DE calculation
      dds_lrt <- DESeq(dds, test = "LRT", reduced = as.formula(design_matrix_reduced), parallel = parallel, BPPARAM = BPPARAM)

      tmp <- results(dds_lrt)
      tmp <- tmp[!is.na(tmp$pvalue),]
      tmp$regulation <- ifelse(tmp$log2FoldChange > 0, "up", "down")
      tmp$significance <- ifelse(tmp$padj < 0.05, "sign", "ns")
      tmp$gene <- rownames(tmp)
      tmp$comparison <- "longitudinal"

      DEG.power[[m]] <- rbind(DEG.power[[m]], as.data.frame(tmp))

    }#close LRT loop

    if(save_DEG){saveRDS(DEG.power[[m]], file = paste0(outdir,"DEG_",m,"_power.rds"))}
  }
  return(DEG.power)
}


#' Detection power calculated from the Power DE call
#'
#' @param DEG.power Output of the Power DE call
#' @param padj Name of the adjusted pvalue column. If NULL the adjusted pvalue column is generated based on regular expression.
#'
#' @return Detection power calculated from the Power DE call
#' @export detection_power_padj
#'
#' @examples \dontrun{detection_power_padj(DEG.power)}
detection_power_padj <- function(DEG.power, padj = NULL){

  fdr <- list()
  # Generate data for plot
  for (m in names(DEG.power)){
    DEG_df <- DEG.power[[m]]

    if(is.null(padj)){
      p <- grep(".adj", names(DEG_df), perl = T, value = T)}
    else{p<- padj[[m]]}

    fdr[[m]] <- DEG_df[,p]
  }
  avg_power <- list()
  for (n in names(fdr)){
    fdr_method <- fdr[[n]]
    for (i in c(0.01, 0.025, 0.05, 0.075, 0.1)){
      tmp <- length(which(fdr_method<i))/length(fdr_method)
      avg_power[[n]] <- c(avg_power[[n]], tmp)}
    names(avg_power[[n]]) <-  c(0.01 , 0.025, 0.05, 0.075, 0.1)}

  avg_power_plot <- reshape2::melt(avg_power)
  names(avg_power_plot) <- c("value", "method")
  avg_power_plot$alpha <- c(0.01, 0.025, 0.05, 0.075, 0.1)


  avg_power_plot$method <- factor(avg_power_plot$method, levels = names(DEG.power))
  return(avg_power_plot)

}

#' Plot detection power
#'
#' @param avg_power Output from detection_power_padj
#' @param color color vector
#'
#' @return Plot detection power
#' @export plot_detection_power
#'
#' @examples \dontrun{plot_detection_power(avg_power)}
plot_detection_power <- function(avg_power, color = color_method){
  ggplot(data= avg_power, aes(x = alpha, y = value, group.by = method, color = method))+
    geom_line()+
    theme_classic()+
    scale_color_manual(values = color)+
    xlab("Alpha")+ ylab("Average detection power")+
    theme(legend.title = element_blank())+ ggtitle("Average detection power")
}
