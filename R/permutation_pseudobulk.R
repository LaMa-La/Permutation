#' Generate Pseudobulk from seurat object
#'
#' @param seurat_analysis seurat object
#' @param gene_filter Filter genes by list
#' @param min.pt Min percent filter
#' @param perm List of permutation names
#' @param assay Name of the Assay in the seurat object
#' @param design_matrix Design matrix for DESeq2
#' @param save_DDS Logical. Save output as rds
#' @param outdir Output location for saved files
#' @param verbose Logical
#'
#' @return DESeq2 pseudobulk object
#'
#' @export generate_pseudobulk_perm
#'
#' @examples \dontrun{dds_perm <- generate_pseudobulk_perm(seurat_analysis = seurat_analysis)}
generate_pseudobulk_perm <- function(seurat_analysis,
                                     gene_filter = seurat_analysis@assays$RNA@var.features,
                                     min.pt = 0,
                                     perm = c("perm1","perm2","perm3", "perm4", "perm5"),
                                     assay = "RNA",
                                     design_matrix = ~ time + donor,
                                     save_DDS = FALSE,
                                     outdir = NULL,
                                     verbose = TRUE){

  pseudobulk_counts_perm <- list()
  sample_table_perm <- list()
  dds_perm <- list()
  norm_anno_perm <- list()
  for(i in perm){
    if(verbose){print(i)}
    pseudobulk_counts_perm[[i]] <- Seurat::AggregateExpression(seurat_analysis, assays = assay , slot = "counts" ,return.seurat = F, group.by = paste0("PID_", i))[[assay]]

    # Filter for NAs
    pseudobulk_counts_perm[[i]] <- pseudobulk_counts_perm[[i]][!is.na(rownames(pseudobulk_counts_perm[[i]])),]

    # Generate Sample table from column names
    tmp <- data.frame(row.names = seq(1, length(as.character(colnames(pseudobulk_counts_perm[[i]])))))
    # mandatory columns

    ## Ensure you generate a variable "ID" which matches the colnames in your pseudobulk counts matrix
    tmp$ID <- as.character(colnames(pseudobulk_counts_perm[[i]]))
    rownames(tmp) <- tmp$ID

    m <- stringr::str_split_fixed(as.character(colnames(pseudobulk_counts_perm[[i]])),pattern =  "_", n=2)
    # Ensure the presence of a column "condition" which will be used for the grouping of your samples
    tmp$condition <- factor(paste(m[,1], m[,2], sep="_"), levels = unique(paste(m[,1], m[,2], sep="_")))



    tmp$sample_name <- paste(m[,1], m[,2], sep="_")
    tmp$donor <- factor(m[,1], levels = sort(unique(m[,1])))


    # addtional columns
    tmp$time <- factor(m[,2], levels = sort(unique(m[,2])))


    # double-check your tmp includes all necessary information
    c("ID", "condition") %in% colnames(tmp)
    sample_table_perm[[i]] <- tmp

    dds_txi <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_counts_perm[[i]], colData = sample_table_perm[[i]] , design = as.formula(design_matrix))


    threshold <- min.pt * ncol(seurat_analysis)
    genes_to_keep <- rowSums(seurat_analysis@assays$RNA@counts > 0) >= threshold & rownames(seurat_analysis) %in% gene_filter
    dds <- dds_txi[genes_to_keep,]

    dds_perm[[i]] <- DESeq(dds)

  }
  if(save_DDS){saveRDS(dds_perm, file = paste0(outdir, "Pseudobulk_dds_perm.rds"))}
  return(dds_perm)
}


generate_norm_anno <- function(dds_object = dds){

  norm_anno <- as.data.frame(DESeq2::counts(dds_object, normalized=T))
  norm_anno$GENE <- row.names(norm_anno)

  return(norm_anno)
}


gernerate_norm_anno_perm <- function(dds_perm, verbose = T){
  for(i in names(dds_perm)){
    if(verbose){print(i)}
    #Does contain the normalized data
    norm_anno_perm[[i]] <- generate_norm_anno(dds_object = dds_perm[[i]])
  }
  return(norm_anno_perm)}

#' Permutation DE call pseudobulk
#'
#' @param dds_perm Output from generate_pseudobulk_perm
#' @param method "Wald" or "LRT"
#' @param comparison.group1 Define pair-wise comparisons
#' @param comparison.group2 Define pair-wise comparison control groups
#' @param time.var Name of the time variable
#' @param design_matrix Design matrix for the DE call
#' @param design_matrix_reduced reduced design matrix
#' @param donor.var Name of donor ID column in the metadata of the seurat object
#' @param save_DEG Logical. Save the output of the DE call as RDS
#' @param outdir Define location where the output files are saved
#' @param verbose Logical
#' @param ncores Numbers of cores used
#' @param alpha Significance threshold
#' @param multiple_testing DE call parameter
#' @param pAdjustMethod Pvalue adjustment method
#' @param independentFiltering DE call parameter
#' @param shrinkage DE call parameter
#' @param shrinkType DE call parameter
#'
#' @return DEG of Permuted Pseudobulk DE call
#'
#'
#' @export Perm_DEcall_pseudobulk
#'
#' @examples \dontrun{
#' DEG.perm <- Perm_pseudobulk_DEcall(dds_perm, method = c("pseudobulk_DESeq2_Wald","pseudobulk_DESeq2_LRT"),
#' comparison.group1 = c("W3", "W4", "W5", "W6","W7"),comparison.group2 = rep("W2", 5),save_DEG = T, outdir = outdir,
#' shrinkType = "apeglm")
#' }

Perm_DEcall_pseudobulk <- function(dds_perm, method = c("pseudobulk_DESeq2_Wald","pseudobulk_DESeq2_LRT"),
                                   comparison.group1 = c("V2", "V3", "V4"),
                                   comparison.group2 = rep("V1", 3),
                                   time.var = "time",
                                   design_matrix = ~ time + donor,
                                   design_matrix_reduced = ~ donor,
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

  DEG.perm <- list()

  for (m in method){
    if(verbose){print(m)}
    DEG.perm[[m]] <- data.frame()

    if(grepl("Wald", m, ignore.case = T)){

      # Generate Comparison table
      comparison_table <- data.frame(comparison = comparison.group1,
                                     control = comparison.group2)

      for (ds in names(dds_perm)){
        if(verbose){print(ds)}
        dds_dea <- dds_perm[[ds]]

        DESeq2::design(dds_dea) <- as.formula(design_matrix)

        norm_anno <- generate_norm_anno(dds_object = dds_dea)

        # DE calculation
        DEresults <- data.frame()

        for (i in unique(comparison_table$control)){
          dds_dea$time <- stats::relevel(dds_dea$time, i)
          dds_dea <- DESeq2::nbinomWaldTest(object = dds_dea)

        for(x in 1:nrow(comparison_table)){
          if(verbose){print(paste("Doing DE call for", comparison_table$comparison[x] , "vs", comparison_table$control[x]))}
          if (multiple_testing=="IHW") {
            tmp <- DESeq2::results(dds_dea,
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
            tmp <- DESeq2::results(dds_dea,
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

              tmp <- DESeq2::lfcShrink(dds_dea,
                                                 contrast = c(time.var,
                                                              paste(comparison_table$comparison[x]),
                                                              paste(comparison_table$control[x])),
                                                 res= tmp,
                                                 type = shrinkType)

            }else if(shrinkType == "apeglm"){

              tmp <- DESeq2::lfcShrink(dds_dea,
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
          tmp$permutation_round <- ds
          tmp$regulation <- ifelse(tmp$log2FoldChange > 0, "up", "down")
          tmp$significance <- ifelse(tmp$padj < 0.05, "sign", "ns")
          tmp$comparison_regulation <- paste(tmp$comparison, tmp$regulation, sep="_")
          tmp$gene <- rownames(tmp)

          DEresults <- rbind(DEresults, tmp)
        }}



        DEG.perm[[m]] <- rbind(DEG.perm[[m]],
                               DEresults)
      }} #close Wald loop

    if(grepl("LRT", m, ignore.case = T)){

      for (ds in names(dds_perm)){
        if(verbose){print(ds)}
        dds_dea <- dds_perm[[ds]]

        design(dds_dea) <- stats::as.formula(design_matrix)

        # DE calculation
        dds_lrt <- DESeq2::DESeq(dds_dea, test = "LRT", reduced = as.formula(design_matrix_reduced), parallel = parallel, BPPARAM = BPPARAM)

        tmp <- DESeq2::results(dds_lrt)
        tmp <- tmp[!is.na(tmp$pvalue),]
        tmp$permutation_round <- ds
        tmp$regulation <- ifelse(tmp$log2FoldChange > 0, "up", "down")
        tmp$significance <- ifelse(tmp$padj < 0.05, "sign", "ns")
        tmp$gene <- rownames(tmp)
        tmp$comparison <- "longitudinal"

        DEG.perm[[m]] <- rbind(DEG.perm[[m]], as.data.frame(tmp))

      }}#close LRT loop

    if(save_DEG){saveRDS(DEG.perm[[m]], file = paste0(outdir,"DEG_",m,"_perm.rds"))}
  }
  return(DEG.perm)
}
