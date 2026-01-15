#' DE anaylsis of permuted time labels
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
#' @param perm  Name of all permutation columns in the metadata of the seurat object
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
#' @export Perm_DEcall
#'
#' @examples
#' \dontrun{
#' DEG.perm <- Perm_DEcall(seurat_analysis,
#' method = c("seurat_wilcox","linear_TDE", "mixed_TDE","seurat_MAST","seurat_MAST_latent",
#' "MAST_longitudinal","MAST_longitudinal_RE","MAST_longitudinal_fixed"),
#' comparison.group1 = c("W3", "W4", "W5", "W6","W7"),comparison.group2 = rep("W2", 5),save_DEG = T, ncores = 5,
#' gene_filter = NULL, tde_method = "pseudocell", outdir = outdir)
#' }
#'
Perm_DEcall <- function(seurat_analysis,
                        method = c("seurat_wilcox","seurat_MAST", "seurat_MAST_latent"), #Does support Seurat, MAST, TDESeq and Tradeseq DE call
                        donor.var = "PID",
                        perm = c("perm1","perm2","perm3", "perm4", "perm5"),
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

    sca <- MAST::SceToSingleCellAssay(Seurat::as.SingleCellExperiment(seurat_analysis))
    if(!is.null(gene_filter)){sca <- sca[gene_filter,]}
    sca <- MAST::filterLowExpressedGenes(sca, min.pt)
    options("mc.cores" = ncores)

    # if(length(seurat_analysis@meta.data[,donor.var])> 60000 & any(grepl("longitudinal", method, ignore.case = T))){
    #   print("Warning: Data with more than 60.000 cells causes errors in longitudinal MAST model fitting with ebayes. Setting ebayes to FALSE")
    #   ebayes <- FALSE
    # }else{ebayes <- TRUE}
  }

  ## For TDESeq
  if(any(grepl("TDE", method, ignore.case = T))){
    if(!is.null(gene_filter)){ tde_input <- suppressMessages(subset(seurat_analysis, features = gene_filter))}
    else {tde_input <- seurat_analysis}
    counts<- Seurat::GetAssayData(tde_input,'counts')
    norm.data<- Seurat::GetAssayData(tde_input,'data')
    meta.data<- tde_input@meta.data

    tde_perm <- TDEseq::CreateTDEseqObject(counts = counts, data=norm.data, meta.data = meta.data)
    rm(tde_input)
  }

  ## TradeSeq
  if(any(grepl("tradeseq", method, ignore.case = T))){
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM$workers <- ncores
  }

  DEG.perm <- list()

  for (m in method){
    if(verbose){print(m)}
    DEG.perm[[m]] <- data.frame()

    # Seurat based tests
    if(grepl("seurat", m, ignore.case = T)){
      #Determine the test to use based on the method string
      test <- ifelse(grepl("wilcox", m, ignore.case = T),  "wilcox",
                            ifelse(grepl("MAST", m, ignore.case = T),"MAST", seurat.test.use[[m]]))

      # Include latent variable based on method string
      if (grepl("latent", m, ignore.case = TRUE)){
        latent.var <- donor.var} else {
          latent.var <- NULL }



      # Loop through the permutations
      for (p in perm){
        if(verbose){print(p)}
        #Loop through comparisons
        for(c in comparisons.df$comparison){
          if(verbose){print(paste0("Performing DE call for ", c))}
          Seurat::Idents(seurat_analysis) <- p
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
            tmp$permutation_round <- p
            tmp$time <- str_split(tmp$comparison, "_", simplify = T)[,1]

            DEG.perm[[m]] <- rbind(DEG.perm[[m]], tmp)
          }}}} #close seurat loop


    # longitudinal MAST based tests
    if(!grepl("seurat", m, ignore.case = T) & grepl("MAST", m, ignore.case = T)& grepl("longitudinal", m, ignore.case = T)){

      zlmTime_perm <- list()
      if(load_model){zlmTime_perm <- readRDS(model[[m]])}
      else{
      for (p in perm){


        ifelse(!is.numeric(colData(sca)[,p]), x <- paste0(p, "_num"), x <- p) # numerical input is required
        # Create formula based on method string
        formula <- as.formula(ifelse(grepl("fixed", m, ignore.case = T)|grepl("FE", m, ignore.case = T),  paste("~", x, "+", donor.var),
                                     ifelse(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T),paste("~", x, "+ (1|", donor.var, ")"),
                                            paste("~", x))))

        if(verbose){print(paste("Fitting model for", p, "with formula:", deparse(formula)))}

        # Fit the MAST model

        zlmTime_perm[[p]] <- if(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T)){
          MAST::zlm(formula, sca, parallel = T, method = "glmer", ebayes = F)}
        else {
          MAST::zlm(formula, sca = sca, parallel = T, ebayes = ebayes)}



      }
        if(save_model){saveRDS(zlmTime_perm, file = paste0(outdir, "zlm_",m,"_perm.rds"))}
        }

      for (p in perm){
        if(verbose){print(paste("Doing DE call for", p))}
        x <- paste0(p, "_num")
        summaryCond_perm <- MAST::summary(zlmTime_perm[[p]], doLRT= x, parallel = T)


        summaryDt_perm <- summaryCond_perm$datatable
        tmp <- data.frame(merge(summaryDt_perm[contrast==x & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                summaryDt_perm[contrast==x & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')) #logFC coefficients

        colnames(tmp) <- c("gene", "p_val", "avg_log2FC", "ci.hi", "ci.lo")
        tmp <- tmp[!is.na(tmp$p_val),]
        tmp$regulation <- ifelse(tmp$avg_log2FC > 0, "up", "down")
        rownames(tmp) <- tmp$gene
        tmp$p_val_adj <- p.adjust(tmp$p_val, 'fdr')
        tmp$permutation_round <- p
        tmp$comparison <- "longitudinal"


        DEG.perm[[m]] <- rbind(DEG.perm[[m]], tmp)
      }}# closing MAST loop


    # MAST pairwise
    if(!grepl("seurat", m, ignore.case = T) & grepl("MAST", m, ignore.case = T)& !grepl("longitudinal", m, ignore.case = T)){

      zlmTime_perm <- list()
      if(load_model){zlmTime_perm <- readRDS(model[[m]])}
      else{
      for (p in perm){
        for(control in unique(comparisons.df$group2)){
          x <- p # categorical input is required
          colData(sca)[,x] <- factor(colData(sca)[,x])
          colData(sca)[,x] <- relevel(colData(sca)[,x],control)
          # Create formula based on method string
          formula <- as.formula(ifelse(grepl("fixed", m, ignore.case = T)|grepl("FE", m, ignore.case = T),  paste("~", x, "+", donor.var),
                                       ifelse(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T),paste("~", x, "+ (1|", donor.var, ")"),
                                              paste("~", x))))

          if(verbose){print(paste("Fitting model for reference level:",control, "and", p, "with formula:", deparse(formula)))}

          #Fit the MAST model
          zlmTime_perm[[paste0(p,control)]] <- if(grepl("random", m, ignore.case = T)|grepl("RE", m, ignore.case = T)){
            MAST::zlm(formula, sca, parallel = T, method = "glmer", ebayes = F)}
          else {
            MAST::zlm(formula, sca = sca, parallel = T, ebayes = ebayes)}


          if(save_model){saveRDS(zlmTime_perm, file = paste0(outdir, "zlm_",m,"_perm.rds"))}
        }}}

      for (p in perm){
        for (c in unique(comparisons.df$comparison)){
          if(verbose){print(paste("Doing DE call for",c, "and" ,p))}
          a <- paste0(p,comparisons.df[comparisons.df$comparison == c,]$group1)
          summaryCond_perm <- MAST::summary(zlmTime_perm[[paste0(p,comparisons.df[comparisons.df$comparison == c,]$group2)]], doLRT= a , parallel = T)


          summaryDt_perm <- summaryCond_perm$datatable
          tmp <- data.frame(merge(summaryDt_perm[contrast==a & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                                  summaryDt_perm[contrast==a & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')) #logFC coefficients

          colnames(tmp) <- c("gene", "p_val", "avg_log2FC", "ci.hi", "ci.lo")
          tmp <- tmp[!is.na(tmp$p_val),]
          tmp$regulation <- ifelse(tmp$avg_log2FC > 0, "up", "down")
          rownames(tmp) <- tmp$gene
          tmp$p_val_adj <- p.adjust(tmp$p_val, 'fdr')
          tmp$permutation_round <- p
          tmp$comparison <- c


          DEG.perm[[m]] <- rbind(DEG.perm[[m]], tmp)
        }}}# closing MAST pairwise loop

    # TDESeq tests
    if(grepl("TDE", m, ignore.case = T)){
      ifelse(grepl("linear", m, ignore.case = T), tde_method_use <- tde_method["linear"], tde_method_use <- tde_method["mixed"] )

      model <- ifelse(grepl("mixed", m, ignore.case = T), "lmm", "lm")

      #Loop through the permutations
      for (p in perm){
        if(verbose){print(p)}
        tde_param <- list(sample.var = donor.var,
                          stage.var = p,
                          fit.model = model, # This defines the model to be linear or mixed
                          tde.thr = 0.05,
                          pct = min.pt,
                          mod = 'FastLMM')

        tde_perm <- TDEseq::tdeseq(object = tde_perm, tde.method=tde_method_use, tde.param=tde_param, num.core= ncores)
        tmp <- TDEseq::GetTDEseqAssayData(tde_perm, slot='tde')
        tmp <- tmp[!is.na(tmp$pvalue),]
        tmp$permutation_round <- p
        tmp$regulation <- ifelse(tmp$patter %in% c("Recession","Trough"),"down",
                                 ifelse(tmp$patter %in% c("None"), "None",
                                        "up"))
        tmp$comparison <- "longitudinal"
        DEG.perm[[m]] <- rbind(DEG.perm[[m]], tmp )
      }}#closing TDESeq loop


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


      gam_perm <- list()
      if(load_model){gam_perm <- readRDS(model[[m]])}
      else{
      #Loop through perm to fit model
      for (p in perm){
        set.seed(42)
        ifelse(!is.numeric(seurat_use@meta.data[,p]), x <- paste0(p, "_num"), x <- p) # numerical input is required
        formula <- ifelse(is.null(U), paste("~", x), paste("~", x, "+ PID"))
        if(verbose){print(paste("Fitting model with formula:", formula))}
        ifelse(is.null(cellWeights),cellWeights <- data.frame("lineage" = rep(1, length(seurat_analysis@meta.data[,x])), row.names = names(seurat_analysis@meta.data[,x])), cellWeights <- cellWeights) # To define the lineages

        gam_perm[[p]] <- tradeSeq::fitGAM(counts = seurat_use@assays$RNA@counts,
                                          pseudotime = seurat_use@meta.data[,x],
                                          cellWeights = cellWeights,
                                          U = U,
                                          nknots = nknots,
                                          parallel = TRUE,
                                          BPPARAM = BPPARAM,
                                          verbose = verbose)
      }
      if(save_model){saveRDS(gam_perm, file = paste0(outdir,"GAM_", m, "_perm.rds"))} }



      #Loop through perm for DE testing
      for (p in perm){
        if(verbose){print(paste("Doing DE call for", p))}
        tmp <- tradeSeq::associationTest(gam_perm[[p]])
        tmp <- tmp[!(is.na(tmp$pvalue)), ]
        tmp$padj <- p.adjust(tmp$pvalue, method = "fdr")
        tmp <- tmp[order(tmp$padj), ]
        tmp$gene <- rownames(tmp)
        tmp$permutation_round <- p
        tmp$regulation <- NA
        tmp$comparison <- "longitudinal"

        DEG.perm[[m]] <- rbind(DEG.perm[[m]], tmp)

      }}# close tradeseq loop


    if(save_DEG){saveRDS(DEG.perm[[m]], file = paste0(outdir, "DEG_",m,"_perm.rds"))}
    invisible(gc())}#closing methods loop
  return(DEG.perm)
}




#' Plot number of DEG genes per permutation
#'
#' @param DEG.perm Output of Permutation DE call or load DEG
#' @param alpha Significance cut off
#'
#' @return Plot number of DEG genes per permutation
#' @export plot_DEG
#'
#' @examples \dontrun{plot_DEG(DEG.perm = DEG.perm)}
plot_DEG <- function(DEG.perm = DEG.perm,
                     alpha = 0.05){
  for (n in names(DEG.perm)){
    plt_list <- list()
    data_plot <- DEG.perm[[n]]
    pvalue <- grep("^p.?val(ue)*((?!.adj))", names(data_plot), perl = T, value = T)
    data_plot.sign <- data_plot[data_plot[,pvalue] <= alpha, ]

    plt_list[[paste0(n, "_p1")]]  <- ggplot(data_plot, aes(comparison, fill=regulation)) +
      geom_bar(position = position_dodge()) +
      scale_fill_manual(values = c(up = "firebrick4", down = "dodgerblue3")) +
      theme_bw() + xlab("")+ylab("Number of DE genes")+
      geom_text(stat='count', aes(label=after_stat(count)), position=position_dodge2(width=0.8), size=5, vjust=0.5, hjust=0, angle= 90) +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0, size =14), axis.text.y = element_text(size =14), axis.title=element_text(size=14), strip.text.x = element_text(size = 14)) +
      ggtitle("All genes tested") +
      scale_y_continuous(expand = expansion(mult = 0.1))+
      facet_wrap(.~ permutation_round, ncol = 5, scales = "free_x")

    plt_list[[paste0(n, "_p2")]]<- ggplot(data_plot.sign, aes(comparison, fill=regulation))+
      geom_bar(position = position_dodge())+
      scale_fill_manual(values = c(up = "firebrick4", down = "dodgerblue3"))+
      theme_bw() + xlab("")+ylab("Number of DE genes")+
      geom_text(stat='count', aes(label=after_stat(count)), position=position_dodge2(width=0.8), size=5, vjust=0.5, hjust=0, angle= 90)+
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0, size =14), axis.text.y = element_text(size =14), axis.title=element_text(size=14), strip.text.x = element_text(size = 14)) +
      ggtitle(paste0("Permutation DEGs (Pvalue <= ", alpha, ")"))+
      scale_y_continuous(expand = expansion(mult = 0.1))+
      facet_wrap(.~ permutation_round, ncol = 5, scales = "free_x")

    p <- ggpubr::ggarrange(plotlist = plt_list, ncol = 2, common.legend = T, legend = "right")
    plot(annotate_figure(p, top = text_grob(n,
                                            color = "black", size = 16), fig.lab.pos = c("top.left")))
  }}

#' Plot number of DEG genes  for power DE call
#'
#' @param DEG.power Output of Permutation DE call or load DEG
#' @param alpha Significance cut off
#'
#' @return Plot number of DEG genes  for power DE call
#' @export plot_DEG_power
#'
#' @examples \dontrun{plot_DEG_power(DEG.perm = DEG.perm)}
plot_DEG_power <- function(DEG.power = DEG.power,
                     alpha = 0.05){
  for (n in names(DEG.power)){
    plt_list <- list()
    data_plot <- DEG.power[[n]]
    pvalue <- grep("^p.?val(ue)*((?!.adj))", names(data_plot), perl = T, value = T)
    data_plot.sign <- data_plot[data_plot[,pvalue] <= alpha, ]

    plt_list[[paste0(n, "_p1")]]  <- ggplot(data_plot, aes(comparison, fill=regulation)) +
      geom_bar(position = position_dodge()) +
      scale_fill_manual(values = c(up = "firebrick4", down = "dodgerblue3")) +
      theme_bw() +
      geom_text(stat='count', aes(label=after_stat(count)), position=position_dodge2(width=0.8), size=5, vjust=0.5, hjust=0, angle= 90) +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0, size =14), axis.text.y = element_text(size =14), axis.title=element_text(size=14), strip.text.x = element_text(size = 14)) +
      ggtitle("All genes tested") +
      scale_y_continuous(expand = expansion(mult = 0.25))+
      xlab("")+ylab("Number of DE genes")

    plt_list[[paste0(n, "_p2")]]<- ggplot(data_plot.sign, aes(comparison, fill=regulation))+
      geom_bar(position = position_dodge())+
      scale_fill_manual(values = c(up = "firebrick4", down = "dodgerblue3"))+
      theme_bw()+
      geom_text(stat='count', aes(label=after_stat(count)), position=position_dodge2(width=0.8), size=5, vjust=0.5, hjust=0, angle= 90) +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0, size =14), axis.text.y = element_text(size =14), axis.title=element_text(size=14), strip.text.x = element_text(size = 14)) +
      ggtitle(paste0("DEGs (Pvalue <= ", alpha, ")"))+
      scale_y_continuous(expand = expansion(mult = 0.25))+
      xlab("")+ylab("Number of DE genes")

    p <- ggarrange(plotlist = plt_list, ncol = 2, common.legend = T, legend = "right")
    plot(annotate_figure(p, top = text_grob(n,
                                            color = "black", size = 16), fig.lab.pos = c("top.left")))
  }}

#' QQ plot of pvalues
#'
#' @param DEG.perm Output of Permutation DE call or load DEG
#' @param pvalue Name of the pvalue column. If NULL the pvalue column is generated based on regular expression.
#' @param color color vector
#'
#' @return QQ plot of pvalues
#' @import ggplot2
#'
#' @export Perm_QQplot
#'
#' @examples \dontrun{Perm_QQplot(DEG.perm = DEG.perm)}
Perm_QQplot <- function(DEG.perm = list("seurat_wilcox" = DEG.seurat.perm.null,"pseudobulk_DESeq2_Wald"= DEG.pseudobulk_Wald.perm.null),
                        pvalue = NULL,
                        plot_anno = T,
                        color = c("tradeseq_time"= "#c6dbef", "tradeseq_slingshot" = "#2171b5", "tradeseq_psupertime"= "#032f5f",
                                  "pseudobulk_DESeq2_Wald" = "#A65628","pseudobulk_DESeq2_LRT" = "#ef7c00","seurat_wilcox" = "#e7bd00",
                                  "seurat_MAST" = "#fb6a4a","MAST_RE" = "#cb181d","seurat_MAST_latent" ="#a50f15" ,
                                  "MAST_longitudinal"= "#a1d99b",  "MAST_longitudinal_RE"= "#006d2c" ,"MAST_longitudinal_PID"= "#00441b",
                                  "linear_TDE" = "#ffb3ff", "mixed_TDE" =  "#984EA3")
){

  data_plot <- data.frame()

  if(is.null(pvalue)){
    create_pvalue <- TRUE }
  else{create_pvalue <- FALSE}

  # Generate data for plot
  for (m in names(DEG.perm)){
    DEG_df <- DEG.perm[[m]]

    if(create_pvalue){
      pvalue <- grep("^p.?val(ue)*((?!.adj))", names(DEG_df), perl = T, value = T)}
    else{ pvalue <- pvalue }


    data_plot<- rbind(data_plot,
                      data.frame("method" = m,
                                 "p_val_null" = -log10(sort(ppoints(length(DEG_df[,pvalue])))),
                                 "p_val_test"  = -log10(sort(DEG_df[,pvalue]))
                      ))
  }

  data_plot$method <- factor(data_plot$method, levels = names(DEG.perm))

  p <- ggscatter(data = data_plot, x = "p_val_null",y = "p_val_test", group.by = "method", color = "method")+
    geom_abline(intercept = 0, slope = 1, show.legend = T)+
    ylim(0,5)+
    scale_color_manual(values = color)+
    xlab("-log10 expected pvalues")+ ylab("-log10 observed pvalues")+
    theme(legend.position = "right", legend.title = element_blank(), axis.title=element_text(size=14))+
    ggtitle("Pvalues under the null hypothesis")

  if(plot_anno){
    p <- p +
      annotate("label", label = "conservative", x = 3, y=0.5)+
      annotate("label", label = "lenient", x = 0.35, y=3)
    }

  return(p)
}
