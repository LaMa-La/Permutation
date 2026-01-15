#' FPR data frame
#'
#' @param DEG.perm Output of Permutation DE call or load DEG
#' @param permutations Name of all permutations in the metadata of the seurat object
#' @param permutation_col Name of the permutation column in the DEG data frame
#' @param pvalue Name of the pvalue column. If NULL the pvalue column is generated based on regular expression.
#' @param padj Name of the adjusted pvalue column. If NULL the adjusted pvalue column is generated based on regular expression.
#'
#' @return Data frame of FPRs
#' @export perm_FPR_df
#'
#' @examples \dontrun{fpr_df <- perm_FPR_df(DEG.perm = DEG.perm)}
perm_FPR_df <- function(DEG.perm = list("seurat_wilcox" = DEG.seurat.perm.null,"pseudobulk_DESeq2_Wald"= DEG.pseudobulk_Wald.perm.null),
                        permutations = c("perm1","perm2","perm3", "perm4", "perm5"),
                        permutation_col = "permutation_round",
                        pvalue = NULL,
                        padj = NULL
){

  fpr_df <- data.frame()

  if(is.null(pvalue)){
    create_pvalue <- TRUE
  }
  else{create_pvalue <- FALSE}
  if(is.null(padj)){
    create_padj <- TRUE
  }else{create_padj <- FALSE}


  for (m in names(DEG.perm)){
    DEG_df <- DEG.perm[[m]]

    if(create_pvalue){
      pvalue <- grep("^p.?val(ue)*((?!.adj))", names(DEG_df), perl = T, value = T)}
    else{
      pvalue <- pvalue
    }
    if(create_padj){
      padj <- grep(".adj", names(DEG_df), perl = T, value = T)}
    else{
      padj <- padj
    }

    for (p in permutations){
      fpr_df<- rbind(fpr_df,
                     data.frame("method" = m, "perm" = p,"method_perm" = paste(m, p, sep = "_"),
                                "fpr" = round(length(DEG_df[grepl(p, DEG_df[,permutation_col]) & DEG_df[,pvalue] < 0.05,"gene"])/length(DEG_df[grepl(p, DEG_df[,permutation_col]) , "gene"]),3), #fpr = significant genes/all genes tested
                                "false_positive_padj" = length(DEG_df[grepl(p, DEG_df[,permutation_col]) & DEG_df[,padj] < 0.05, "gene"])
                     ))
    }}
  fpr_df$method <- factor(fpr_df$method, levels = names(DEG.perm))
  return(fpr_df)
}

#' FPR boxplot
#'
#' @param fpr_df Output from perm_fpr_df
#' @param color color vector
#'
#' @return Boxplot of FPR
#' @export plot_FPR_box
#'
#' @examples \dontrun{plot_FPR_box(fpr_df)}
plot_FPR_box <- function(fpr_df, color = color_method){
  p <- ggplot(fpr_df, aes(method, fpr, color = method))+
    geom_boxplot()+
    scale_color_manual(values= color)+
    geom_hline(yintercept = 0.05, linetype = "dashed")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1), legend.title = element_blank())+
    ylab("FPR")+ xlab("")
  return(p)
}


#' FPR scatterplot
#'
#' @param fpr_df Output from perm_fpr_df
#' @param color color vector
#'
#' @return Scatter plot of FPR values
#' @export plot_FPR_point
#'
#' @examples \dontrun{plot_FPR_point(fpr_df)}
plot_FPR_point <- function(fpr_df, color = color_method){
  p <- ggplot(fpr_df, aes(method, fpr, color = perm))+
    geom_point()+
    scale_color_manual("Permutation round", values= color_clusters)+
    geom_hline(yintercept = 0.05, linetype = "dashed")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))+
    ylab("FPR")+ xlab("")
  return(p)
}

#' Plot false positives
#'
#' @param fpr_df Output from perm_fpr_df
#' @param color color vector
#'
#' @return Plot of false positives
#' @export plot_false_positives
#'
#' @examples \dontrun{plot_false_positives(fpr_df)}
plot_false_positives <- function(fpr_df, color = color_method){
  aggregated_df <- fpr_df %>%
    group_by(method) %>%
    summarise(
      mean_fpr = mean(fpr),
      sum_false_positive_padj = sum(false_positive_padj)
    )

  p <- ggplot(aggregated_df, aes(method, sum_false_positive_padj, fill = method))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values= color)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1),legend.title = element_blank())+
    geom_text(stat='identity', aes(label=sum_false_positive_padj), position=position_dodge2(width=0.8), size=5, vjust=0, hjust=0.5)+
    ylab("DEG (padj < 0.05)")+ xlab("")+
    ggtitle("False positive DEG (padj < 0.05)")
  return(p)
}
