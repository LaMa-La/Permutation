#' Correltation matrix between original and permuted data
#'
#' @param seurat_analysis seurat object
#' @param time.var Name of time label column in the metadata of the seurat object
#' @param perm.list Name of all permutation columns in the metadata of the seurat object to be plotted
#' @param numeric.function function to generate numeric version of the time and permutation labels
#'
#' @return Correlation matrix plot
#'
#'
#' @export corr.plot.permutation
#'
#' @examples \dontrun{corr.plot.permutation(seurat_analysis, time.var = "time",
#' numeric.function = function(x){as.numeric(str_split(x, pattern =  "W", simplify = T)[,2])} )}
corr.plot.permutation <- function(seurat_analysis,
                                  time.var = "time_num",
                                  perm.list = c("perm1", "perm2", "perm3", "perm4", "perm5"),
                                  numeric.function = function(x){as.numeric(x)}
){

  #Check if there is residual correlation between the original time label and the permutations

  # create dataframe
  cor <- data.frame(seurat_analysis[[c(time.var, perm.list)]])

  colnames(cor) <- c("original", paste0("perm", 1:length(perm.list)))
  cor <- data.frame(lapply(cor, numeric.function #data needs to be numeric
  ))

  p <- wrap_elements(~corrplot::corrplot.mixed(cor(cor), tl.col ="black", lower.col ="black"))


  return(p)
}


#' Plot cell numbers per donor
#'
#' @param seurat_analysis seurat object
#' @param donor.var Name of donor ID column in the metadata of the seurat object
#' @param time.var Name of time label column in the metadata of the seurat object
#' @param color color vector
#' @param plots Choose between "relative" for relative cell numbers per donor, "absolute" for absolute cell numbers per donor
#'
#' @return Plot of cell numbers per donor
#'
#'
#' @export plot.donor.permutation
#'
#'
#'
#' @examples \dontrun{plot.donor.permutation(seurat_analysis, color = color_time)}
plot.donor.permutation <- function(seurat_analysis,
                                   donor.var = "PID",
                                   time.var = c("time", "perm1", "perm2", "perm3", "perm4", "perm5"),
                                   color = brewer.pal(n = 9,name = "Set1"),
                                   plots = c("relative", "absolute")){

  plt_list <- list()

  if("absolute" %in% plots){
    for (p in time.var){
      plt_list[[paste0("absolute_",p)]] <- ggplot(seurat_analysis@meta.data, aes(x=.data[[donor.var]], fill=.data[[p]])) +
        geom_bar() +
        scale_fill_manual(values=color) +
        theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), legend.position = "right") +
        xlab("Donor") + ylab("Number of cells")+ ggtitle(p)

    }}

  if("relative" %in% plots){
    for (p in time.var){
      plt_list[[paste0("relative_",p)]]<- ggplot(seurat_analysis@meta.data, aes(x= .data[[donor.var]], fill=.data[[p]])) +
        geom_bar(position = "fill") +
        scale_fill_manual(values=color) +
        theme_bw() + theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), legend.position = "right") +
        xlab("Donor") + ylab("Percentage of cells")+ ggtitle(p)
    }}

  ggpubr::ggarrange(plotlist = plt_list, ncol = 3, nrow = length(plt_list)/3, common.legend = T, legend = "right")

}



#' Plot time and donor variable over pseudotime per permutation
#'
#'
#' @param seurat_analysis seurat object
#' @param donor.var Name of donor ID column in the metadata of the seurat object
#' @param time.var Name of original time label column in the metadata of the seurat object
#' @param perm Name of pseudotime and permutation columns to plot
#' @param color color vector for time and donor
#' @param plots Choose between "time" and "donor" for time or donor labels, respectively
#'
#' @return
#' @export plot.pseudotime.permutation
#'
#' @examples
plot.pseudotime.permutation <- function(seurat_analysis,
                                        donor.var = "PID",
                                        time.var = "time",
                                        perm = c("pseudotime", "pseudotime_perm1", "pseudotime_perm2", "pseudotime_perm3", "pseudotime_perm4", "pseudotime_perm5"),
                                        color = list("time" = brewer.pal(n = 9,name = "Set1"), "donor" = brewer.pal(n = 9,name = "Set2")),
                                        plots = c("time", "donor")){

data_plot <- data.frame()

for (p in perm){
  tmp <- data.frame("perm" = seurat_analysis@meta.data[,p], "permutation_round" = rep(p, length(seurat_analysis@meta.data[,p])), factor(seurat_analysis@meta.data[,time.var]), seurat_analysis@meta.data[,donor.var])
  colnames(tmp) <- c("perm","permutation_round", "time","PID")
  data_plot <- rbind(data_plot, tmp)
}

plt <- list()
if(any(grepl("time", plots, ignore.case = T))){
plt[["p1"]] <- ggplot(as.data.frame(data_plot), aes(x=perm, y = permutation_round, group.by = time,color= time, fill=time)) +
  geom_density_ridges(alpha = 0.5)+
  scale_fill_manual(values=color[["time"]])+
  scale_color_manual(values=color[["time"]])+
  theme_bw()+
  xlab("pseudotime")+ ylab("")+
  ggtitle("Time lables over pseudotime")

}

if(any(grepl("donor", plots, ignore.case = T))){
plt[["p2"]] <- ggplot(as.data.frame(data_plot), aes(x=perm,y = permutation_round, group.by = PID ,color= PID, fill=PID)) +
  geom_density_ridges(alpha = 0.5)+
  scale_fill_manual(values=color[["donor"]])+
  scale_color_manual(values=color[["donor"]])+
  theme_bw()+
  xlab("pseudotime")+ ylab("")+ labs(fill = "donor ID", color = "donor ID")+
  ggtitle("Donor lables over pseudotime")

}
return(plt)
}
