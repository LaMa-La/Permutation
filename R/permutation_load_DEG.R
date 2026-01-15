#' Load to DEG perm function
#'
#' @param files List of files to load
#'
#' @return List of DE call data frames
#' @export load_DEG
#'
#' @examples \dontrun{
#' DEG.perm <- load_DEG(c("seurat_wilcox" = "/data/Femmunity/analysis/laura/external_data/PALMO/v3/241212_DEG_seurat_wilcox_perm.rds"))
#' }

load_DEG <- function(files = c("seurat_wilcox" = "DEG_seurat_wilcox_perm.rds")){
  DEG <- list()

  for (f in names(files)){
    DEG[[f]] <- readRDS(files[[f]])
  }

  return(DEG)
}
