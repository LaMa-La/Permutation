# Permutation workflow

The permutation package offers functions to determine the FDR and p-value calibration in multiple-time point scRNA seq data. It takes the data in the commonly used Seurat object format and carries out permutation of time labels. Further, FDR can be calculated after DE analysis with Seurat (Wilcoxon rank sum), MAST, DESeq2, TDESeq and Tradeseq. 


# Install package 
```
# install.packages("devtools") 
devtools::install_git("LaMa-La/Permutation", upgrade = "never")
```
The package is optimized for the use with of the docker container from Jonas Schulte-Schrepping (https://github.com/jsschrepping/r_docker/tree/r_docker_4.3.0#)

# Example workflow 

```{r}
library(future)
library(future.apply)
library(parallel)

library(Permutation)
library(Seurat)
library(ggplot2)
library(tidydr)
library(patchwork)
library(DESeq2)
library(psupertime)
library(SingleCellExperiment)
```


```{r}
color_clusters <- c(RColorBrewer::brewer.pal(n = 9,name = "Set1"),
                    RColorBrewer::brewer.pal(n = 8,name = "Set2"),
                    RColorBrewer::brewer.pal(n = 12,name = "Set3"),
                    RColorBrewer::brewer.pal(n = 12,name = "Paired"),
                    RColorBrewer::brewer.pal(n = 9,name = "Pastel1"),
                    RColorBrewer::brewer.pal(n = 8,name = "Pastel2"),
                    RColorBrewer::brewer.pal(n = 8,name = "Accent"))

color_method <- c("tradeseq_time"= "#c6dbef", "tradeseq_slingshot" = "#2171b5", "tradeseq_psupertime"= "#032f5f", 
                  "pseudobulk_DESeq2_Wald" = "#A65628","pseudobulk_DESeq2_LRT" = "#ef7c00","seurat_wilcox" = "#e7bd00",  
                  "seurat_MAST" = "#fb6a4a","MAST_RE" = "#cb181d","seurat_MAST_latent" ="#a50f15" , 
                  "MAST_longitudinal"= "#a1d99b",  "MAST_longitudinal_RE"= "#006d2c" ,"MAST_longitudinal_PID"= "#00441b","MAST_longitudinal_fixed"= "#00441b",
                  "linear_TDE" = "#ffb3ff", "mixed_TDE" =  "#984EA3")

```



```{r}
cl <- parallel::makeCluster(8)
registerDoParallel(cl)
plan("multisession", workers = 5)
options(future.globals.maxSize= 16000 * 1024^2)
```

------
------

# Prep 

```{r}
# Download example data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190992
seurat_analysis <- readRDS("/input/PALMO_seurat.rds")
```

```{r}
color_time <- viridis::viridis(6)
names(color_time) <- c("W2", "W3", "W4", "W5", "W6", "W7")
```

```{r fig.width=12, fig.height=5}
p1 <- Seurat::DimPlot(seurat_analysis, 
              reduction = "umap", 
              group.by = "celltype", 
              cols = color_clusters,
              shuffle = TRUE)+ ggplot2::coord_fixed() + theme_dr() + ggplot2::theme(panel.grid = ggplot2::element_blank())

p2 <- Seurat::DimPlot(seurat_analysis, 
              reduction = "umap", 
              group.by = "time", 
              cols = color_time,
              shuffle = TRUE)+ ggplot2::coord_fixed() + theme_dr() + ggplot2::theme(panel.grid = ggplot2::element_blank())


p3 <- Seurat::DimPlot(seurat_analysis, 
              reduction = "umap", 
              group.by = "PID", 
              cols = color_clusters,
              shuffle = TRUE)+ ggplot2::coord_fixed() + theme_dr() + ggplot2::theme(panel.grid = ggplot2::element_blank())

p1 + p2 + p3 + plot_layout(ncol = 3)
```

```{r}
# Path to your output directory
outdir <- "/data/output/"
```

```{r}
seurat_analysis <- subset(seurat_analysis, subset = celltype %in% c("CD14 Mono"))

seurat_analysis$time <- seurat_analysis$time #assign your time variable to this slot
seurat_analysis$time_num <- as.numeric(stringr::str_split(seurat_analysis$time,"W", simplify = T)[,2])  #put a numeric version of your time variable in this slot, 0 are not possible 
seurat_analysis$PID <- seurat_analysis$PID #assign your donor variable to this slot

#Create Donor time combined variable
seurat_analysis$PID_time <- paste(seurat_analysis$PID, seurat_analysis$time, sep = "_")
seurat_analysis$PID_time_num <- paste(seurat_analysis$PID, seurat_analysis$time_num, sep = "_")

## assign your pseudotime variable if needed
#seurat_analysis$pseudotime <- seurat_analysis$slingshot 
```

```{r, message=FALSE}
Seurat::Idents(seurat_analysis) <- "time"

## variable genes
seurat_analysis <- Seurat::FindVariableFeatures(object = seurat_analysis, assay="RNA", selection.method = 'vst', nfeatures=5000)
```

# Permutation workflow

```{r}
seurat_analysis <- permutation(seurat_analysis, method = "longitudinal", numeric = "custom", numeric.function = function(x){as.numeric(str_split(x, pattern =  "W", simplify = T)[,2])})
```


```{r}
corr.plot.permutation(seurat_analysis, time.var = "time", numeric.function = function(x){as.numeric(str_split(x, pattern =  "W", simplify = T)[,2])} )
```

```{r, fig.width= 10}
plot.donor.permutation(seurat_analysis, color = color_time)
```

## Permutation DE call

```{r}
dds_perm <- generate_pseudobulk_perm(seurat_analysis = seurat_analysis)
```

```{r}
DEG.perm <- Perm_DEcall_pseudobulk(dds_perm, method = c("pseudobulk_DESeq2_Wald","pseudobulk_DESeq2_LRT"),
comparison.group1 = c("W3", "W4", "W5", "W6","W7"),
comparison.group2 = rep("W2", 5),
 save_DEG = F, outdir = outdir, shrinkType = "apeglm")
```

```{r}
DEG.perm <- Perm_DEcall(seurat_analysis = seurat_analysis, method = c(
"seurat_wilcox",
"linear_TDE",
"tradeseq_latent"
 ),
comparison.group1 = c("W3", "W4", "W5", "W6","W7"),
comparison.group2 = rep("W2", 5),
 save_DEG = F, ncores = 5, gene_filter = NULL, tde_method = c("linear" = "cell", "mixed"="pseudocell"), outdir = outdir)
```

## Permutation pseudotime

```{r}
Idents(seurat_PALMO_mono) <- "orig.ident"
seurat_analysis_small <- subset(seurat_analysis, downsample = 1000, seed = 42)
sce <- as.SingleCellExperiment(seurat_analysis_small)
sce@colData$time <- factor(sce@colData$time, levels = sort(unique(sce@colData$time)))
```

```{r}
library(psupertime)
```

```{r}
y <- sce@colData$time
psuper_palmo <-psupertime(sce@assays@data$logcounts, y, sel_genes='list',gene_list = seurat_PALMO_mono@assays$RNA@var.features, seed=1234)
psuper_palmo
```

```{r}
seurat_analysis_small <- AddMetaData(object = seurat_analysis_small, metadata = psuper_palmo$proj_dt$psuper, col.name = "psupertime")
```

```{r}
seurat_analysis_small <- permutation(seurat_analysis_small, method = "longitudinal", time.var = "psupertime", perm_prefix = "psupertime_")
```

```{r}
plot.pseudotime.permutation(seurat_analysis_small, perm = c("psupertime", "psupertime_perm1","psupertime_perm2","psupertime_perm3","psupertime_perm4","psupertime_perm5"), color = list("time" = color_time, "donor" = color_clusters))
```


```{r}
corr.plot.permutation(seurat_analysis_small, time.var = "psupertime", perm.list = c("psupertime_perm1", "psupertime_perm2", "psupertime_perm3", "psupertime_perm4", "psupertime_perm5"),)
```

```{r}
DEG.perm <- Perm_DEcall(seurat_analysis = seurat_analysis_small, method = c(
"tradeseq_psupertime_latent"
 ),
perm = c("psupertime_perm1", "psupertime_perm2", "psupertime_perm3", "psupertime_perm4", "psupertime_perm5"),
 save_DEG = F, ncores = 5, gene_filter = seurat_analysis_small@assays$RNA@var.features[1:50], outdir = outdir)
```

## Plots

```{r, fig.width=15, fig.height=8}
plot_DEG(DEG.perm = DEG.perm)
```

```{r, fig.width= 10, fig.height=5}
Perm_QQplot(DEG.perm = DEG.perm  , color = color_method)
```

## FPR df

```{r}
fpr_df <- perm_FPR_df(DEG.perm = DEG.perm)
```

```{r}
p1 <- plot_FPR_box(fpr_df)
p1  
```

```{r}
p2 <- plot_false_positives(fpr_df)
p2
```

```{r, fig.height= 10, fig.width= 10}
p1 + p2 + plot_layout(nrow= 2)
```

```{r}
plot_FPR_point(fpr_df)
```

## Detection power

```{r}
outdir <- "/data/outdir/detection_power/"
```

```{r}
dds <- generate_pseudobulk(seurat_analysis)
```

```{r}
DEG.power <- Power_DEcall_pseudobulk(dds, method = c("pseudobulk_DESeq2_Wald", "pseudobulk_DESeq2_LRT"),
comparison.group1 = c("W3", "W4", "W5", "W6","W7"),
comparison.group2 = rep("W2", 5),
 save_DEG = F, outdir = outdir, shrinkType = "apeglm")
```

```{r}
DEG.power <- Power_DEcall(seurat_analysis = seurat_analysis, method = c(
"seurat_wilcox",
"linear_TDE",
"mixed_TDE",
"seurat_MAST",
"seurat_MAST_latent",
"MAST_longitudinal",
"MAST_longitudinal_RE",
"MAST_longitudinal_fixed",
"tradeseq_latent"
 ),
comparison.group1 = c("W3", "W4", "W5", "W6","W7"),
comparison.group2 = rep("W2", 5),
gene_filter = seurat_analysis@assays$RNA@var.features[1:50], tde_method = c("linear" = "cell", "mixed"="pseudocell"),
save_DEG = F, ncores = 5,  outdir = outdir)
```

```{r}
avg_power <- detection_power_padj(DEG.power)
```

```{r, fig.height= 5}
plot_detection_power(avg_power)
```

------
------

```{r}
sessionInfo()
```




