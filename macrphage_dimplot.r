require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)

grid.col<-colorRampPalette(brewer.pal(9, "Set1"))(8)[1:8]
rds<-readRDS('Mac_without_UK.rds')
print(head(rds@meta.data))
DefaultAssay(object = rds) <- "RNA"
Idents(rds)<-rds@meta.data$cell.type
DimPlot(rds, reduction = "umap", group.by = "cell.type",pt.size=2,cols=grid.col,raster=TRUE)
ggsave('Macrophage_dimplot.pdf', dpi=900,width=5, height=4, units="in")
DimPlot(rds, reduction = "umap", group.by = "cell.type",split.by='sampleID',pt.size =2,cols=grid.col,raster=TRUE,ncol=4)
ggsave('Macrophage_dimplot_group.pdf', dpi=900,width=17, height=4, units="in")
