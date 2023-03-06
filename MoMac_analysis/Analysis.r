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
rds<-readRDS('MonMac.rds')
rds@meta.data$cell_type<-rds@meta.data$seurat_clusters
rds@meta.data$cell_type<-plyr::mapvalues(x=rds@meta.data$cell_type,from =as.vector(c(0,1,2,3,4)),to=as.vector(c('Cx3cr1+','CD36+','C5ar1+','F13a1+','Isg15+Ly6a+')))
saveRDS(rds,'MonMac.anno.rds')
rds<-readRDS('MonMac.anno.rds')
grid.col<-colorRampPalette(brewer.pal(8, "Set1"))(5)[1:5]
DefaultAssay(object = rds) <- "RNA"
#Idents(rds)<-rds@meta.data$cell_type
DimPlot(rds, reduction = "umap", group.by = "cell_type",pt.size=2,cols=grid.col,raster=TRUE)
ggsave('MonMac_dimplot.pdf', dpi=900,width=6, height=4, units="in")
DimPlot(rds, reduction = "umap", group.by = "cell_type",split.by='sampleID',pt.size =2,cols=grid.col,raster=TRUE,ncol=2)
ggsave('MonMac_dimplot_group.pdf', dpi=900,width=10, height=8, units="in")
