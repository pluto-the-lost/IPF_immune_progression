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

grid.col<-colorRampPalette(brewer.pal(9, "Set2"))(13)[1:13]
rds<-readRDS('/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/TET_PUBLIC/zhaoyue/project/B_TET_074/IPF_dbfree.rds')
print(head(rds@meta.data))
DefaultAssay(object = rds) <- "RNA"
Idents(rds)<-rds@meta.data$cell.type
DimPlot(rds, reduction = "umap", group.by = "cell.type",pt.size=2,cols=grid.col,raster=TRUE)
ggsave('IPF_dimplot.pdf', dpi=900,width=5, height=4, units="in")
DimPlot(rds, reduction = "umap", group.by = "cell.type",split.by='orig.ident',pt.size =2,cols=grid.col,raster=TRUE,ncol=4)
ggsave('IPF_dimplot_group.pdf', dpi=900,width=17, height=4, units="in")
