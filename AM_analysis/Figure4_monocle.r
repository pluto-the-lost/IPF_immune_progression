require(monocle)
require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(RColorBrewer)

rds<-readRDS('AM_cds.rds')
colnames(pData(rds))[1]<-'cell_type'
pData(rds)$cell_type<-plyr::mapvalues(x=pData(rds)$cell_type,from =as.vector(c(0,1,2,3,4,5,6)),to=as.vector(c('Fabp5+Gpnmd+','Lung resident Ams','CCL17-producing','Nr4a1+','Fabp1-Car4-Krt79+Cidec+','Fabp1+Car4+Krt79+Cidec+','Isg15+Ly6e+')))
pdf('AM_monocle2.pdf',w=6,h=5)
p1<-plot_cell_trajectory(rds, color_by = "cell_type",cell_size = 2.5)+ scale_color_manual(breaks = waiver(),values=colorRampPalette(brewer.pal(9, "Set3"))(7)[1:7])
print (p1)
dev.off()

pdf('AM_time.pdf',w=6,h=5)
p1<-plot_cell_trajectory(rds, color_by = "sampleID",cell_size = 2.5)+ scale_color_manual(breaks = waiver(),values=c('#1F3563','#26629D','#4E95C1','#93C1D9'))
print (p1)
dev.off() 

pdf('AM_ptime.pdf',w=6,h=5)
p1<-plot_cell_trajectory(rds, color_by = "Pseudotime",cell_size = 2.5)
print (p1)
dev.off() 



