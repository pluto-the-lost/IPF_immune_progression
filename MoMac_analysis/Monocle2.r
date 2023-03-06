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

rds<-readRDS('MonMac_cds.rds')
colnames(pData(rds))[1]<-'cell_type'
pData(rds)$cell_type<-plyr::mapvalues(x=pData(rds)$cell_type,from =as.vector(c(0,1,2,3,4)),to=as.vector(c('Cx3cr1+','CD36+','C5ar1+','F13a1+','Isg15+Ly6a+')))
pdf('MonMac_monocle2.pdf',w=6,h=5)
p1<-plot_cell_trajectory(rds, color_by = "cell_type",cell_size = 2.5)+ scale_color_manual(breaks = waiver(),values=colorRampPalette(brewer.pal(8, "Set1"))(5)[1:5])
print (p1)
dev.off()

pdf('MonMac_time.pdf',w=6,h=5)
p1<-plot_cell_trajectory(rds, color_by = "sampleID",cell_size = 2.5)+ scale_color_manual(breaks = waiver(),values=c('#1F3563','#26629D','#4E95C1','#93C1D9'))
print (p1)
dev.off() 

pdf('MonMac_ptime.pdf',w=6,h=5)
p1<-plot_cell_trajectory(rds, color_by = "Pseudotime",cell_size = 2.5)
print (p1)
dev.off() 



