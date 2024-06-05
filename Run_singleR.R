library(Seurat)
library(ggsci)
library(tidyverse)
library(purrr)
library(dplyr)
library(plyr)
library(RColorBrewer)
require(scales)
library(clustree)
library(ggplot2)
library(sctransform)
library(data.table)
library(future)
library(future.apply)
library(biomaRt)
library(lattice)
library('org.Hs.eg.db')
library(DT)
library(singleR)
library(celldex)

set.seed(42)

macs.integrated <- readRDS("macsIntegrationObject_10_10_BPH.RDS")
T.integrated <- readRDS("../Tcells/T_IntegrationObject_10_10_BPH.RDS")

mac_colors <- c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499","#44AA99","#999933","#882255","#661100","#6699CC","#888888","deeppink3","darkslategrey")
my_color_palette.t <- colorRampPalette(brewer.pal(n=8, name="Dark2"))(14)


ref <- celldex::HumanPrimaryCellAtlasData()

macs.data <- macs.integrated@assays$RNA@data 
T.data <- T.integrated@assays$RNA@data 

#singleR on individual cells
singler <- SingleR(test = macs.data, ref = ref, assay.type.test=1,  labels = ref$label.main)
macs.integrated$celltype <- singler$pruned.labels
saveRDS(macs.integrated,"../macs/macsIntegrationObject_10_10_BPH_v2_singleR.RDS")
pdf("mac_Cell_identification_singleR.pdf",height=6,width=8)
DimPlot(macs.integrated,group.by = "celltype")
dev.off()
pdf("mac_Cell_identification_singleR_labeled.pdf",height=6,width=8)
DimPlot(macs.integrated,group.by = "celltype",label=TRUE,repel=TRUE,label.size = 6)
dev.off()

singler.T <- SingleR(test = T.data, ref = ref, assay.type.test=1,  labels = ref$label.main)
T.integrated$celltype <- singler.T$pruned.labels
saveRDS(T.integrated,"../Tcells/T_IntegrationObject_10_10_BPH.RDS")
pdf("T_Cell_identification_singleR.pdf",height=6,width=8)
DimPlot(T.integrated,group.by = "celltype")
dev.off()
pdf("T_Cell_identification_singleR_labeled.pdf",height=6,width=8)
DimPlot(T.integrated,group.by = "celltype",label=TRUE,repel=TRUE,label.size = 6)
dev.off()

#singleR on clusters
T.clusters <- Idents(T.integrated)
mac.clusters <- Idents(macs.integrated)
singler.macs.c <- SingleR(test = macs.data, ref = ref, assay.type.test=1,  labels = ref$label.main, clusters=mac.clusters)
singler.T.c <- SingleR(test = T.data, ref = ref, assay.type.test=1,  labels = ref$label.main, clusters=T.clusters)


#hacky way of doing this but it works.improve in future
T.integrated <- RenameIdents(T.integrated, `0`=singler.T.c$pruned.labels[1],`1`=singler.T.c$pruned.labels[2],
                             `2`=singler.T.c$pruned.labels[3],`3`=singler.T.c$pruned.labels[4],
                             `4`=singler.T.c$pruned.labels[5],`5`=singler.T.c$pruned.labels[6],
                             `6`=singler.T.c$pruned.labels[7],`7`=singler.T.c$pruned.labels[8],
                             `8`=singler.T.c$pruned.labels[9],`9`=singler.T.c$labels[10],
                             `10`=singler.T.c$pruned.labels[11],`11`=singler.T.c$pruned.labels[12],
                             `12`=singler.T.c$pruned.labels[13],`13`=singler.T.c$pruned.labels[14])
macs.integrated <- RenameIdents(macs.integrated, `0`=singler.macs.c$pruned.labels[1],`1`=singler.macs.c$pruned.labels[2],
                             `2`=singler.macs.c$pruned.labels[3],`3`=singler.macs.c$pruned.labels[4],
                             `4`=singler.macs.c$pruned.labels[5],`5`=singler.macs.c$pruned.labels[6],
                             `6`=singler.macs.c$pruned.labels[7],`7`=singler.macs.c$pruned.labels[8],
                             `8`=singler.macs.c$pruned.labels[9],`9`=singler.macs.c$labels[10],
                             `10`=singler.macs.c$pruned.labels[11],`11`=singler.macs.c$pruned.labels[12],
                             `12`=singler.macs.c$pruned.labels[13],`13`=singler.macs.c$pruned.labels[14])

