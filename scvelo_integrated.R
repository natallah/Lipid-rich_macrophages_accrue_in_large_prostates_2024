library(SeuratWrappers)
library("Seurat")
library(Matrix)
library(dplyr)
library(stringr)
library(SeuratDisk)
library(remotes)


t_cells <- readRDS("T.integrated_ADTnormalized.RDS")
macs <- readRDS("macs.integratedADTnorm.RDS")

t_cells$cluster.type <- paste(Idents(t_cells),t_cells$sample.Type,sep="_")
macs$cluster.type <- paste(Idents(macs),macs$sample.Type,sep="_")

#rename the cells
unique(t_cells@meta.data$orig.ident)

#save metadata
t_cells$barcode <- colnames(t_cells)
t_cells$UMAP_1 <- t_cells@reductions$umap@cell.embeddings[,1]
t_cells$UMAP_2 <- t_cells@reductions$umap@cell.embeddings[,2]
t.meta <- t_cells@meta.data
unique(str_extract(t.meta$barcode,"-\\d+[_\\d+]*|_\\d+[_\\d+]*"))
bar.t.temp <- str_replace(t.meta$barcode,"-\\d+[_\\d+]*|_\\d+[_\\d+]*","")
bar.t.temp2 <- paste(bar.t.temp,t.meta$orig.ident,sep="_")
t.meta$barcode <- bar.t.temp2
dup.t <- t.meta[duplicated(t.meta$barcode),]
t.meta2 <- t.meta[!duplicated(t.meta$barcode),]
keep_cells <- rownames(t.meta2)
rownames(t.meta2) <- t.meta2$barcode

write.csv(t.meta2, file='t_cells_metadata_renamed.csv', quote=F, row.names=F)

macs$barcode <- colnames(macs)
macs$UMAP_1 <- macs@reductions$umap@cell.embeddings[,1]
macs$UMAP_2 <- macs@reductions$umap@cell.embeddings[,2]
macs.meta <- macs@meta.data
unique(str_extract(macs.meta$barcode,"-\\d+[_\\d+]*|_\\d+[_\\d+]*"))
bar.macs.temp <- str_replace(macs.meta$barcode,"-\\d+[_\\d+]*|_\\d+[_\\d+]*","")
bar.macs.temp2 <- paste(bar.macs.temp,macs.meta$orig.ident,sep="_")
macs.meta$barcode <- bar.macs.temp2
rownames(macs.meta) <- macs.meta$barcode
write.csv(macs.meta, file='macs_metadata_renamed.csv', quote=F, row.names=F)

#write expression counts matrix
t_counts_matrix <- GetAssayData(t_cells, assay = "RNA", slot='counts')
keep.t_counts_matrix <- t_counts_matrix[,colnames(t_counts_matrix) %in% keep_cells]
colnames(keep.t_counts_matrix) <- rownames(t.meta2)
writeMM(keep.t_counts_matrix, file="t_cell_counts.mtx")

mac_counts_matrix <- GetAssayData(macs, assay = "RNA", slot='counts')
colnames(mac_counts_matrix) <- rownames(macs.meta)
writeMM(mac_counts_matrix, file="macs_counts.mtx")

#write PCA matrix
t.embeddings <- t_cells@reductions$pca@cell.embeddings
keep.t.embeddings <- t.embeddings[rownames(t.embeddings) %in% keep_cells,]
rownames(keep.t.embeddings) <- rownames(t.meta2)
write.csv(keep.t.embeddings, file='t_cells_pca.csv', quote = F, row.names=F)
macs.embeddings <- macs@reductions$pca@cell.embeddings
rownames(macs.embeddings) <- rownames(macs.meta)
write.csv(macs.embeddings, file='macs_pca.csv', quote = F, row.names=F)

#write gene names
write.table(data.frame('gene'=rownames(keep.t_counts_matrix)),file='t_cell_gene_names.csv',
                       quote=F, row.names=F, col.names=F)
write.table(data.frame('gene'=rownames(mac_counts_matrix)),file='macs_gene_names.csv',
            quote=F, row.names=F, col.names=F)


######################Try another way to make filtering simpler################

#get Cell IDs into one file
test <- Cells(t_cells)
write.csv(rownames(t.meta2),file="T_cellID_obs.csv", row.names = FALSE)
write.csv(rownames(macs.meta),file="macs_cellID_obs.csv", row.names = FALSE)

#get Cell Embeddings
t.embeddings2 <- Embeddings(t_cells,reduction="umap")
keep.t.embeddings2 <- t.embeddings2[rownames(t.embeddings2) %in% keep_cells,]
rownames(keep.t.embeddings2) <- rownames(t.meta2)
write.csv(keep.t.embeddings2, file="t_cell_embeddings.csv")
macs.embeddings2 <- Embeddings(macs,reduction="umap")
rownames(macs.embeddings2) <- rownames(macs.meta)
write.csv(macs.embeddings2, file="macs_embeddings.csv")

##### separate out CD4+ and CD8+ T cells
my_color_palette.t <- pal_futurama("planetexpress")(12)

#use ridge plots
png(filename= "RNA_CD4expression_Tcells_RidgePlot.png", height=4, width=8, units='in', res=300)
RidgePlot(t_cells,features="CD4",assay = "SCT",cols = my_color_palette.t)
dev.off()

png(filename= "RNA_CD8Aexpression_Tcells_RidgePlot.png", height=4, width=8, units='in', res=300)
RidgePlot(t_cells,features="CD8A",assay = "SCT",cols = my_color_palette.t)
dev.off()

DefaultAssay(t_cells) <- "ADT"
t_cells <- NormalizeData(t_cells, normalization.method = "CLR", margin=2)

png(filename= "Protein_CD4expression_Tcells_RidgePlot.png", height=4, width=8, units='in', res=300)
RidgePlot(t_cells,features="CD4",assay = "ADT",cols = my_color_palette.t)
dev.off()

png(filename= "Protein_CD8aexpression_Tcells_RidgePlot.png", height=4, width=8, units='in', res=300)
RidgePlot(t_cells,features="CD8a",assay = "ADT",cols = my_color_palette.t)
dev.off()

#cluster 13 coexpressed neither CD8 and CD4 and wasnt identified by projectTILS
#so it will be removed from the analyses.
# CD4+ clusters: 3, 5, 6, 8, 9
# CD8+ clusters: 0, 1, 2, 4, 7, 10, 11,12
cd4 <- subset(
  t_cells,
  idents = c("3", "5","6","8","9"))
#saveRDS(cd4,"cd4.RDS")
cd8 <- subset(
  t_cells,
  idents = c("0","1","2","4","7","10","11", "12"))
#saveRDS(cd8,"cd8.RDS")


#save metadata
#annoyingly, barcodes are renamed differently.  Some are
#original barcode appended with -\\d+* and others are
#original barcode appended with _\\d+*
#examples: -1_4_4_1 and _8_8_1 and -1_17_1

cd4$barcode <- colnames(cd4)
cd4$UMAP_1 <- cd4@reductions$umap@cell.embeddings[,1]
cd4$UMAP_2 <- cd4@reductions$umap@cell.embeddings[,2]
cd4.meta <- cd4@meta.data
unique(str_extract(cd4.meta$barcode,"-\\d+[_\\d+]*|_\\d+[_\\d+]*"))
bar.cd4.temp <- str_replace(cd4.meta$barcode,"-\\d+[_\\d+]*|_\\d+[_\\d+]*","")
bar.cd4.temp2 <- paste(bar.cd4.temp,cd4.meta$orig.ident,sep="_")
cd4.meta$barcode <- bar.cd4.temp2
dup.cd4 <- cd4.meta[duplicated(cd4.meta$barcode),]
cd4.meta2 <- cd4.meta[!duplicated(cd4.meta$barcode),]
cd4keep_cells <- rownames(cd4.meta2)
rownames(cd4.meta2) <- cd4.meta2$barcode

write.csv(cd4.meta2, file='cd4_t_cells_metadata_renamed.csv', quote=F, row.names=F)

#write expression counts matrix
cd4_counts_matrix <- GetAssayData(cd4, assay = "RNA", slot='counts')
keep.cd4_counts_matrix <- cd4_counts_matrix[,colnames(cd4_counts_matrix) %in% cd4keep_cells]
colnames(keep.cd4_counts_matrix) <- rownames(cd4.meta2)
writeMM(keep.cd4_counts_matrix, file="cd4_t_cell_counts.mtx")


cd4.embeddings <- cd4@reductions$pca@cell.embeddings
keep.cd4.embeddings <- cd4.embeddings[rownames(cd4.embeddings) %in% cd4keep_cells,]
rownames(keep.cd4.embeddings) <- rownames(cd4.meta2)
write.csv(keep.cd4.embeddings, file='cd4_t_cells_pca.csv', quote = F, row.names=F)

#write gene names
write.table(data.frame('gene'=rownames(keep.cd4_counts_matrix)),file='cd4_t_cell_gene_names.csv',
            quote=F, row.names=F, col.names=F)


