
# Read in merged DGE matrix
combined <- data.frame(fread('data/DROP_combined_m400.tsv.gz',sep='\t'),row.names=1)

# Extract cell identity
cells <- strsplit(colnames(combined),split='_')
day <- factor(sapply(cells, FUN=function(x){x[2]}),levels = c('0','1','3','7','15'))
names(day) = colnames(combined)

batch <-  as.factor(sapply(cells, FUN=function(x){x[3]}))
names(batch) = colnames(combined)

batch_day <- factor(sapply(cells, FUN=function(x){paste0('Day ',x[2],' ',x[3])}), levels=c('Day 0 Rep1','Day 0 Rep2','Day 1 Rep1','Day 1 Rep2','Day 3 Rep1','Day 3 Rep2',
                                                                                          'Day 7 Rep1','Day 7 Rep2','Day 15 Rep2'))
names(batch_day) = colnames(combined)

# Seurat

ipsc <- CreateSeuratObject(combined)
rm(combined)
gc()

ipsc <- NormalizeData(ipsc)

ipsc <- FindVariableGenes(object = ipsc, mean.function = ExpMean, dispersion.function = LogVMR, 
                       x.low.cutoff = 0.15, x.high.cutoff = 3, y.cutoff = 1.5)
ipsc <- ScaleData(ipsc)

ipsc <- RunPCA(ipsc, pc.genes=ipsc@var.genes, do.print=F)

PCElbowPlot(ipsc)

ipsc <- FindClusters(ipsc, reduction.type = 'pca', dims.use = 1:7, resolution = 0.13, print.output = F, force.recalc = TRUE)

ipsc <- RunUMAP(ipsc, reduction.use = "pca", dims.use = 1:7)

DimPlot(ipsc, reduction.use = 'umap')

curr.idents <- as.factor(c(0,1,2,3,4,5))
new.idents <- c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Alt. Lineage 1','Cardiomyocyte 2','Alt. Lineage 2')
cluster.ids <- plyr::mapvalues(ipsc@ident, from = curr.idents, to = new.idents)
cluster.ids <- factor(cluster.ids, levels = c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2'))
ipsc@ident <- cluster.ids

save(cluster.ids,file='drop_cluster_ids')

markers <- FindAllMarkers(ipsc,
                         only.pos = TRUE, 
                         min.pc = 0.25,
                         test.use = 'negbinom')
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(ipsc, genes.use=top10$gene, slim.col.label= TRUE, remove.key = TRUE, group.label.rot = T)

save(markers,file = 'drop_markers')
save(ipsc, file = 'drop_seurat_object.Robj')

# Plotting

png('figures/drop_gene_markers.png',res=200,width=2500,height=2500)
DoHeatmap(ipsc, genes.use=top10$gene, slim.col.label= TRUE, remove.key = TRUE,group.label.rot = T)
dev.off()

png('figures/drop_UMAP_clusters.png',res=200, width = 1500, height=800)
ipsc@ident <- cluster.ids
DimPlot(ipsc, reduction.use = 'umap', pt.size = 0.5, cols.use = c('orange','darkgreen','cyan','blue','magenta','purple')) + 
  theme(legend.position = "right",axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18))
dev.off()

png('figures/drop_umap_batch.png',res=200, width = 1300, height=800)
ipsc@ident <- batch_day
DimPlot(ipsc, reduction.use = 'umap',pt.size = 0.5, cols.use=c('red','orange','darkgreen','green','blue','cyan','purple','pink','magenta')) + 
  theme(legend.position = "right",legend.text = element_text(size=14), axis.text = element_text(size=18), axis.title = element_text(size=14))
dev.off()

png('figures/drop_FeaturePlot_exon.png',res=200, width = 1400, height=1500)
GENES  <- c('DPPA4','EOMES','APLNR','FOXA2','TTR','MYH6','TNNT2','MYL7','AFP','SERPINA1','FLT1')
p <- FeaturePlot(ipsc, features.plot = GENES, cols.use = c('yellow','red'),pt.size = 0.1,no.axes=TRUE,do.return = T,reduction.use = 'umap')
p + theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18))
dev.off()


# Calculate proportion of celltypes per day

days <- unique(day)
frequency.mat <- matrix(0, ncol=6, nrow=length(days))
for(i in 1:length(days)){
  frequency.mat[i,] = table(cluster.ids[day==days[i]])
}

frequency.mat <- data.frame(frequency.mat/rowSums(frequency.mat))
colnames(frequency.mat) = factor(c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2'))
row.names(frequency.mat) = days
frequency.mat$day = days

map_col <- c('orange','darkgreen','cyan','blue','purple','magenta')
names(map_col) <- c('iPSC','Cardiac Progenitor', 'Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2')

frequency.mat <- melt(frequency.mat, id=c("day"))

png('figures/drop_day_prop.png',res=200,width=1200,height=800)
ggplot(frequency.mat, aes(x=day, y=value, fill=variable)) + geom_bar(stat="identity") + 
  scale_color_manual(values = map_col, aesthetics = "fill") + 
  theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18)) +
  xlab('Day') + ylab('Cell Type Proportion') + ggtitle('Drop-seq') + 
  theme(legend.title=element_blank())
dev.off()




