# Read in merged DGE matrix

combined <- data.frame(fread('data/DRONC_combined_allcells.tsv.gz',sep='\t'),row.names=1)

# Basic filtering

combined <- filter(data = combined, 
                  min.genes = 300, 
                  max.genes = 2000)

cells <- strsplit(colnames(combined),split='_')

day <- factor(sapply(cells, FUN=function(x){x[2]}),levels = c('0','1','3','7','15'))
names(day) <- colnames(combined)

batch <- as.factor(sapply(cells, FUN=function(x){x[3]}))
names(batch) <- colnames(combined)

batch_day <- factor(sapply(cells, FUN=function(x){paste0('Day ',x[2],' ',x[3])}), 
                   levels=c('Day 0 Rep1','Day 0 Rep2','Day 1 Rep2','Day 3 Rep1','Day 3 Rep2','Day 7 Rep1','Day 7 Rep2','Day 15 Rep2'))
names(batch_day) <- colnames(combined)

# Seurat

ipsc=CreateSeuratObject(combined)
rm(combined)
gc()

ipsc <- NormalizeData(ipsc)

ipsc <- FindVariableGenes(object = ipsc, mean.function = ExpMean, dispersion.function = LogVMR, 
                       x.low.cutoff = 0.15, x.high.cutoff = 3, y.cutoff = 1.5)
ipsc <- ScaleData(ipsc)

ipsc <- RunPCA(ipsc, pc.genes=ipsc@var.genes, do.print=F,pcs.compute = 20)

PCElbowPlot(ipsc)

ipsc <- FindClusters(ipsc, reduction.type = 'pca', dims.use = 1:7, resolution = 0.13, print.output = F, force.recalc = TRUE)

ipsc <- RunUMAP(ipsc, reduction.use = "pca", dims.use = 1:7, min_dist = 0.5)

curr.idents <- as.factor(c(0,1,2,3,4))
new.idents <- c('iPSC','Cardiac Progenitor','Cardiomyocyte 2','Hepatocyte-like','Cardiomyocyte 1')
cluster.ids <- plyr::mapvalues(ipsc@ident, from = curr.idents, to = new.idents)
cluster.ids <- factor(cluster.ids, levels= c('iPSC','Cardiac Progenitor', 'Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 2'))
ipsc@ident <- cluster.ids

save(cluster.ids, file='dronc_cluster_ids.Robj')

markers <- FindAllMarkers(ipsc,
                         only.pos = TRUE, 
                         min.pc = 0.25,
                         test.use = 'negbinom')

top <- markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

save(markers, file='dronc_markers.Robj')
save(ipsc, file='dronc_seurat_object.Robj')

# Plotting

png('figures/dronc_gene_markers.png',res=200,width=2500,height=2500)
DoHeatmap(ipsc, genes.use=top10$gene, slim.col.label= TRUE, remove.key = TRUE,group.label.rot = T)
dev.off()

png('figures/dronc_umap_cluster.png',res=200, width = 1500, height=800)
ipsc@ident <- cluster.ids
DimPlot(ipsc, reduction.use = 'umap', pt.size = 0.5, cols.use = c('orange','darkgreen','cyan','blue','magenta')) + 
  theme(legend.position = "right",axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18))
dev.off()

png('figures/dronc_umap_batch.png',res=200, width = 1300, height=800)
ipsc@ident <- batch_day
DimPlot(ipsc, reduction.use = 'umap',pt.size = 0.5, cols.use=c('red','orange','darkgreen','blue','cyan','purple','pink','magenta')) + 
  theme(legend.position = "right",legend.text = element_text(size=14), axis.text = element_text(size=18), axis.title = element_text(size=14))
dev.off()

png('figures/dronc_FeaturePlots.png',res=200, width = 1500, height=1500)
GENES <- c('DPPA4','EOMES','APLNR','FOXA2','TTR','MYH6','TNNT2','MYL7','AFP','SERPINA1','FLT1')
p <- FeaturePlot(ipsc, features.plot = GENES, cols.use = c('yellow','red'),pt.size = 0.1,no.axes=TRUE,do.return = T,reduction.use = 'umap')
p + theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18))
dev.off()


# Calculate celltype proportion per day

days <- unique(day)
frequency.mat <- matrix(0, ncol=5, nrow=length(days))
for(i in 1:length(days)){
  frequency.mat[i,] <- table(cluster.ids[day==days[i]])
}
frequency.mat <- data.frame(frequency.mat/rowSums(frequency.mat))
colnames(frequency.mat) <- factor(c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Cardiomyocyte 2','Hepatocyte-like'))
row.names(frequency.mat) <- days
frequency.mat$day <- days

map_col <- c('orange','darkgreen','cyan','blue','magenta')
names(map_col) <- c('iPSC','Cardiac Progenitor', 'Cardiomyocyte 1','Cardiomyocyte 2','Hepatocyte-like')

frequency.mat <- melt(frequency.mat, id=c("day"))

# Plot
png('figures/dronc_day_prop.png',res=200,width=1200,height=800)
ggplot(frequency.mat, aes(x=day, y=value, fill=variable)) + geom_bar(stat="identity") + 
  scale_color_manual(values = map_col, aesthetics = "fill") + 
  theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18)) +
  xlab('Day') + ylab('Cell Type Proportion') + ggtitle('DroNc-seq') + 
  theme(legend.title=element_blank())
dev.off()

