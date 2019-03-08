# Read in merged DGE matrix
combined = data.frame(fread('data/DROP_combined_m400.tsv.gz',sep='\t'),row.names=1)

# Extract cell identity
cells = strsplit(colnames(combined),split='_')
day = factor(sapply(cells, FUN=function(x){x[2]}),levels = c('0','1','3','7','15'))
names(day) = colnames(combined)

batch =  as.factor(sapply(cells, FUN=function(x){x[3]}))
names(batch) = colnames(combined)

batch_day = factor(sapply(cells, FUN=function(x){paste0('Day ',x[2],' ',x[3])}), levels=c('Day 0 Rep1','Day 0 Rep2','Day 1 Rep1','Day 1 Rep2','Day 3 Rep1','Day 3 Rep2',
                                                                                          'Day 7 Rep1','Day 7 Rep2','Day 15 Rep2'))
names(batch_day) = colnames(combined)

# Seurat
ipsc=CreateSeuratObject(combined)
rm(combined)
gc()

ipsc=NormalizeData(ipsc)
ipsc=FindVariableGenes(object = ipsc, mean.function = ExpMean, dispersion.function = LogVMR, 
                       x.low.cutoff = 0.15, x.high.cutoff = 3, y.cutoff = 1.5)
ipsc = ScaleData(ipsc)

ipsc=RunPCA(ipsc, pc.genes=ipsc@var.genes, do.print=F)

PCElbowPlot(ipsc)

ipsc = FindClusters(ipsc, reduction.type = 'pca', dims.use = 1:7, resolution = 0.13, print.output = F, force.recalc = TRUE)

ipsc = RunUMAP(ipsc, reduction.use = "pca", dims.use = 1:7)

DimPlot(ipsc, reduction.use = 'umap')

curr.idents = as.factor(c(0,1,2,3,4,5))
new.idents = c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Alt. Lineage 1','Cardiomyocyte 2','Alt. Lineage 2')
cluster.ids = plyr::mapvalues(ipsc@ident, from = curr.idents, to = new.idents)
cluster.ids = factor(cluster.ids, levels = c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2'))
ipsc@ident = cluster.ids

save(cluster.ids,file='drop_cluster_ids')

markers = FindAllMarkers(ipsc,
                         only.pos = TRUE, 
                         min.pc = 0.25,
                         test.use = 'negbinom')
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(ipsc, genes.use=top10$gene, slim.col.label= TRUE, remove.key = TRUE, group.label.rot = T)

save(markers,file='drop_markers')
save(ipsc, file='drop_seurat_object.Robj')

# Plotting

png('figures/drop_gene_markers.png',res=200,width=2500,height=2500)
DoHeatmap(ipsc, genes.use=top10$gene, slim.col.label= TRUE, remove.key = TRUE,group.label.rot = T)
dev.off()

png('figures/drop_UMAP_clusters.png',res=200, width = 1500, height=800)
ipsc@ident = cluster.ids
DimPlot(ipsc, reduction.use = 'umap', pt.size = 0.5, cols.use = c('orange','darkgreen','cyan','blue','magenta','purple')) + 
  theme(legend.position = "right",axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18))
dev.off()

png('figures/drop_umap_batch.png',res=200, width = 1300, height=800)
ipsc@ident = batch_day
DimPlot(ipsc, reduction.use = 'umap',pt.size = 0.5, cols.use=c('red','orange','darkgreen','green','blue','cyan','purple','pink','magenta')) + 
  theme(legend.position = "right",legend.text = element_text(size=14), axis.text = element_text(size=18), axis.title = element_text(size=14))
dev.off()

png('figures/drop_FeaturePlot_exon.png',res=200, width = 1400, height=1500)
GENES = c('DPPA4','EOMES','APLNR','FOXA2','TTR','MYH6','TNNT2','MYL7','AFP','SERPINA1','FLT1')
p = FeaturePlot(ipsc, features.plot = GENES, cols.use = c('yellow','red'),pt.size = 0.1,no.axes=TRUE,do.return = T,reduction.use = 'umap')
p + theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18))
dev.off()


# Calculate proportion of celltypes per day

days = unique(day)
frequency.mat = matrix(0, ncol=6, nrow=length(days))
for(i in 1:length(days)){
  frequency.mat[i,] = table(cluster.ids[day==days[i]])
}
frequency.mat = data.frame(frequency.mat/rowSums(frequency.mat))
colnames(frequency.mat) = factor(c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2'))
row.names(frequency.mat) = days
frequency.mat$day = days

map_col = c('orange','darkgreen','cyan','blue','purple','magenta')
names(map_col) = c('iPSC','Cardiac Progenitor', 'Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2')

frequency.mat = melt(frequency.mat, id=c("day"))

png('figures/drop_day_prop.png',res=200,width=1200,height=800)
ggplot(frequency.mat, aes(x=day, y=value, fill=variable)) + geom_bar(stat="identity") + 
  scale_color_manual(values = map_col, aesthetics = "fill") + 
  theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18), text=element_text(size=18)) +
  xlab('Day') + ylab('Cell Type Proportion') + ggtitle('Drop-seq') + 
  theme(legend.title=element_blank())
dev.off()

# Monocle

load('drop_markers')
load('drop_cluster_ids')

combined <- data.frame(fread('data/DROP_combined_m400.tsv.gz',sep='\t'),row.names=1)

day <- getCellTag(colnames(combined), 'day')
lookup = data.frame(row.names=unique(day), cols= c('orange','darkgreen','blue','purple','magenta'))

sub.cells = c()
for(d in unique(day)){
  curr_cells = names(day[which(day==d)])
  top = names(sort(colSums(combined[,curr_cells]),decreasing = TRUE))
  sub.cells = c(sub.cells, top[1:700])
}

log.cpm <- computeLogCPM(combined)
rm(combined)
gc()

sub_combined = log.cpm[,sub.cells]
sub_ids = cluster.ids[sub.cells]
sub_days = day[sub.cells]

top = markers %>% group_by(cluster) %>% top_n(60, avg_logFC)
var_genes = unique(top$gene)

samples.df = data.frame(row.names = colnames(sub_combined), Day = sub_days, CellType = sub_ids)
features.df = data.frame(row.names = row.names(sub_combined), gene_short_name = row.names(sub_combined))
pd = new("AnnotatedDataFrame",data=samples.df)
fd = new("AnnotatedDataFrame",data=features.df)

ipsc_m = newCellDataSet(sub_combined,phenoData = pd, featureData = fd,expressionFamily = VGAM::inv.gaussianff())
ipsc_m = setOrderingFilter(ipsc_m, var_genes)
ipsc_m = reduceDimension(ipsc_m, max_components = 2, method = 'DDRTree',norm_method = 'none')
ipsc_m = orderCells(ipsc_m,reverse = T)

png('figures/drop_monocle_days.png',res = 200, width = 1000, height = 700)
p = plot_cell_trajectory(ipsc_m,color_by = "Day",cell_size = 0.3)
map_col = as.character(lookup$cols)
names(map_col) = row.names(lookup)
p + scale_color_manual(values = map_col, name = "Day",labels = c("0","1","3","7","15")) + 
  theme(legend.position = "right",legend.text = element_text(size=12)) + 
  guides(color = guide_legend(override.aes = list(size=2)),text = element_text(size = 20))
dev.off()

png('figures/drop_monocle_celltype.png',res = 200, width = 1000, height = 800)
p = plot_cell_trajectory(ipsc_m,color_by = "CellType",cell_size = 0.5)
map_col = c('orange','darkgreen','cyan','blue','purple','magenta')
names(map_col) = c('iPSC','Cardiac Progenitor','Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1','Alt. Lineage 2')
p + scale_color_manual(values = map_col, name = "Cell Type") +
  theme(legend.position="none",text = element_text(size = 20))
dev.off()

png('figures/drop_monocle_time.png',res = 200, width = 800, height = 700)
plot_cell_trajectory(ipsc_m,color_by = "Pseudotime",cell_size = 0.5) + theme(legend.text = element_text(size=6)) + 
  theme(legend.position="none",text = element_text(size = 14))
dev.off()



