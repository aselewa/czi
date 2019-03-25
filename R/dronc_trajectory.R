# Monocle

load('dronc_markers.Robj')
load('dronc_cluster_ids.Robj')

combined <- data.frame(fread('data/DRONC_combined_m300.tsv.gz',sep='\t'),row.names=1)

day <- getCellTag(colnames(combined), 'day')

sub.cells <- c()
for(d in unique(day)){
  curr_cells <- names(day[which(day==d)])
  top <- names(sort(colSums(combined[,curr_cells]),decreasing = TRUE))
  sub.cells <- c(sub.cells, top[1:700])
}

log.cpm <- computeLogCPM(combined)
rm(combined)
gc()

sub_combined <- log.cpm[,sub.cells]
sub_ids <- cluster.ids[sub.cells]
sub_days <- day[sub.cells]

top <- markers %>% group_by(cluster) %>% top_n(60, avg_logFC)
var_genes <- unique(top$gene)

samples.df <- data.frame(row.names = colnames(sub_combined), Day = sub_days, CellType = sub_ids)
features.df <- data.frame(row.names = row.names(sub_combined), gene_short_name = row.names(sub_combined))
pd <- new("AnnotatedDataFrame",data=samples.df)
fd <- new("AnnotatedDataFrame",data=features.df)

ipsc_m <- newCellDataSet(sub_combined,phenoData = pd, featureData = fd,expressionFamily = VGAM::inv.gaussianff())
ipsc_m <- setOrderingFilter(ipsc_m, var_genes)
ipsc_m <- reduceDimension(ipsc_m, max_components = 2, method = 'DDRTree',norm_method = 'none')
ipsc_m <- orderCells(ipsc_m, reverse = T)

# Plot
lookup <- data.frame(row.names=unique(day), cols= c('orange','darkgreen','blue','purple','magenta'))

png('figures/dronc_monocle_days.png',res = 200, width = 1000, height = 700)
p <- plot_cell_trajectory(ipsc_m,color_by = "Day",cell_size = 0.5)
map_col <- as.character(lookup$cols)
names(map_col) <- row.names(lookup)
p + scale_color_manual(values = map_col, name = "Day",labels = c("0","1","3","7","15")) + 
  theme(legend.position = "right",legend.text = element_text(size=12)) + 
  guides(color = guide_legend(override.aes = list(size=2)))
dev.off()

png('figures/dronc_monocle_celltype.png',res = 200, width = 1000, height = 800)
p <- plot_cell_trajectory(ipsc_m,color_by = "CellType",cell_size = 0.3)
map_col <- c('orange','darkgreen','cyan','blue','magenta')
names(map_col) <- c('iPSC','Cardiac Progenitor', 'Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 2')
p + scale_color_manual(values = map_col, name = "Cell Type") +
  theme(legend.position="none",text = element_text(size = 20))
dev.off()

png('figures/dronc_monocle_time.png',res = 200, width = 800, height = 700)
plot_cell_trajectory(ipsc_m,color_by = "Pseudotime",cell_size = 0.5) + theme(legend.text = element_text(size=6)) + 
  theme(text = element_text(size = 14),legend.position="none")
dev.off()


