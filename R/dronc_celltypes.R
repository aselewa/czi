
# Read in merged DGE matrix
combined <- data.frame(fread('data/DRONC_combined_m300_exons.tsv.gz',sep='\t'),row.names=1)

# top expressed genes
tot <- colSums(combined)
comb.top <- combined[order(rowSums(combined),decreasing = T)[1:15],]
comb.percent <- t(t(comb.top)/tot)*100
comb.melt <- melt(comb.percent)
comb.melt <- comb.melt[comb.melt$value < 10,]

comb.melt$Var1 <- factor(comb.melt$Var1, levels = rev(row.names(comb.top)))
ggplot(comb.melt, aes(x=Var1, y=value,fill=Var1)) + geom_boxplot() + ylab('% of Total Counts') + xlab('Genes') + 
  theme(legend.position = 'none',text = element_text(size = 18)) + coord_flip()

# Extract metadata from column names
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

ipsc <- FindClusters(ipsc, reduction.type = 'pca', dims.use = 1:7, resolution = 0.15, print.output = F, force.recalc = TRUE)

ipsc <- RunUMAP(ipsc, reduction.use = "pca", dims.use = 1:7, min_dist = 0.5)

curr.idents <- as.factor(c(0,1,2,3,4))
new.idents <- c('iPSC','Cardiac Progenitor','Cardiomyocyte 2','Alt. Lineage 1','Cardiomyocyte 1')
cluster.ids <- plyr::mapvalues(ipsc@ident, from = curr.idents, to = new.idents)
cluster.ids <- factor(cluster.ids, levels= c('iPSC','Cardiac Progenitor', 'Cardiomyocyte 1','Cardiomyocyte 2','Alt. Lineage 1'))
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
p <- FeaturePlot(ipsc, features.plot = genes[1:10], cols.use = c('yellow','red'),pt.size = 0.1,no.axes=TRUE,do.return = T,reduction.use = 'umap')
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

# DroNc on human heart

load('dronc_heart_seurat.Robj')

object <- FindClusters(object, reduction.type = 'pca', dims.use = 1:5, resolution = 0.1, print.output = F, force.recalc = TRUE)

object <- RunUMAP(object, reduction.use = 'pca', dims.use = 1:5)

GENES <- c('TNNT2','MYH6','MYBPC3','MYL7','COL6A3','COL5A2','FLT1','VWF','ENPEP','CPE')
VlnPlot(object, features.plot = GENES, point.size.use = 0.2, nCol = 5, size.x.use = 0, cols.use = c('deeppink','darkgreen','grey48','black'))

og = c(0,1,2,3)
new.ids = c('Cardiomyocytes','Connective Tissue','Endothelial','Other')
celltypes = plyr::mapvalues(object@ident,from = og, to = new.ids)
object@ident = celltypes

png('../figures/dronc_heart_umap.png',res=200, width = 1600, height=800)
DimPlot(object,pt.size = 0.5, reduction.use = 'umap', cols.use = c('deeppink','darkgreen','grey48','black')) + 
  theme(legend.position = "right",text=element_text(size=24))
dev.off()

#### Correlation with bulk

bulk.heart.tpm = log(read.delim('GSE110471_RNA_tpms.txt',sep='\t',row.names=1)+1)

sample.types = sapply(strsplit(colnames(bulk.heart.tpm),split='_'),function(x){x[3]})
sample.types[sample.types=='1'] = 'iPSC'
sample.types[sample.types=='2' | sample.types=='2A' | sample.types=='2B' | sample.types=='2C'] = 'iPSC-Cardiomyocytes'
sample.types[sample.types=='3'] = 'Primary Heart Tissue'

hrt.counts = object@raw.data
same_genes = intersect(row.names(bulk.heart.tpm), row.names(hrt.counts))
hrt.counts = hrt.counts[same_genes,]
bulk.heart.tpm = bulk.heart.tpm[same_genes,]

#cm.counts = rowSums(hrt.counts[,object@ident=='Cardiomyocytes'])
#cm.pseudo = log((cm.counts/sum(cm.counts))*1e4 + 1)
#cm.pseudo = cm.pseudo[cm.pseudo>1.25]

tot.counts = rowSums(hrt.counts)
tot.pseudo = log((tot.counts/sum(tot.counts))*1e4 + 1)
tot.pseudo = tot.pseudo[tot.pseudo>1]

#cm.scores = cor(bulk.heart.tpm[names(cm.pseudo),], cm.pseudo, method='pearson')
tot.scores = cor(bulk.heart.tpm[names(tot.pseudo),], tot.pseudo, method='pearson')

#n <- length(cm.scores)
#type <- factor(c(rep('Heart Pseudobulk',n),rep('CM Pseudobulk',n)),levels = c('Heart Pseudobulk','CM Pseudobulk'))

scores.df = data.frame(scores = tot.scores, type=sample.types)

ggplot(scores.df, aes(y=scores, x=type, fill=type)) + geom_boxplot() +
  ylab('Pearson') + xlab('Heart Pseudo-bulk') +
  theme_bw() +
  theme(text=element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=18), axis.text.x = element_blank(), axis.ticks.x = element_blank())


#### Correlation with iPSC-CMs 

ipsc_cts = ipsc@raw.data
hrt_cts = hrt@raw.data
cluster.ids = ipsc@ident

hrt_cts = hrt_cts[rowSums(hrt_cts>0)>100,]
same_genes = intersect(row.names(ipsc_cts), row.names(hrt_cts))

ipsc_exp = as.matrix(ipsc@data)[same_genes,]
hrt_exp = as.matrix(hrt@data)[same_genes,]

ids = c('Cardiomyocyte 1','Cardiomyocyte 2', 'iPSC')
k = 100
top.cells = c()
for(i in 1:length(ids)){
  curr.ids = names(cluster.ids[which(cluster.ids==ids[i])])
  if(ids[i] == 'iPSC'){
    topk = curr.ids[order(colSums(ipsc_cts[,curr.ids]),decreasing = T)][1:50]
    top.cells = c(top.cells, topk)
  }
  else{
    topk = curr.ids[order(colSums(ipsc_cts[,curr.ids]),decreasing = T)][1:k]
    top.cells = c(top.cells, topk)
  }
}

scores = cor(ipsc_exp[,top.cells],hrt_exp, method='pearson')

og = c('Cardiomyocytes','Connective Tissue','Endothelial','Other')
clust.cols =c('deeppink','darkgreen','grey48','black')
names(clust.cols) = og
colAnnots = data.frame(Type=hrt@ident, row.names = names(hrt@ident))

ipsc.cols =  c('cyan','blue','orange')
names(ipsc.cols) = ids
rowAnnots = data.frame(iPSC_Derived=cluster.ids[top.cells], row.names = top.cells)

library(pheatmap)
library(RColorBrewer)
base = brewer.pal(8, "RdBu")
pallete = colorRampPalette(base)(50)

pheatmap(scores, 
         show_rownames = F, 
         show_colnames = F,
         annotation_col = colAnnots, 
         annotation_row = rowAnnots,
         annotation_colors = list(Type=clust.cols, iPSC_Derived=ipsc.cols))





