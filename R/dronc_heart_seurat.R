
# Compare Human Heart against DroNc-seq
load('/project2/onibasu/data/public_html/Bing_analysis/HeartTissue/humanHeart/intron_exon_sum/top_30/top_30_sum.Robj')
load('dronc_cluster_ids.txt')

object = FindClusters(object, reduction.type = 'pca', dims.use = 1:5, resolution = 0.1, print.output = F, force.recalc = TRUE)

object = RunTSNE(object, reduction.use = 'pca', dims.use = 1:5, check_duplicates = FALSE, perplexity=100)

markers = FindAllMarkers(object, only.pos = TRUE, min.pc = 0.25, logfc.threshold = 0.25,test.use = 'negbinom')
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

GENES = c('TNNT2','MYH6','MYBPC3','MYL7','COL6A3','COL5A2','FLT1','VWF','ENPEP','CPE')
VlnPlot(object, features.plot = GENES,point.size.use = 0.3, nCol = 5,size.x.use = 0,cols.use = c('blue','red','orange','black'))

og = c(0,1,2,3)
new.ids = c('Cardiomyocytes','Connective Tissue','Endothelial','Other')
celltypes = plyr::mapvalues(object@ident,from = og, to = new.ids)
object@ident = celltypes

png('figures/dronc_heart_tsne.png',res=200, width = 1600, height=800)
TSNEPlot(object,pt.size = 0.5, colors.use = c('blue','red','orange','black')) + 
  theme(legend.position = "right",text=element_text(size=24))
dev.off()

ipsc_cts = combined
hrt_cts = object@raw.data

same_genes = intersect(rownames(hrt_cts),rownames(ipsc_cts))
ngenes = length(same_genes)

ipsc_cts = ipsc_cts[same_genes,]
hrt_cts = hrt_cts[same_genes,]

ipsc_exp = t(log((t(ipsc_cts)/colSums(ipsc_cts))*1e4 + 1))
hrt_exp = t(log((t(hrt_cts)/colSums(hrt_cts))*1e4 + 1))

ipsc_cm = ipsc_cts[,cluster.ids=='Cardiomyocyte']
ipsc_cm = rowMeans(ipsc_cm)
ipsc_cm = log((ipsc_cm/sum(ipsc_cm))*10^4 + 1)


x= cor(ipsc_cm, hrt_exp[,sidecols=='blue'],method='pearson')
y = cor(ipsc_cm, hrt_exp[,sidecols=='red'],method='pearson')
z = cor(ipsc_cm, hrt_exp[,sidecols=='orange'],method='pearson')
w = cor(ipsc_cm, hrt_exp[,sidecols=='black'],method='pearson')
type=c(rep('Cardiomyocytes',length(x)),rep('Connective Tissue',length(y)),rep('Endothelial',length(z)),rep('Other',length(w)))
cors_CM = c(x,y,z,w)

ipsc_cm = ipsc_cts[,cluster.ids=='Hepatocyte/Cardiomyocyte']
ipsc_cm = rowMeans(ipsc_cm)
ipsc_cm = log((ipsc_cm/sum(ipsc_cm))*10^4 + 1)

x= cor(ipsc_cm, hrt_exp[,sidecols=='blue'],method='pearson')
y = cor(ipsc_cm, hrt_exp[,sidecols=='red'],method='pearson')
z = cor(ipsc_cm, hrt_exp[,sidecols=='orange'],method='pearson')
w = cor(ipsc_cm, hrt_exp[,sidecols=='black'],method='pearson')
cors_CM2= c(x,y,z,w)

data.df = data.frame(cors = c(cors_CM,cors_CM2), type=c(type,type),from=c(rep('Cardiomyocyte',length(cors_CM)),rep('Hepatocyte/Cardiomyocyte',length(cors_CM2))))

cols = c('purple','red','orange','black')
names(cols) = c('Cardiomyocytes','Connective Tissue','Endothelial','Other')

ggplot(data.df, aes(y=cors, x=from, fill=type)) + geom_boxplot() +  theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18),text=element_text(size = 20),legend.position = 'none') +
  scale_color_manual(values = cols)

ids = c('Cardiomyocyte','Hepatocyte/Cardiomyocyte')
k=150
top.cells = c()
for(i in 1:length(ids)){
  
  curr.ids = names(cluster.ids[which(cluster.ids==ids[i])])
  topk = curr.ids[order(colSums(ipsc_cts[,curr.ids]),decreasing = T)][1:k]
  top.cells = c(top.cells, topk)
  
}

scores = cor(ipsc_exp[,top.cells],hrt_exp, method='pearson')

og = c('Cardiomyocytes','Connective Tissue','Endothelial','Other')
new.ids = c('blue','red','orange','black')
sidecols = as.character(plyr::mapvalues(object@ident,from = og, to = new.ids))

newcols = c('grey34','grey88')
rowcols = as.character(plyr::mapvalues(cluster.ids[top.cells],from=ids, to=newcols))

base = brewer.pal(8, "RdBu")
pallete = colorRampPalette(base)(50)

heatmap.2(scores, 
          trace='none',
          dendrogram = 'column',
          col = pallete,
          ColSideColors = sidecols,
          RowSideColors = rowcols,
          labCol = NA,
          labRow = NA
)
