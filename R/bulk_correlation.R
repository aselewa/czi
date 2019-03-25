
# Load data
drop = data.frame(fread('data/DROP_combined_m400.tsv.gz',sep='\t'),row.names=1)
drop = drop[rowSums(drop>0)>1000,]
load('~/drop_cluster_ids.Robj')
drop.ids = cluster.ids

dronc = data.frame(fread('data/DRONC_combined_m300.tsv.gz',sep='\t'),row.names=1)
dronc = dronc[rowSums(dronc>0)>1000,]
load('~/dronc_cluster_ids.Robj')
dronc.ids = cluster.ids

# Load bulk RNA-seq
bulk.heart.tpm = log(read.delim('GSE110471_RNA_tpms.txt',sep='\t',row.names=1)+1)
bulk.heart.mean = rowMeans(bulk.heart.tpm)

# Compute psuedo-bulk for iPSCs and CM clusters in Drop and DroNc-seq
drop.pseudo.ipsc = rowSums(drop[,names(drop.ids)[drop.ids=='iPSC']])
drop.pseudo.cm = rowSums(drop[,names(drop.ids)[drop.ids=='Cardiomyocyte 1' | drop.ids=='Cardiomyocyte 2']])
drop_pseudo.ipsc = log((drop.pseudo.ipsc/sum(drop.pseudo.ipsc))*1e4 + 1)
drop_pseudo.cm = log((drop.pseudo.cm/sum(drop.pseudo.cm))*1e4 + 1)

dronc.pseudo.ipsc = rowSums(dronc[,names(dronc.ids)[dronc.ids=='iPSC']])
dronc.pseudo.cm = rowSums(dronc[,names(dronc.ids)[(dronc.ids=='Cardiomyocyte 1' | dronc.ids=='Cardiomyocyte 2')]])
dronc_pseudo.ipsc = log((dronc.pseudo.ipsc/sum(dronc.pseudo.ipsc))*1e4 + 1)
dronc_pseudo.cm = log((dronc.pseudo.cm/sum(dronc.pseudo.cm))*1e4 + 1)

same_genes = intersect(names(drop_pseudo.cm),names(dronc_pseudo.cm))
same_genes_w_bulk = intersect(same_genes, names(bulk.heart.mean))

# Define bulk RNA-seq sample types
sample.types = sapply(strsplit(colnames(bulk.heart.tpm),split='_'),function(x){x[3]})
sample.types[sample.types=='1'] = 'iPSC'
sample.types[sample.types=='2' | sample.types=='2A' | sample.types=='2B' | sample.types=='2C'] = 'iPSC-Cardiomyocytes'
sample.types[sample.types=='3'] = 'Primary Heart Tissue'

# Plotting
theme_set(theme_cowplot())

drop.ipsc.cor = cor(drop_pseudo.ipsc[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')
dronc.ipsc.cor = cor(dronc_pseudo.ipsc[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')

drop.cm.cor = cor(drop_pseudo.cm[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')
dronc.cm.cor = cor(dronc_pseudo.cm[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')

n <- length(sample.types)
type <- factor(c(rep('Drop iPSC',n),rep('Drop CM',n),rep('DroNc iPSC',n),rep('DroNc CM',n)),levels = c('Drop iPSC','Drop CM','DroNc iPSC','DroNc CM'))
data.df = data.frame(cors=c(drop.ipsc.cor,drop.cm.cor,dronc.ipsc.cor,dronc.cm.cor),
                     type=type,
                     samples=rep(sample.types,4))

ggplot(data.df,aes(y=cors,x=type,fill=samples)) + geom_boxplot() + theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=20),text=element_text(size = 24, angle=90), legend.position = 'none') +
  ylab('Pearson') + xlab('') 
