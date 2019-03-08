

## Correlation analysis

drop = data.frame(fread('data/DROP_combined_m400.tsv.gz',sep='\t'),row.names=1)
drop = filter(data = drop, 
                  min.genes = 400, 
                  max.genes = 2000)
load('drop_cluster_ids.txt')
drop.ids = cluster.ids

dronc = data.frame(fread('data/DRONC_combined_m300.tsv.gz',sep='\t'),row.names=1)
dronc = filter(data = dronc, 
              min.genes = 300, 
              max.genes = 2000)
load('dronc_cluster_ids.txt')
dronc.ids = cluster.ids

bulk.heart.tpm = log2(read.delim('metadata/GSE110471_RNA_tpms.txt',sep='\t',row.names=1)+1)
bulk.heart.mean = rowMeans(bulk.heart.tpm)

drop.pseudo.ipsc = rowSums(drop[,drop.ids=='iPSC'])
drop.pseudo.cm = rowSums(drop[,drop.ids=='Cardiomyocyte'|drop.ids=='Hepatocyte/Cardiomyocyte'])
drop_pseudo.ipsc = log2((drop.pseudo.ipsc/sum(drop.pseudo.ipsc))*1e4 + 1)
drop_pseudo.cm = log2((drop.pseudo.cm/sum(drop.pseudo.cm))*1e4 + 1)

dronc.pseudo.ipsc = rowSums(dronc[,dronc.ids=='iPSC'])
dronc.pseudo.cm = rowSums(dronc[,dronc.ids=='Cardiomyocyte'|dronc.ids=='Hepatocyte/Cardiomyocyte'])
dronc_pseudo.ipsc = log2((dronc.pseudo.ipsc/sum(dronc.pseudo.ipsc))*1e4 + 1)
dronc_pseudo.cm = log2((dronc.pseudo.cm/sum(dronc.pseudo.cm))*1e4 + 1)

same_genes = intersect(rownames(drop),rownames(dronc))
same_genes_w_bulk = intersect(same_genes, names(bulk.heart.mean))


sample.types = sapply(strsplit(colnames(bulk.heart.tpm),split='_'),function(x){x[3]})
sample.types[sample.types=='1'] = 'iPSC'
sample.types[sample.types=='2' | sample.types=='2A' | sample.types=='2B' | sample.types=='2C'] = 'iPSC-Cardiomyocytes'
sample.types[sample.types=='3'] = 'Primary Heart Tissue'

theme_set(theme_cowplot())

drop.ipsc.cor = cor(drop_pseudo.ipsc[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')
dronc.ipsc.cor = cor(dronc_pseudo.ipsc[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')

data.df = data.frame(cors=c(drop.cor,dronc.cor),type=c(rep('Drop',length(drop.cor)),rep('DroNc',length(dronc.cor))),samples=c(sample.types,sample.types))

drop.cm.cor = cor(drop_pseudo.cm[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')
dronc.cm.cor = cor(dronc_pseudo.cm[same_genes_w_bulk], bulk.heart.tpm[same_genes_w_bulk,],method = 'pearson')

n <- length(sample.types)
type <- factor(c(rep('Drop iPSC',n),rep('Drop CM',n),rep('DroNc iPSC',n),rep('DroNc CM',n)),levels = c('Drop iPSC','Drop CM','DroNc iPSC','DroNc CM'))
data.df = data.frame(cors=c(drop.ipsc.cor,drop.cm.cor,dronc.ipsc.cor,dronc.cm.cor),
                     type=type,
                     samples=rep(sample.types,4))

ggplot(data.df,aes(y=cors,x=type,fill=samples)) + 
  geom_boxplot() +  coord_flip() + 
  theme(axis.text.y = element_text(size=18), axis.text.x = element_text(size=18),text=element_text(size = 20)) +
  ylab('Pearson') + xlab('')

