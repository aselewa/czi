
# Load data

drop <- data.frame(fread('DROP_combined_allcells.tsv.gz',sep='\t'),row.names=1)
drop <- filter(data = drop, 
               min.genes = 400, 
               max.genes = 2000)


species_assign = data.frame(fread('metadata/species_assign_matrix.txt',sep='\t'))

chimp.cells = species_assign$CB[species_assign$hg_specificity_score<0.5]

cells = strsplit(colnames(drop),split='_')
barcodes = sapply(cells, function(x){x[1]})
drop = drop[,!(barcodes %in% chimp.cells)]


dronc <- data.frame(fread('DRONC_combined_allcells.tsv.gz',sep='\t'),row.names=1)
dronc <- filter(data = dronc, 
                min.genes = 300, 
                max.genes = 2000)


colnames(drop) <- paste0(colnames(drop),'_DROP')
colnames(dronc) <- paste0(colnames(dronc),'_DRONC')

same_genes <- intersect(
  row.names(dronc),
  row.names(drop)
)

comb <- as(as.matrix(cbind(
  drop[same_genes,],
  dronc[same_genes,])
),"dgCMatrix")

rm(dronc,drop)
gc()

cellnames <- colnames(comb)

cells <- strsplit(cellnames,split='_')

day <- factor(sapply(cells, FUN=function(x){x[2]}),levels=c('0','1','3','7','15'))
batch <- factor(sapply(cells, FUN=function(x){if(x[3]=='Rep1'){return(0)} else{return(1)}}),levels=c('0','1'))
tech <- sapply(cells, FUN=function(x){if(x[4]=='DROP'){return(0)} else{return(1)}})

## DROP-seq

## Differential expression between days

scomb <- CreateSeuratObject(comb)

cellnames <- colnames(comb)

cells <- strsplit(cellnames,split='_')
ident <- factor(sapply(cells, FUN=function(x){paste0(x[2],'_',x[3],'_',x[4])}))
names(ident) <- cellnames
scomb@ident <- ident

res <- list()

days <- c('0','1','3','7','15')
batch_tech <- c('Rep1_DROP','Rep2_DROP','Rep1_DRONC','Rep2_DRONC')

for(k in batch_tech){
  for(i in 1:(length(days)-1)){
    for(j in (i+1):length(days)){
      ident1 <- paste0(days[i],'_',k)
      ident2 <- paste0(days[j],'_',k)
      tag <- paste0(ident1,',',ident2)
      if(ident1 != '1_Rep1_DRONC' & ident2 != '1_Rep1_DRONC' & ident2 != '15_Rep1_DROP' & ident2 != '15_Rep1_DRONC'){
        curr <- FindMarkers(scomb, ident.1=ident1, ident.2 =  ident2, min.pct = 0.25, test.use='MAST', check_sanity = F)
        curr$gene <- row.names(curr)
        curr$id <- tag
        res[[tag]] <- curr
      }
    }
  }
}

markers_days <- Reduce(rbind, res)

save(markers_days, file='metadata/markers_days')

## Differential between Reps 

res <- list()
test <- 'negbinom'

res[[1]] <- FindMarkers(scomb, ident.1 = '0_Rep1_DROP',ident.2 = '0_Rep2_DROP', min.pct = 0.25, test.use = test)
res[[2]] <- FindMarkers(scomb, ident.1 = '0_Rep1_DROP',ident.2 = '0_Rep2_DROP', min.pct = 0.25, test.use = test)

res[[3]] <- FindMarkers(scomb, ident.1 = '1_Rep1_DROP',ident.2 = '1_Rep2_DROP', min.pct = 0.25, test.use = test)

res[[4]] <- FindMarkers(scomb, ident.1 = '3_Rep1_DROP',ident.2 = '3_Rep2_DROP', min.pct = 0.25, test.use = test)
res[[5]] <- FindMarkers(scomb, ident.1 = '3_Rep1_DRONC',ident.2 = '3_Rep2_DRONC', min.pct = 0.25, test.use = test)

res[[6]] <- FindMarkers(scomb, ident.1 = '7_Rep1_DROP',ident.2 = '7_Rep2_DROP', min.pct = 0.25, test.use = test)
res[[7]] <- FindMarkers(scomb, ident.1 = '7_Rep1_DRONC',ident.2 = '7_Rep2_DRONC', min.pct = 0.25, test.use = test)

for(i in 1:length(res)){
  res[[i]]$gene <- row.names(res[[i]])
  res[[i]]$id <- i
}

markers_batch <- Reduce(rbind, res)


save(markers_batch, file='metadata/markers_batch')

## Between DroNc and Drop

res <- list()

res[[1]] <- FindMarkers(scomb, ident.1 = '0_Rep1_DROP',ident.2 = '0_Rep1_DRONC', min.pct = 0.25, test.use = test)
res[[2]] <- FindMarkers(scomb, ident.1 = '0_Rep2_DROP',ident.2 = '0_Rep2_DRONC', min.pct = 0.25, test.use = test)

res[[3]] <- FindMarkers(scomb, ident.1 = '1_Rep2_DROP',ident.2 = '1_Rep2_DRONC', min.pct = 0.25, test.use = test)

res[[4]] <- FindMarkers(scomb, ident.1 = '3_Rep1_DROP',ident.2 = '3_Rep1_DRONC', min.pct = 0.25, test.use = test)
res[[5]] <- FindMarkers(scomb, ident.1 = '3_Rep2_DROP',ident.2 = '3_Rep2_DRONC', min.pct = 0.25, test.use = test)

res[[6]] <- FindMarkers(scomb, ident.1 = '7_Rep1_DROP',ident.2 = '7_Rep1_DRONC', min.pct = 0.25, test.use = test)
res[[7]] <- FindMarkers(scomb, ident.1 = '7_Rep2_DROP',ident.2 = '7_Rep2_DRONC', min.pct = 0.25, test.use = test)

res[[8]] <- FindMarkers(scomb, ident.1 = '15_Rep2_DROP',ident.2 = '15_Rep2_DRONC', min.pct = 0.25, test.use = test)

for(i in 1:length(res)){
  res[[i]]$gene <- row.names(res[[i]])
  res[[i]]$id <- i
}

markers_tech <- Reduce(rbind, res)

save(markers_tech, file='metadata/markers_tech')


## PLOTTING
### Start here if the above files are already generated
load('markers_days.Robj')
load('markers_batch.Robj')
load('markers_tech.Robj')

pval.thresh <- 0.05
logfc.thresh <- 4

DEGs_per_day <- as.vector(table(markers_days$id[markers_days$p_val_adj < pval.thresh & abs(markers_days$avg_logFC) > logfc.thresh]))
DEGs_per_rep <- as.vector(table(markers_batch$id[markers_batch$p_val_adj < pval.thresh & abs(markers_batch$avg_logFC) > logfc.thresh]))
DEGs_per_tech <- as.vector(table(markers_tech$id[markers_tech$p_val_adj < pval.thresh & abs(markers_tech$avg_logFC) > logfc.thresh]))

DEGs.df <- data.frame(DEGs=c(DEGs_per_day,DEGs_per_rep,DEGs_per_tech), 
                      type=c(rep('Between Days',length(DEGs_per_day)),rep('Between Cell Lines',length(DEGs_per_rep)),rep('Between Methods',length(DEGs_per_tech))))

ggplot(DEGs.df, aes(x=type, y=DEGs,fill=type)) + 
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none',text=element_text(size = 18), axis.text = element_text(size=16)) +
  xlab('') + 
  ylab('Differentially Expressed Genes') + 
  scale_y_continuous(breaks = pretty(DEGs.df$DEGs, n = 5))

tech.genes <- unique(markers_tech[markers_tech$p_val_adj < pval.thresh & abs(markers_tech$avg_logFC) > logfc.thresh,'gene'])
n <- length(tech.genes)

drop.genes <- unique(markers_tech[markers_tech$p_val_adj < pval.thresh & markers_tech$avg_logFC > logfc.thresh,'gene'])
linc <- readLines('metadata/lincRNAs.txt')

mito.genes <- length(grep('^MT',drop.genes,value=T))
ribo.genes <- length(grep('^RP',drop.genes,value=T))
RNA.genes <- length(intersect(drop.genes,linc))
drop.frac <- c(mito.genes,ribo.genes,RNA.genes)/length(drop.genes)

dronc.genes <- unique(markers_tech[markers_tech$p_val_adj < pval.thresh & markers_tech$avg_logFC < -logfc.thresh,'gene'])

mito.genes <- length(grep('^MT',dronc.genes,value=T))
ribo.genes <- length(grep('^RP',dronc.genes,value=T))
RNA.genes <- length(intersect(dronc.genes,linc))
dronc.frac <- c(mito.genes,ribo.genes,RNA.genes)/length(dronc.genes)

drdnc.compare <- data.frame(frac=c(drop.frac, dronc.frac),
                            type=rep(c('Mitochondrial','Ribosomal','lncRNAs'),2),
                            tech=rep(c(rep('Drop',3),rep('DroNc',3))))

ggplot(drdnc.compare, aes(x=type, y=frac,fill=tech)) + 
  geom_bar(stat = 'identity',position="dodge",color="black") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none',text=element_text(size = 18), axis.text = element_text(size=16)) +
  xlab('') + 
  ylab('Fraction of DEGs') 


gene <- markers_tech$gene[markers_tech$avg_logFC>0]
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)

ego_drop <- enrichGO(gene     = gene.df$ENSEMBL,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

barplot(ego_drop, drop=TRUE, showCategory=20, title = 'Upregulated Dropseq GO Analysis')


gene <- row.names(markers_tech[markers_tech$avg_logFC<0,])
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)

ego_dronc <- enrichGO(gene          = gene.df$ENTREZID,
                      universe      = names(geneList),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

barplot(ego_dronc, drop=TRUE, showCategory=20, title = 'Upregulated DroNC GO Analysis')



