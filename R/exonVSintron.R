
source('R/required_funs_libs.R')

## Plot detection rate with/without introns

drop.introns <- data.frame(fread('data/DROP_combined_allcells.tsv.gz',sep='\t'),row.names=1)
drop.exons <- data.frame(fread('data/DROP_combined_allcells_exons.tsv',sep='\t'),row.names=1)

same_cells <- intersect(colnames(drop.introns),colnames(drop.exons))
same_genes <- intersect(row.names(drop.introns), row.names(drop.exons))
drop.intron.rate <- colMeans(drop.introns[same_genes,same_cells]>0)
drop.exon.rate <- colMeans(drop.exons[same_genes,same_cells]>0)

dronc.introns <- data.frame(fread('data/DRONC_combined_allcells.tsv',sep='\t'),row.names=1)
dronc.exons <- data.frame(fread('data/DRONC_combined_allcells_exons.tsv',sep='\t'),row.names=1)

same_cells <- intersect(colnames(dronc.introns),colnames(dronc.exons))
same_genes <- intersect(row.names(dronc.introns), row.names(dronc.exons))
dronc.intron.rate <- colMeans(dronc.introns[same_genes,same_cells]>0)
dronc.exon.rate <- colMeans(dronc.exons[same_genes,same_cells]>0)

det.rate.df <- data.frame(intron = c(drop.intron.rate,dronc.intron.rate), exon = c(drop.exon.rate,dronc.exon.rate),
                          type = c(rep('Drop',length(drop.intron.rate)), rep('DroNc',length(dronc.intron.rate))))
                          
ggplot(det.rate.df, aes(x=exon, y=intron, color=type)) + geom_point() + geom_abline(slope=1) + theme_bw() +
  xlab('Detection Rate (Exons)') + ylab('Detection Rate (Introns+Exons)') + 
  theme(text=element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Exon vs intron for a few genes

dronc.introns <- data.frame(fread('data/DRONC_combined_allcells.tsv',sep='\t'),row.names=1)
dronc.exons <- data.frame(fread('data/DRONC_combined_allcells_exons.tsv',sep='\t'),row.names=1)

same_cells <- intersect(colnames(dronc.introns),colnames(dronc.exons))
intron <- CreateSeuratObject(dronc.introns[,same_cells],min.genes = 0)
exons <- CreateSeuratObject(dronc.exons[,same_cells], min.genes = 0)

intron <- NormalizeData(intron)
exons <- NormalizeData(exons)


bmp.intron <- intron@data['BMPER',colnames(exons@data)]
bmp.exon <- exons@data['BMPER',]
tnnt2.intron <- intron@data['TNNT2',colnames(exons@data)]
tnnt2.exon <- exons@data['TNNT2',]
myh.intron <- intron@data['MYH6',colnames(exons@data)]
myh.exon <- exons@data['MYH6',]

exp.df <- data.frame(tnnt2.exon,tnnt2.intron,bmp.intron,bmp.exon,myh.intron,myh.exon)

ggplot(exp.df, aes(x=tnnt2.exon,y=tnnt2.intron)) + geom_point() + theme_bw() +
  xlab('') + ylab('') + ggtitle('TNNT2') +
  theme(text=element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_rect(aes(xmin=-0.1,xmax=0.1,ymin=0.8,ymax=8.5), alpha=0, color="red") + 
  annotate("text", x = 1.25, y = 6, label="N[recovered] == 188", parse=T)

ggplot(exp.df, aes(x=myh.exon,y=myh.intron)) + geom_point() + theme_bw() +
  xlab('') + ylab('') + ggtitle('MYH6') + 
  theme(text=element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_rect(aes(xmin=-0.1,xmax=0.1,ymin=0.5,ymax=7.5), alpha=0, color="red") + 
  annotate("text", x = 1.5, y = 6, label="N[recovered] == 555", parse=T)

ggplot(exp.df, aes(x=bmp.exon,y=bmp.intron)) + geom_point() + theme_bw() +
  xlab('') + ylab('') + ggtitle('BMPER') +
  theme(text=element_text(size=20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_rect(aes(xmin=-0.1,xmax=0.1,ymin=0.1,ymax=6), alpha=0, color="red") + 
  annotate("text", x = 1, y = 4, label="N[recovered] == 3207", parse=T)
