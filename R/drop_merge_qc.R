
### Drop-seq
############################################ Rep 1 ############################################
###############################################################################################

# Day 0
day0_Rep1 = combine_mats(data.frame(fread('data/Rep1/AB-HE1120A-Drop-Day0_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep1/AB-HE1120A-Drop-Day0_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day0_Rep1) = paste0(colnames(day0_Rep1),'_0_Rep1')

# Day 1
day1_Rep1 = combine_mats(data.frame(fread('data/Rep1/AB-HE0105C-Drop-Day1_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep1/AB-HE0105C-Drop-Day1_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day1_Rep1) = paste0(colnames(day1_Rep1),'_1_Rep1')

# Day 3
day3_Rep1 = combine_mats(data.frame(fread('data/Rep1/AB-HE1120B-Drop-Day3_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep1/AB-HE1120B-Drop-Day3_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day3_Rep1) = paste0(colnames(day3_Rep1),'_3_Rep1')

# Day 7
day7_Rep1 = combine_mats(data.frame(fread('data/Rep1/AB-EH1129-Drop-Day7_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep1/AB-EH1129-Drop-Day7_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day7_Rep1) = paste0(colnames(day7_Rep1),'_7_Rep1')

############################################ Rep 2 ############################################
###############################################################################################

# Day 0
day0_Rep2 = combine_mats(data.frame(fread('data/Rep2/AB-HE0202A-CZI-Drop-Day0_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep2/AB-HE0202A-CZI-Drop-Day0_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day0_Rep2) = paste0(colnames(day0_Rep2),'_0_Rep2')

# Day 1

day1_Rep2 = combine_mats(data.frame(fread('data/Rep2/AB-HE0202B-CZI-Drop-Day1_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep2/AB-HE0202B-CZI-Drop-Day1_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day1_Rep2) = paste0(colnames(day1_Rep2),'_1_Rep2')

# Day 3

day3_Rep2 = combine_mats(data.frame(fread('data/Rep2/AB-HE0202B-CZI-Drop-Day3_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep2/AB-HE0202B-CZI-Drop-Day3_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day3_Rep2) = paste0(colnames(day3_Rep2),'_3_Rep2')

# Day 7

day7_Rep2 = combine_mats(data.frame(fread('data/Rep2/AB-HE0202-C-CZI-Drop-Day7_counts_exon.tsv.gz',sep='\t'),row.names=1),
                         data.frame(fread('data/Rep2/AB-HE0202-C-CZI-Drop-Day7_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day7_Rep2) = paste0(colnames(day7_Rep2),'_7_Rep2')

# Day 15
day15_Rep2 = combine_mats(data.frame(fread('data/Rep2/AB-HE0307-CZI-Drop-Day15_counts_exon.tsv.gz',sep='\t'),row.names=1),
                          data.frame(fread('data/Rep2/AB-HE0307-CZI-Drop-Day15_counts_intron.tsv.gz',sep='\t'),row.names=1))

colnames(day15_Rep2) = paste0(colnames(day15_Rep2),'_15_Rep2')


############################################ Merge into one matrix ############################
###############################################################################################

genes = c(row.names(day0_Rep1),row.names(day0_Rep2),row.names(day1_Rep1),row.names(day1_Rep2),row.names(day3_Rep1),row.names(day3_Rep2),
          row.names(day7_Rep1),row.names(day7_Rep2),row.names(day15_Rep2))
genes = unique(genes)
cells = c(colnames(day0_Rep1),colnames(day0_Rep2),colnames(day1_Rep1),colnames(day1_Rep2),colnames(day3_Rep1),colnames(day3_Rep2),
          colnames(day7_Rep1),colnames(day7_Rep2),colnames(day15_Rep2))

comb = as.data.frame(matrix(0, nrow=length(genes), ncol=length(cells)))
row.names(comb) = genes
colnames(comb) = cells

comb[row.names(day0_Rep1),colnames(day0_Rep1)] = day0_Rep1
comb[row.names(day0_Rep2),colnames(day0_Rep2)] = day0_Rep2
comb[row.names(day1_Rep1),colnames(day1_Rep1)] = day1_Rep1
comb[row.names(day1_Rep2),colnames(day1_Rep2)] = day1_Rep2
comb[row.names(day3_Rep1),colnames(day3_Rep1)] = day3_Rep1
comb[row.names(day3_Rep2),colnames(day3_Rep2)] = day3_Rep2
comb[row.names(day7_Rep1),colnames(day7_Rep1)] = day7_Rep1
comb[row.names(day7_Rep2),colnames(day7_Rep2)] = day7_Rep2
comb[row.names(day15_Rep2),colnames(day15_Rep2)] = day15_Rep2

rm(day0_Rep1,day0_Rep2, day1_Rep1,day1_Rep2,day3_Rep1,day3_Rep2,day7_Rep1,day7_Rep2,day15_Rep2)
gc()

save(comb, file='DROP_combined_all')

comb = comb[rowSums(comb>0)>10,]

genenames = read.delim(file='metadata/gencode.v27.annotation.ID_to_Symbol.txt', header=F, sep="\t", row.names=1)
comb = matchIDToSymbol(comb, genenames)

# Remove low quality cells
comb = filter(data = comb, 
                  min.genes = 400, 
                  max.genes = 2000)
# Remove chimp cells
species_assign = data.frame(fread('metadata/species_assign_matrix.txt',sep='\t'))
chimp.cells = species_assign$CB[species_assign$hg_specificity_score<0.5]
cells = strsplit(colnames(comb),split='_')
barcodes = sapply(cells, function(x){x[1]})
comb = comb[,!(barcodes %in% chimp.cells)]

png("figures/species_mixing.png", res=200, width = 1000, height = 800)
ggplot(species_assign, aes(x=hg_specificity_score)) +
  geom_histogram(bins = 30) +
  theme(text=element_text(size=18)) +
  xlab('hg38 Specificity Score') + 
  ylab('Number of Cells')
dev.off()

fwrite(comb, file = "data/DROP_combined_m400.tsv",quote = F, sep='\t',row.names = T, col.names = T)
system('gzip data/DROP_combined_m400.tsv')


############################################ Plot QC metrics ##################################
###############################################################################################

cells = strsplit(colnames(comb),split='_')
uniq = c('Day 0_Rep1','Day 0_Rep2','Day 1_Rep1','Day 1_Rep2','Day 3_Rep1','Day 3_Rep2','Day 7_Rep1','Day 7_Rep2','Day 15_Rep2')
batch_day = factor(sapply(cells, FUN=function(x){paste0('Day ',x[2],'_',x[3])}), levels=uniq)

total.genes <- c()
total.nUMI <- c()
total.cells <- c()
experiment <- c()
for(i in 1:length(uniq)){
  bol <- batch_day == uniq[i]
  ngenes <- colSums(comb[,bol]>0)
  nUMI <- colSums(comb[,bol])
  total.nUMI <- c(total.nUMI,nUMI)
  total.genes <- c(total.genes,ngenes)
  
  ncells <- sum(bol)
  total.cells <- c(total.cells, ncells)
  experiment <- c(experiment, rep(uniq[i], ncells))
}

temp <- strsplit(as.character(experiment), split = '_')
days <- factor(sapply(temp, function(x){x[1]}), levels=c('Day 0','Day 1','Day 3','Day 7','Day 15'))
batch <- factor(sapply(temp, function(x){x[2]}), levels=c('Rep1','Rep2'))

qc.df <- data.frame(nGenes = total.genes, nUMI = total.nUMI ,Day = days, Batch = batch)
qc2.df <- data.frame(nCells = total.cells, 
                     Day = factor(c('Day 0','Day 1','Day 3','Day 7','Day 0','Day 1','Day 3','Day 7','Day 15'),levels = c('Day 0','Day 1','Day 3','Day 7','Day 15')),
                     Batch = c(rep('Rep1',4),rep('Rep2',5)))


p1 <- ggplot(data = qc.df, aes(x = Day, y=nUMI, fill=Batch)) + geom_boxplot() + ylab('Number of UMI') +  
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), text=element_text(size=15))

p2 <- ggplot(data = qc.df, aes(x = Day, y=nGenes, fill=Batch)) + geom_boxplot() + ylab('Number of Genes') +  
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), text=element_text(size=15))

p3 <- ggplot(data = qc2.df, aes(x=Day, y=nCells, fill=Batch)) + geom_bar(stat = "identity") + ylab('Number of Cells') +
  theme(axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), text=element_text(size=15))


png('figures/drop_QC.png',res=200, width = 1000, height=1500)
multiplot(p1,p2,p3)
dev.off()


fracs <- scale(as.matrix(comb),center=FALSE, scale=colSums(comb))
fracMeds <- rowMedians(fracs)
ord <- order(fracMeds,decreasing=T)[1:15]
fracsOrd <- fracs[ord,]*100
fracsOrdMelt <- melt(fracsOrd)
fracsOrdMelt$Var1 <- factor(fracsOrdMelt$Var1, levels = rev(row.names(fracs)[ord]))

fracsOrdMelt <- fracsOrdMelt[fracsOrdMelt$value < 7.5,]
png('figures/PercentCounts_Drop.png',res=200,width=1500,height=1000)
ggplot(fracsOrdMelt, aes(x=Var1, y=value,fill=Var1)) + geom_boxplot() + coord_flip() + theme(legend.position = 'none',text=element_text(size=18)) + 
  ylab('% of counts') + xlab('Gene Symbol')
dev.off()

