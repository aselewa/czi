
# DEPENDENCIES

suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

#' Combines intron/exon expression matrices into a single matrix
#' 
#' @param exon_mat - exon-based count matrix
#' @param intron_mat - intron-based count matrix
#' 
#' @return combined_mat
combine_mats = function(exon_mat, intron_mat){
  
  cells = unique(c(colnames(exon_mat),colnames(intron_mat)))
  genes = unique(c(row.names(exon_mat),row.names(intron_mat)))
  
  combined_mat = data.frame(matrix(0,nrow=length(genes),ncol=length(cells)), row.names=genes)
  colnames(combined_mat) = cells
  
  combined_mat[row.names(exon_mat),colnames(exon_mat)] = exon_mat
  combined_mat[row.names(intron_mat),colnames(intron_mat)] = combined_mat[row.names(intron_mat),colnames(intron_mat)] + intron_mat
  
  return(combined_mat)
  
}

#' Filters cells with number of genes outside of (min.genes, max.genes)
#' 
#' @param data - expression matrix
#' @param min.genes - minimum number of genes (atleast 1 UMI)
#' @param max.genes - maximum number of genes
#' 
#' @return data
filter = function(data, min.genes, max.genes){
  
  nGenes = colSums(data>0)

  keep.cells = nGenes>min.genes & nGenes<max.genes
  data = data[,keep.cells]
  
  return(data)
  
}
#' Converts Ensembl IDs to Short Gene Names
#' 
#' @param DGEmatrix - expression matrix
#' @param genenames - lookup table for Ensembl ID/short gene name
#' 
#' @return DGEmatrix 
matchIDToSymbol = function(DGEmatrix, genenames){
  
  g = genenames[row.names(DGEmatrix),1]
  bol = !duplicated(g)
  DGEmatrix = DGEmatrix[bol,]
  row.names(DGEmatrix) = g[bol]
  
  return(DGEmatrix)
}

#' Get metadata based on cell names from DGE
#' 
#' @param cell.ids Cell (column) names in DGE
#' @param type - which information to return: day, batch
#' 
#' @return cell tag
getCellTag <- function(cell.ids,type){
  
  cells <- strsplit(cell.ids,split='_')
  tag <- c()
  if(type=='day'){
    tag <- factor(sapply(cells, FUN=function(x){x[2]}),levels = c('0','1','3','7','15'))
    names(tag) <- cell.ids
  } else if(type=='batch'){
    tag <- as.factor(sapply(cells, FUN=function(x){x[3]}))
    names(tag) <- cell.ids
  } else{
    print('type not found')
  }
  return(tag)
}

#' computes log(CPM) on DGE matrix
#' 
#' @param DGEmatrix digital expression matrix
#' 
#' @return logCPM
computeLogCPM <- function(DGEmatrix){
  
  log.cpm <- log(t((t(DGEmatrix)/colSums(DGEmatrix))*10^4)+1)
  
  return(log.cpm)
  120808
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


