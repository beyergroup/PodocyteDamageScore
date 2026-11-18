#### function for generating a Damage Signature ####

# Output: Table of genes by average rank, in ascending order. 
# Contains the following columns:
# gene_symbol - conventional Gene names, mean_rank - average rank  and 
# direction_foldchange - common direction of fold change, where '-1' and '1'
# represent negative and positive fold changes, respectively.

damage_signature.func <- function( DElist , thresh= 0.75 ,
                                   geneFilter )
  {
  require(matrixStats)
  ### *file_DElist* argument specifies an input: a .rda file, containing a list of DE tables, 
  # each table represents an experiment, and necessarily contains:
  # a column "Gene.Symbol" with conventional gene names
  # a column "log2FoldChange" with log2 fold changes, !!! experimental samples are compared against controls
  # a column "pvalue" with raw p-values
  
  ### *thresh* argument sets a threshold on proportion of studies in which 
  # a) a gene should be detected, and then
  # b) the gene should have the same direction of expression changes upon damage,
  # to be considered as a damage marker. 
  
  ### geneFilter are genes that should be considered, 
  # e.g. only genes detected in sc studies
  
  # extract gene names from all tables 
  allGenes <- Reduce( union, lapply( DElist, function(X) X$Gene.Symbol ) )
  
  # pass through a specific gene filter
  allGenes <- allGenes[ allGenes %in% geneFilter ]
  
  # create table p-value per gene for each study
  allStud_pval <- Reduce( cbind, lapply(DElist, function(X) 
    X$pvalue[match(allGenes, X$Gene.Symbol)]))
  
  # create table lfc per gene for each study
  allStud_lfc <- Reduce(cbind, lapply(DElist, function(X)
    X$log2FoldChange[match(allGenes, X$Gene.Symbol)]))
  
  rownames(allStud_pval) <- rownames(allStud_lfc ) <- allGenes
  
  
  #first sub-setting step: gene should be present in at least 75% of the studies 
  allStud_pval <- allStud_pval[rowSums( !is.na(allStud_pval)) >= thresh*ncol(allStud_pval),]
  allStud_lfc <- allStud_lfc[rownames( allStud_lfc ) %in% rownames(allStud_pval),]
  #second subs-setting step
  allStud_lfc <- allStud_lfc[rowSums(allStud_lfc > 0, na.rm = T) >= thresh*ncol(allStud_lfc) |
                               rowSums(allStud_lfc <0, na.rm = T )>= thresh*ncol(allStud_lfc) ,]
  allStud_pval <- allStud_pval[rownames(allStud_pval) %in% 
                                 rownames(allStud_lfc ),]
  
  # generate ranking 
  # be aware, lowest p-value should be on top of list (rank 1)
  
  allStud_ranks <- apply( allStud_pval,2, rank )
  
  # output table: gene symbol / averaged rank / positive or negative fold change 
  final_table <- data.frame('gene_symbol' = rownames(allStud_ranks), 
                            'mean_rank' = rowMeans(allStud_ranks, na.rm = T),
                            'mean_rankSC' = rowMeans(allStud_ranks[,
                                                                   grep("sc",names(DElist),ignore.case = T)], na.rm = T),
                            'mean_rankBulk' = rowMeans(allStud_ranks[,
                                                                   -grep("sc",names(DElist),ignore.case = T)], na.rm = T),
                            'direction_foldchange' = sign(rowMeans(apply(allStud_lfc,2,sign),
                                                                   na.rm = T)))
  
  final_table <- final_table[ order( final_table$mean_rank, decreasing=F,
                                   na.last=T), ]
  
  return(final_table)
}

