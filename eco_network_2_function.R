eco_network_2<-function(sparcc.boot,phylo, threshold, corThreshold=0, sign="positive"){
  
  #do a first run
  sparcc <- sparcc(otu_table(phylo))
  #calculate p-values
  sparcc.pvals <- pval.sparccboot(sparcc.boot)
  pvals.dat <- data.frame(sparcc.pvals$cors, sparcc.pvals$pvals) 
  #reorganize the p-value matrix
  cors <- sparcc.pvals$cors
  pvals <- sparcc.pvals$pvals
  sparCCpcors <- diag(0.5, nrow = dim(sparcc$Cor)[1], ncol = dim(sparcc$Cor)[1])
  sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
  sparCCpcors <- sparCCpcors + t(sparCCpcors)
  
  sparCCpval <- diag(0.5, nrow = dim(sparcc$Cor)[1], ncol = dim(sparcc$Cor)[1])
  sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
  sparCCpval <- sparCCpval + t(sparCCpval)
  
  rownames(sparCCpcors) <- colnames(otu_table(phylo))
  colnames(sparCCpcors) <- colnames(otu_table(phylo))
  rownames(sparCCpval) <- colnames(otu_table(phylo))
  colnames(sparCCpval) <- colnames(otu_table(phylo))
  
  sparCCpcors[1:5, 1:5]
  sparCCpval[1:5, 1:5]
  
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  reorder_cor_and_p <- function(cormat, pmat){
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
    pmat <- pmat[hc$order, hc$order]
    list(r = cormat, p = pmat)
  }
  
  reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
  reordered_sparccCor <- reordered_all_sparcc$r
  reordered_sparccP<- reordered_all_sparcc$p
  
  
  
  sparccCor_processed <- reordered_sparccCor  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(cor = value)
  sparccP_processed <- reordered_sparccP  %>% get_upper_tri() %>% reshape2::melt() %>% na.omit() %>% rename(p = value)
  
  SparccP <- left_join(sparccCor_processed, sparccP_processed, by = c("Var1", "Var2")) %>%
    
    mutate(fdr = p.adjust(p, method = "BH"))
  
  #remove NAs
  SparccP$p[is.na(SparccP$p)] <- 1
  SparccP$fdr[is.na(SparccP$fdr)] <- 1
  
  
  sparcc.high.corr <- SparccP%>% filter(abs(cor) > corThreshold)
  sparccOkP <- sparcc.high.corr%>% filter(fdr < threshold)
  #keep only positive interactions
  sparcc.pos <- sparcc.high.corr%>% filter(cor > 0.00)
  sparcc.neg <- sparcc.high.corr%>% filter(cor < 0.00)
  
  print(paste(dim(sparccOkP)[1]," correlations with adjusted p-values under your", " X ", "threshold were found!", sep=""))
  
  #generate a matrix of fdr-adjusted p-values for igraph
  if (sign == "both"){
    mat.sparccP <- acast(sparcc.high.corr, Var1~Var2 )
    sparcc.pos$sign <- "pos"
    sparcc.neg$sign <- "neg"
    list <- rbind(sparcc.pos,sparcc.neg)
  }
  else if (sign == "positive"){
    mat.sparccP <- acast(sparcc.pos, Var1~Var2 )
    list <- sparcc.pos
  }
  else if (sign == "negative"){
    mat.sparccP <- acast(sparcc.neg, Var1~Var2 )
    list <- sparcc.neg
  }
  else {
    return('Wrong sign argument, please choose between "positive", "negative" and "both"')
  }
  
  mat.sparccP[lower.tri(mat.sparccP)] <- t(mat.sparccP)[lower.tri(mat.sparccP)]
  mat.sparccP[is.na(mat.sparccP)] <- 1
  
  SparccP_plot <- sparcc.high.corr %>% ggplot(aes(x = Var2, y = Var1, fill = cor)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = 1)
  
  sparcc.graph <- abs(mat.sparccP) < threshold
  print(paste("So a matrix with ",table(sparcc.graph > 0)[2]/2," non-zero relationships was made",sep=""))
  diag(sparcc.graph) <- 0
  #keep a vector of which taxa have more than one relationship over the threshhold
  posit.vector <- apply(sparcc.graph,2,sum) > 0
  table(posit.vector)
  
  print(paste("Out of the ",length(posit.vector)," taxa in the analysis, ",length(posit.vector)-table(posit.vector)["TRUE"]," did not have any significant relationship with another taxa",sep=""))
  #keep a sparcc.graph version of non-zero relationships
  sparcc.graph.filt <- sparcc.graph[posit.vector,posit.vector]
  #edit our phyloseq object to keep ASVs with a least one significant relationship
  vect.filt <- rownames(tax_table(phylo)) %in% rownames(sparcc.graph.filt)
  
  #make a new phyloseq object
  OTU <- otu_table(phylo)[,vect.filt]
  TAX <- tax_table(phylo)[vect.filt,]
  DATA <- sample_data(phylo)
  phylo.filt <- phyloseq(OTU,TAX,DATA)
  print(paste("So ",length(rownames(tax_table(phylo.filt)))," out of the ",length(rownames(tax_table(phylo)))," taxa were kept in the final phyloseq object",sep=""))
  sparcc.graph.filt <- Matrix(sparcc.graph.filt, sparse=TRUE)
  
  
  
  return(list(SparccP_plot, sparcc.graph.filt,phylo.filt,list))
  
}
