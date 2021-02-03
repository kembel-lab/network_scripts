library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(picante); packageVersion("picante")
library(reshape2)
library(seqTools)
library(adespatial)
library(ade4)
library(adegraphics)
library(spdep)
library(maptools)
library(akima)
library(gridExtra)
library(bipartite)
library(cooccur)
library(SpiecEasi)
library(Matrix)
library(igraph)
library(compositions)
library(tidyverse)
library(ggrepel)
library(decontam)
##INPUTS
meta <- read.table(file="organized2018_extraction_metadata.txt", sep="\t", header=TRUE, row.names="ID",stringsAsFactors = FALSE)
taxo <- read.csv(file="taxonomy_DADA2.csv", header=TRUE, row.names=1)
comm <- read.csv(file="community_DADA2.csv", header=TRUE, row.names=1)

##PREPROCESSING OF THE TABLES
par(mfrow=c(1,1))
netw.comm <- comm
netw.meta <- meta
netw.taxo <- taxo

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(netw.taxo)){
  levels(netw.taxo[,i])<-c(levels(netw.taxo[,i]), "unassigned", NA)
  for(j in rownames(netw.taxo)){
    if(is.na(netw.taxo[j,i])){
      netw.taxo[j,i] <- "unassigned"
    }
  }
}


netw.trash <- netw.taxo[c(netw.taxo$Kingdom == "unassigned" |netw.taxo$Kingdom == "Eukaryota" | netw.taxo$Family == "Mitochondria"| netw.taxo$Class == "Chloroplast"  ),]
netw.taxo <- subset(netw.taxo, !(rownames(netw.taxo) %in% rownames(netw.trash)))



#match the metadata and community tables
netw.comm <- netw.comm[rownames(netw.meta),]
netw.comm <- netw.comm[,rownames(netw.taxo)]


netw.taxo <- data.frame(lapply(netw.taxo, as.character), stringsAsFactors=FALSE)

#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(netw.taxo)))){
  for(j in rownames(netw.taxo)){
    if(netw.taxo[j,i]=="unassigned"){
      if(substr(netw.taxo[j,i-1], start=1,stop=2)=="UA"){
        netw.taxo[j,i] <- netw.taxo[j,i-1]
      }
      else{
        netw.taxo[j,i] <- paste0("UA_",netw.taxo[j,i-1])
      }
    }
  }
}

#remove ASVs with less than 0,1% reads
#netw.comm <- netw.comm[,apply(netw.comm,2,sum)> (sum(apply(netw.comm,2,sum))/1000)]
#netw.taxo <- netw.taxo[colnames(netw.comm),]

#change the names of our ASVs to something more practical
n <- 0
sequences <- character()
for(i in colnames(netw.comm)){
  n <- n +1
  sequences[n] <- colnames(netw.comm)[n]
  colnames(netw.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(netw.comm)[n]
}
rm(n,i,j)
rownames(netw.taxo) <- colnames(netw.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(netw.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(netw.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(netw.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation_level"){
    colnames(netw.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#remove the positive controls
netw.meta <- netw.meta[netw.meta$sample_type!="pcr_pos",]
netw.comm <- netw.comm[rownames(netw.meta),]

#decontam with extraction negatives
ps<-phyloseq(otu_table(netw.comm,taxa_are_rows = FALSE), tax_table(as.matrix(netw.taxo)), sample_data(netw.meta))
sample_data(ps)$is.neg <- sample_data(ps)$Site == "extraction_neg"
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

notcontam <- rownames(contamdf.prev05)[contamdf.prev05$contaminant==FALSE]
netw.nc.comm <- netw.comm[,notcontam]
netw.nc.taxo <- netw.taxo[notcontam,]
netw.nc.meta <- netw.meta

#remove negative controls
netw.nc.meta <- netw.nc.meta[netw.nc.meta$sample_type!="pcr_neg"&netw.nc.meta$sample_type!="extraction_neg",]
netw.nc.comm <- netw.nc.comm[rownames(netw.nc.meta),]

#add more metadata to our sites
sitesSpecs <- data.frame(read.csv(file="site_specs_2.csv", header=TRUE, row.names="ID1", stringsAsFactors = FALSE))
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(netw.nc.meta)){
  for (siteName in rownames(sitesSpecs)){
    if (siteName == netw.nc.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,sitesSpecs[siteName,c(5:11)])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
netw.nc.meta <- cbind(netw.nc.meta[rownames(tempMeta),], tempMeta)
rm(n,sitesSpecs,tempMeta,netw.nc.trash, colname, rowname, siteName)
netw.nc.meta$site_specie <- paste(netw.nc.meta$Site, netw.nc.meta$tree_species)

#make sure all the numbers are numeric
list_of_numerics <- c("host_tree_damage","Latitude","Longitude","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity")
for (number in list_of_numerics){
  netw.nc.meta[,number] <- as.numeric(netw.nc.meta[,number])
}
rm(number,list_of_numerics)

#keep only larvae with all the metadata
netw.nc.meta <- netw.nc.meta[netw.nc.meta$sample_type == "L6"&!is.na(netw.nc.meta$host_tree_damage),]
netw.nc.comm <- netw.nc.comm[rownames(netw.nc.meta),]


#try occurence filtration (transform into presence/absence then look at distribution of columns totals) before rarefaction
pa <- decostand(netw.nc.comm,method = "pa")
hist(apply(pa,2,sum))

#Do first steps for 5 occ

#try occurence in 5 samples
netw.nc.comm <- netw.nc.comm[,apply(pa,2,sum) > 5]
##rarefy our dataset
quantile(apply(netw.nc.comm, 1,sum))
netw.nc.ab.comm <- rrarefy(netw.nc.comm, 2000)
netw.nc.ab.comm <- netw.nc.ab.comm[apply(netw.nc.ab.comm,1,sum)>= 2000,]
rarecurve(netw.nc.ab.comm, step=100, label=FALSE)
netw.nc.ab.meta <- netw.nc.meta[rownames(netw.nc.ab.comm),]
netw.nc.ab.taxo <- netw.nc.taxo[colnames(netw.nc.ab.comm),]


#make a phyloseq object
OTU <- otu_table(netw.nc.ab.comm, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(netw.nc.ab.taxo))
samples <- sample_data(netw.nc.ab.meta)
phylo <- phyloseq(OTU, TAX, samples)
phylo

#calculate bootstraps
tp0 <- proc.time()
sparcc.boot <- sparccboot(otu_table(phylo), R=1000, ncpus=4)
tp1 <- proc.time()
tp1 - tp0



#do a first cooccurence table
sparcc <- sparcc(otu_table(phylo))
# sparcc$Cor[1:5, 1:5]
# sparcc$Cov[1:5, 1:5]
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



sparcc.high.corr <- SparccP%>% filter(abs(cor) > 0.0)
sparccOkP <- sparcc.high.corr%>% filter(fdr < 0.01)

print(paste(dim(sparccOkP)[1]," correlations with adjusted p-values under your", " threshold were found!", sep=""))

#generate a matrix of fdr-adjusted p-values for igraph
mat.sparccP <- acast(sparcc.high.corr, Var1~Var2 )
mat.sparccP[lower.tri(mat.sparccP)] <- t(mat.sparccP)[lower.tri(mat.sparccP)]
mat.sparccP[is.na(mat.sparccP)] <- 1

SparccP_plot <- sparcc.high.corr %>% ggplot(aes(x = Var2, y = Var1, fill = cor)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = 1)
SparccP_plot

sparcc.graph <- abs(mat.sparccP) < 0.01
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
phylo.filt
print(paste("So ",length(rownames(tax_table(phylo.filt)))," out of the ",length(rownames(tax_table(phylo)))," taxa were kept in the final phyloseq object",sep=""))
sparcc.graph.filt <- Matrix(sparcc.graph.filt, sparse=TRUE)



ig.sparcc <- adj2igraph(sparcc.graph.filt,  vertex.attr=list(name=taxa_names(phylo.filt)))


vsize    <- colMeans(clr(otu_table(phylo.filt),1))*5 + 5
#am.coord <- layout.fruchterman.reingold(ig.sparcc)

network.plot <- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Phylum", label="Genus") 
network.plot$layers <- network.plot$layers[-2]
network.plot <-  network.plot +
  geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
  theme(legend.title = element_text(size=10,face="bold"),
        legend.text = element_text(size=10))
network.plot

write.graph(ig.sparcc,file = "network_test.graphml",format = "graphml")
