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
library(econullnetr)
library(intergraph)
library(ggnetwork)

##INPUTS
meta <- read.table(file="organized2018_extraction_metadata.txt", sep="\t", header=TRUE, row.names="ID",stringsAsFactors = FALSE)
taxo <- read.csv(file="taxonomy_DADA2.csv", header=TRUE, row.names=1)
comm <- read.csv(file="community_DADA2.csv", header=TRUE, row.names=1)

##PREPROCESSING OF THE TABLES
par(mfrow=c(1,1))
net.comm <- comm
net.meta <- meta
net.taxo <- taxo

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(net.taxo)){
  levels(net.taxo[,i])<-c(levels(net.taxo[,i]), "unassigned", NA)
  for(j in rownames(net.taxo)){
    if(is.na(net.taxo[j,i])){
      net.taxo[j,i] <- "unassigned"
    }
  }
}
net.trash <- net.taxo[c(net.taxo$Kingdom == "unassigned" |net.taxo$Kingdom == "Eukaryota" | net.taxo$Family == "Mitochondria"| net.taxo$Class == "Chloroplast"  ),]
net.taxo <- subset(net.taxo, !(rownames(net.taxo) %in% rownames(net.trash)))

#match the metadata and community tables
net.comm <- net.comm[rownames(net.meta),]
net.comm <- net.comm[,rownames(net.taxo)]

net.taxo <- data.frame(lapply(net.taxo, as.character), stringsAsFactors=FALSE)

#replace unassigned by the lowest  assigned taxonomic level
for(i in c(1:length(colnames(net.taxo)))){
  for(j in rownames(net.taxo)){
    if(net.taxo[j,i]=="unassigned"){
      if(substr(net.taxo[j,i-1], start=1,stop=2)=="UA"){
        net.taxo[j,i] <- net.taxo[j,i-1]
      }
      else{
        net.taxo[j,i] <- paste0("UA_",net.taxo[j,i-1])
      }
    }
  }
}

#remove ASVs with less than 0,1% reads
#net.comm <- net.comm[,apply(net.comm,2,sum)> (sum(apply(net.comm,2,sum))/1000)]
#net.taxo <- net.taxo[colnames(net.comm),]

#change the names of our ASVs to something more practical
n <- 0
sequences <- character()
for(i in colnames(net.comm)){
  n <- n +1
  sequences[n] <- colnames(net.comm)[n]
  colnames(net.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(net.comm)[n]
}
rm(n,i,j)
rownames(net.taxo) <- colnames(net.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(net.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(net.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(net.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation_level"){
    colnames(net.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#decontam
#remove the positive controls
net.meta <- net.meta[net.meta$sample_type!="pcr_pos",]
net.comm <- net.comm[rownames(net.meta),]

#decontam for pcr

ps<-phyloseq(otu_table(net.comm,taxa_are_rows = FALSE), tax_table(as.matrix(net.taxo)), sample_data(net.meta))
sample_data(ps)$is.neg <-  sample_data(ps)$sample_type == "pcr_neg"

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

notcontam <- rownames(contamdf.prev05)[contamdf.prev05$contaminant==FALSE]
net.nc.comm <- net.comm[,notcontam]
net.nc.taxo <- net.taxo[notcontam,]
net.nc.meta <- net.meta

#decontam for DNA extraction
#remove the controls
net.nc.meta <- net.nc.meta[net.nc.meta$sample_type != "pcr_neg",]
net.nc.comm <- net.nc.comm[rownames(net.nc.meta),]

ps<-phyloseq(otu_table(net.nc.comm,taxa_are_rows = FALSE), tax_table(as.matrix(net.nc.taxo)), sample_data(net.nc.meta))
sample_data(ps)$is.neg <-  sample_data(ps)$sample_type == "extraction_neg"

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

notcontam <- rownames(contamdf.prev05)[contamdf.prev05$contaminant==FALSE]
net.nc.comm <- net.nc.comm[,notcontam]
net.nc.taxo <- net.nc.taxo[notcontam,]
net.nc.meta <- net.nc.meta

apply(net.nc.comm,1,sum)
#remove the controls
net.nc.meta <- net.nc.meta[net.nc.meta$sample_type != "extraction_neg",]
net.nc.comm <- net.nc.comm[rownames(net.nc.meta),]

#add more metadata to our sites
sitesSpecs <- data.frame(read.csv(file="site_specs_2.csv", header=TRUE, row.names="ID1", stringsAsFactors = FALSE))
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(net.nc.meta)){
  for (siteName in rownames(sitesSpecs)){
    if (siteName == net.nc.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,sitesSpecs[siteName,c(5:11)])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
net.nc.meta <- cbind(net.nc.meta[rownames(tempMeta),], tempMeta)
rm(n,sitesSpecs,tempMeta, colname, rowname, siteName)
net.nc.meta$site_specie <- paste(net.nc.meta$Site, net.nc.meta$tree_species)

#make sure all the numbers are numeric
list_of_numerics <- c("host_tree_damage","Latitude","Longitude","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity")
for (number in list_of_numerics){
  net.nc.meta[,number] <- as.numeric(net.nc.meta[,number])
}
rm(number,list_of_numerics)

#keep only larvae with all the metadata
net.la.meta <- net.nc.meta[net.nc.meta$sample_type == "L6"&!is.na(net.nc.meta$host_tree_damage),]
net.la.comm <- net.nc.comm[rownames(net.la.meta),]
net.la.taxo <- net.nc.taxo
#keep only foliage with all the metadata
net.fo.meta <- net.nc.meta[net.nc.meta$sample_type == "Foliage"&!is.na(net.nc.meta$host_tree_damage),]
net.fo.comm <- net.nc.comm[rownames(net.fo.meta),]
net.fo.taxo <- net.nc.taxo

#aggregate by genus
gen.la.comm <- aggregate(t(net.la.comm), by=list(net.la.taxo$Genus), sum)
rownames(gen.la.comm) <- gen.la.comm[,1]
gen.la.comm <- t(gen.la.comm[,-1])

gen.la.taxo <- data.frame()
done <- c()
#n<-0
#make a taxonomy table for the aggregated community file
for (ASV in rownames(net.nc.taxo)){
  
  taxa <- net.nc.taxo[ASV,"Genus"]
  
  if (!(taxa%in%done)){
    #n <- n + 1
    print(net.nc.taxo[ASV,c(1:5)])
    gen.la.taxo <- rbind(gen.la.taxo,net.nc.taxo[ASV,c(1:6)])
    
    done <- c(done,taxa)
  }
}
rm(taxa,done,ASV)
rownames(gen.la.taxo) <- gen.la.taxo$Genus
gen.la.taxo <- gen.la.taxo[colnames(gen.la.comm),]

gen.la.meta <- net.la.meta[rownames(gen.la.comm),]



#try occurence filtration (transform into presence/absence then look at distribution of columns totals) before rarefaction
pa <- decostand(gen.la.comm,method = "pa")
hist(apply(pa,2,sum))


#try occurence in 5 samples
gen.la.comm <- gen.la.comm[,apply(pa,2,sum) > 5]
gen.la.taxo <- gen.la.taxo[colnames(gen.la.comm),]
##rarefy our dataset
quantile(apply(gen.la.comm, 1,sum))
gen.la.comm <- rrarefy(gen.la.comm, 3000)
gen.la.comm <- gen.la.comm[apply(gen.la.comm,1,sum)>= 3000,]
rarecurve(gen.la.comm, step=100, label=FALSE)
gen.la.meta <- gen.la.meta[rownames(gen.la.comm),]
gen.la.taxo <- gen.la.taxo[colnames(gen.la.comm),]

#make a phyloseq object
OTU <- otu_table(gen.la.comm, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(gen.la.taxo))
samples <- sample_data(gen.la.meta)
phylo <- phyloseq(OTU, TAX, samples)
phylo

?cooccur
test <- cooccur(t(gen.la.comm),spp_names=TRUE)
test$co_occurrences

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
#keep only positive interactions
sparcc.pos <- sparcc.high.corr%>% filter(cor > 0.00)
sparcc.neg <- sparcc.high.corr%>% filter(cor < 0.00)

print(paste(dim(sparccOkP)[1]," correlations with adjusted p-values under your", " threshold were found!", sep=""))
###only positive relationships
#generate a matrix of fdr-adjusted p-values for igraph
mat.sparccP <- acast(sparcc.pos, Var1~Var2 )
mat.sparccP[lower.tri(mat.sparccP)] <- t(mat.sparccP)[lower.tri(mat.sparccP)]
mat.sparccP[is.na(mat.sparccP)] <- 1

SparccP_plot <- sparcc.high.corr %>% ggplot(aes(x = Var2, y = Var1, fill = cor)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(data = sparccOkP, shape = 1)
SparccP_plot




sparcc.graph <- mat.sparccP < 0.01 ######### all relationships are kept, negative ones



print(paste("So a matrix with ",table(sparcc.graph > 0)[2]/2," positive non-zero relationships was made",sep=""))
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



ig.sparcc <- adj2igraph(sparcc.graph.filt,  vertex.attr=list(name=rownames(sparcc.graph.filt)))

phylo_comm <- otu_table(phylo.filt, taxa_are_rows = FALSE)@.Data
phylo_comm <- phylo_comm[,rownames(sparcc.graph.filt)]
vsize    <- colMeans(clr(phylo_comm,1))*5 + 5
#am.coord <- layout.fruchterman.reingold(ig.sparcc)
library(GGally)
library(intergraph)
library(ggnetwork)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


network.plot <- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Phylum", label="Genus") 
network.plot$layers <- network.plot$layers[-2]
network.plot <-  network.plot +
  geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
  theme(legend.title = element_text(size=10,face="bold"),
        legend.text = element_text(size=10))
network.plot

write.graph(ig.sparcc,file = "network_test.graphml",format = "graphml")

#generate null model
cons <- sparcc$Cor
rownames(cons) <- colnames(otu_table(phylo))
colnames(cons) <- colnames(otu_table(phylo))
res <- otu_table(phylo)
generate_null_net(consumers=cons,resources = res)



###both positive and negative###
#generate a matrix of fdr-adjusted p-values for igraph
mat.sparccP <- acast(sparcc.high.corr, Var1~Var2 )
mat.sparccP[lower.tri(mat.sparccP)] <- t(mat.sparccP)[lower.tri(mat.sparccP)]
mat.sparccP[is.na(mat.sparccP)] <- 1

sparccOkP.pos <- sparccOkP%>% filter(cor > 0.00)
sparccOkP.neg <- sparccOkP%>% filter(cor < 0.00)

sparcc.graph <- mat.sparccP < 0.01 



print(paste("So a matrix with ",table(sparcc.graph > 0)[2]/2,"positive non-zero relationships was made",sep=""))
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



ig.sparcc <- adj2igraph(sparcc.graph.filt,  vertex.attr=list(name=rownames(sparcc.graph.filt)))
edge_list <- data.frame(as_edgelist(ig.sparcc)) #Sometimes it is useful to work with a standard representation of a graph, like an edge list
sparccOkP.pos$pos <- "pos"
sparccOkP.neg$pos <- "neg"
list <- rbind(sparccOkP.pos,sparccOkP.neg)


#put the correlation table in the same order as the list of  edges
rownames(list) <- paste(list$Var1,list$Var2)

rownames(edge_list) <- paste(edge_list$X1,edge_list$X2)
rownames(list) %in% rownames(edge_list)

list <- list[rownames(edge_list),]

rownames(list) == rownames(edge_list)




phylo_comm <- otu_table(phylo.filt, taxa_are_rows = FALSE)@.Data
phylo_comm <- phylo_comm[,rownames(sparcc.graph.filt)]
phylo_taxa <- tax_table(phylo.filt)@.Data
phylo_taxa <- data.frame(phylo_taxa[rownames(sparcc.graph.filt),])
vsize    <- colMeans(clr(phylo_comm,1))*15 + 5
#am.coord <- layout.fruchterman.reingold(ig.sparcc)
#add the sign as a edge attribute 
edge_attr(ig.sparcc, "sign" ) <- list$pos
vertex_attr(ig.sparcc, "Phylum" ) <- as.character(phylo_taxa$Phylum)
vertex_attr(ig.sparcc, "Genus" ) <- as.character(phylo_taxa$Genus)
vertex_attr(ig.sparcc, "vsize" ) <- vsize
#add attributes to vertices for phyla, genus and abundance




#convert to a network object
ne.sparcc <- asNetwork(x = ig.sparcc)

#ggnet2(ne.sparcc,node.size = vsize, color = phylo_taxa$Phylum,edge.color = list$pos) + guides( size = FALSE) + geom_text_repel(label=rownames(sparcc.graph.filt)) + scale_color_manual(values=gg_color_hue(8)) + labs(color="Phylum")

ggnet <- ggnetwork(ne.sparcc, layout = "fruchtermanreingold")
#with red segments to labels
ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend)) + geom_edges(aes(linetype = sign), color = "red",size=0.75) + theme_blank() + geom_nodes(aes(color = Phylum, size = vsize)) +geom_nodetext_repel (aes(label = Genus),fontface = "bold", box.padding = unit(1, "lines"), segment.size = 0.5, segment.colour = "gray",cex = 3) + guides( size = FALSE) + scale_linetype_manual(labels=c("Negative","Positive"), values=c(2,1)) + labs(linetype="Type of correlation") + scale_size_continuous(range=c(3,18))

#only the labels 
ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend)) + geom_edges(aes(linetype = sign), color = "red",size=0.75) + theme_blank() + geom_nodes(aes(color = Phylum, size = vsize)) +geom_nodetext_repel (aes(label = Genus),fontface = "bold", box.padding = unit(1, "lines"), segment.size = 0, segment.colour = "gray",cex = 3) + guides( size = FALSE) + scale_linetype_manual(labels=c("Negative","Positive"), values=c(2,1)) + labs(linetype="Type of correlation") + scale_size_continuous(range=c(3,18))

#no labels
ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend)) + geom_edges(aes(linetype = sign), color = "red",size=0.75) + theme_blank() + geom_nodes(aes(color = Phylum, size = vsize))  + guides( size = FALSE) + scale_linetype_manual(labels=c("Negative","Positive"), values=c(2,1)) + labs(linetype="Type of correlation") + scale_size_continuous(range=c(3,18))

network.plot <- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Phylum", label="Genus") 
network.plot$layers <- network.plot$layers[-2]
network.plot <-  network.plot +
  geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
  theme(legend.title = element_text(size=10,face="bold"),
        legend.text = element_text(size=10))
network.plot

write.graph(ig.sparcc,file = "network_test.graphml",format = "graphml")







.



####Vertex (taxa or nodes) metrics####

alpha_centrality(ig.sparcc) #alpha_centrality calculates the alpha centrality of some (or all) vertices in a graph
authority_score(ig.sparcc)$vector #The authority scores of the vertices are defined as the principal eigenvector of AT A, where A is the adjacency matrix of the graph
closeness(ig.sparcc) #Closeness centrality measures how many steps is required to access every other vertex from a given vertex
coreness(ig.sparcc) #The k-core of graph is a maximal subgraph in which each vertex has at least degree k. The coreness of a vertex is k if it belongs to the k-core but not to the (k+1)-core
degree(ig.sparcc) #The degree of a vertex is its most basic structural property, the number of its adjacent edges
estimate_betweenness(ig.sparcc,cutoff = 0) #The vertex and edge betweenness are (roughly) defined by the number of geodesics (shortest paths) going through a vertex or an edge
knn(ig.sparcc)$knn #Calculate the average nearest neighbor degree of the given vertices and the same quantity in the function of vertex degree
strength(ig.sparcc) #Summing up the edge weights of the adjacent edges for each vertex

vertex.metrics <- data.frame(authority_score(ig.sparcc)$vector,closeness(ig.sparcc),coreness(ig.sparcc),degree(ig.sparcc),estimate_betweenness(ig.sparcc,cutoff = 0),knn(ig.sparcc)$knn,strength(ig.sparcc))
View(vertex.metrics)

####Graph metrics####
articulation_points(ig.sparcc) #Articulation points or cut vertices are vertices whose removal increases the number of connected components in a graph.
centralize() #Centralization is a method for creating a graph level centralization measure from the centrality scores of the vertices
cliques(ig.sparcc) #These functions find all, the largest or all the maximal cliques in an undirected graph. The size of the largest clique can also be calculated.
component_distribution(ig.sparcc) #Calculate the maximal (weakly or strongly) connected components of a graph
count_motifs(ig.sparcc) #Graph motifs are small connected subgraphs with a well-defined structure. These functions search a graph for various motifs
diameter(ig.sparcc) #The diameter of a graph is the length of the longest geodesic
edge_connectivity(ig.sparcc) #The edge connectivity of a graph or two vertices, this is recently also called group adhesion
edge_density(ig.sparcc) #The density of a graph is the ratio of the number of edges and the number of possible edges
girth(ig.sparcc) #The girth of a graph is the length of the shortest circle in it
gorder(ig.sparcc) #Order (number of vertices) of a graph
gsize(ig.sparcc) #The size of the graph (number of edges)
motifs(ig.sparcc) #Graph motifs are small connected subgraphs with a well-defined structure. These functions search a graph for various motifs, and returns a vector with different motif types
radius(ig.sparcc) #The eccentricity of a vertex is its shortest path distance from the farthest other node in the graph.The smallest eccentricity in a graph is called its radius
transitivity(ig.sparcc) #Transitivity measures the probability that the adjacent vertices of a vertex are connected. This is sometimes also called the clustering coefficient
vertex_connectivity(ig.sparcc) #The vertex connectivity of a graph or two vertices, this is recently also called group cohesion

####Null model generation####

hrg-methods() #Fitting and sampling hierarchical random graph models
permute() #Create a new graph, by permuting vertex ids
stochastic_matrix() #Retrieves the stochastic matrix of a graph of class igraph


####Graph manipulation####
as_edgelist() #Sometimes it is useful to work with a standard representation of a graph, like an edge list
as_adj_list() #Create adjacency lists from a graph, either for adjacent edges or for neighboring vertices
as_adjacency_matrix() #Sometimes it is useful to work with a standard representation of a graph, like an adjacency matrix
complementer () #A complementer graph contains all edges that were not present in the input graph
V() #Create a vertex sequence (vs) containing all vertices of a graph

####Graph comparison####
intersection() #The intersection of two or more graphs are created. The graphs may have identical or overlapping vertex sets



