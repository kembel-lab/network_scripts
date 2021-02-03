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

pdf("HTS_plots.pdf")

for (n in c(1:10)){
  
  print(paste("doing", n))
  ##rarefy our dataset
  raref.comm <- rrarefy(netw.nc.comm, 2000)
  raref.comm <- raref.comm[apply(raref.comm,1,sum)>= 2000,]
  
  raref.meta <- netw.nc.meta[rownames(raref.comm),]
  raref.taxo <- netw.nc.taxo[colnames(raref.comm),]
  
  #make a phyloseq object
  OTU <- otu_table(raref.comm, taxa_are_rows = FALSE)
  TAX <- tax_table(as.matrix(raref.taxo))
  samples <- sample_data(raref.meta)
  phylo <- phyloseq(OTU, TAX, samples)
  phylo
  assign(paste0("n",n,".phylo" ),phylo)
  
  #do on different host-trees
  host.meta <- sample_data(phylo)
  phylo.fir <- subset_samples(phylo, tree_species=="Balsam_fir")
  
  phylo.spruce <- subset_samples(phylo, tree_species=="Spruce")
  fir_sites <- sample_data(phylo.fir)$Site
  spruce_sites <- sample_data(phylo.spruce)$Site
  
  common_sites <- fir_sites[fir_sites%in%spruce_sites]
  phylo.spruce <- subset_samples(phylo.spruce, Site%in%common_sites)
  phylo.fir <- subset_samples(phylo.fir, Site%in%common_sites)
  
  #remove taxa that have less than two reads in that group
  phylo.spruce <- subset_taxa(phylo.spruce,apply(otu_table(phylo.spruce),2,sum) > 2)
  phylo.fir <- subset_taxa(phylo.fir,apply(otu_table(phylo.fir),2,sum) > 2)
  
  assign(paste0("n",n,".spruce.phylo" ),phylo.spruce)
  assign(paste0("n",n,".fir.phylo" ),phylo.fir)
  
  
  
  sparcc.boot.fir <- eco_network_1(phylo.fir,1000,8)
  sparcc.boot.spruce <- eco_network_1(phylo.spruce,1000,8)
  
  
  for (occ in c("fir","spruce")){
    print(paste0("doing ", n,occ))
    sparcc.boot <- get(paste0("sparcc.boot.",occ))
    phylo <- get(paste0("phylo.",occ))
    results.1 <- eco_network_2(sparcc.boot, phylo,0.01)
    sparcc.graph.filt <- results.1[[1]]
    phylo.filt <- results.1[[2]]
    phylo.filt
    ig.sparcc <- adj2igraph(sparcc.graph.filt,  vertex.attr=list(name=taxa_names(phylo.filt)))
    vsize    <- colMeans(clr(otu_table(phylo.filt),1))*10 + 5
    am.coord <- layout.fruchterman.reingold(ig.sparcc)
    par(mfrow=c(1,1))
    
    
    phylo.plot<- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Class",label="Genus", title=paste0("n",n,".",occ))
    #remove labels and re-do them
    phylo.plot$layers <- phylo.plot$layers[-2]
    phylo.plot <-  phylo.plot +
      geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
      theme(legend.title = element_text(size=10,face="bold"),
            legend.text = element_text(size=10))
    
    assign(paste0("n",n,".",occ,".phylo.plot"),phylo.plot)
    grid.arrange(phylo.plot)
    filename <- paste0("n",n,".",occ,".graphml")
    write.graph(ig.sparcc,file = filename,format = "graphml")
  }
}

#use this if the pdf did not work
pdf("HTS_plots.pdf")
for (n in c(1:10)){
  for (occ in c("fir","spruce")){
    phylo.plot <- get(paste0("n",n,".",occ,"phylo.plot"))
    grid.arrange(phylo.plot)

  }
}










