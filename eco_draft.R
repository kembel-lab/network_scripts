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

##INPUTS
meta <- read.table(file="organized2018_extraction_metadata.txt", sep="\t", header=TRUE, row.names="ID",stringsAsFactors = FALSE)
taxo <- read.csv(file="taxonomy_DADA2.csv", header=TRUE, row.names=1)
comm <- read.csv(file="community_DADA2.csv", header=TRUE, row.names=1)

##PREPROCESSING OF THE TABLES
par(mfrow=c(1,1))
all.comm <- comm
all.meta <- meta
all.taxo <- taxo

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(all.taxo)){
  levels(all.taxo[,i])<-c(levels(all.taxo[,i]), "unassigned", NA)
  for(j in rownames(all.taxo)){
    if(is.na(all.taxo[j,i])){
      all.taxo[j,i] <- "unassigned"
    }
  }
}
all.trash <- all.taxo[c(all.taxo$Kingdom == "unassigned" |all.taxo$Kingdom == "Eukaryota" | all.taxo$Family == "Mitochondria"| all.taxo$Class == "Chloroplast"  ),]
all.taxo <- subset(all.taxo, !(rownames(all.taxo) %in% rownames(all.trash)))

#match the metadata and community tables
all.comm <- all.comm[rownames(all.meta),]
all.comm <- all.comm[,rownames(all.taxo)]

#remove ASVs with less than 0,1% reads
all.comm <- all.comm[,apply(all.comm,2,sum)> (sum(apply(all.comm,2,sum))/1000)]
all.taxo <- all.taxo[colnames(all.comm),]

#change the names of our ASVs to something more practical
n <- 0
sequences <- character()
for(i in colnames(all.comm)){
  n <- n +1
  sequences[n] <- colnames(all.comm)[n]
  colnames(all.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(all.comm)[n]
}
rm(n,i,j)
rownames(all.taxo) <- colnames(all.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(all.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(all.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(all.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation_level"){
    colnames(all.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#remove the controls
all.meta <- all.meta[!is.na(all.meta$sample_type),]

#add more metadata to our sites
sitesSpecs <- data.frame(read.csv(file="site_specs_2.csv", header=TRUE, row.names="ID1", stringsAsFactors = FALSE))
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(all.meta)){
  for (siteName in rownames(sitesSpecs)){
    if (siteName == all.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,sitesSpecs[siteName,c(5:11)])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
all.meta <- cbind(all.meta[rownames(tempMeta),], tempMeta)
rm(n,sitesSpecs,tempMeta,all.trash, colname, rowname, siteName)
all.meta$site_specie <- paste(all.meta$Site, all.meta$tree_species)

#make sure all the numbers are numeric
list_of_numerics <- c("host_tree_damage","Latitude","Longitude","Elevation","Mean_temperature","Total_precipitations","Solar_radiation","Reference_evaporation","Climatic_moisture_deficit","Relative_humidity")
for (number in list_of_numerics){
  all.meta[,number] <- as.numeric(all.meta[,number])
}
rm(number,list_of_numerics)

#keep only larvae with all the metadata
all.meta <- all.meta[all.meta$sample_type == "L6"&!is.na(all.meta$host_tree_damage),]
all.comm <- all.comm[rownames(all.meta),]








#try occurence filtration (transform into presence/absence then look at distribution of columns totals) before rarefaction
# all.pa <- decostand(all.comm,method = "pa")
# hist(apply(all.pa,2,sum))
# apply(all.pa,2,sum) > 5

#Do first steps for 5 occ

#try occurence in 5 samples
all.5.comm <- all.comm[,apply(all.pa,2,sum) > 5]
##rarefy our dataset
quantile(apply(all.5.comm, 1,sum))
all.5.ab.comm <- rrarefy(all.5.comm, 4000)
all.5.ab.comm <- all.5.ab.comm[apply(all.5.ab.comm,1,sum)>= 4000,]
rarecurve(all.5.ab.comm, step=100, label=FALSE)
all.5.ab.meta <- all.meta[rownames(all.5.ab.comm),]
all.5.ab.taxo <- all.taxo[colnames(all.5.ab.comm),]


#make a phyloseq object
OTU <- otu_table(all.5.ab.comm, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(all.5.ab.taxo))
samples <- sample_data(all.5.ab.meta)
phylo.5 <- phyloseq(OTU, TAX, samples)
phylo.5

sparcc.boot.5 <- eco_network_1(phylo.5,1000,4)
#and export it

for (occ in c(2,5,10,15,20,25,30,35,40)){
  print(paste("doing", occ))
  
  all.occ.comm <- all.comm[,apply(all.pa,2,sum) > occ]
  ##rarefy our dataset
  quantile(apply(all.occ.comm, 1,sum))
  all.occ.ab.comm <- rrarefy(all.occ.comm, 4000)
  all.occ.ab.comm <- all.occ.ab.comm[apply(all.occ.ab.comm,1,sum)>= 4000,]
  rarecurve(all.occ.ab.comm, step=100, label=FALSE)
  all.occ.ab.meta <- all.meta[rownames(all.occ.ab.comm),]
  all.occ.ab.taxo <- all.taxo[colnames(all.occ.ab.comm),]
        
        
  #make a phyloseq object
  OTU <- otu_table(all.occ.ab.comm, taxa_are_rows = FALSE)
  TAX <- tax_table(as.matrix(all.occ.ab.taxo))
  samples <- sample_data(all.occ.ab.meta)
  phylo.occ <- phyloseq(OTU, TAX, samples)
  print(phylo.occ)
   
  sparcc.boot.occ <- eco_network_1(phylo.occ,1000,4)
  assign(paste0("phylo.", occ),phylo.occ)
  assign(paste0("sparcc.boot.", occ),sparcc.boot.occ)
  
}
par(mfrow=c(3,3))
pdf("occurence_plots.pdf")
for (occ in c(2,5,10,15,20,25,30,35,40)){
  print(paste("doing", occ))
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
  
  
  phylo.plot<- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Class",label="Genus", title=occ)
  #remove labels and re-do them
  phylo.plot$layers <- phylo.plot$layers[-2]
  phylo.plot <-  phylo.plot +
    geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
    theme(legend.title = element_text(size=10,face="bold"),
          legend.text = element_text(size=10))
  
  assign(paste0("phylo.plot.", occ),phylo.plot)
  grid.arrange(get(paste0("phylo.plot.", occ)))
  filename <- paste0(occ,".graphml")
  write.graph(ig.sparcc,file = filename,format = "graphml")
}
  #function to graph
dev.off()








#do on different host-trees
host.meta <- sample_data(phylo.5)
phylo.fir <- subset_samples(phylo.5, tree_species=="Balsam_fir")

phylo.spruce <- subset_samples(phylo.5, tree_species=="Spruce")
fir_sites <- sample_data(phylo.fir)$Site
spruce_sites <- sample_data(phylo.spruce)$Site

common_sites <- fir_sites[fir_sites%in%spruce_sites]
phylo.spruce <- subset_samples(phylo.spruce, Site%in%common_sites)
phylo.fir <- subset_samples(phylo.fir, Site%in%common_sites)

sparcc.boot.fir <- eco_network_1(phylo.fir,1000,4)
sparcc.boot.spruce <- eco_network_1(phylo.spruce,1000,4)


pdf("HTS_plots.pdf")
for (occ in c("fir","spruce")){
  print(paste("doing", occ))
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
  
  
  phylo.plot<- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Class",label="Genus", title=occ)
  #remove labels and re-do them
  phylo.plot$layers <- phylo.plot$layers[-2]
  phylo.plot <-  phylo.plot +
    geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
    theme(legend.title = element_text(size=10,face="bold"),
          legend.text = element_text(size=10))
  
  assign(paste0("phylo.plot.", occ),phylo.plot)
  grid.arrange(get(paste0("phylo.plot.", occ)))
  filename <- paste0(occ,".graphml")
  write.graph(ig.sparcc,file = filename,format = "graphml")
}
dev.off()

#look at anuran results
anuran.HST <- read.graph("demo_HTS_1_intersection.graphml",format = "graphml" )
plot(anuran.HST)














#do guts vs whole insects
##INPUTS
comm <- read.table(file="community_table_ch2.csv", sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE)
taxo <- read.table(file="taxonomy_table_ch2.csv", sep=",", header=TRUE, row.names=1,stringsAsFactors = FALSE)
meta <- read.table(file="metadata_ch2.csv", sep=",", row.names=1, header=TRUE,stringsAsFactors = FALSE)

##PREPROCESSING OF THE TABLES
par(mfrow=c(1,1))
sexing.comm <- as.data.frame(comm)
sexing.meta <- as.data.frame(meta)
sexing.taxo <- as.data.frame(taxo)

sexing.meta$SexFact <- as.factor(sexing.meta$Sex)
sexing.meta$SiteFact <- as.factor(sexing.meta$Site)

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(sexing.taxo)){
  levels(sexing.taxo[,i])<-c(levels(sexing.taxo[,i]), "unassigned", NA)
  for(j in rownames(sexing.taxo)){
    if(is.na(sexing.taxo[j,i])){
      sexing.taxo[j,i] <- "unassigned"
    }
  }
}
sexing.trash <- sexing.taxo[c(sexing.taxo$Kingdom == "unassigned" |sexing.taxo$Kingdom == "Eukaryota" | sexing.taxo$Family == "Mitochondria"| sexing.taxo$Class == "Chloroplast"  ),]
sexing.taxo <- subset(sexing.taxo, !(rownames(sexing.taxo) %in% rownames(sexing.trash)))

#match the metadata and community tables
sexing.comm <- sexing.comm[rownames(sexing.meta),]
sexing.comm <- sexing.comm[,rownames(sexing.taxo)]

#remove ASVs with less than 0,1% reads
sexing.comm <- sexing.comm[,apply(sexing.comm,2,sum)> (sum(apply(sexing.comm,2,sum))/10000)]
sexing.taxo <- sexing.taxo[colnames(sexing.comm),]

#change the names of our ASVs to something more practical
n <- 0
sequences <- character()
for(i in colnames(sexing.comm)){
  n <- n +1
  sequences[n] <- colnames(sexing.comm)[n]
  colnames(sexing.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(sexing.comm)[n]
}
rm(n,i,j)
rownames(sexing.taxo) <- colnames(sexing.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(sexing.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(sexing.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(sexing.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation"){
    colnames(sexing.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#remove the controls
#sexing.meta <- sexing.meta[sexing.meta$sample_type!="control",]
#try decontam
library(ggplot2)
library(phyloseq)
library(decontam)

#remove the positive controls
sexing.meta <- sexing.meta[sexing.meta$sample_type!="poscontrol"&sexing.meta$sample_type!="pcrnegcontrol",]
sexing.comm <- sexing.comm[rownames(sexing.meta),]

ps<-phyloseq(otu_table(sexing.comm,taxa_are_rows = FALSE), tax_table(as.matrix(sexing.taxo)), sample_data(sexing.meta))
sample_data(ps)$is.neg <- sample_data(ps)$Site == "negcontrol"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

notcontam <- rownames(contamdf.prev05)[contamdf.prev05$contaminant==FALSE]
sexing.nc.comm <- sexing.comm[,notcontam]
sexing.nc.taxo <- sexing.taxo[notcontam,]
sexing.nc.meta <- sexing.meta

#do decontam again for only N3
ps<-phyloseq(otu_table(sexing.nc.comm,taxa_are_rows = FALSE), tax_table(as.matrix(sexing.nc.taxo)), sample_data(sexing.nc.meta))
sample_data(ps)$is.neg <- rownames(sample_data(ps)) == "N3"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)

notcontam <- rownames(contamdf.prev05)[contamdf.prev05$contaminant==FALSE]
sexing.nc.comm <- sexing.nc.comm[,notcontam]
sexing.nc.taxo <- sexing.nc.taxo[notcontam,]
sexing.nc.meta <- sexing.nc.meta

#remove negative controls
sexing.nc.meta <- sexing.nc.meta[sexing.meta$Site!="negcontrol",]
sexing.nc.comm <- sexing.nc.comm[rownames(sexing.nc.meta),]

sexing.pa <- decostand(sexing.nc.comm,method = "pa")
hist(apply(sexing.pa,2,sum))
# apply(all.pa,2,sum) > 5

#try occurence in 5
all.gut.comm <- sexing.nc.comm[,apply(sexing.pa,2,sum) > 5]
##rarefy our dataset
quantile(apply(all.gut.comm, 1,sum))
all.gut.ab.comm <- rrarefy(all.gut.comm, 2000)
all.gut.ab.comm <- all.gut.ab.comm[apply(all.gut.ab.comm,1,sum)>= 2000,]
rarecurve(all.gut.ab.comm, step=100, label=FALSE)
all.gut.ab.meta <- sexing.nc.meta[rownames(all.gut.ab.comm),]
all.gut.ab.taxo <- sexing.nc.taxo[colnames(all.gut.ab.comm),]


#make a phyloseq object
OTU <- otu_table(all.gut.ab.comm, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(all.gut.ab.taxo))
samples <- sample_data(all.gut.ab.meta)
phylo.gut <- phyloseq(OTU, TAX, samples)
phylo.gut
phylo.whole<-phylo.5 


#try setting a correlation threshold also
### work on figure text
sparcc.boot.whole <- eco_network_1(phylo.5,1000,4)
sparcc.boot.gut <- eco_network_1(phylo.gut,1000,4)


pdf("insect_plots.pdf")
for (occ in c("gut","whole")){
  print(paste("doing", occ))
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
  
  
  phylo.plot<- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Class",label="Genus", title=occ)
  #remove labels and re-do them
  phylo.plot$layers <- phylo.plot$layers[-2]
  phylo.plot <-  phylo.plot +
    geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
    theme(legend.title = element_text(size=10,face="bold"),
          legend.text = element_text(size=10))
  
  assign(paste0("phylo.plot.", occ),phylo.plot)
  grid.arrange(get(paste0("phylo.plot.", occ)))
  filename <- paste0(occ,".graphml")
  write.graph(ig.sparcc,file = filename,format = "graphml")
}
#function to graph
dev.off()


#check sex difference
phylo.M <- subset_samples(phylo.gut, Sex=="M")
phylo.M
phylo.F <- subset_samples(phylo.gut, Sex=="F")
phylo.F

sparcc.boot.F <- eco_network_1(phylo.F,1000,4)
sparcc.boot.M <- eco_network_1(phylo.M,1000,4)

pdf("insect_sex.pdf")
for (occ in c("F","M")){
  print(paste("doing", occ))
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
  
  
  phylo.plot<- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Class",label="Genus", title=occ)
  #remove labels and re-do them
  phylo.plot$layers <- phylo.plot$layers[-2]
  phylo.plot <-  phylo.plot +
    geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
    theme(legend.title = element_text(size=10,face="bold"),
          legend.text = element_text(size=10))
  
  assign(paste0("phylo.plot.", occ),phylo.plot)
  grid.arrange(get(paste0("phylo.plot.", occ)))
  filename <- paste0(occ,".graphml")
  write.graph(ig.sparcc,file = filename,format = "graphml")
}
#function to graph
dev.off()










#do larvae vs foliage
###LARVAE VS FOLIAGE
##INPUTS
meta <- read.table(file="organized2018_extraction_metadata.txt", sep="\t", header=TRUE, row.names="ID",stringsAsFactors = FALSE)
taxo <- read.csv(file="taxonomy_DADA2.csv", header=TRUE, row.names=1)
comm <- read.csv(file="community_DADA2.csv", header=TRUE, row.names=1)

##PREPROCESSING OF THE TABLES
compare.comm <- comm
compare.meta <- meta
compare.taxo <- taxo

#replace all NAs by unassigned if we want the following part to work
for(i in colnames(compare.taxo)){
  levels(compare.taxo[,i])<-c(levels(compare.taxo[,i]), "unassigned", NA)
  for(j in rownames(compare.taxo)){
    if(is.na(compare.taxo[j,i])){
      compare.taxo[j,i] <- "unassigned"
    }
  }
}
compare.trash <- compare.taxo[c(compare.taxo$Kingdom == "unassigned" |compare.taxo$Kingdom == "Eukaryota" | compare.taxo$Family == "Mitochondria"| compare.taxo$Class == "Chloroplast"  ),]
compare.taxo <- subset(compare.taxo, !(rownames(compare.taxo) %in% rownames(compare.trash)))

#match the metadata and community tables
compare.comm <- compare.comm[rownames(compare.meta),]
compare.comm <- compare.comm[,rownames(compare.taxo)]

#remove ASVs with less than 0,1% reads
compare.comm <- compare.comm[,apply(compare.comm,2,sum)> (sum(apply(compare.comm,2,sum))/1000)]
compare.taxo <- compare.taxo[colnames(compare.comm),]

#change the names of our ASVs to something more practical
n <- 0
compare.sequences <- character()
for(i in colnames(compare.comm)){
  n <- n +1
  sequences[n] <- colnames(compare.comm)[n]
  colnames(compare.comm)[n] <- paste("ASV",n, sep="")
  names(sequences)[n] <- colnames(compare.comm)[n]
}
rm(n,i,j,compare.trash)
rownames(compare.taxo) <- colnames(compare.comm)

#change headers in the metadata
n <- 0
for(colname in colnames(compare.meta)){
  n <- n + 1
  if(colname == "Dev._Stage"){
    colnames(compare.meta)[n] <-"sample_type"
  }
  else if(colname == "Tree_spec."){
    colnames(compare.meta)[n] <-"tree_species"
  }
  if(colname == "Defoliation_level"){
    colnames(compare.meta)[n] <-"host_tree_damage"
  }
}
rm(n)

#remove the controls
compare.meta <- compare.meta[!is.na(compare.meta$sample_type),]
compare.meta <- compare.meta[!is.na(compare.meta$host_tree_damage),]
compare.comm <- compare.comm[rownames(compare.meta),]

#add more metadata to our sites
sitesSpecs <- data.frame(read.csv(file="site_specs_2.csv", header=TRUE, row.names="ID1", stringsAsFactors = FALSE))
tempMeta <- data.frame()
n <- 0
for (rowname in rownames(compare.meta)){
  for (siteName in rownames(sitesSpecs)){
    if (siteName == compare.meta[rowname,"Site"]){
      tempMeta <- rbind(tempMeta,sitesSpecs[siteName,c(5:11)])
      n <- n + 1
      rownames(tempMeta)[n] <- rowname
    }
  }
}
compare.meta <- cbind(compare.meta[rownames(tempMeta),], tempMeta)
rm(n,sitesSpecs,tempMeta,compare.trash, colname, rowname, siteName)
compare.meta$site_specie <- paste(compare.meta$Site, compare.meta$tree_species)








#function to identify samples for which foliage and larvae are present
pairs <- table(compare.meta$Tree_ID) == 2
compareList <- c()
for (element in names(pairs)){
  if (pairs[element] == TRUE){
    print(element)
    compareList <- c(compareList,element) 
  }
}
rm(element)
rm(pairs)
compareListOut <- c()
for (tree_ID in compareList){
  for (rowname in rownames(compare.comm)){
    if (tree_ID == compare.meta[rowname,"Tree_ID"]){
      compareListOut <- c(compareListOut,rowname)
    }
  }
}
rm(rowname,tree_ID,compareList)

#make tables
compare.comm <- compare.comm[compareListOut,]
compare.meta <- compare.meta[rownames(compare.comm),]
compare.taxo <- compare.taxo[colnames(compare.comm),]

#try occurence in 5
comp.pa <- decostand(compare.comm,method = "pa")
hist(apply(comp.pa,2,sum))
all.comp.comm <- compare.comm[,apply(comp.pa,2,sum) > 5]
##rarefy our dataset
quantile(apply(all.comp.comm, 1,sum))
all.comp.ab.comm <- rrarefy(all.comp.comm, 2000)
all.comp.ab.comm <- all.comp.ab.comm[apply(all.comp.ab.comm,1,sum)>= 2000,]
rarecurve(all.comp.ab.comm, step=100, label=FALSE)
all.comp.ab.meta <- compare.meta[rownames(all.comp.ab.comm),]
all.comp.ab.taxo <- compare.taxo[colnames(all.comp.ab.comm),]


#make a phyloseq object
OTU <- otu_table(all.comp.ab.comm, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(all.comp.ab.taxo))
samples <- sample_data(all.comp.ab.meta)
phylo.comp <- phyloseq(OTU, TAX, samples)
phylo.comp


phylo.larv <- subset_samples(phylo.comp,sample_type=="L6")
phylo.larv
phylo.foliage <- subset_samples(phylo.comp,sample_type=="Foliage")
phylo.foliage

#try setting a correlation threshold also
### work on figure text
sparcc.boot.larv <- eco_network_1(phylo.larv,1000,4)
sparcc.boot.foliage <- eco_network_1(phylo.foliage,1000,4)


pdf("diet_plots.pdf")
for (occ in c("larv","foliage")){
  print(paste("doing", occ))
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
  
  
  phylo.plot<- plot_network(ig.sparcc, point_size = vsize, phylo.filt, type='taxa', color="Class",label="Genus", title=occ)
  #remove labels and re-do them
  phylo.plot$layers <- phylo.plot$layers[-2]
  phylo.plot <-  phylo.plot +
    geom_text_repel(aes(label=Genus), size = 3, color = "dimgray")+
    theme(legend.title = element_text(size=10,face="bold"),
          legend.text = element_text(size=10))
  
  assign(paste0("phylo.plot.", occ),phylo.plot)
  grid.arrange(get(paste0("phylo.plot.", occ)))
  filename <- paste0(occ,".graphml")
  write.graph(ig.sparcc,file = filename,format = "graphml")
}
#function to graph
dev.off()




