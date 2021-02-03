eco_network_1<-function(phylo,bootstraps,ncpus){

#bootstrap for p-values
tp0 <- proc.time()
sparcc.boot <- sparccboot(otu_table(phylo), R=bootstraps, ncpus = ncpus)
tp1 <- proc.time()
print(paste("It took ", (tp1 - tp0)["elapsed"][[1]], " seconds for ", bootstraps, " bootstraps to perform with the help of ", ncpus, " processors",sep=""))


return(sparcc.boot)

}
