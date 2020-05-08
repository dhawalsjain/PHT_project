
rm(list=ls())
WORKDIR="A:/work/Lucy_Parathyroid/"


library(Rtsne)
library(fpc)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(igraph)
library(FNN)
myplotdf <- function(res,e){
  
  z <- read.delim(paste0(WORKDIR,"GTEx_v7_Annotations_SampleAttributesDS.txt"),header=T,stringsAsFactors = F)
  z <- z[,c("SAMPID","SMTS","SMTSD")] %>% unique
  pl <- res$Y %>% as.data.frame
  pl$SAMPID <- rownames(e)
  k = 100
  knn.real = get.knn(as.matrix(res$Y), k = k)
  knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index),k), 
                        to = as.vector(knn.real$nn.index), weight = 1/(1 +as.vector(knn.real$nn.dist)))
  nw.real = graph_from_data_frame(knn.real, directed = FALSE)
  nw.real = simplify(nw.real)
  lc.real = cluster_louvain(nw.real)
  pl$louvain = as.factor(membership(lc.real))
  pl <- merge(pl,z,by="SAMPID",all.x=T)
  pl$SMTS <- ifelse(is.na(pl$SMTS),"Parathyroid",as.character(pl$SMTS))
  pl$SMTSD <- ifelse(is.na(pl$SMTSD),"Parathyroid",as.character(pl$SMTSD))
  return(pl)
  
}


