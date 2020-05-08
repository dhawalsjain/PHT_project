setwd("A:/work/scripts/Lucy_Parathyroid")
source("vars.R")

#################################################################################################33
############ dataframes
#################################################################################################
if(F){
  
  d <- read.delim(paste0(WORKDIR,"table.TPM.all.genes.8.PTG.txt"),comment.char = "#",stringsAsFactors = F,header = T)
  d$transcript_id.s. <- NULL
  d$Description <- NULL
  d <- d[,c("gene_id", "Symbol", "ETR862_TPM", "ETR868_TPM", "ETR892_TPM", "ETR902_TPM", 
            "ETR904_TPM", "ETR918_TPM", "ETR929_TPM", "PTH.3518_TPM")] %>% unique
  
  e <- read.delim(paste0(WORKDIR,"GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct"),
                  stringsAsFactors = F,header = F,skip = 1)
  cl <- unlist(e[2,])
  e <- e[3:nrow(e),]
  e <- e[e$V2%in%d$Symbol,]
  e <- e[!is.na(e$V2),]
  d <- d[match(e$V2,d$Symbol),]
  sum(e$V2==d$Symbol)
  e <- cbind(e,d[,3:10])
  rm(d)
  rownames(e) <- e$V1
  e$V1 <- e$V2 <- NULL
  rn <- rownames(e)
  e <- as.data.frame(sapply(e, as.numeric))
  cl <- c(cl,colnames(e)[11689:11696])
  cl <- cl[3:length(cl)]
  colnames(e) <- cl
  rownames(e) <- rn
  e <- e[rowSums(e)>0,]
  save(e,file=paste0(WORKDIR,"GTEX_counts.RData"))
  
}

#################################################################################################33
############ feature selection
#################################################################################################
 ## features ranking
if(F){
  load(paste0(WORKDIR,"GTEX_counts.RData"))
  
  rq <- data.frame(id=rownames(e),mean = rowMeans(e))
  rq$variance <- apply(e,1,var)
  rq$variance.expected <- 0
  rq$variance.standardized <- 0
  not.const <- rq$variance > 0
  fit <- loess(formula = log10(x = variance) ~ log10(x = mean), 
               data = rq[not.const, ], span = 0.3)
  rq$variance.expected[not.const] <- 10^fit$fitted
  z <- ((e - rq$mean)/sqrt(rq$variance.expected))
  rq$variance.standardized <- apply(z,1,var)
  rm(z,fit)
  rq <- rq[order(rq$variance.standardized,decreasing = T),]
  
  pdf(paste0(WORKDIR,"MeanVarRelationShip.pdf"),width = 6,height = 4)
  plot(log(rq$mean),log(rq$variance),col=rgb(0,0,0,0.2),xlab = "mean, log10", ylab = "var, log10",main="mean variance relationship")
  points(log(rq$mean),log(rq$variance.expected),col=rgb(1,0,0,0.2))
  points(log(rq$mean),log(rq$variance.standardized),col=rgb(0,0,1,0.2))
  legend("topleft", legend=c("observed", "expected","corrected"),
         col=c("black", "red","blue"), lty=1:2, cex=0.8)
  dev.off()
  
  
  cq <- rq[1:2000,]$id
  e <- e[rownames(e)%in%cq,]
  save(e,rq,file=paste0(WORKDIR,"SelectedFeatures.RData"))
  
}
## 2k and 5k features
if(F){
  load(paste0(WORKDIR,"SelectedFeatures.RData"))
  load(paste0(WORKDIR,"GTEX_counts.RData"))
  cq <- rq[1:2000,]
  e <- e[rownames(e)%in%cq$id,]
  save(e,cq,file=paste0(WORKDIR,"2k_selectedFeatures.RData"))
  
  
  load(paste0(WORKDIR,"SelectedFeatures.RData"))
  load(paste0(WORKDIR,"GTEX_counts.RData"))
  cq <- rq[1:5000,]
  e <- e[rownames(e)%in%cq$id,]
  save(e,cq,file=paste0(WORKDIR,"5k_selectedFeatures.RData"))
  
  load(paste0(WORKDIR,"SelectedFeatures.RData"))
  load(paste0(WORKDIR,"GTEX_counts.RData"))
  cq <- rq[1:10000,]
  e <- e[rownames(e)%in%cq$id,]
  save(e,cq,file=paste0(WORKDIR,"10k_selectedFeatures.RData"))
  
}


#################################################################################################33
############ t-SNE on selected features
#################################################################################################
## 5k rowMeanScaled features (perplexity 30, default)
if(T){
  source("vars.R")
  
  load(paste0(WORKDIR,"5k_rowMeanScaledCounts.RData"))
  res <- Rtsne(t(f), dims = 2, initial_dims = 50,
               perplexity = 30, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  
  pl <- myplotdf(res,f)
  
  pdf(paste0(WORKDIR,"5kFeatures_rowMeanScaledCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_rowMeanScaledCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
}

## 5k rowMeanScaled features (perplexity 5)
if(T){
  source("vars.R")
  load(paste0(WORKDIR,"5k_rowMeanScaledCounts.RData"))
  res <- Rtsne(t(f), dims = 2, initial_dims = 50,
               perplexity = 5, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  
  pl <- myplotdf(res,f)
  
  pdf(paste0(WORKDIR,"5kFeatures_p5rowMeanScaledCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_p5rowMeanScaledCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
}

## 5k rowMeanScaled features (perplexity 10)
if(T){
  source("vars.R")
  load(paste0(WORKDIR,"5k_rowMeanScaledCounts.RData"))
  res <- Rtsne(t(f), dims = 2, initial_dims = 50,
               perplexity = 10, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  
  pl <- myplotdf(res,f)
  
  pdf(paste0(WORKDIR,"5kFeatures_p10rowMeanScaledCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_p10rowMeanScaledCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
}

## 5k rowMeanScaled features (perplexity 20)
if(T){
  source("vars.R")
  load(paste0(WORKDIR,"5k_rowMeanScaledCounts.RData"))
  res <- Rtsne(t(f), dims = 2, initial_dims = 50,
               perplexity = 20, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  
  pl <- myplotdf(res,f)
  
  pdf(paste0(WORKDIR,"5kFeatures_p20rowMeanScaledCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_p20rowMeanScaledCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
}

## 5k rowMeanScaled features (perplexity 40)
if(T){
  source("vars.R")
  load(paste0(WORKDIR,"5k_rowMeanScaledCounts.RData"))
  res <- Rtsne(t(f), dims = 2, initial_dims = 50,
               perplexity = 40, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  
  pl <- myplotdf(res,f)
  
  pdf(paste0(WORKDIR,"5kFeatures_p40rowMeanScaledCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_p40rowMeanScaledCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
}


## 5k Zscores
if(T){
  source("vars.R")
  load(paste0(WORKDIR,"5k_selectedFeatures.RData"))
  
  e <- apply(e,1,function(x){
    round((x -mean(x))/sd(x),2)
  }) %>% t %>% as.data.frame
  
  res <- Rtsne(t(e), dims = 2, initial_dims = 50,
               perplexity = 30, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  pl <- myplotdf(res,e)
  
  pdf(paste0(WORKDIR,"5kFeatures_ZscoreCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_ZscoreCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
}

## 5k standardized Zscores
if(T){
  source("vars.R")
  load(paste0(WORKDIR,"5k_selectedFeatures.RData"))
  
  cq <- cq[match(rownames(e),cq$id),]
  g <- round((e - cq$mean)/cq$variance.standardized,2)
  g <- t(g) %>% as.data.frame
  
  res <- Rtsne(g, dims = 2, initial_dims = 50,
               perplexity = 30, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  pl <- myplotdf(res,t(g))
  
  pdf(paste0(WORKDIR,"5kFeatures_stdZscoreCounts_SMTS.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTS),size=3, data = lc.cent) 
  dev.off()
  
  pdf(paste0(WORKDIR,"5kFeatures_stdZscoreCounts_SMTSD.pdf"),width = 12,height = 10)
  lc.cent = pl %>% group_by(louvain, SMTSD) %>% select(V1,V2) %>% summarize_all(mean)
  ggplot(pl, aes(x = V1, y = V2, colour = louvain)) + 
    geom_point(alpha = 0.8) + theme_bw()+ guides(colour = FALSE)  +
    geom_label_repel(aes(label = SMTSD),size=3, data = lc.cent) 
  dev.off()
  
  
}


#################################################################################################33
############ SCVIS write down the data 
#################################################################################################
if(F){
  source("vars.R")
  load(file=paste0(WORKDIR,"5k_rowMeanScaledCounts.RData"))
  f <- t(f) %>% as.data.frame
  write.table(f[1:11688,],file=paste0(WORKDIR,"GTEX_t5kRowMeanScaled.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  write.table(f[11689:11696,],file=paste0(WORKDIR,"parathyroid_t5kRowMeanScaled.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  
  source("vars.R")
  load(paste0(WORKDIR,"GTEX_counts.RData"))
  e <- log2(e+1)
  e <- apply(e,1,function(x){
    round((x-mean(x))/sd(x),2)
  }) %>% t %>% as.data.frame
  write.table(e[1:11688,],file=paste0(WORKDIR,"GTEX_t.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  write.table(e[11689:11696,],file=paste0(WORKDIR,"parathyroid_t.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  
  
  source("vars.R")
  load(paste0(WORKDIR,"5k_selectedFeatures.RData"))
  e <- apply(e,1,function(x){
    round((x -mean(x))/sd(x),2)
  }) %>% t %>% as.data.frame
  e <- t(e)
  e <- as.data.frame(e)
  write.table(e[1:11688,],file=paste0(WORKDIR,"GTEX_t5kZscore.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  write.table(e[11689:11696,],file=paste0(WORKDIR,"parathyroid_t5kZscore.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  
  
  source("vars.R")
  load(paste0(WORKDIR,"10k_selectedFeatures.RData"))
  e <- apply(e,1,function(x){
    round((x -mean(x))/sd(x),2)
  }) %>% t %>% as.data.frame
  e <- t(e)
  e <- as.data.frame(e)
  write.table(e[1:11688,],file=paste0(WORKDIR,"GTEX_t10kZscore.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  write.table(e[11689:11696,],file=paste0(WORKDIR,"parathyroid_t10kZscore.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
  
   
}
 ## 






#################################################################################################33
############ Trash 
#################################################################################################
## 3. perform t-sne (2000 features)
if(F){
  library(Rtsne)
  library(fpc)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  load(paste0(WORKDIR,"SelectedFeatures.RData"))
  e <- unique(e)
  f <- apply(e,2,function(x){
    round((x-mean(x))/sd(x),2)
  }) %>% as.data.frame
  
  rq1 <- rq[match(rownames(e),rq$id),]
  g <- round((e - rq1$mean)/rq1$variance.standardized,2)
  
  
  res <- Rtsne(t(g), dims = 2, initial_dims = 50,
               perplexity = 30, theta = 0.5, check_duplicates = TRUE,
               pca = TRUE, partial_pca = FALSE, max_iter = 1000)
  
  ds = dbscan(res$Y, 2)
  z <- read.delim(paste0(WORKDIR,"GTEx_v7_Annotations_SampleAttributesDS.txt"),header=T,stringsAsFactors = F)
  z <- z[,c("SAMPID","SMTS")] %>% unique
  pl <- res$Y %>% as.data.frame
  pl$SAMPID <- colnames(e)
  pl$density = ds$cluster
  
  pl <- merge(pl,z,by="SAMPID",all.x=T)
  pl$SMTS <- ifelse(is.na(pl$SMTS),"Parathyroid",as.character(pl$SMTS))
  
  ds.cent = pl %>% group_by(SMTS) %>% select(V1,V2) %>% summarize_all(mean)
  
  ggplot(pl, aes(x = V1, y = V2, colour = SMTS)) + 
    geom_point(alpha = 0.3) + theme_bw()  + guides(colour = FALSE)+
    geom_label_repel(aes(label = SMTS), data = ds.cent) 
  
  
  
  pl.tpm <- pl
  pl.zscore <- pl  
  pl.stdvar <- pl
  save(pl.tpm,pl.zscore,pl.stdvar,file=paste0(WORKDIR,"tSNE_tests.RData"))
  
  
}


## 
