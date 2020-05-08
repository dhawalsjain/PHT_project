setwd("A:/work/scripts/Lucy_Parathyroid")
source("vars.R")

if(F){
  load(paste0(WORKDIR,"5k_selectedFeatures.RData"))
  labs <- data.frame(SAMPID=colnames(e))
  write.table(labs, file=paste0(WORKDIR,"scvis_5krowMeanScaled/labs.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(labs, file=paste0(WORKDIR,"scvis_5krowZscaled/labs.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
  write.table(labs, file=paste0(WORKDIR,"scvis_10krowZscaled/labs.txt"),quote = F,sep = "\t",row.names = F,col.names = F)
}

## 30 distinct colors
if(F){
library(randomcoloR)
   cols <- unique(pl$SMTS)
   cols <- distinctColorPalette(k = length(cols), altCol = FALSE, runTsne = T)
   names(cols) <- unique(pl$SMTS)
}
cols <- c(`Adipose Tissue` = "#E3E7A7", Muscle = "#C0B68C", `Blood Vessel` = "#65D9D7", 
          Heart = "#D8B589", Ovary = "#8475DC", Uterus = "#D050D7", Vagina = "#78EBA2", 
          Breast = "#DEC2ED", Skin = "#CFDFEA", `Salivary Gland` = "#8585B4", 
          Brain = "#87C8D4", `Adrenal Gland` = "#F8F0D5", Thyroid = "#D8E141", 
          Lung = "#A27867", Spleen = "#A8D074", Pancreas = "#64A9D4", Esophagus = "#E69556", 
          Stomach = "#DA7ED1", Colon = "#F4D5DD", `Small Intestine` = "#E5C555", 
          Prostate = "#56DABD", Testis = "#9ECAA0", Nerve = "#A2F975", 
          Blood = "#8450F0", Pituitary = "#DA8BA7", Liver = "#B9E0D4", 
          Kidney = "#F26063", `Fallopian Tube` = "#F7F08D", Bladder = "#CC5497", 
          `Cervix Uteri` = "#CC95D8")


gg_aes <- theme(axis.text = element_text(size=15,color="black"),
                axis.title = element_text(size=15,color="black"),
                legend.text = element_text(size=15,colour = "black"),
                strip.text = element_blank())

myf <- function(p0,p1,p2,p3,p4){
  labs <- read.delim(paste0(WORKDIR,p0),header=F,stringsAsFactors = F)
  z <- read.delim(paste0(WORKDIR,"GTEx_v7_Annotations_SampleAttributesDS.txt"),header=T,stringsAsFactors = F)
  z <- z[,c("SAMPID","SMTS","SMTSD")] %>% unique
  gtex <- read.delim(paste0(WORKDIR,p1),header = T,stringsAsFactors = F)
  pth <- read.delim(paste0(WORKDIR,p2),header = T,stringsAsFactors = F)
  gtex1 <- read.delim(paste0(WORKDIR,p3),header = T,stringsAsFactors = F)
  pth1 <- read.delim(paste0(WORKDIR,p4),header = T,stringsAsFactors = F)
  gtex$log_likelihood <- gtex1$log_likelihood
  pth$log_likelihood <- pth1$log_likelihood
  rm(gtex1,pth1)
  
  pl <- data.frame(SAMPID=labs$V1[1:nrow(gtex)])
  pl <- cbind(pl,gtex[,2:4])
  names(pl)[2:3] <- c("V1","V2")
  gl <- data.frame(SAMPID=labs$V1[11689:11696])
  gl <- cbind(gl,pth[,2:4])
  names(gl)[2:3] <- c("V1","V2")
  pl$data <- "GTEx"
  gl$data <- "PTH"
  pl <- rbind(pl,gl)
  rm(gl,pth)
  
  k = 100
  knn.real = get.knn(as.matrix(pl[,2:3]), k = k)
  knn.real = data.frame(from = rep(1:nrow(knn.real$nn.index),k), 
                        to = as.vector(knn.real$nn.index), weight = 1/(1 +as.vector(knn.real$nn.dist)))
  nw.real = graph_from_data_frame(knn.real, directed = FALSE)
  nw.real = simplify(nw.real)
  lc.real = cluster_louvain(nw.real)
  pl$louvain = as.factor(membership(lc.real))
  pl <- merge(pl,z,by="SAMPID",all.x=T)
  pl$SMTS <- ifelse(is.na(pl$SMTS),"Parathyroid",as.character(pl$SMTS))
  pl$SMTSD <- ifelse(is.na(pl$SMTSD),"Parathyroid",as.character(pl$SMTSD))
  rm(knn.real,lc.real,nw.real)
  pl
}

###############################################################################################3
## 1. rowMeanScaled 

if(F){
  
  
 pl <- myf(p0 = "scvis_5krowMeanScaled/labs.txt",
      p1="scvis_5krowMeanScaled/gtex.tsv",
     p2="scvis_5krowMeanScaled/pth.tsv",
     p3="scvis_5krowMeanScaled/gtex_log_likelihood.tsv",
     p4="scvis_5krowMeanScaled/pth_log_likelihood.tsv")
 
  xrange <- range(pl$V1)
  yrange <- range(pl$V2)
  
  lc = pl %>% group_by(SMTS,louvain) %>% select(V1,V2) %>% summarize_all(c(mean,length))
  lc = data.table(lc)
  lc = lc[,mx:=max(V1_fn2),by=list(louvain)]
  lc = data.frame(lc)
  lc = lc[lc$V2_fn2==lc$mx,]
  lc$mx <- NULL 
  gl <- pl[pl$SMTS=="Parathyroid",]
  gl$SMTS <- NULL
  gl <- merge(gl,unique(lc[,c("louvain","SMTS")]),by="louvain")
  lc = data.table(lc)
  lc = lc[,mx:=max(V1_fn2),by=list(SMTS)]
  #lc$SMTS <- ifelse(lc$V2_fn2==lc$mx,as.character(lc$SMTS),as.character(""))
  names(lc) <- c("SMTS","lopuvain","V1","V2","m1","m2","m3")
  pl <- pl[pl$SMTS!="Parathyroid",]
  
  pllist <- list()
  pllist[[1]] <- ggplot(pl,aes(x = V1, y = V2, colour = SMTS)) + geom_point(alpha = 0.95) + 
    xlim(xrange)+ylim(yrange)+
    xlab("z coordinate one") + ylab("z coordinate two")+
    guides(colour = FALSE) +theme_bw()+
    scale_color_manual(values = cols)+
    geom_label_repel(aes(label = SMTS),size=5,data=lc)+gg_aes
  
  pllist[[2]]  <- ggplot(gl,aes(x = V1, y = V2, colour = SMTS)) + geom_point(alpha = 0.95) + 
  xlim(xrange)+ylim(yrange)+
    xlab("z coordinate one") + ylab("z coordinate two")+
    geom_label_repel(aes(label = SAMPID),size=5)+ 
    scale_color_manual(values = cols)+
    guides(colour = FALSE) +theme_bw()+gg_aes
  
  cl <- rbind(pl,gl)
  pllist[[3]]  <- ggplot(cl,aes(x = V1, y = V2, colour = (log_likelihood))) + 
    geom_point() + 
    xlim(xrange)+ylim(yrange)+theme_bw()+
    xlab("z coordinate one") + ylab("z coordinate two")+
    facet_wrap(~data,nrow=1)+gg_aes+
    scale_colour_gradient(low = "#ffeda0", high = "#f03b20", na.value = NA,space = "Lab")
  
  lay <- rbind(c(1,2),
               c(3,3))
  
  pdf(file = paste0(WORKDIR,"5k_rowMeanScaled_SCVIS.pdf"),width = 15,height = 15)
  grid.arrange(grobs = pllist, layout_matrix = lay)
  dev.off()


    
}

## 2. rowZScaled _5k
if(T){
  pl <- myf("scvis_5krowZscaled/labs.txt",
            "scvis_5krowZscaled/gtex.tsv",
            "scvis_5krowZscaled/pth.tsv",
            "scvis_5krowZscaled/gtex_log_likelihood.tsv",
            "scvis_5krowZscaled/pth_log_likelihood.tsv")
  
  xrange <- range(pl$V1)
  yrange <- range(pl$V2)
  
  lc = pl %>% group_by(SMTS,louvain) %>% select(V1,V2) %>% summarize_all(c(mean,length))
  lc = data.table(lc)
  lc = lc[,mx:=max(V1_fn2),by=list(louvain)]
  lc = data.frame(lc)
  lc = lc[lc$V2_fn2==lc$mx,]
  lc$mx <- NULL 
  gl <- pl[pl$SMTS=="Parathyroid",]
  gl$SMTS <- NULL
  gl <- merge(gl,unique(lc[,c("louvain","SMTS")]),by="louvain")
  lc = data.table(lc)
  lc = lc[,mx:=max(V1_fn2),by=list(SMTS)]
  #lc$SMTS <- ifelse(lc$V2_fn2==lc$mx,as.character(lc$SMTS),as.character(""))
  names(lc) <- c("SMTS","lopuvain","V1","V2","m1","m2","m3")
  pl <- pl[pl$SMTS!="Parathyroid",]
  
  pllist <- list()
  pllist[[1]] <- ggplot(pl,aes(x = V1, y = V2, colour = SMTS)) + geom_point(alpha = 0.95) + 
    xlim(xrange)+ylim(yrange)+
    xlab("z coordinate one") + ylab("z coordinate two")+
    guides(colour = FALSE) +theme_bw()+
    scale_color_manual(values = cols)+
    geom_label_repel(aes(label = SMTS),size=5,data=lc)+gg_aes
  
  pllist[[2]]  <- ggplot(gl,aes(x = V1, y = V2, colour = SMTS)) + geom_point(alpha = 0.95) + 
    xlim(xrange)+ylim(yrange)+
    xlab("z coordinate one") + ylab("z coordinate two")+
    geom_label_repel(aes(label = SAMPID),size=5)+ 
    scale_color_manual(values = cols)+
    guides(colour = FALSE) +theme_bw()+gg_aes
  
  cl <- rbind(pl,gl)
  pllist[[3]]  <- ggplot(cl,aes(x = V1, y = V2, colour = (log_likelihood))) + 
    geom_point() + 
    xlim(xrange)+ylim(yrange)+theme_bw()+
    xlab("z coordinate one") + ylab("z coordinate two")+
    facet_wrap(~data,nrow=1)+gg_aes+
    scale_colour_gradient(low = "#ffeda0", high = "#f03b20", na.value = NA,space = "Lab")
  
  lay <- rbind(c(1,2),
               c(3,3))
  
  pdf(file = paste0(WORKDIR,"5k_rowZScaled_SCVIS.pdf"),width = 15,height = 15)
  grid.arrange(grobs = pllist, layout_matrix = lay)
  dev.off()
  
  
  
}
#

## 3. rowZScaled _5k
if(T){
  pl <- myf("scvis_10krowZscaled/labs.txt",
            "scvis_10krowZscaled/gtex.tsv",
            "scvis_10krowZscaled/pth.tsv",
            "scvis_10krowZscaled/gtex_log_likelihood.tsv",
            "scvis_10krowZscaled/pth_log_likelihood.tsv")
  
  xrange <- range(pl$V1)
  yrange <- range(pl$V2)
  
  lc = pl %>% group_by(SMTS,louvain) %>% select(V1,V2) %>% summarize_all(c(mean,length))
  lc = data.table(lc)
  lc = lc[,mx:=max(V1_fn2),by=list(louvain)]
  lc = data.frame(lc)
  lc = lc[lc$V2_fn2==lc$mx,]
  lc$mx <- NULL 
  gl <- pl[pl$SMTS=="Parathyroid",]
  gl$SMTS <- NULL
  gl <- merge(gl,unique(lc[,c("louvain","SMTS")]),by="louvain")
  lc = data.table(lc)
  lc = lc[,mx:=max(V1_fn2),by=list(SMTS)]
  #lc$SMTS <- ifelse(lc$V2_fn2==lc$mx,as.character(lc$SMTS),as.character(""))
  names(lc) <- c("SMTS","lopuvain","V1","V2","m1","m2","m3")
  pl <- pl[pl$SMTS!="Parathyroid",]
  
  pllist <- list()
  pllist[[1]] <- ggplot(pl,aes(x = V1, y = V2, colour = SMTS)) + geom_point(alpha = 0.95) + 
    xlim(xrange)+ylim(yrange)+
    xlab("z coordinate one") + ylab("z coordinate two")+
    guides(colour = FALSE) +theme_bw()+
    scale_color_manual(values = cols)+
    geom_label_repel(aes(label = SMTS),size=5,data=lc)+gg_aes
  
  pllist[[2]]  <- ggplot(gl,aes(x = V1, y = V2, colour = SMTS)) + geom_point(alpha = 0.95) + 
    xlim(xrange)+ylim(yrange)+
    xlab("z coordinate one") + ylab("z coordinate two")+
    geom_label_repel(aes(label = SAMPID),size=5)+ 
    scale_color_manual(values = cols)+
    guides(colour = FALSE) +theme_bw()+gg_aes
  
  cl <- rbind(pl,gl)
  pllist[[3]]  <- ggplot(cl,aes(x = V1, y = V2, colour = (log_likelihood))) + 
    geom_point() + 
    xlim(xrange)+ylim(yrange)+theme_bw()+
    xlab("z coordinate one") + ylab("z coordinate two")+
    facet_wrap(~data,nrow=1)+gg_aes+
    scale_colour_gradient(low = "#ffeda0", high = "#f03b20", na.value = NA,space = "Lab")
  
  lay <- rbind(c(1,2),
               c(3,3))
  
  pdf(file = paste0(WORKDIR,"10k_rowZScaled_SCVIS.pdf"),width = 15,height = 15)
  grid.arrange(grobs = pllist, layout_matrix = lay)
  dev.off()
  
  
  
}
#
#
#