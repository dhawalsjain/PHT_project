#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Arguments are not provided!\n\nUsage:\n\n SNV_seqExtract.R vcf_file up-window down-window", call.=FALSE)
}else if(length(args)==2){
  args[3] = args[2]
}
#Command line usage:
## SNV_seqExtract.R 'path2VCF.vcf' 200 200
# Inline use
##args <- c("path to vcf",200,200)

suppressPackageStartupMessages(library(GenomicRanges,quietly = T))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19,quietly = T))
suppressPackageStartupMessages(library(Biostrings,quietly = T))

genome <- BSgenome.Hsapiens.UCSC.hg19
WIN=c(args[2],args[3])
vcf <- read.delim(args[1],header=F,comment.char = "#",stringsAsFactors = F)
vcf$end <- vcf$V2 + width(DNAStringSet(vcf$V4))-1

if(!grep("chr",vcf$V1[1])){
  vcf$V1 <- paste0("chr",vcf$V1)
}

vcf <- with(vcf,GRanges(V1,IRanges(V2,end),"*",REF=V4,ALT=V5))
vcf$id <- paste0(seqnames(vcf),":",start(vcf),"|",vcf$REF,">",vcf$ALT)

sqb <- getSeq(genome,vcf)
if(sum(sqb==vcf$REF)!=length(sqb)){
  stop("the genome build installed here does not match that from the vcf file\n")
}

up <- BSgenome::getSeq(genome,flank(vcf,WIN[1],start = T))
dn <- BSgenome::getSeq(genome,flank(vcf,WIN[2],start = F))
ref = Biostrings::xscat(up,DNAStringSet(as.character(vcf$REF),start=1) ,dn)
alt = Biostrings::xscat(up,DNAStringSet(as.character(vcf$ALT),start=1),dn)
names(ref) = vcf$id
names(alt) = paste0(vcf$id,"_Alt")

writeXStringSet(ref,file="refAllele.fa")
writeXStringSet(alt,file="AltAllele.fa")



