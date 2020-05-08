# PHT_project
The scripts and data files used for generating figures for Lucy's parathyroid project



# Script for extracting genomic sequences around the variant (small indels)

```
SNV_seqExtract.R vcf_file up-window down-window

```


# Embedding of the PTH samples into GTEx data

Feature selection:
- RNA-seq data shows dependance of expression mean on the variance for a gene
- In order to correct for it, we used the strategies applied in the field to correct for the mean-variance relationship (e.g. by DESeq2, Seurat)
- Briefly, we calculated mean and variance for each gene across the combined GTEx and PTH dataset
- We fit a polynomial curve of 2nd degree to estimate the dependance
- Next, we transformed the count table using the mean and predicted variance of the gene to Z-scores 
- Based on the variance of the Z-score, top 5000 and 10000 genes were identified as highly variable features


Low dimesion visualization of the data:
- We used SCVIS method to visualize the high dimension data into a 2D space.    
- SCVIS method uses neural network based approach to derive probabilistic generative model from the high dimension data.
- The method performs better than routinely used approxiate BH t-SNE alogirthm for dimensionality reduction. 
- Furthermore, it allows for embendding other dataset on to the original 2D projection based on the learned model.
- We used the default setting with 100 iterations to calculate optimized parameters for the SCVIS model based on GTEx data and then embeded PTH data to it using the same learned model
- We used 5k and 10k features respectively where the count tables were normalized by Z-transformation of the gene signal across cell-types
- In addition, we explored the 5k features such that prior to Z-scaling of the gene signal across celltypes, we first log transformed the TPM counts by addition of a pseudocount 1.


Clustering of the data points:
- We used modularity optimization algorithm (Louvain) to cluster the 2D representation into groups of cells. 
- We used 100 iterations to group the 2D projection of the cell into nearest group of cells. 
