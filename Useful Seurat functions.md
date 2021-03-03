# Useful scRNAseq Functions

## scRNAseq Analysis notes

Standard workflows expects:

1. **Filtering**

   1. Filter cells with low counts and expecially low number different of genes expressed (automatically performed by STAR or Cell Ranger, but you can be more inclusive doing it manually). The general approach is to plot the total number of counts and features on a log_y scale and cells ranked from highest to lowest and identify the elbow of the plot
   2. Filter genes with at leas 1 counts in more the a certain number of cells (just 1 if you want to be super conservative). It makes no sense to perform downstream analysis on genes with 0 count an almost all the cells

2. Evaluate % of mitochondrial DNA. Usually a common treshold is below 5%. But for cells sorted for long time the % can reach up to 20-25%

3. **Normalize** each cell for the total number of counts (even better for total counts excluding mt-DNA). Than, to avoid to deal with small numbers, multiply back for a scaling factor. Somebody multiplies for 10^6 to obtain CountsPerMillion. This measure is not really meaningful for scRNAseq since individual cells don't express millions of transcripts (extimated 350.000 mRNA molecules per cell). The best approach is to scale for the median number of reads for the sequencing experiments (10-40.000 reads). This gives back a dimension of values close to the original number of counts retrieved for each cell.

   *Some people now log1p transforms the counts and scales them (z-score) to identify highly variable genes. Then they retain only highly variable genes and renormalize the reads*

4. **Log transform**: The variance is proportional to the average expression of a gene and is not normally distributed. Thus to stabilize the variance across all the genes we log transform the data (this is automatically performed by default in the seurat function ` NormalizeData()` ). In order avoid to get -infinity values for 0 reads, we add a *pseudocount* of 1 so that they get log scaled to 0. Unfortunately this approach overextimates the variance of lovly epxressed genes. There are other approaches as rlog but they are mainly used in bulk-rnaseq and are not available in routinely workflows for scRNAseq. Note: in bulk RNAseq, log transformed (and than z-scored reads) are used uniquely for heatmap, PCA, clustering etc but not for differential expression becuase edgeR and DEseq2 packages requires raw counts and uses custom scaling methods

5.  **Scale**: for PCA, clustering etc, you need to center (mean = 0) and scale (sd =1) all the genes actually converting them in z-scores. You can also use scaled counts for plotting on scatters and violin plots, but people are more used to see log transformed expression values. Note: you should not scale non log-transformed expression values since they are not normally distributed.



## Matrix

Filter low quality cells

https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html

```R
counts_per_cell <- Matrix::colSums(counts)
counts_per_gene <- Matrix::rowSums(counts)
genes_per_cell <- Matrix::colSums(counts>0) # count gene only if it has non-zero reads mapped.
cells_per_gene <- Matrix::rowSums(counts>0) # only count cells where the gene is expressed
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
title('counts vs genes per cell')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
```



## Seurat

### Plot Mt DNA pct

```R
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Use "^Mt-" for mice
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

### Evaluate Variance explained by PCs to select the right amount for UMAP

```R
ElbowPlot(object = obj)
# The same as 
plot(obj@reductions$pca@stdev)
```

### Find markers for a specified clustering meta_field

```r
# the same as 
Idents(tongue) <- tongue@meta.data$free_annotation
markers <- FindAllMarkers(tongue, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")
```

### Annotate cluster based on example clustering with SingleR

```r
library(SingleR)
tongue.ref <- as.SingleCellExperiment(tongue, assay = "RNA")
lfile.ann <- as.SingleCellExperiment(lfile.filtered, assay = "RNA")

predictions <- SingleR(test=lfile.ann, assay.type.test=1, 
                       ref=tongue.ref, labels=tongue.ref$free_annotation)

plotScoreHeatmap(predictions)

# Import annotations
lfile.filtered <- AddMetaData(lfile.filtered, metadata = predictions$labels, col.name = "free_annotation")
```

### Identify relevant number of PCAs to use for tSNE and UMAP dimensionality reduction

```r
#Jackstraw 
ElbowPlot(epithelial, ndims = 80)
epithelial <- JackStraw(epithelial, num.replicate = 100,dims=80)
epithelial <- ScoreJackStraw(epithelial, dims = 1:80, do.plot = TRUE)
# JackStrawPlot(epithelial, dims = 1:80)
dev.off()
#From JackStraw we see 66 dim are significant
```

## Magic imputation of low expressed genes

https://github.com/satijalab/seurat/issues/612

in `genes` you have to select highly variable genes + your genes of interest. If the number of total genes selected it's too slow, results will become too sensitive to the input.

```r
pbmc <- CreateSeuratObject(...)
pbmc <- ScaleData(pbmc)
MAGIC_pbmc <- magic(pbmc, genes=pbmc@var.genes)
# or MAGIC_pbmc <- magic(pbmc, genes="all_genes")
MAGIC_pbmc <- RunPCA(MAGIC_pbmc )
MAGIC_pbmc <- RunTSNE(MAGIC_pbmc, dims.use = 1:10)
TSNEPlot(MAGIC_pbmc)
```

Some people recommend to keep clustering and running PCA on original data and not on the MAGIC-imputed one (it also takes longer).

Or in the developmental version

```r
devtools::install_github("KrishnaswamyLab/MAGIC", subdir='Rmagic')
# Run
seurat_obj <- magic(seurat_obj, genes = append(c("GENES", "TO", "IMPUTE"),
                                               VariableFeatures(seurat_obj)))
# and you'll add a new assay named MAGIC_RNA

# Visualize it
VlnPlot(seurat_obj, features = c("CDK2AP1", "PTEN", "PDCD4"), group.by = "clusters", assay = "MAGIC_RNA")
```

To visualize imputed genes in DimPlot a simple trick is to copy that gene expression level in the metadata

## Phate dimensionality reduction

```r
library("phateR")
normalized_data <- GetAssayData(obj)
# Don't run it in multicore (at leas on macOS)
# Transpose the sparse matrix, phate expects cells as row, genes as columns
data_phate <- phate(t(normalized_data), npca = 12, gamma = 1)

# Preview the embedding
quickplot(x=data_phate$embedding[,1], y=data_phate$embedding[,2 ])

colnames(data_phate$embedding) <- paste0("PHATE_", 1:2)
# We will now store this as a custom dimensional reduction called 'mds'
obj[["phate"]] <- CreateDimReducObject(embeddings = data_phate$embedding, key = "PHATE_", assay = DefaultAssay(obj))
```

PhateR runs in reticulate virtual env, thus you have to install python phate in that virtual env trough reticulate `reticulate::py_install("phate", pip=TRUE)` and eventual other missing python dependencies as `packaging`

## Nebulosa density gene expression plot

I don't linke this kind of plot because it is more representative of the cell nubmer than of the real expression level

```r
# install from github
library(Nebulosa)
plot_density(seur_obj, c("Actb", "Gapdh"))
```

## Mathijs recomandation for Gene Signature Scoring

1. Compute Z score of your genes/magic imputed genes
2. Take the average of the genes in the signature (made only of highly expressed genes)

