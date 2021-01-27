# Useful scRNAseq Functions

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

