# Seurat Cheat Sheet

```R
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```



## Statistics



## Plotting

#### Correlation of features

```R
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
```

