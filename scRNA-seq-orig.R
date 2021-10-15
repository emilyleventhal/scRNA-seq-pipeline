install.packages('Seurat')
# devtools::install_github("satijalab/seurat-data")
# install.packages('SeuratData')
# install.packages('BiocManager')
# BiocManager::install('multtest')
# BiocManager::install('limma')
# install.packages('metap')
# remotes::install_github("ChristophH/sctransform@develop")


library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)

# (1) find files that contain "SampleX" + barcodes.tsv.gz, genes.tsv.gz, and matrix.mtx.gz at end
#     change names to barcodes.tsv.gz, genes.tsv.gz, matrix.mtx
# (2) unzip all 
#     gzip -d barcodes.tsv.gz
#     must unzip them 

## move features.tsv to genes.tsv

# deleted third column of genes.tsv for sample_01 but that didn't work 

#setwd("~/Desktop/-/School/Scottlab")
#project_name = "sample1"
#data_dir = "GSE124952_RAW_ALL/Sample2"
scRNAseq <- function(data_dir, project_name) {
  dir.name = paste(project_name, "_data", sep = "")
  dir.create(dir.name)
  
  sample.data = Read10X(data.dir = data_dir) # expression matrix with genes 
  ?Read10x
  sample <- CreateSeuratObject(counts=sample.data, project = project_name, min.cells=3, min.features=20) 
  # 54 features for > 100 
  # 16k features for > 200 in the sample 
  
  ### Standard pre-processing workflow ### 
  # [[ operator adds columns to object metadata 
  
  # stash QC stats, MT is set of mitochondrial genes
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
  
  # view QC metrics 
  name = paste(dir.name, "/QC_metrics.txt", sep = "")
  sink(name)
  print(sample@meta.data) # percent.mt is 0 tho
  sink()
  
  
  # visualize QC metrics as violin plot 
  name = paste(dir.name, "/QC_metrics_vln.pdf", sep = "")
  pdf(name)
  print(VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3))
  dev.off()
  # FeatureScatter to visualize feature-feature relationships 
  # can be used for anything calculated by object (columns in object metadata, PC scores, etc.)
  name = paste(dir.name, "/feature-feature-relationships.pdf", sep = "")
  pdf(name, width=10)
  plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
  
  # filter out cells that have unique feature counts over 2500 or < 200
  # filter out cells that have > 5 % mitochondrial counts 
  sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  ### normalizing data ### 
  # global-scaling normalization method "LogNormalize" - normalizes feature expression measurements for each cell by total expression, multiplies by scale factor (10,000) + log-transforms resutl
  # normalized values stored in sample[["RNA"]]@data
  sample <- NormalizeData(sample) # sample1 <- NormalizeData(sample1, normalization.method = "LogNormalize", scale.factor = 10000) # with default values ]
  
  ### identification of highly variable features (feature selection) ### 
  # calculate subset of features that exhibit high cell-to-cell variation in dataset
  # helps highlight biological signal in single-cell datasets (procedure here: https://linkinghub.elsevier.com/retrieve/pii/S0092867419305598)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  
  # identify 10 most highly variable genes
  top10 <- head(VariableFeatures(sample), 10)
  
  # top variable genes in general: 
  name = paste(dir.name, "/top_variable_genes.txt", sep = "")
  sink(name)
  print("Top highly variable genes:")
  print(cat(VariableFeatures(sample), sep="\n"))
  sink()
  
  # plot variable features with and without labels 
  plot1 <- VariableFeaturePlot(sample)
  plot2 <- LabelPoints(plot=plot1, points = top10, repel = TRUE)
  name = paste(dir.name, "/variable_features.pdf", sep = "")
  pdf(name)
  print(plot2)
  dev.off()
  
  ### scaling the data ### 
  # apply linear transformation 'scaling' 
  #     -> shifts expression so mean expression is 0
  #     -> scales expression so that variance is 1 
  #           -> gives equal weight in downstream analyses so that highly-expressed genes do not dominate
  # stored in sample[["RNA"]]@scale.data
  
  all.genes <- rownames(sample)
  sample <- ScaleData(sample, features = all.genes)
  name = paste(dir.name, "/scaled_data.txt", sep = "")
  sink(name)
  options(max.print = .Machine$integer.max)
  print(sample[["RNA"]]@scale.data)
  sink()
  
  ### perform linear dimensional reduction ###
  # perform PCA on scaled data - only variable features are used as input by defulat 
  
  sample <- RunPCA(sample) # sample1 <-RunPCA(sample1, features = VariableFeatures(object = sample1))
  # visualize cells and features that define the PCA 
  
  # examine and visualize PCA results in diff ways
  name = paste(dir.name, "/PCA_results.txt", sep = "")
  sink(name)
  print(sample[["pca"]], dims = 1:50, nfeatures = 5) # how many dimensions should this be? 
  sink()
  pca.dir.name = paste(dir.name, "/PCA_visualizations", sep = "")
  dir.create(pca.dir.name)
  
  # how many clusters do we want to save? 
  # save pc_o1
  answer = readline(prompt = "Do you want to visualize a cluster? Y/N: ")
  while(answer == "Y"){
    cluster <- as.integer(readline(prompt = "Which cluster number? "))
    name = paste(pca.dir.name, "/PC_", cluster, ".pdf", sep = "")
    pdf(name)
    print(VizDimLoadings(sample, dim = cluster, reduction = "pca")) # dims = 1:2,
    dev.off()
    answer = readline(prompt = "Do you want to visualize another cluster? Y/N: ")
  }
  
  answer = readline(prompt = "Do you want to compare two clusters? Y/N: ")
  while (answer == "Y"){
    print("Type the two clusters for the comparison below, and press enter twice when finished: ")
    clusters <- scan()
    clusters_str = paste(clusters, sep = "", collapse = "_")
    clusters_str
    name = paste(pca.dir.name, "/comparison_PCA_", clusters_str, ".pdf", sep = "")
    pdf(name)
    print(DimPlot(sample, reduction = "pca", dims = clusters))
    dev.off()
    answer = readline(prompt = "Do you want to compare more clusters? Y/N: ")
  }
  
  #print(DimPlot(sample, reduction = "pca", dims = 1:2))
  # explore primary sources of heterogeneity in a dataset 
  # DimHeatmap() helpful when trying to decide which PCs to include for further downstream anlyses 
  #     -> cells and features are ordered according to their PCA scores
  # DimHeatmap(sample, dims = 1, cells = 500, balanced = TRUE)
  
  # how many for the heat maps? 
  answer = readline(prompt = "Do you want to create heatmaps for clusters? Y/N: ")
  while (answer == "Y"){
    oneBool = readline(prompt = "Do you want to create one per page or multiple? one/multiple: ")
    if(oneBool == "multiple"){
      print("Type the clusters, and press enter twice when finished: ")
      clusters <- scan()
      clusters_str = paste(clusters, sep = "", collapse = "_")
      name = paste(pca.dir.name, "/heatmap_comp_PC_", clusters_str, ".pdf", sep = "")
      pdf(name)
      print(DimHeatmap(sample, dims = clusters, cells = 500, balanced = TRUE))
      dev.off()
    }
    else{
      cluster = as.integer(readline(prompt = "Which cluster? "))
      name = paste(pca.dir.name, "/heatmap_PC", cluster, ".pdf", sep = "")
      pdf(name)
      print(DimHeatmap(sample, dims = cluster, cells = 500, balanced = TRUE))
      dev.off()
    }
    answer = readline(prompt = "Do you want to create another heatmap? Y/N: ")
    
  }

  
  ### determine the 'dimensionality' of the dataset ###
  # overcome technical noise in any single feature for scRNA-seq data 
  # cluster cells based on PCA scores (each PC represents a 'metafeature' that combines info across a correlated feature set)
  # resampling test to permute subset of data (1%) and rerun PCA, constructing 'null distribution' of feature scores
  # identify 'significant' PCs as those with strong enrichment of low p-value features
  sample <- JackStraw(sample, num.replicate = 100)
  sample <- ScoreJackStraw(sample, dims = 1:20)
  # visualization to compare distribution of p-values for each PC with a uniform distribution (dashed line) 
  # significant PCs will show strong encrichment of features with low p-values (solid curve above dashed line)
  name = paste(dir.name, "/JackStrawPlot-p-values-PCs.pdf", sep ="")
  pdf(name)
  print(JackStrawPlot(sample, dims = 1:20))
  dev.off()
  
  # Elbow plot is a ranking of principle components based on % of variance explained by each one 
  name = paste(dir.name, "/ElbowPlot-ranking-by-variance.pdf", sep = "")
  pdf(name)
  print(ElbowPlot(sample))
  dev.off()
  # where the elbow is = majority of true signal is captured in first X before that
  # finding 'optimal number of clusters' --> it's 10 in this case
  
  ### IDENTIFYING true dimensionality of a dataset - 3 approaches ###
  # (1) more supervised, explore PCs to determine sources of heterogenity (could be used in conjunction with GSEA)
  # (2) statistical test based on random null model (time-consuming for large datasets), may not return a clear PC cutoff
  # (3) heuristic commonly used, can be calculated instantly 
  
  # chose 10, because that is what we found in the elbow plot --> HOW WILL WE DO THIS IN THE FUNCTION?
  # "We advise users to err on the higher side when choosing this parameter."
  opt_cluster_num=as.integer(readline(prompt="What is the optimal number of clusters? Check the 'elbow' in the elbow plot: "))
  
   ### cluster the cells ###
  sample <- FindNeighbors(sample, dims = 1:opt_cluster_num)
  sample <- FindClusters(sample, resolution = 0.5) # 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells.--> larger for larger datasets
  # look at cluster IDs of first 5 cells
  name = paste(dir.name, "/ClusterIDs.txt", sep = "")
  sink(name)
  print("Cluster IDs of first 100 cells:")
  print(head(Idents(sample), 100))  
  sink()

  ### run non-linear dimensional reduction (UMAP/tSNE) ###
  # non-linear dimensional reduction techniques to visualize datasets
  # place similar cells together in low-dimensional space 
  # cells within graph-based clusters determined above should co-localize on these dimension reduction plots
  # use the same PCs as the input to clustering analysis 
  sample <- RunUMAP(sample, dims = 1:opt_cluster_num)
  name = paste(dir.name, "/umap.pdf", sep = "")
  pdf(name)
  print(DimPlot(sample, reduction = "umap"))
  dev.off()
  
  ### Finding differentially expressed features (cluster biomarkers) ###
  # find markers that define clusters via differential expression 
  # identifies positive and negative markers of a single cluster (idnet.1) 
  # FindAllMarkers() automates this process, but you can also test groups of clusters vs. each other 
  # min.pct requires feature to be detected at a minimum % in either of the two groups of cells 
  # thresh.test requires feature to be differentially epxressed by some amount 
  # max.cells.per.ident will downsample each identity class to have no more cells than whatever
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  sample.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # name = paste(dir.name, "/all_markers.txt", sep = "")
  sample.markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = avg_log2FC)
  
  answer=readline(prompt="Do you want to compare two cluster biomarkers? Y/N: ")
  while(answer == "Y"){
    print("Type the FIRST cluster(s) for the comparison below, and press enter twice when finished: " )
    x <- scan()
    clusters1 = paste(x, sep = "", collapse = "_")
    print("Type the SECOND cluster(s) for the comparison below, and press enter twice when finished: " )
    y <- scan()
    clusters2 = paste(y, sep = "", collapse = "_")
    print_statement = paste("Comparing cluster(s) ", clusters1, " vs ", clusters2, "...", sep = "")
    print(print_statement)
    # cluster_1 = as.integer(readline(prompt = "What is the first cluster(s) you want to compare? e.g. 1, c(1,5) "))
    cluster.markers <- FindMarkers(sample, ident.1 = x, ident.2 = y, min.pct = 0.25)
    name = paste(dir.name, "/biomarkers_PC_", clusters1, "_vs_", clusters2, ".txt", sep = "")
    # =here, 252
    sink(name)
    print(head(cluster.markers, n = 100))
    sink()
    answer=readline(prompt="Do you want to compare two more clusters? Y/N: ")
  }
  

  # tests for differential expression : https://satijalab.org/seurat/articles/de_vignette.html
  # ^ set with test.use parameter test.use = "roc" in FindMarkers
  
  # visualize marker expression 
  answer = readline(prompt = "Do you want to visualize marker expression? Y/N: ")
  while(answer == "Y"){
    print("What genes do you want to visualize the expression of? ")
    x <- scan(what= "String")
    genes = paste(x, sep = "", collapse="_")
    name = paste(dir.name, "/expression_", genes, ".pdf", sep = "")
    rawBool <- readline("Do you want to plot raw counts? Y/N: ")
    pdf(name)
    if(rawBool == "Y"){
      print(VlnPlot(sample, features = x, slot = "counts", log = TRUE))
    }
    else{
      print(VlnPlot(sample, features = x))
    }
    dev.off()
    answer = readline(prompt = "Do you want to visualize another marker expression? Y/N: ")
  }
  print("What genes do you want a feature plot of? Press enter after each gene:")
  x <- scan(what = "String") #what = "String"
  print("Scanned")
  name = paste(dir.name, "/FeaturePlot.pdf", sep = "")
  pdf(name)
  print(FeaturePlot(sample, features = x))
  dev.off()
  # generate expression heatmap for cells and features 
  # plot top 10 markers for each cluster 
  name = paste(dir.name, "/Top10Markers.pdf", sep = "")
  top10 <- sample.markers %>% 
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  pdf(name)
  print(DoHeatmap(sample, features = top10$gene) + NoLegend())
  #DoHeatmap(sample, features = top10$gene) + NoLegend()
  dev.off()

}

working_d <- readline(prompt = "What is the working directory? ")
# ~/Desktop/main/School/scottlab
setwd(working_d)

data_dir <- readline(prompt = "What is the directory of the data? ")
# GSE124952_RAW_ALL/Sample1
# sample_01/outs/filtered_feature_bc_matrix
# combined/outs/filtered_feature_bc_matrix
project_name <- readline(prompt = "What is the name of the project? ")
#sample1 

scRNAseq(data_dir, project_name)


#sample1.data <- Read10X(data.dir = "GSE124952_RAW_ALL/Sample2")

#sample1 <- CreateSeuratObject(counts=sample1.data, project = "Sample1", min.cells=3, min.features=200)
##pbmc <- CreateSeuratObject(counts = pbmc.data, project="pbmc3k", min.cells = 3, min.features = 200)

#sample1



##### integrate multiple 
# https://satijalab.org/seurat/articles/integration_introduction.html
