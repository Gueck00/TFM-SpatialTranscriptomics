library(SpatialExperiment)
library(spatialLIBD)
library(Seurat)
library(BayesSpace)
library(DR.SC)
library(SC.MEB)
library(scran)
library(BASS)
library(ggpubr)

# Visium data (DLPFC)

dlpfc_fulldata <- fetch_data(type = "spe")

colData(dlpfc_fulldata) <- colData(dlpfc_fulldata)[,c("sample_id","cell_count","subject","layer_guess","in_tissue","array_row","array_col","key")] 

colData(dlpfc_fulldata)$row <- colData(dlpfc_fulldata)$array_row
colData(dlpfc_fulldata)$col <- colData(dlpfc_fulldata)$array_col

dlpfc_fulldata <- dlpfc_fulldata[rowSums(assay(dlpfc_fulldata)) != 0,] # Quita todos los genes no presentes en las muestras

for (k in 1:length(reducedDimNames(dlpfc_fulldata))){
  reducedDim(dlpfc_fulldata) <- NULL
}

sample_names <- as.character(unique(dlpfc_fulldata$sample_id))
dlpfc_list <- lapply(sample_names, function(x) dlpfc_fulldata[, dlpfc_fulldata$sample_id == x])

rm("dlpfc_fulldata")

## Spatially-aware clustering

### BayesSpace

dlpfc_list <- lapply(dlpfc_list, function(x) spatialPreprocess(x, platform="Visium",  n.PCs = 50, log.normalize = F))

set.seed(08042024)
dlpfc_list <-  lapply(dlpfc_list, function(x) spatialCluster(x, use.dimred = "PCA", q = length(unique(x$layer_guess))-1, nrep = 100000, burn.in = 20000)) 
# Por la obvia ingente cantidad de iteraciones y datos, se ejecutó en el cluster del CIPF, de forma que cargo los resultados en el chunk posterior

for (i in 1:length(dlpfc_list)){
  dlpfc_list[[i]]$spatial.cluster <- as.factor(dlpfc_list[[i]]$spatial.cluster)
}

### DR-SC

dlpfc_list.seu <- dlpfc_list

dlpfc_list.seu <- lapply(dlpfc_list.seu, function(x) CreateSeuratObject(counts = assay(x, "counts"), data = assay(x, "logcounts"), assay = "Spatial", meta.data = as.data.frame(colData(x))))

dlpfc_list.seu <- lapply(dlpfc_list.seu, function(x) FindVariableFeatures(object = x, nfeatures = 500, verbose = F))

dlpfc_list.seu <- lapply(dlpfc_list.seu, function(x) DR.SC(x, K = length(unique(x$layer_guess))-1, platform = 'Visium', verbose=F))

for (i in 1:length(dlpfc_list)){
  dlpfc_list[[i]]$spatial.drsc.cluster <- as.factor(dlpfc_list.seu[[i]]$spatial.drsc.cluster)
}

rm("dlpfc_list.seu")

### SC.MEB

platform = "Visium"
beta_grid = seq(0,4,0.2)
K_set= 2:10
parallel=TRUE
num_core = 3
PX = TRUE
maxIter_ICM = 10
maxIter = 50

.find_neighbors2 <- function (sce, platform) 
{
  library(purrr)
  library(Matrix)
  if (platform == "Visium") {
    offsets <- data.frame(x.offset = c(-2, 2, -1, 1, -1, 
                                       1), y.offset = c(0, 0, -1, -1, 1, 1))
  }
  else if (platform == "ST") {
    offsets <- data.frame(x.offset = c(0, 1, 0, -1), y.offset = c(-1, 
                                                                  0, 1, 0))
  }
  else {
    stop(".find_neighbors: Unsupported platform \"", platform, 
         "\".")
  }
  spot.positions <- SingleCellExperiment::colData(sce)[, c("col", "row")]
  spot.positions$spot.idx <- seq_len(nrow(spot.positions))
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
  neighbors <- merge(as.data.frame(neighbor.positions), as.data.frame(spot.positions), 
                     by.x = c("x.pos", "y.pos"), by.y = c("col", "row"), suffixes = c(".primary", 
                                                                                      ".neighbor"), all.x = TRUE)
  neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- unname(df_j)
  df_j <- purrr::map(df_j, function(nbrs) purrr::discard(nbrs, function(x) is.na(x))) 
  n_with_neighbors <- length(purrr::keep(df_j, function(nbrs) length(nbrs) > 
                                           0))
  message("Neighbors were identified for ", n_with_neighbors, 
          " out of ", ncol(sce), " spots.")
  n <- length(df_j)
  D <- sparseMatrix(i = 1:n, j = 1:n, x = 0)
  for (i in 1:n) {
    if (length(df_j[[i]]) != 0) 
      D[i, df_j[[i]]] <- 1
  }
  D
}

Adj_sp <- lapply(dlpfc_list, function(x) .find_neighbors2(x, platform = "Visium"))

fit <- list()

set.seed(08042024)
for(i in 1:length(dlpfc_list)){
  fit[[i]] <- SC.MEB(reducedDim(dlpfc_list[[i]], type = "PCA")[,1:15], Adj_sp = Adj_sp[[i]], beta_grid = beta_grid, K_set= length(unique(dlpfc_list[[i]]$layer_guess))-1 , parallel=parallel, num_core = num_core, PX = PX, maxIter_ICM=maxIter_ICM, maxIter=maxIter)
}

for (i in 1:length(dlpfc_list)){
  dlpfc_list[[i]]$spatialcluster.sc_meb <- as.factor(fit[[i]][,1]$x)
}

rm("Adj_sp","fit")

### BASS

set.seed(03052024)
# Set up BASS object
Bass_list <- lapply(dlpfc_list, function(x) createBASSObject(
  X = list(as.matrix(assay(x,"counts"))), 
  xy = list(as.data.frame(colData(x)[,c("row","col")])), 
  C = 20, 
  R = length(unique(x$layer_guess))-1, 
  beta_method = "SW", 
  init_method = "mclust", 
  nsample = 10000
))

Bass_list <- lapply(Bass_list, function(x) BASS.preprocess(x, doLogNormalize = TRUE,
                                                           geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
                                                           scaleFeature = FALSE, nPC = 20, doBatchCorrect = F))

Bass_list <- lapply(Bass_list, function(x) BASS.run(x))

Bass_list <- lapply(Bass_list, function(x) BASS.postprocess(x))

for (i in 1:length(dlpfc_list)){
  dlpfc_list[[i]]$bass.spatialcluster <- as.factor(Bass_list[[i]]@results$z[[1]])
}

## Non-spacially aware clustering

### k-means (base R)

set.seed(02052024)
camelias <- lapply(dlpfc_list, function(x) kmeans(x = reducedDim(x, type = "PCA"), centers = length(unique(x$layer_guess))-1, nstart = 50))

for (i in 1:length(dlpfc_list)){
  dlpfc_list[[i]]$kmeans.cluster <- as.factor(camelias[[i]]$cluster)
}

### Louvain

set.seed(02052024)
k <- 15

for (i in 1:length(dlpfc_list)){
  g <- buildSNNGraph(dlpfc_list[[i]], k = k, use.dimred = "PCA")
  g_walk <- igraph::cluster_walktrap(g)
  dlpfc_list[[i]]$louvain.cluster <- as.factor(g_walk$membership)
}

# Checkpoint Visium
save(dlpfc_list, file = "./datos/definitivos/dlpfc_list_definitive.Rdata")
