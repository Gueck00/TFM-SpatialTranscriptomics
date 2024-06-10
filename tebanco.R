library(SpatialExperiment)
library(spatialLIBD)
library(Seurat)
library(infotheo)
library(clevr)
library(ggplot2)
library(cluster)
library(scran)
library(scater)
library(pheatmap)

load("./datos/definitivos/dlpfc_list_definitive.Rdata")

for(i in 1:length(dlpfc_list)){
  dlpfc_list[[i]] <- dlpfc_list[[i]][,!is.na(dlpfc_list[[i]]$layer_guess)]
}

Clustemios_list <- lapply(dlpfc_list, function(x) as.data.frame(colData(x)[,c("layer_guess","spatial.cluster",
                                                        "spatial.drsc.cluster","spatialcluster.sc_meb",
                                                        "bass.spatialcluster","kmeans.cluster",
                                                        "louvain.cluster")]))

# Accuracy

Acc <- function(dt){
  if(!is.data.frame(dt)){
    stop("dt is not a data.frame")
  }
  else{
    GT <- dt[,1]
    Acc_values <- data.frame(row.names = c("NMI","ARI","HOM","COM"))
    for (i in 1:dim(dt)[2]){
      nmi <- infotheo::mutinformation(GT, dt[,i]) / sqrt(infotheo::entropy(GT) * infotheo::entropy(dt[,i]))
      ari <- mclust::adjustedRandIndex(GT,dt[,i])
      hom <- clevr::homogeneity(GT,dt[,i])
      com <- clevr::completeness(GT,dt[,i])
      aceituna <- c(nmi,ari,hom,com)
      Acc_values[,i] <- aceituna
    }
    
    colnames(Acc_values) <- colnames(dt)
    return(Acc_values)
  }
}

accuracy_dlpfc <- lapply(Clustemios_list, function(x) Acc(x)) # It works!

# Continuity

pasotismo <- function(vector) {
  
  primer_elemento <- vector[1]
  iguales <- sum(vector[-1] == primer_elemento) >= 4
  
  return(iguales)
}

compute_pas <- function(se, clusters){
  location <- as.matrix(SpatialExperiment::spatialCoords(se))
  dist_matrix <- as.matrix(dist(location, method = "euclidean"))
  diag(dist_matrix) <- Inf
  Pas_value <- data.frame(row.names = c("PAS"))
  
  for(i in 1:dim(clusters)[2]){
    vecinos_iguales <- apply(X = dist_matrix, MARGIN = 1, function(x) 
      pasotismo(c(clusters[which(x== Inf),i],clusters[names(head(x[order(x)], 6)),i])))
    PAS <- 1 - sum(vecinos_iguales)/dim(clusters)[1]
    Pas_value[,i] <- PAS
  }
  colnames(Pas_value) <- colnames(clusters)
  return(Pas_value)
}

PAS_list <- mapply(compute_pas,dlpfc_list,Clustemios_list, SIMPLIFY = F) #It works!

compute_asw <- function(se, clusters){
  location <- as.matrix(SpatialExperiment::spatialCoords(se))
  dist_matrix <- dist(location, method = "euclidean")
  Asw_value <- data.frame(row.names = c("ASW"))
  for (i in 1:dim(clusters)[2]){
    silk <- cluster::silhouette(x = as.numeric(clusters[,i]), dist = dist_matrix)
    silk <- as.data.frame(silk)
    Asw_value[,i] <- mean((silk$sil_width + 1)/2)
  }
  colnames(Asw_value) <- colnames(clusters)
  return(Asw_value)
}

ASW_list <- mapply(compute_asw,dlpfc_list,Clustemios_list, SIMPLIFY = F) #It works!

# Genes marcadores

MoranI <- function(se, clusters){
  location <- as.matrix(SpatialExperiment::spatialCoords(se))
  
  dist_matrix <- as.matrix(dist(location, method = "euclidean"))
  
  diag(dist_matrix) <- Inf
  
  for(i in 1:dim(dist_matrix)[2]){
    dist_matrix[order(dist_matrix[,i])[1:6],i] <- rep(1, 6)
    dist_matrix[-order(dist_matrix[,i])[1:6],i] <- rep(0, dim(dist_matrix)[1]-6)
  }
  
  rownames(se) <- rowData(se)$gene_name
  Moran_value <- data.frame(row.names = c("MoranI"))
  
  for (i in 1:dim(clusters)[2]){
    colLabels(se) <- as.factor(as.character(clusters[,i]))
    markers <- scran::findMarkers(se, test = "binom", direction = "up")
    marc_gen <- c()
    
    for (j in 1:length(markers)){
      interesting <- markers[[j]]
      best_set <- interesting[interesting$Top <= 1, ]
      marc_gen <- unique(c(marc_gen,rownames(best_set)))
    }
    
    trimed_gen <- as.matrix(assay(se, "logcounts")[marc_gen,])
    
    average_gen <- rowSums(trimed_gen)/dim(trimed_gen)[2]
    
    spots <- colnames(trimed_gen)
    
    MoranIgen <- c()
    for (j in marc_gen){
      
      xi <- trimed_gen[j,]
      xmedia <- average_gen[j]
      
      xj <- data.frame()
      for (k in 1:length(spots)){
        vecinos_gen <- trimed_gen[j,names(dist_matrix[which(dist_matrix[,spots[k]] == 1),spots[k]])]
        
        xj <- rbind(xj,vecinos_gen)
      }
      
      rownames(xj) <- spots
      xij <- c()
      for (k in spots){
        suma <- sum((xi[k]-xmedia)*(xj[k,]-xmedia))
        xij[k] <- suma
      }
      
      fullsumatorios <- sum(xij)/sum((xi-xmedia)^2)
      MoranI <- (fullsumatorios*length(spots))/(6*length(spots))
      MoranIgen <- c(MoranIgen,MoranI)
    }
    Moran_value[,i] <- mean(MoranIgen)
  }
  colnames(Moran_value) <- colnames(clusters)
  return(Moran_value)
}

MoranI_list <- mapply(MoranI,dlpfc_list,Clustemios_list, SIMPLIFY = F) #It kinda works, not sure though

print("hola1")

GearyC <- function(se, clusters){
  location <- as.matrix(SpatialExperiment::spatialCoords(se))
  
  dist_matrix <- as.matrix(dist(location, method = "euclidean"))
  
  diag(dist_matrix) <- Inf
  
  for(i in 1:dim(dist_matrix)[2]){
    dist_matrix[order(dist_matrix[,i])[1:6],i] <- rep(1, 6)
    dist_matrix[-order(dist_matrix[,i])[1:6],i] <- rep(0, dim(dist_matrix)[1]-6)
  }
  
  rownames(se) <- rowData(se)$gene_name
  Geary_value <- data.frame(row.names = c("GearyC"))
  
  for (i in 1:dim(clusters)[2]){
    colLabels(se) <- as.factor(as.character(clusters[,i]))
    markers <- scran::findMarkers(se, test = "binom", direction = "up")
    marc_gen <- c()
    
    for (j in 1:length(markers)){
      interesting <- markers[[j]]
      best_set <- interesting[interesting$Top <= 1, ]
      marc_gen <- unique(c(marc_gen,rownames(best_set)))
    }
    
    trimed_gen <- as.matrix(assay(se, "logcounts")[marc_gen,])
    
    average_gen <- rowSums(trimed_gen)/dim(trimed_gen)[2]
    
    spots <- colnames(trimed_gen)
    
    GearyCgen <- c()
    for (j in marc_gen){
      
      xi <- trimed_gen[j,]
      xmedia <- average_gen[j]
      
      xj <- data.frame()
      for (k in 1:length(spots)){
        vecinos_gen <- trimed_gen[j,names(dist_matrix[which(dist_matrix[,spots[k]] == 1),spots[k]])]
        
        xj <- rbind(xj,vecinos_gen)
      }
      
      rownames(xj) <- spots
      xij <- c()
      for (k in spots){
        suma <- sum((xi[k]-xj[k,])^2)
        xij[k] <- suma
      }
      
      fullsumatorios <- sum(xij)/sum((xi-xmedia)^2)
      GearyC <- (fullsumatorios*length(spots))/(2*6*length(spots))
      GearyCgen <- c(GearyCgen,GearyC)
    }
    Geary_value[,i] <- mean(GearyCgen)
  }
  colnames(Geary_value) <- colnames(clusters)
  return(Geary_value)
}

GearyC_list <- mapply(GearyC,dlpfc_list,Clustemios_list, SIMPLIFY = F)

print("hola2")

## Vamos a juntarlo todo en un dt

bench_list <- list()
for (i in 1:length(dlpfc_list)){
  bench_list[[i]] <- rbind(accuracy_dlpfc[[i]],PAS_list[[i]],ASW_list[[i]],MoranI_list[[i]],GearyC_list[[i]])
}

save(bench_list, file = "./datos/definitivos/benchmarking/bench_list.Rdata")
