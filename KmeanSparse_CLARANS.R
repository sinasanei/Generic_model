#### R spatial data packages
library(raster)
library(rgdal)
library(sp)
library(sf)

#########
library(tidyverse)
library(lattice)
library(USAboundaries)
library(spData)
library(spDataLarge)
library(prism)
library(rasterVis)
#####
library(sparcl)
library(qtcat)
clarans(ras_cmpl_DF, k=2, maxNeigbours = 100, nLocal = 10, mc.cores = 4)
km.out <- KMeansSparseCluster(ras_cmpl_DF, K = 3, wbounds = 5, nstart = 20, 
                              silent = FALSE, maxiter = 6)

res(ras_stack)
ras_stack_aggr <- aggregate(ras_stack, fact = 3)
res(ras_stack_aggr)
#extract all values from the raster into a data frame
ras_DF_aggr <- values(ras_stack_aggr)
ras_DF_aggr[is.nan(ras_DF_aggr)] <- NA
colnames(ras_DF_aggr)<- rasname
#check for NA's in the data
idxa <- complete.cases(ras_DF_aggr)
ras_cmpl_DF_aggr <- ras_DF_aggr[idxa,]
apply(ras_cmpl_DF_aggr,2, function(x) sum(is.na(x))) #number of missing values in each column
apply(ras_cmpl_DF_aggr, 2,mean)
apply(ras_cmpl_DF_aggr, 2,sd)
ras_cmpl_DF_aggr <- scale(ras_cmpl_DF_aggr, scale = TRUE)
#########
km.out <- KMeansSparseCluster(ras_cmpl_DF_aggr, K = 3, wbounds = 5, nstart = 20, 
                              silent = FALSE, maxiter = 6)
print(km.out)
plot(km.out)
km_Clust <- vector(mode = "integer", length = dim(ras_DF_aggr)[[1]])
km_Clust[!idxa] <- NA
c <- as.vector(unlist(km.out[[1]][2]))
km_Clust[idxa] <- c
KM_rst <- raster(ras_stack_aggr[[1]])
values(KM_rst) <- km_Clust
plot(KM_rst)
#####Tunning the kmeans algorithm 

km.perm <- KMeansSparseCluster.permute( ras_cmpl_DF_aggr , K=2 ,
                  wbounds = seq(1.1,10, len=10), nperms = 5)
saveRDS(km.perm, file="km_perm.Rdata")
plot(km.perm)
km.perm$nnonzerows

######
km.bestw <- rep(NA,6)
km.res <- list()
for(i in 1:6){
  cat(">>>iteration",i,":")
  km.perm <- KMeansSparseCluster.permute( ras_cmpl_DF_aggr , K=i+1 ,
                      wbounds = seq(2,10, len=10), nperms = 3)
  km.bestw[i] <- km.perm$bestw
  km.res[[i]] <- km.perm
}

saveRDS(km.res, file="km_result.Rdata")
saveRDS(km.bestw, file="km_bestw.Rdata")
######
km.out.list <- list()
for(nClust in 1:6){
  cat(">>> Clustering data for nClust =",nClust+1,"......")
  km.out <- KMeansSparseCluster(ras_cmpl_DF_aggr, K = nClust+1 , wbounds = km.bestw[[i]],
                                nstart = 25, silent = FALSE, maxiter = 10)
  km.out.list[[nClust]] <- km.out
  
}
  
saveRDS(km.out.list, file="km_out.Rdata")

KM_Sparse_rst <- raster(ras_stack_aggr[[1]])

for(i in 1:6){ 
  kmSparse <- vector(mode = "integer", length = dim(ras_DF_aggr)[[1]])
  kmSparse[!idxa] <- NA
  c <- as.vector(unlist(km.out.list[[i]][[1]][2]))
  kmSparse[idxa] <- c
  Tmp_Sparse_rst <- raster(ras_stack_aggr[[1]])
  values(Tmp_Sparse_rst) <- kmSparse
  if(i==1){
    KM_Sparse_rst <- Tmp_Sparse_rst
  }else{
    KM_Sparse_rst <- stack(KM_Sparse_rst, Tmp_Sparse_rst)
  }
}
KM_Sparse_rst
saveRDS(KM_Sparse_rst, file="KM_Sparse_rst.Rdata")

####
km1.1<- levelplot(KM_Sparse_rst[[4]], scales = list(draw = FALSE),  sub="k=5"  ,
                  maxpixels=ncell(KM_Sparse_rst[[4]]),col.regions = rev(terrain.colors(255)))  # "scales..." removes axes

km1.2 <- levelplot(ffp_smi_KM[[4]], scales = list(draw = FALSE), sub="k=5" ,
                   maxpixels=ncell(ffp_smi_KM[[4]]),col.regions = rev(terrain.colors(255)))

km2.1 <- levelplot(KM_Sparse_rst[[5]], scales = list(draw = FALSE), sub="k=6" ,
                   maxpixels=ncell(KM_Sparse_rst[[5]]),col.regions = rev(terrain.colors(255)))

km2.2 <- levelplot(ffp_smi_KM[[5]], scales = list(draw = FALSE), sub="k=6" ,
                   maxpixels=ncell(ffp_smi_KM[[5]]),col.regions = rev(terrain.colors(255)))

km3.1 <- levelplot(KM_Sparse_rst[[6]], scales = list(draw = FALSE), sub="k=7" ,
                   maxpixels=ncell(KM_Sparse_rst[[6]]),col.regions = rev(terrain.colors(255)))

km3.2 <- levelplot(ffp_smi_KM[[6]], scales = list(draw = FALSE), sub="k=7" ,
                   maxpixels=ncell(ffp_smi_KM[[6]]),col.regions = rev(terrain.colors(255)))

# Plot prints
print(km1.1, split = c(1, 1, 3, 2), more = TRUE)
print(km1.2, split = c(1, 2, 3, 2), more = TRUE)
print(km2.1, split = c(2, 1, 3, 2), more = TRUE)
print(km2.2, split = c(2, 2, 3, 2), more = TRUE)  
print(km3.1, split = c(3, 1, 3, 2), more = TRUE)
print(km3.2, split = c(3, 2, 3, 2))

#########
clustPerf_kmSparse <- data.frame(nClust = 2:7, SI_KM = NA, SDbw_KM = NA)
library(clusterCrit)
for(i in 1:6){ 
  cat(">>> Iteration : ",i)
  clust.vec <- as.integer(unlist(km.out.list[[i]][[1]][2]))
  crit <- intCriteria(traj = ras_cmpl_DF_aggr, 
              part = clust.vec, 
              crit = c("Silhouette","S_Dbw"))
  clustPerf_kmSparse[i, "SI_KM"]    <- unlist(crit)[[1]]
  clustPerf_kmSparse[i, "SDbw_KM"] <- unlist(crit)[[2]]
}
########
########DBSCAN
library(dbscan)  
snn <- sNNclust(ras_cmpl_DF_aggr, k=4,  minPts= 5, eps=5, borderPoints = TRUE)
nn <- sNN(ras_cmpl_DF_aggr, k = 5)

knn <- kNN(ras_cmpl_DF_aggr, k=15, sort = TRUE, search = "kdtree", bucketSize = 10,
    splitRule = "suggest", approx = 0)
snn <- sNNclust(ras_cmpl_DF_aggr, k=5,  minPts= 20, eps=0.4, borderPoints = TRUE)


cl <- dbscan(ras_cmpl_DF_aggr, eps = 0.55, minPts = 25)
pairs(ras_cmpl_DF_aggr[,c(1,5,18)], col = snn$cluster+1L)

cl_rst <- raster(ras_stack_aggr[[1]])
clSparse <- vector(mode = "integer", length = dim(ras_DF_aggr)[[1]])
clSparse[!idxa] <- NA
clSparse[idxa] <- cl$cluster 
values(cl_rst) <- clSparse
cl_cl<- levelplot(cl_rst, scales = list(draw = FALSE),  sub="DBCSAN"  ,
                  maxpixels=ncell(KM_Sparse_rst[[4]]),col.regions = rev(terrain.colors(255)))

print(cl_cl)



snn <- sNNclust(ras_cmpl_DF_aggr, k=5,  minPts= 30, eps=0.4, borderPoints = TRUE)


cl <- dbscan(ras_cmpl_DF_aggr, eps = 0.55, minPts = 20)
cl_rst <- raster(ras_stack_aggr[[1]])
clSparse <- vector(mode = "integer", length = dim(ras_DF_aggr)[[1]])
clSparse[!idxa] <- NA
clSparse[idxa] <- cl$cluster 
values(cl_rst) <- clSparse
cl_cl<- levelplot(cl_rst, scales = list(draw = FALSE),  sub="DBCSAN"  ,
                  maxpixels=ncell(KM_Sparse_rst[[4]]),col.regions = rev(terrain.colors(255)))
print(cl_cl)

cl <- dbscan(ras_cmpl_DF_aggr, eps = 0.60, minPts = 40)
cl_rst <- raster(ras_stack_aggr[[1]])
clSparse <- vector(mode = "integer", length = dim(ras_DF_aggr)[[1]])
clSparse[!idxa] <- NA
clSparse[idxa] <- cl$cluster 
values(cl_rst) <- clSparse
cl_cl<- levelplot(cl_rst, scales = list(draw = FALSE),  sub="DBCSAN"  ,
                  maxpixels=ncell(KM_Sparse_rst[[4]]),col.regions = rev(terrain.colors(255)))
print(cl_cl)




