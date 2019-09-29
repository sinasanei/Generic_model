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

#Reading in raster clips of PA
# get all file names with .tif extension
files <- list.files(path = "PA_Clips/",pattern = "*.tif$")
rasname <- rep(NA,length(files)) # names of raster objects
for (i in 1:length(files)){
  assign(strsplit(files[i], split = "_" )[[1]][1], raster(paste("PA_Clips/",files[i],sep = "")))
  rasname[i] <- strsplit(files[i], split = "_" )[[1]][1]
}
rasname <- c(rasname,"smi")
smi<- gsdd5/ gsp
# Plot assignments
pd100 <- levelplot(d100, scales = list(draw = FALSE),  sub="d100"  ,
                   maxpixels=ncell(mat),col.regions = rev(terrain.colors(255)))  # "scales..." removes axes
pffp <- levelplot(ffp, scales = list(draw = FALSE), sub="ffp" ,
                  maxpixels=ncell(mat),col.regions = rev(terrain.colors(255)))
pmat <- levelplot(mat, scales = list(draw = FALSE), sub="mat" ,
                  maxpixels=ncell(mat),col.regions = rev(terrain.colors(255)))
pmap <- levelplot(map, scales = list(draw = FALSE), sub="map" ,
                  maxpixels=ncell(mat),col.regions = rev(terrain.colors(255)))
pmtcm <- levelplot(mtcm, scales = list(draw = FALSE), sub="mtcm" ,
                   maxpixels=ncell(mat),col.regions = rev(terrain.colors(255)))
psmrsprpb <- levelplot(smrsprpb, scales = list(draw = FALSE), sub="smrsprpb" ,
                       maxpixels=ncell(mat),col.regions = rev(terrain.colors(255)))

# Plot prints
print(pd100, split = c(1, 1, 3, 2), more = TRUE)
print(pffp, split = c(2, 1, 3, 2), more = TRUE)
print(pmat, split = c(3, 1, 3, 2), more = TRUE)
print(pmap, split = c(1, 2, 3, 2), more = TRUE)  
print(pmtcm, split = c(2, 2, 3, 2), more = TRUE)
print(psmrsprpb, split = c(3, 2, 3, 2))

sum(is.na(values(d100)))
sum(is.na(values(map)))

########
########
#######
#######

layers_pa <- stack(d100, dd0)
ras_mat <- as.data.frame(matrix(data = rep(NA,length(rasname)*length(rasname)), nrow= length(rasname)))
for (i in 1:length(rasname)){
  for (j in 1:length(rasname)){
    ras_mat[i,j] <- compareRaster(get(rasname[i]),get(rasname[j]),stopiffalse=FALSE)
  }
}


#ex <- extent(dd0)
#extent(d100)<- ex
#d100 <- projectRaster(d100,ffp,method='bilinear')

#put all climate variables in a Raster stack 
ras_stack <- stack(get(rasname[1]),get(rasname[2]))

for( i in 3:length(rasname)){
  ras_stack <- stack(ras_stack,get(rasname[i]))
}
names(ras_stack) <- rasname
#extract all values from the raster into a data frame
ras_DF <- values(ras_stack)
colnames(ras_DF)<- rasname
#check for NA's in the data
idx <- complete.cases(ras_DF)
ras_cmpl_DF <- ras_DF[idx,]
apply(ras_cmpl_DF,2, function(x) sum(is.na(x))) #number of missing values in each column
apply(ras_cmpl_DF, 2,mean)
apply(ras_cmpl_DF, 2,sd)
ras_cmpl_DF <- scale(ras_cmpl_DF, scale = TRUE)
######################################################
######################################################
write.csv(ras_DF, file = "raster_df.csv")
writeRaster(ras_stack,  filename='raster_stack.tif',overwrite=TRUE)
           
######################################################
######################################################
library(FactoMineR)
library(factoextra)
ras_pc <- PCA(ras_cmpl_DF  , scale.unit = FALSE, ncp = 5, graph = FALSE)
# Contributions of variables to PC1
fviz_contrib(ras_pc, choice = "var", axes = 1, top = 15)
# Contributions of variables to PC2
fviz_contrib(ras_pc, choice = "var", axes = 2, top = 10)
# Contributions of variables to PC3
fviz_contrib(ras_pc, choice = "var", axes = 3, top = 10)
fviz_contrib(ras_pc, choice = "var", axes = 4, top = 10)


summary(ras_pc)

ras_pc <- PCA(ras_cmpl_DF  , scale.unit = FALSE, ncp = 2, graph = FALSE)
df <- ras_pc$ind[[1]]
apply(df, 2,mean)
apply(df, 2,sd)

colnames(df) <- c("PC1", "PC2"  )


write.csv(df, file = "df.csv")
PC_DF <- matrix(data = NA, nrow = dim(ras_DF)[[1]], ncol = 2)
PC_DF[,1]<- ifelse(idx, df[,1], NA)
PC_DF[,2]<- ifelse(idx, df[,2], NA)

######################################################
############Principle Component Clustering #############
######################################################

library(cluster)
library(clusterCrit)
PC_KM <- raster(ras_stack[[1]])
PC_CLARA <- raster(ras_stack[[1]])

for(nClust in 2:12){
  
  cat("-> Clustering data for nClust =",nClust,"......")
  # Perform K-means clustering
  km <- kmeans(df, centers = nClust, iter.max = 50 , nstart = 25)
  
  # Perform CLARA's clustering (using manhattan distance)
  cla <- clara(df, k = nClust,stand = FALSE, samples = 25, metric = "manhattan")
  
  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = dim(ras_DF)[[1]])
  claClust <- vector(mode = "integer", length = dim(ras_DF)[[1]])
  
  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster
  
  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering
  
  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(ras_stack[[1]])
  # CLARA
  tmpRstCLARA <- raster(ras_stack[[1]])
  
  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust
  
  # Stack the temporary rasters onto the final ones
  if(nClust==2){
    PC_KM    <- tmpRstKM
    PC_CLARA <- tmpRstCLARA
  }else{
    PC_KM    <- stack(PC_KM, tmpRstKM)
    PC_CLARA <- stack(PC_CLARA, tmpRstCLARA)
  }
  
  cat(" done!\n\n")
}

# Write the clustering solutions for each algorithm
writeRaster(PC_KM,"PC_KM", overwrite=TRUE)
writeRaster(PC_CLARA,"PC_CLARA", overwrite=TRUE)

#The silhouette index ranges from ???1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. If most objects have a high value, then the clustering configuration is considered appropriate. If many points have a low or negative value, then the clustering configuration may have too many or too few clusters.

# Start a data frame that will store all silhouette values
# for k-means and CLARA   
clustPerfSI <- data.frame(nClust = 2:12, SI_KM = NA, SI_CLARA = NA)

set.seed(17)
for(i in 1:nlayers(PC_KM)){ # Iterate through each layer
  
  cat("-> Evaluating clustering performance for nClust =",(2:12)[i],"......")
  
  # Extract random cell samples stratified by cluster
  cellIdx_PC_KM <- sampleStratified(PC_KM[[i]], size = 5000)
  cellIdx_PC_CLARA <- sampleStratified(PC_CLARA[[i]], size = 5000)
  
  # Get cell values from the Stratified Random Sample from the raster 
  # data frame object (PC_DF)
  DFStRS_KM <- PC_DF[cellIdx_PC_KM[,1], ]
  DFStRS_CLARA <- PC_DF[cellIdx_PC_CLARA[,1], ]
  
  # Make sure all columns are numeric (intCriteria function is picky on this)
  DFStRS_KM[] <- sapply(DFStRS_KM, as.numeric)
  DFStRS_CLARA[] <- sapply(DFStRS_CLARA, as.numeric)
  
  # Compute the sample-based Silhouette index for: 
  #    
  # K-means
  clCritKM <- intCriteria(traj = DFStRS_KM, 
                          part = as.integer(cellIdx_PC_KM[,2]), 
                          crit = "S_Dbw")
  # and CLARA
  clCritCLARA <- intCriteria(traj = DFStRS_CLARA, 
                             part = as.integer(cellIdx_PC_CLARA[,2]), 
                             crit = "S_Dbw")
  
  # Write the silhouette index value to clustPerfSI data frame holding 
  # all results
  clustPerfSI[i, "SI_KM"]    <- clCritKM[[1]][1]
  clustPerfSI[i, "SI_CLARA"] <- clCritCLARA[[1]][1]
  
  cat(" done!\n\n")
  
}

write.csv(clustPerfSI, file = "cluster_info.csv", row.names = FALSE)

# print a table with teh silhouttee index results to compare each clustering solution
knitr::kable(clustPerfSI, digits = 3, align = "c", 
             col.names = c("#clusters","Avg. Silhouette (k-means)","Avg. Silhouette (CLARA)"))

# make a plot to comare the two algorithms
plot(clustPerfSI[,1], clustPerfSI[,2], 
     xlim = c(1,13), ylim = range(clustPerfSI[,2:3]), type = "n", 
     ylab="Avg. Silhouette Index", xlab="# of clusters",
     main="Silhouette index by # of clusters")

# Plot Avg Silhouette values across # of clusters for K-means
lines(clustPerfSI[,1], clustPerfSI[,2], col="red")
# Plot Avg Silhouette values across # of clusters for CLARA
lines(clustPerfSI[,1], clustPerfSI[,3], col="blue")

# Grid lines
abline(v = 1:13, lty=2, col="light grey")
abline(h = seq(0.52,0.62,0.02), lty=2, col="light grey")

legend("topright", legend=c("K-means","CLARA"), col=c("red","blue"), lty=1, lwd=1) # you want to select the number of clusters that has the highest Avg. silhouette index

#plot the stacked raster with the recommended number of clusters
plot(PC_KM[[1]])
levelplot(PC_KM[[3]], scales = list(draw = FALSE),
          col.regions = rev(terrain.colors(26)))

plot(ffp_smi_KM[[1]])         

plot(PC_KM[[8]])
plot(ffp_smi_KM[[8]])    

####
####
km1.1<- levelplot(PC_KM[[4]], scales = list(draw = FALSE),  sub="km5"  ,
                   maxpixels=ncell(PC_KM[[4]]),col.regions = rev(terrain.colors(255)))  # "scales..." removes axes

km1.2 <- levelplot(ffp_smi_KM[[4]], scales = list(draw = FALSE), sub="km5" ,
                  maxpixels=ncell(ffp_smi_KM[[4]]),col.regions = rev(terrain.colors(255)))

km2.1 <- levelplot(PC_KM[[5]], scales = list(draw = FALSE), sub="km6" ,
                  maxpixels=ncell(PC_KM[[5]]),col.regions = rev(terrain.colors(255)))

km2.2 <- levelplot(ffp_smi_KM[[5]], scales = list(draw = FALSE), sub="km6" ,
                  maxpixels=ncell(ffp_smi_KM[[5]]),col.regions = rev(terrain.colors(255)))

km3.1 <- levelplot(PC_KM[[6]], scales = list(draw = FALSE), sub="km7" ,
                   maxpixels=ncell(PC_KM[[6]]),col.regions = rev(terrain.colors(255)))

km3.2 <- levelplot(ffp_smi_KM[[6]], scales = list(draw = FALSE), sub="km7" ,
                       maxpixels=ncell(ffp_smi_KM[[6]]),col.regions = rev(terrain.colors(255)))

# Plot prints
print(km1.1, split = c(1, 1, 3, 2), more = TRUE)
print(km1.2, split = c(1, 2, 3, 2), more = TRUE)
print(km2.1, split = c(2, 1, 3, 2), more = TRUE)
print(km2.2, split = c(2, 2, 3, 2), more = TRUE)  
print(km3.1, split = c(3, 1, 3, 2), more = TRUE)
print(km3.2, split = c(3, 2, 3, 2))

######################################################
############ mtcm-gsp Clustering #############
######################################################
#extract all values from the raster into a data frame
library(dbscan)  
ras_DF <- read.csv("raster_df.csv")
ras_DF <- ras_DF[,-1]
ras_stack <- raster("raster_stack.tif")

ffp_smi <- ras_DF[,c(5,21)]
idx <- complete.cases(ffp_smi)
ffp_smi_comp <- ffp_smi[idx,]


cl <- hdbscan(ffp_smi_comp, minPts = 5)
cl
snn <- dbscan(ffp_smi_comp, eps=0.99,  minPts= 45, borderPoints = TRUE)
cl_DB <- vector(mode = "integer", length = dim(ffp_smiDF)[[1]])
cl_DB[!idx2] <- NA
cl_DB[idx2] <- snn$cluster 
DB_rst <- raster(ras_stack[[1]])
values(DB_rst) <- cl_DB
cl_cl<- levelplot(DB_rst, scales = list(draw = FALSE),  sub="DBCSAN"  ,
                  maxpixels=ncell(ras_stack[[1]]),col.regions = rev(terrain.colors(255)))
print(cl_cl)

kNNdistplot(ffp_smi_BD, k = 5000)

idex3 <- sample(1:dim(ffp_smi_BD)[[1]], 5000 , replace = FALSE)
ffp_smiBD <- ffp_smi_BD[idex3,]

kNNdistplot(ffp_smiBD, k = 500)

snn2 <- dbscan(ffp_smiBD, eps=1.001,  minPts= 50, borderPoints = TRUE)
unique(snn2$cluster)

cl_DB[idx2] <- predict(snn2, ffp_smi_BD, data = ffp_smiBD)

