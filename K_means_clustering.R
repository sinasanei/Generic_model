# Trial to conduct a k-means analysis of the different combinations of the climate data layers

#ffp_smi

#set working directory
##setwd("~/Climate Maps/Combination Trials")

#library
library(raster)
library(rgdal)
library(RStoolbox)
library(cluster)
library(clusterCrit)


#input raster layers
#ffp <- raster("ffp.tif")
#mmax <- raster("mmax.tif")
#mtcm <- raster("mtcm.tif")
#smi <- raster ("smi.tif")

#view and compare the structure of the rasters
ffp
mmax
mtcm
smi<- gsdd5/ gsp

#plot a few to see what they look like and check to make sure that the resolutions are the same
plot(ffp)
plot(mmax)
crs(mmax)
crs(ffp)

#view a histogram of the data
hist(ffp, maxpixels=ncell(ffp)) # normally distributed
hist(mmax, maxpixels=ncell(mmax)) # skewed to the right
hist(mtcm, maxpixels=ncell(mtcm)) # kind of normally distributed
hist(smi, maxpixels=ncell(smi)) # normally distributed

ffp_smi <- stack(ffp, smi)

#extract all values from the raster into a data frame
ffp_smiDF <- values(ffp_smi)
apply(ffp_smiDF,2, function(x) sum(is.na(x))) #number of missing values in each column
identical(is.na(ffp_smiDF[,1]),is.na(ffp_smiDF[,2]))
#check for NA's in the data
idx2 <- complete.cases(ffp_smiDF)

# Initiate the raster datasets that will hold all clustering solutions 
# from 2 groups/clusters up to 12
ffp_smi_KM <- raster(ffp_smi[[1]])
ffp_smi_CLARA <- raster(ffp_smi[[1]])

for(nClust in 2:12){
  
  cat("-> Clustering data for nClust =",nClust,"......")
  # Perform K-means clustering
  km <- kmeans(ffp_smiDF[idx2,], centers = nClust, iter.max = 50)
  
  # Perform CLARA's clustering (using manhattan distance)
  cla <- clara(ffp_smiDF[idx2, ], k = nClust, metric = "manhattan")
  
  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(ffp_smi))
  claClust <- vector(mode = "integer", length = ncell(ffp_smi))
  
  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx2] <- NA
  kmClust[idx2] <- km$cluster
  
  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx2] <- NA
  claClust[idx2] <- cla$clustering
  
  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(ffp_smi[[1]])
  # CLARA
  tmpRstCLARA <- raster(ffp_smi[[1]])
  
  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust
  
  # Stack the temporary rasters onto the final ones
  if(nClust==2){
    ffp_smi_KM    <- tmpRstKM
    ffp_smi_CLARA <- tmpRstCLARA
  }else{
    ffp_smi_KM    <- stack(ffp_smi_KM, tmpRstKM)
    ffp_smi_CLARA <- stack(ffp_smi_CLARA, tmpRstCLARA)
  }
  
  cat(" done!\n\n")
}

# Write the clustering solutions for each algorithm
writeRaster(ffp_smi_KM,"ffp_smi_KM", overwrite=TRUE)
writeRaster(ffp_smi_CLARA,"ffp_smi_CLARA", overwrite=TRUE)

# evaluate the unsupervised classification/clustering performance
#Use the silhouette index to evaluate the performance of each clustering solution and select the "best" number of clusters for partitioning the sample data
# this will give you a measure of how similar an object is to its own cluster (intra-cluster cohesion) compared to other clusters (between-cluster separation)

#The silhouette index ranges from ???1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. If most objects have a high value, then the clustering configuration is considered appropriate. If many points have a low or negative value, then the clustering configuration may have too many or too few clusters.

# Start a data frame that will store all silhouette values
# for k-means and CLARA   
clustPerfSI2 <- data.frame(nClust = 2:12, SI_KM = NA, SI_CLARA = NA)


for(i in 1:nlayers(ffp_smi_KM)){ # Iterate through each layer
  
  cat("-> Evaluating clustering performance for nClust =",(2:12)[i],"......")
  
  # Extract random cell samples stratified by cluster
  cellIdx_ffp_smi_KM <- sampleStratified(ffp_smi_KM[[i]], size = 2000)
  cellIdx_ffp_smi_CLARA <- sampleStratified(ffp_smi_CLARA[[i]], size = 2000)
  
  # Get cell values from the Stratified Random Sample from the raster 
  # data frame object (ffp_smiDF)
  ffp_smi_DFStRS_KM <- ffp_smiDF[cellIdx_ffp_smi_KM[,1], ]
  ffp_smi_DFStRS_CLARA <- ffp_smiDF[cellIdx_ffp_smi_CLARA[,1], ]
  
  # Make sure all columns are numeric (intCriteria function is picky on this)
  ffp_smi_DFStRS_KM[] <- sapply(ffp_smi_DFStRS_KM, as.numeric)
  ffp_smi_DFStRS_CLARA[] <- sapply(ffp_smi_DFStRS_CLARA, as.numeric)
  
  # Compute the sample-based Silhouette index for: 
  #    
  # K-means
  clCritKM <- intCriteria(traj = ffp_smi_DFStRS_KM, 
                          part = as.integer(cellIdx_ffp_smi_KM[,2]), 
                          crit = "Silhouette")
  # and CLARA
  clCritCLARA <- intCriteria(traj = ffp_smi_DFStRS_CLARA, 
                             part = as.integer(cellIdx_ffp_smi_CLARA[,2]), 
                             crit = "Silhouette")
  
  # Write the silhouette index value to clustPerfSI data frame holding 
  # all results
  clustPerfSI2[i, "SI_KM"]    <- clCritKM[[1]][1]
  clustPerfSI2[i, "SI_CLARA"] <- clCritCLARA[[1]][1]
  
  cat(" done!\n\n")
  
}

write.csv(clustPerfSI, file = "cluster_info.csv", row.names = FALSE)

# print a table with teh silhouttee index results to compare each clustering solution
knitr::kable(clustPerfSI2, digits = 3, align = "c", 
             col.names = c("#clusters","Avg. Silhouette (k-means)","Avg. Silhouette (CLARA)"))

# make a plot to comare the two algorithms
plot(clustPerfSI2[,1], clustPerfSI2[,2], 
     xlim = c(1,13), ylim = range(clustPerfSI2[,2:3]), type = "n", 
     ylab="Avg. Silhouette Index", xlab="# of clusters",
     main="Silhouette index by # of clusters")

# Plot Avg Silhouette values across # of clusters for K-means
lines(clustPerfSI2[,1], clustPerfSI2[,2], col="red")
# Plot Avg Silhouette values across # of clusters for CLARA
lines(clustPerfSI2[,1], clustPerfSI2[,3], col="blue")

# Grid lines
abline(v = 1:13, lty=2, col="light grey")
abline(h = seq(0.52,0.62,0.02), lty=2, col="light grey")

legend("topright", legend=c("K-means","CLARA"), col=c("red","blue"), lty=1, lwd=1) # you want to select the number of clusters that has the highest Avg. silhouette index

#plot the stacked raster with the recommended number of clusters
plot(ffp_smi_KM[[3]])
levelplot(ffp_smi_KM[[3]], scales = list(draw = FALSE),
          col.regions = rev(terrain.colors(26)))
          
plot(ffp_smi_KM[[1]])         
     
     
          
           