#################################################################################
# Script:  Broughton_main_DFO.R - DFO CLASSIFICATION VERSION
# Created: February 2024. EJG
# 
# This script sources the necessary libraries and functions, coordinates the 
# analysis, and creates data structures that are then 'knitted' together (I think).
# So the idea is to run bits in here, and then 'Render' the RMD script. 
# Seems straightforward. :)
#
# Updates: 
# 2024/04/29: Steady development past few weeks; alll necessary pieces now developed. 
# 2024/04/29: Git repo created and working version pushed prior to RMarkdown development.
# 2024/05/02: Completed smoothing pass thru code; sorted raster plotting. Pushed.
# 2024/05/07: Another pass thru, adding some controls. Ready for RMD work.Pushed.
# 2024/05/29: After some exploratory work with Romina's tifs from SPECTRAL lab, forked the 
#   code so this version focuses on DFO data/objectives and Broughton2 attends to the LSSM needs.
# 2024/05/29: Back after summer. Nothing this WILL NOT use Romina's data.
# 2024/08/30: Back on this. Renamed the 2 projects for clarity, and updated here with LSSM progress. 
# 2024/09/10: Working now. Have spent days looking at distributions, outliers, and skew. 
#   Substrate has now joined bathy as a necessary characteristic. Almost ready to to start running 
#   some RMD reports and comparing results.
# 2024/09/11: A few minor(ish) changes: consolidate all changes to data (ie, transforming, centering, 
#   scaling) in one place. Bathy and Substrate applied as exclusions. 
# 2024/10/02: Update RMD based on findings from LSSM work. 

# To Do:
#   Find where we document what happened to FW index
#   Include the updated names (or maybe not?
#   See about gap statistics and fix the transformation table)

#################################################################################

print('Starting Classification - DFO Version ...')
rm(list=ls(all=T))  # Erase environment.

# Load necessary packages and functions ... 
source( "classification_functions.R" )
# source( "Plot_Functions.R" )

# Directories ...
#-- Source and output directories. Will be created if doesn't exist, overwritten if it does.
raster_dir <- 'C:/Data/SpaceData/Classification/MSEA'
data_dir   <- 'C:/Data/Git/classification-DFO/Data'
results_dir<- 'C:/Data/Git/classification-DFO/Results' 

# Processing FLAGS...
loadtifs <- F # If true the data will be re-loaded from TIFs, else it will be loaded from rData.
clipdata <- T # If true a spatial subset of the data will be taken based on a polygon shape file. 
reclust  <- T # If true, re-cluster full data set prior to mapping, else predict to unclassified pixels.
addKmDat <- T

#---- Part 1 of 3: Load, clean, and prepare predictor data.  ----
# If loadtifs == TRUE then run all this else load the processed data.

tif_stack <- stack()
today <- format(Sys.Date(), "%Y-%m-%d")

if (loadtifs) {
  print( "Loading predictors ... ")
  src_stack <- LoadPredictors( raster_dir, addKmDat )
  print( "Data loaded.")

  tif_stack <- src_stack
  
  # bathymetry is trimmed for landside and unsuitable depths for kelps.
  # Land and deep elevations are removed from the MSEA bathymetry
  # From Substrate produce hard only, marking kelp suitable areas. 
  # These restrictions manifest when completeCases() are selected.
  tif_stack <- DropNonHabitat( tif_stack, -5, 40 )

  if (clipdata) {
    print( "clipping TIFs ... ")
    amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_region.shp")
    #amask <- shapefile("C:\\Data\\SpaceData\\Broughton\\broughton_smaller.shp")
    tif_stack <- ClipPredictors( tif_stack, amask )
    print('Rasters clipped.')
  }

  save( tif_stack, file = paste0( data_dir, '/tif_stack_', today, '.rData' ))

} else {
  print( 'Loading project data ... ')
  # Ideally meaningfully named and tested so no source is required.
#  load( paste0( data_dir, '/tifs_DFO_scaled_QCS_2024-09-05.rData' ))
#  load( paste0( data_dir, '/tifs_DFO_centred_QCS_2024-09-05.rData' ))

 load( paste0( data_dir, '/tif_stack_2024-10-02.rData' ))
#  load( paste0( data_dir, '/t_stack_data_2024-09-14.rData' ))
}


#---- Final data preparation ----

# Quick correlation across data layers
#-- Move to matrix space from raster space 
x <- getValues( tif_stack )
#-- Use only pixels with all data
x_clean <- x[ complete.cases(x), ]

cor_table <- cor( x_clean )
cor_table[lower.tri(cor_table)] <- NA
(cor_table >= 0.6) & (cor_table != 1)
#--> cor_table prepared for printing in the RMD script.

#-- Remove correlated layers from raster stack
selected_stack <- dropLayer( tif_stack, c("bathymetry", "SUBSTRATE", "qcs_freshwater_index", "salinity_range", "temp_mean_summer", "circ_mean_summer",
                                    "sst_sentinel_20m_bi_mean", "sst_sentinel_20m_bi_sd") )
stack_data <- getValues( selected_stack )

#-- Visualize unmodified source raster data
dim( selected_stack )
names( selected_stack )
plot( selected_stack )
histogram(selected_stack, nclass=50)

# Prepare the data for classification.
if (scaledat) {
  
  print( "Transforming data  ... ")
  # REMOVE fix some (hard-coded) distributions by adding ceilings and root transforms.
  t_stack_data <- MakeMoreNormal( stack_data )
  
  print( "Centering and scaling  ... ")
  tmp_stack <- scale( t_stack_data, center=T,  scale=T )
  t_stack_data <- tmp_stack
  print('Data prepped.')
  save( t_stack_data, file = paste0( data_dir, '/t_stack_data', today, '.rData' ))
  print('Scaled data saved.')
}

### Histograms of unscaled and scaled vars in the RMD. 

# Compare pre-post skew? Put it on the histograms? :)
# RENAME variables after selection for plot prettiness.

new_names <- c("bathymetry", "substrate", "standard_dev_slope, arc-chord rugosity",
               "circ_mean_summer, tidal_mean_summer", "freshwater_index, salt_mean_summer, 
               salt_range", "rei", "sentinel_max, sentinel_mean, sentinel_sd, temp_mean_summer, temp_range")


# remove any rows with an NA
# NB: This decouples the data from the RasterStack and requires re-assembly for mapping
# THESE are the two key data structures used in subsequent steps
clean_idx <- complete.cases(t_stack_data)
stack_data_clean <- t_stack_data[ clean_idx, ]

dim( stack_data )
dim( stack_data_clean )


#---- Part 2 of 3: Cluster number selection ----

set.seed <- 42 # Seed for reproducibility
randomz  <- 20 # the number of randomizations for kmeans to do.
imax     <- 25 # maximum iterations to try for convergence

#---- Part 2a: Explore number of clusters using Within-sum-of-squares scree plot ----
# Runs kmeans with increasing number of clusters

nclust   <- 18 # number of clusters for scree plot
nsample  <- 50000 # Scree plot needs a subsample to run reasonably. 
plotme <- MakeScreePlot( stack_data_clean, nclust, randomz, imax, nsample )
plotme

#---- Create a working set of N clusters (N based on scree plot) to further assess cluster number. ----

nclust  <- 5 # the number of clusters based on scree plot, above.
nsample <- 500000 # a larger sample for more robust classification

sidx <- sample( 1:length( stack_data_clean[ , 1] ), nsample )
samp <- stack_data_clean[ sidx, ]
cluster_result <- kmeans(samp, centers=nclust, nstart=randomz, iter.max=imax) 

#---- Part 2b: Create heat map of within-cluster standard deviations ----

# Define color palette
pal_heat <- rev( brewer.pal(n = nclust, name = "RdYlBu")) # heat map palette

profile_data <- as.data.frame( cbind(cluster = cluster_result$cluster, samp ) )

cluster_sd <- profile_data %>%
  group_by(cluster) %>%
  summarise_all(sd)

x <- as.data.frame( cluster_sd )
head(x)
xm <- melt( x, id.var = "cluster" )

z_heat <- ggplot(xm, aes(x=cluster, y=variable, fill=value) ) +
  geom_tile() +
  scale_fill_gradientn(colours = pal_heat) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "Within-cluster Standard Deviation", x = "Clusters", y = "Attributes", fill = "Value")
z_heat

#---- Part 2c: Examine silhouette plot of the WORKING clusters  ----
# Uses the predictor values and the corresponding assigned cluster
# Need to subsample from the cluster result above as distance matrix take long time.

# Take a subsample of the clusters and the predictors for the silhouette plot. 
sil_n <- 10000
silx <- sample( 1:length( samp[ , 1] ), sil_n )

cs <- cluster_result$cluster[ silx ]
ss <- samp[ silx, ]

# Calculate a distance matrix between the predictors and plot assigned to each cluster.
# both these steps are time consuming, hence a smaller sample.
c_dist <- dist(ss)
sk <- silhouette(cs, c_dist)
#mean( sk[,"sil_width"] )

par(mfrow = c(1, 1))
plot(sk, col = 1:nclust, border=NA, main = "Hi World" )


#---- Part 3 of 3: Detailed examination of N clusters  ----
#---- Part 3a: Show cluster groupings using PCA ----

#-- Can take some time so it makes its own cluster 
pca_n <- 25000

#Returns a list of results (loadings table, and 2 plots)
pca_results <- ClusterPCA( pca_n, nclust ) # uses global variable stack_data_clean
names(pca_results) <- c("loadings", "plot1","plot2")

#Percentage of variance explained by dimensions
#eigenvalue <- round(get_eigenvalue(res_pca), 1)
#var_percent <- eigenvalue$variance.percent


#---- Part 3b: Violins of predictor contributions to WORKING clusters ----

x <- as.data.frame( samp )
x$cluster <- as.factor( cluster_result$cluster )

y <- x %>%
  pivot_longer(cols = -cluster, names_to = "predictor", values_to = "value")

# Create violin plot
vplots <- 
  ggplot(y, aes(x = cluster, y = value, fill = cluster)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ predictor, scales = "free_y") +
  theme_minimal() +
  labs(title = "Violin Plots of Predictors Across k-means Clusters",
       x = "Cluster",
       y = "Value")
vplots


#---- Part 4: Spatialize the WORKING clusters ----
# NB: To show a comprehensive map, can either:
#     a) re-cluster the entire data set (using imax and randomz from above) or
#     b) Predict to the unsampled portion of the raster. 

# Building a predictor raster takes time. Not useful to run this until 
# a comparison RMD report is being generated. 

out_tif_fname <- paste0( '/MSEA_gdata_5cluster_', today, '.tif' )

# initialize target data structure 
cluster_raster <- selected_stack[[1]]
dataType(cluster_raster) <- "INT1U"
cluster_raster[] <- NA

if (reclust == T) {
  # Re-cluster using all the clean data  ... 
  # less than 1 min with iter.max = 20, nstart = 20 for smallest region
  cluster_result <- kmeans(stack_data_clean, centers = nclust, nstart = randomz, iter.max = imax)
  
  # Assign the clustered values ... 
  # extract values from the target cluster
  new_values <- values( cluster_raster )
  # replace non-NA values with the cluster results
  new_values[ clean_idx ] <- cluster_result$cluster
  # put the updated values back on the target cluster
  values( cluster_raster ) <- new_values  
} else {
  # Predict values for unclustered cells. Can be more time-consuming than re-classifying everything. 
  values( cluster_raster ) <- transferCluster( values(cluster_raster), cluster_result )
}

#--- Display the results, first as histogram then as map.  
raster::hist( values(cluster_raster ))

# Define color palette
pal_clust <- brewer.pal(n = nclust, "Accent") # Max for Accent is 8

ckey <- list( at=0:nclust, 
              title="Clusters", 
              col = pal_clust)
myTheme <- rasterTheme( region = pal_clust )
z_map <- levelplot( cluster_raster, margin = F, 
           colorkey = ckey,
           par.settings = myTheme,
           main = "K-means Clustering Result - Local extents" )
z_map

writeRaster( cluster_raster, paste0( results_dir, out_tif_fname ), overwrite=TRUE)


#---- Knit and render Markdown file -----
# 2024/04/29: It looks like this has gotten easier in the last 2 years ... version up!


### Process file 
# To HTML ... 
# rmarkdown::render( "Broughton_DFO.Rmd",   
#                    output_format = 'html_document',
#                    output_file = paste0( "C:/Data/Git/Broughton/Results/DFO_Class_Report_", today ) )

# To PDF:
# the tinytex library is necessary for compiling the .tex file to be rendered.
# then run >tinytex::install_tinytex()
rmarkdown::render( "Classification_DFO_PDF.Rmd",   
                   output_format = 'pdf_document',
                   output_dir ="C:/Data/Git/Classification-DFO/Results",
                   output_file = paste0( "MSEA_gdata_5cluster_", today ))

#---- Some details on correlation analysis ... ----
#---- Correlation across UN-scaled datalayers ----
# foo <- getValues( scaled_layers )
# foo_clean <- na.omit(stack_data)
# pick <- sample( 1:length( foo_clean[ , 1] ), 10000 )
# plot( foo_clean[ pick,'rugosity'] ~ foo_clean[ pick,'standard_deviation_slope'] )

# RUGOSITY is a bit of a problem distribution
#cellStats( data_layers$rugosity, stat="range" )
#raster::hist(log( data_layers$rugosity+10 ), nclass=50)
# Look at bottom roughness relationships  
#pick <- sample( 1:length( stack_data_clean[ , 1] ), 10000 )
#plot( stack_data_clean[ pick,'rugosity'] ~ stack_data_clean[ pick,'standard_deviation_slope'] )

#---- 

#---- Plot any classified outliers in classified space --
# create a df w the centers for each cluster to facilitate distance calculation.
centers <- cluster_result$centers[cluster_result$cluster, ] 
head(centers)
# calculate distance
distances <- sqrt(rowSums((stack_data_clean - centers )^2))
head(x)
outlier_idx <- order(distances, decreasing=T)[1:1000000]
# a subsample to reduce the plot size ... 
ssidx <- sample( 1:length( stack_data_clean[ , 1] ), 10000 )

plot(stack_data_clean[ ssidx, c("rei_qcs", "qcs_freshwater_index") ], pch=19, col=cluster_result$cluster[ ssidx ], cex=1)
points(cluster_result$centers[ ,c("rei_qcs", "qcs_freshwater_index") ], col=1:3, pch=15, cex=2)
#----


# FIN.






