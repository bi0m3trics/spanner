#' Segment a terrestrial point cloud using graph theory.
#'
#' `segment_graph` returns a .las object with a new column "treeID".
#'
#' Preforms Individual tree segmentation following ecological principles for “growing” trees 
#' based on these input locations in a graph-theory approach inspired by work of Tao and 
#' others (2015). Point coordinates are linked together based on proximity and turned into 
#' a connected graph object, using the estimated tree bole locations as origin points, connecting individual
#' points back to those tree bole origins based on shortest paths within the graph network, and finally
#' assigning those points a unique tree identification based on the bole coordinate for which 
#' they are connected. Input point cloud is subsampled to a lower resolution before processing to 
#' increase processing efficiency. 
#'
#' @param las LAS normalized las object.
#' @param tree.locations data.frame A data.frame contained the following seed information: `TreeID`,
#' `X`, `Y`, `Z`, and `Radius` in the same units as the .las, output from the raster_eigen function.
#' @param k integer Number of nearest neighbors to be used in processing (k >= 50 suggested)
#' @param distance.threshold numeric Maximum distance (in the same units as the .las) under which 
#' two points are connected in the graph object (greater than point spacing). Two points with a greater
#' distance than this threshold are not connected in the graph for processing. 
#' @param plot_radius numeric Radius (in the same units as the .las) used to define a region of 
#' interest for processing
#' @param use.metabolic.scale bool Use of weights in the assignment of points to a given treeID. Useful
#' when interlocking crowns are present and trees are of different sizes. 
#' @param subsample.graph numeric The subsampled point spacing to use during processing. Note: processing 
#' time increases quickly with smaller point spacing with negligible returns in accuracy. 
#' @param return.dense bool Decision to return the subsampled point cloud or assign treeIDs back to 
#' points in the input dense point cloud. 
#' @param output_location character Where to save processing outputs (rasters). 
#' 
#' @return a .las with the column treeID added.
#' 
#' @examples
#' \dontrun{
#' # set the number of threads to use in lidR
#' set_lidr_threads(8)
#' 
#' # read the laz
#' LASfile <- system.file("extdata", "Pine_Example.laz", package="spanner")
#' las = readLAS(LASfile, select = "xyzc")
#' 
#' # plot(las, color="Z", backend="lidRviewer", trim=30)
#' 
#' # find tree locations and attribute data
#' myTreeLocs = get_raster_eigen_treelocs(las = las, res = 0.05, pt_spacing = 0.02, 
#'                                        dens_threshold = 0.25, 
#'                                        neigh_sizes=c(0.1, 0.1, 0.25), 
#'                                        eigen_threshold = 0.25, 
#'                                        grid_slice_min = 0.666, 
#'                                        grid_slice_max = 2,
#'                                        minimum_polygon_area = 0.01, 
#'                                       cylinder_fit_type = "ransac", 
#'                                        output_location = "E:/GithubRepos/", 
#'                                        max_dia = 0.5, 
#'                                        SDvert = 0.25)
#' 
#' # segment the point cloud 
#' myTreeGraph <- segment_graph(las = las, tree.locations = myTreeLocs, k = 50, 
#'                              distance.threshold = 0.5,
#'                              use.metabolic.scale = FALSE, subsample.graph = 0.1, 
#'                              return.dense = FALSE,
#'                              output_location = "E:/GithubRepos/")
#' 
#' plot(myTreeGraph, color="treeID",  backend="lidRviewer")
#' }
#' 
#' @export
segment_graph <- function(las = las, tree.locations = NULL, k = NULL, distance.threshold = 0.38,
                          use.metabolic.scale = FALSE, subsample.graph = 0.1,
                          return.dense = FALSE, output_location = getwd()){

  ## pcd
  segmented <- las
  
  ## if subsampling, do that now
  if(!is.null(subsample.graph)){
    segmented <- tlsSample(segmented, smp.voxelize(spacing = subsample.graph))
  }
  ## get rid of ground points
  message("Filtering out ground points... (step 1 of 7)\n")
  segmented <- lidR::filter_poi(segmented, Classification != 2, Z > (1/3))
  
  ## how many tree objects
  ntree <- nrow(tree.locations)
  
  ## Give tree locations a uniform Z coordinate
  tree.locations$Z <- 1.3
  tree.locations$Radius[is.na(tree.locations$Radius)] = 0.01
  ## Make sure tree location is close to bole points (if in center, distance threshold could cut it off)
  ## Eventually, this should be a value tied to the actual DBH of the tree
  for(d in 1:nrow(tree.locations)){
    if(tree.locations$Radius[d] > distance.threshold){
      if(tree.locations$Y[d] > mean(tree.locations$Y)){
        tree.locations$Y <- tree.locations$Y - ((tree.locations$Radius)/2)
      } else {
        tree.locations$Y <- tree.locations$Y + ((tree.locations$Radius)/2)
      }
    }}
  
  ## ---- Determining metabolic scaling factors for individual trees ----
  ## This is the metabolic scaling factor which will get multiplied by the weights so that smaller trees
  ## don't "suck up" points from nearby larger trees met_scale <- 1 / ( tree.locations$Radius^(2/3) )
  ## it seems like the equation in the original manuscript doesn't work that well for me... this can be modified
  met_scale <- (1 / ( tree.locations$Radius))^(3/2) ## it seems like the equation in the original manuscript doesn't
                                                    ## work that well for me... this can be modified
  
  ##---------------------- Identify neighbors  -------------------------------------
  ## create the network by connecting lidar points within X meters
  ## and add in the points for the tree centers to make connections
  message("Identifing neighbors... (2/7)\n")
  point_density <- density(segmented)
  
  # point_dist_threshold <- round((1/sqrt(point_density))*k, 1) ## in meters, the distance within which points
                                                                ## will be connected into a network
  # point_dist_threshold <- (1/3)
  # point_dist_threshold <- 0.42
  
  ## k is arbitrary... not sure if this should be based on a calculation, e.g. ...
  # k = ((sqrt(point_density))/(1/3))
  
  ## in meters, the distance within which points will be connected into a network
  k_tree <- FNN::get.knn(rbind(segmented@data[ ,c("X","Y","Z")], tree.locations[ ,c("X","Y","Z")]), k = k,  algorithm="cover_tree")
  
  ##---------------------- Build the Graph OBJ  -------------------------------------
  ## matrix of from, to, weight (dist)
  message("Building the graph object... (3/7)\n")
  weight_matrix = cbind(rep(1:dim(k_tree$nn.index)[1],each=k),
                        as.vector(t(k_tree$nn.index)),
                        as.vector(t(k_tree$nn.dist)))
  # which(weight_matrix[,3] > point_dist_threshold)
  # hist(weight_matrix[which(weight_matrix[,3] > point_dist_threshold),3])
  
  ## get rid of rows in the matrix that have a distance greater than the point_dist_threshold
  weight_matrix <- weight_matrix[-(which(weight_matrix[,3] > distance.threshold)), ]
  
  ## actually make the graph using this matrix
  graph_obj <- cppRouting::makegraph(weight_matrix,
                                     directed = FALSE) ## a new package
  
  ## ---- IDs for the origin tree locations ----
  treeloc_ids <- (nrow(segmented@data)+1):(nrow(k_tree$nn.index))
  treeloc_ids <- treeloc_ids[treeloc_ids %in% c(weight_matrix[,1],weight_matrix[,2])]
  
  ## ---- IDs for the origin tree locations ----
  remaining_ids <- unique(c(weight_matrix[,1],weight_matrix[,2])[!(c(weight_matrix[,1],weight_matrix[,2]) %in% treeloc_ids)])
  
  ## ---- IDs not used in graph ----
  unused_ids <- (1:(nrow(segmented@data)))[!(1:(nrow(segmented@data)) %in% c(weight_matrix[,1],weight_matrix[,2]))]
  
  ##-------------- Calculate the shortest paths from all tree origins  -----------------------------
  message("Calculating the shortest path from tree... (4/7)\n")
  shortpath <- cppRouting::get_distance_matrix(graph_obj, from = treeloc_ids, to = remaining_ids, allcores = FALSE)
  # shortpath[ ,1:100] ## testing
  
  # shortpath <- cppRouting::get_multi_paths(graph_obj, from = treeloc_ids, to = remaining_ids, long = TRUE)
  
  ##-------------- Determine which tree ID is linked to the shortest path  -----------------------------
  message("Linking tree IDs... (5/7)\n")
  ## Use the metabolic scaling factors to modify the distances for each treeID
  if(use.metabolic.scale == TRUE){
    shortpath <- shortpath * met_scale
  }
  
  ## add this so you don't drop columns in the colMins function
  shortpath <- rbind(shortpath, 999)
  
  ## Classify the vector of treeIDs (we need all of the point ids back together)
  treeID <- c( Rfast::colMins(shortpath), rep(0, length(unused_ids)) )
  ## give the values the proper point ID name
  names(treeID) <- c(colnames(shortpath), unused_ids)
  ## now get the treeID values back in order to attach to the las object
  treeID <- treeID[order(as.numeric(names(treeID)))]
  ## change the "extra" treeID to 0
  treeID[treeID > ntree] <- 0
  
  ## -------------------------- add the treeID data to the las object -----------------------------------
  message("Assigning treeIDs to the las object... (6/7)\n")
  segmented <- lidR::add_lasattribute(segmented, (treeID), name = "treeID", desc = "tree id")
  # segmented <- add_lasattribute(segmented, (as.numeric(names(treeID))), name = "name", desc = "tree name")
  
  ## based on whether the user wants the dense point cloud returned
  if(return.dense == FALSE){
    
    return(segmented)
    
  } else {
    ## get the indices of nearest sparse points
    message("Processing dense cloud... (6b/7)\n")
    s2d_lookup <- RANN::nn2(data=data.frame(X=segmented@data$X,Y=segmented@data$Y,Z=segmented@data$Z),
                            query=data.frame(X=las@data$X,Y=las@data$Y,Z=las@data$Z),
                            k=1)
    ## grab the treeID of the sparse point
    las@data$treeID<-segmented@data$treeID[s2d_lookup[[1]]]
    las@data$nnd<-s2d_lookup[[2]]
    
    ## assign the point if within subsampled sparse point density, otherwise 0
    las@data$treeID[las@data$nnd>(subsample.graph * 1.25)]<-0
    las@data$treeID[las@data$Classification==2]<-0
    
    # plot(las, color="treeID", backend="pcv")
    
    ## add the treeID values as attribute the dense point cloud
    las <- lidR::add_lasattribute(las, (las@data$treeID), name = "treeID", desc = "tree id")
    
    return(las)
  }
  message("Done! (7/7)\n")
}