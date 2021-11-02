library(spanner)
library(dbscan)
library(lidR)
library(RANN)

# set the number of threads to use in lidR
set_lidr_threads(8)

# download and read an example laz
LASfile = "D:/LEARN_2020_1_1_1.laz"
las = readTLSLAS(LASfile, select = "xyzc", "-filter_with_voxel 0.0254 -keep_circle 0, 0, 30")
# Don't forget to make sure the las object has a projection
# projection(las) = sp::CRS("+init=epsg:26912")

# pre-process the example lidar dataset by classifying the ground points
# using lidR::csf(), normalizing it, and removing outlier points 
# using lidR::ivf()
las = classify_ground(las, csf(sloop_smooth = FALSE, 
                               class_threshold = 0.5,
                               cloth_resolution = 0.1, rigidness = 1L, 
                               iterations = 500L, time_step = 0.65))
# las = normalize_height(las, tin())
dtm = grid_terrain(las, 0.33, kriging(model = gstat::vgm(0.59, "Sph", 874), k = 10L))
las = normalize_height(las, dtm)


las = classify_noise(las, ivf(0.5, 3))
las = filter_poi(las, Classification != LASNOISE)

eigen = spanner::eigen_metrics(las, radius=0.1, ncpu=8)

las@data<-cbind(las@data,eigen)
for (j in colnames(eigen)) las = add_lasattribute(las, eigen[[j]], name = j, desc = j)

message(paste("Taking the initial las slice...\n"))
las_slice <- filter_poi(las, Z >= 0.5 & Z < 2 & Classification !=2)

message(paste("Building the clusters using dbscan...\n"))
minPts = ceiling(log(nrow(las_slice@data)))
eps = mean(RANN::nn2(data=las_slice@data[,1:3], query=las_slice@data[,1:3], k=30)$nn.dists)

las_slice@data = na.omit(las_slice@data)

db_clusters<-dbscan::dbscan(las_slice@data[, c("X","Y","Z","PCA1","PCA2","Verticality")],
                            eps = 0.1, minPts = 500)
# ln(n) is used in Section 4.1 of https://www.sciencedirect.com/science/article/pii/S0169023X06000218
las_slice@data$treeID <- db_clusters$cluster

plot(filter_poi(las_slice, treeID > 0 & Verticality > 0.90), color="treeID")

library(dplyr)
myTrees <- las_slice@data %>% count(treeID)
myTrees <- myTrees[myTrees$n>5000 & myTrees$treeID!=0]

las_dbh <- filter_poi(las_slice, treeID > 0 & Verticality > 0.90)
las_dbh <- filter_poi(las_dbh, Z>=1.27, Z<=1.47)

pb <- utils::txtProgressBar(min = 0, max = length(unique(las_dbh$treeID)), style = 3)
i <- 0
dbh = data.table()
for(d in unique(myTrees$treeID)) {
  temp = filter_poi(las_dbh, treeID == d)
  
  if(nrow(temp@data)>=minPts)
  {
    tryCatch(
      {
        # message(paste("'Trying' to process tree ",d))
        est = TreeLS:::cylinderFit(temp, method = 'ransac', n=13, inliers=.8, conf=.95, max_angle=45)
        est$TreeID = d
      },
      error=function(cond) {
        message(paste("Tree", d, " caused an error!"))
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        est<-c(rep(NA,9),d)
        names(est)<-c("rho","theta","phi","alpha","radius","err","px", "py", "pz","TreeID")
      },
      warning=function(cond) {
        message(paste("Tree",d," caused a warning!"))
        message("Here's the original warning message:")
        message(cond)
        # Choose a return value in case of warning
        est = TreeLS:::cylinderFit(temp, method = 'ransac', n=13, inliers=.8, conf=.95, max_angle=45)
        est$TreeID = d
      },
      finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        if(!is.null(est)) dbh = rbind(dbh,est)
        # message(paste("Processed Tree:", d))
      }
    )}
  i <- i + 1
  utils::setTxtProgressBar(pb, i)
}

dbh <- dbh[dbh$radius*100 < 100 & dbh$radius*100 >= 0.1 & dbh$err*100 < 0.2,]

###################################################################################
circleFun <- function(center = c(0,0), radius = 11.34, npoints = 100){
  r = radius
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
cir_dat=circleFun(c(0,0), radius = 11.34, npoints = 100)

chm = grid_canopy(filter_poi(las, Z>=0.5), res = 0.1, p2r())
test_spdf <- as(chm, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

library(ggplot2)
ggplot(test_df) +
  geom_tile(data = test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  scale_fill_gradientn(colours = height.colors(10), name="Canopy\nHeight") +
  geom_point(data = dbh, aes(x=px, y=py, size = radius*2*100, color=err), alpha=0.7) +
  geom_text(data = dbh, aes(x=px, y=py, label = round(radius*2*100,1)), size = 2, hjust = 0, nudge_x = 0.5, color="black") +
  labs(size = "Diameter (cm)") + scale_size(range = c(0, 3))+
  scale_color_gradientn(colours = c("black", "tomato2", "firebrick"), name="Ransac\nError") +
  # scale_color_gradientn(colours = RColorBrewer::brewer.pal(10, 'RdYlBu'), name="Ransac\nError") +
  geom_path(data=cir_dat, aes(x,y), linetype = 2) +
  ylab("y-coord (m)") +
  xlab("x-coord (m)") +
  # theme(legend.position = "none") +
  coord_fixed(ratio=1, xlim=c(min(test_df$x)-2, max(test_df$x)+2), ylim=c(min(test_df$y)-2, max(test_df$y)+2))


spdf<-SpatialPointsDataFrame(coords=dbh[,c("px","py")], data=dbh)[,c("px","py","pz","radius", "TreeID")]
names(spdf)<-c("X","Y","Z","radius","treeID")
x <- plot(filter_poi(las, Z>=1.27, Z<=1.47), color="Verticality", bg = "white", size = 1)
add_treetops3d(x, spdf, color="black", radius=spdf$radius)

