### Frozen flowSOM ###

# Apply this script on reference merged NBM sample
# Analyse all generated MSTs and chose the best one as FROZEN flowSOM MST

########################################################

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(FlowSOM)

# Output folder
output.folder <- "D:/school/Stage officieel/frozen_flowSOM/"
dir.create(output.folder)
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S")

#########
# Start #
#########
# read FCS file into flowframe
ff <- read.FCS("D:/school/Stage officieel/DATA/normal bone marrows to normalize/csv2fcs/merged/15_NBM.fcs")

# store AnnotatedDataFrame of original fcs
adf <- parameters(ff)

# Logicle transformation
ff.trf <- transform(ff,transformList(colnames(ff)[6:15],logicleTransform(w = 1, t = 1048576)))

# Density plot to check transformation
channels <- c(6:15) # FL channels
markers <- c("CD58", "CD81", "CD34", "CD22", "CD38", "CD10", "CD19", "CD5", "CD20", "CD45") # B-ALL tube markers
windows(title = "Logicle transformation")
par(mfrow= c(4,3), pty="s")
nb.plots.per.windows <- 10
aa <- 0
for (FL in channels) {
  aa <- aa+1
  d1 <- density(exprs(ff.trf[,FL]))
  plot(d1, xlab = colnames(ff.trf[,FL]), xlim=c(min(exprs(ff.trf[,FL])),max(exprs(ff.trf[,FL]))), main= markers[aa])
}

#  normalize with mean = 0, SD = 1
normsd <- function(x) {
  return ((x - mean(x))/sd(x))
}
FL.norm <- normsd(exprs(ff.trf)[,6:15])
# replace the expression matrix with the normalized values
exprs(ff.trf)[,6:15] <- norm

### preprocessing finished ###

# Assign the FL columns to be used for the flowSOM analysis
parameters.to.use <- c(7:11,14:15)
# Assign parameters for flowSOM
x.dim <- 11
y.dim <- 11
nb.clusters.fsom <- 5
dispersion.coef <- 50  #default=30 how big you want the nodes
fl.names <- ff.trf@parameters@data$desc # for legend plots

# Make 24 MST's of NBM, analyse them and chose the best one as FROZEN flowSOM for further analysis with samples
for (seed in 1:24) {
  # Assign name for file with settings for flowSOM
  file.name.fixed.flowSOM <- paste0(output.folder, date.format, "_Frozen_", seed, ".Rdata")
  # set seed
  set.seed(seed)
  # Use wrapper function to perform flowSOM
  fSOM.res <- FlowSOM(ff.trf, 
                  compensate = FALSE,
                  transform = FALSE,
                  scale = FALSE,
                  scaled.center = FALSE,
                  scaled.scale = FALSE,
                  colsToUse = parameters.to.use,
                  xdim = x.dim,
                  ydim = y.dim,
                  nClus = 100,
                  rlen = 10)
  save(fSOM.res, file = file.name.fixed.flowSOM)
  fSOM <- fSOM.res[[1]]
  # Save coördinates of nodes and add it to flowframe
  nodes.mapping <- fSOM$map$mapping[,1] # value between 1 - 100 (--> which node does the event belong to)
  nb.cells <- nrow(ff.trf)
  node.fSOM <- vector()
  node.fSOM.x <- NULL
  node.fSOM.y <- NULL
  layout.fcs <- fSOM$MST$l # x y coordinates per node??
  # plot(layout.fcs, main= "Frozen flowSOM")
  mapping.fcs <- fSOM$map$mapping[,1]
  vertex.size.fcs <- fSOM$MST$size # size of every node
  nb.clusters <- (fSOM$map$xdim) * (fSOM$map$ydim)
  nb.nodes <- nb.clusters
  nb.points.per.line <- nb.cells/nb.nodes
  for (i in 1:nb.clusters) {
    tempo <- mapping.fcs[mapping.fcs==i]
    tempo_indices <- which(mapping.fcs==i)
    
    random.nb <- runif(n=tempo,min = 0, max = 1 )
    
    random.radius <-  random.nb * vertex.size.fcs[i]/2
    random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
    
    node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
    node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
  }
  # plot(node.fSOM.x, node.fSOM.y, main= "Frozen flowSOM + size")
  
  matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
  matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=10, xmax=1014, ymin=10, ymax=1014)
  colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.1","FlowSOM.2")
  
  ff.fsom <- cbind2 (ff, matrix.node.fSOM.scaled)
  
  # get metaclustering data for fcs file
  metacluster <- fSOM.res [[2]]
  metacluster <- metaClustering_consensus(fSOM$map$codes, k = 6)
  data.metacluster <- metacluster[fSOM$map$mapping[,1]] # = metaClustering_perCell <- metaClustering[fSOM$map$mapping[,1]]
  data.metacluster <- as.matrix(data.metacluster)
  colnames(data.metacluster) <- "Metaclustering Consensus"
  
  ff.meta <-  cbind2(ff.fsom, data.metacluster)
  
  write.FCS(ff.meta, filename = paste0(output.folder, date.format, "_Frozen_flowSOM_", seed, ".fcs"))
  
  # Figures
  pdf.name <- paste0(output.folder, date.format, "_frozen_flowSOM_", seed, ".pdf")
  pdf(file = pdf.name,
      width = 8.3, 
      height = 11.7)
  par(mfrow=c(1,1))
  par(cex.main=1.4)
  fSOM$prettyColnames <- as.character(fl.names)
  PlotStars(fSOM)
  # Plot metaclusters
  PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metacluster), 
            backgroundColor = grDevices::colorRampPalette(rainbow(6)),
            main = "MST with metaclusters")
  # Plot metaclusters
  PlotStars(fSOM, view= "grid", backgroundValues = as.factor(metacluster), 
            backgroundColor = grDevices::colorRampPalette(rainbow(6)),
            main = "MST with metaclusters")
  dev.off()
}

####### without wrapper function ##########
# # Read the data into a flowSOM object
# fSOM <- ReadInput(ff.trf, 
#                   compensate = FALSE, 
#                   transform = FALSE, 
#                   scale = FALSE,
#                   scaled.center = FALSE,
#                   scaled.scale = FALSE)
# 
# # Build the self-organizing map --> adds the information of the map in the flowSOM object with $map
# fSOM <- BuildSOM(fSOM, colsToUse = parameters.to.use)
# 
# # Build the minimal spanning tree --> adds it to flowSOM object with $MST
# fSOM <- BuildMST(fSOM, tSNE = TRUE)
############################################

# Now the flowSOM object can be used for visualization. The noded can be plotted in several layouts: 
# MST = minimal spanning tree = default
# grid = grid of the SOM
# tSNE = alternative layout (only possible when tSNE parameter is set to TRUE in BuildMST)
# PlotStars(fSOM, view = "grid")
# PlotStars(fSOM, view = "tSNE")
# PlotStars(fSOM, view = "MST")
# 
# # reset node size if you don't want the size to depend on the number of cells assigned to a node
# fSOM.reset <- UpdateNodeSize(fSOM, reset = TRUE)
# fSOM.reset$MST$size <- fSOM.reset$MST$size/2
# PlotStars(fSOM.reset)
# 
# # Plot specific marker
# print(colnames(fSOM$map$medianValues))
# PlotMarker(fSOM, "FL3 INT") # CD34
# PlotMarker(fSOM, "FL6 INT") # CD10
# # If you need to refer to the nodes it might be useful to number them
# PlotNumbers(UpdateNodeSize(fSOM, reset = TRUE))
# # This number can be used for a 2D scatter plot indicating the node values
# PlotClusters2D(fSOM, "FL1 INT", "FL2 INT", c(31))
# 
# # Perform metaclustering of the data for further analysis. Gives good approximation of manual gating results.
# # If you have background knowledge about the number of cell types you are looking for, 
# # it might be good to provide the number to the algorithm.
# # all 100 nodes are assigned to one of the 7 clusters
# metaClustering <- metaClustering_consensus(fSOM$map$codes, k=7) # k = number of clusters
# cell_types <- c(1:7)
# metaClustering_perCell <- metaClustering[fSOM$map$mapping[,1]]
# # Plot metaclusters
# PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metaClustering), 
#           backgroundColor = grDevices::colorRampPalette(rainbow(7)),
#           main = "MST with metaclusters")
# 
# 
# # If no manual gating to map: possible to query the tree to indicate nodes similar to a specified pattern
# # e.g. look for CD34+, CD10+ cells
# query <- c("FL3 INT" = "high", "FL6 INT" = "high")
# query_res <- QueryStarPlot(UpdateNodeSize(fSOM, reset=TRUE), query, plot = FALSE)      
# cellTypes <- factor(rep("Unknown", 49), levels=c("Unknown", "CD34+ CD10+"))
# cellTypes[query_res$selected] <- "CD34+ CD10+"
# PlotStars(fSOM, backgroundValues = cellTypes, backgroundColor = c("#FFFFFF00", "#0000FF22"))