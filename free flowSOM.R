### Perform free flowSOM on fcs file ###

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(FlowSOM)
library(rstudioapi)

#########
# Start #
#########
# read FCS file into flowframe
ff <- read.FCS(filename = choose.files(caption = "Select fcs file", multi = FALSE))

# store AnnotatedDataFrame of original fcs
adf <- parameters(ff)

# Logicle transformation
# get index for channels to transform (FL)
FL.index <- grep("FL", colnames(ff))
ff.trf <- transform(ff,transformList(colnames(ff)[FL.index],logicleTransform(w = 0.9, t = 1048576)))

### Density plot to check transformation
# save (pretty: CD) channelnames of flowframe in variable
FL.names <- as.character(parameters(ff)$desc[FL.index])
windows(title = "Logicle transformation")
par(mfrow= c(4,3), pty="s")
nb.plots.per.windows <- 10
aa <- 0
for (FL in FL.index) {
  aa <- aa+1
  d1 <- density(exprs(ff.trf[,FL]))
  plot(d1, xlab = colnames(ff.trf[,FL]), xlim=c(-0.5,4.5), main= FL.names[aa])
}

#  normalize with mean = 0, SD = 1
normsd <- function(x) {
  return ((x - mean(x))/sd(x))
}
FL.norm <- normsd(exprs(ff.trf)[,FL.index])
# replace the expression matrix with the normalized values
exprs(ff.trf)[,FL.index] <- FL.norm

### preprocessing finished ###

# Assign the FL columns to be used for the flowSOM analysis
parameters.to.use <- c(7:11,14:15)
# Assign parameters for flowSOM
x.dim <- 11
y.dim <- 11
nb.metacluster <- 5
dispersion.coef <- 50  #default=30 how big you want the nodes
channel.names <- ff.trf@parameters@data$desc # for legend plots

set.seed(1234)

fSOM.res <- FlowSOM(ff.trf, 
                    compensate = FALSE,
                    transform = FALSE,
                    scale = FALSE,
                    scaled.center = FALSE,
                    scaled.scale = FALSE,
                    colsToUse = parameters.to.use,
                    xdim = x.dim,
                    ydim = y.dim,
                    nClus = nb.metacluster <- 5, # use maxMeta?
                    rlen = 10)
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
colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.X","FlowSOM.Y")

ff.fsom <- cbind2 (ff, matrix.node.fSOM.scaled)

# get metaclustering data for fcs file
metacluster <- fSOM.res [[2]]
metacluster <- metaClustering_consensus(fSOM$map$codes, k = nb.metacluster <- 5)
data.metacluster <- metacluster[fSOM$map$mapping[,1]] # = metaClustering_perCell <- metaClustering[fSOM$map$mapping[,1]]
data.metacluster <- as.matrix(data.metacluster)
colnames(data.metacluster) <- "Metaclustering Consensus"

ff.meta <-  cbind2(ff.fsom, data.metacluster)

# Add range keyword to description slot and ranges to annotated datafram
for (i in c((ncol(ff)+1):ncol(ff.meta))){
  # adjust rowname of adf
  r.name <- paste0("$P",i)
  rownames(ff.meta@parameters@data)[i] <- r.name
  # adjust desc of adf
  desc <- ff.meta@parameters@data$name[i]
  ff.meta@parameters@data$desc[i] <- desc
  # Add $PnN to the desc slot
  name.kw <- paste0("$P",i,"N")
  ff.meta@description[name.kw] <- ff.meta@parameters@data$name[i]
  desc.kw <- paste0("$P",i,"S")
  ff.meta@description[desc.kw] <- ff.meta@parameters@data$desc[i]
  # Adjust range of adf and add $PnR keyword to the desc slot
  range.kw <- paste0("$P",i,"R")
  if (ff.meta@parameters@data$name[i]=="Metaclustering Consensus"){
    ff.meta@parameters@data$range[i] <- ff.meta@parameters@data$maxRange[i] + 1
    ff.meta@parameters@data$maxRange[i] <- ff.meta@parameters@data$maxRange[i] + 1
    ff.meta@parameters@data$minRange[i] <- 0
    ff.meta@description[range.kw] <- ff.meta@parameters@data$range[i]
  } else {
    ff.meta@parameters@data$range[i] <- 1100
    ff.meta@parameters@data$maxRange[i] <- 1100
    ff.meta@parameters@data$minRange[i] <- 0
    ff.meta@description[range.kw] <- ff.meta@parameters@data$range[i]
  }
}

output.file <- selectFile(caption = "Save file as ...", label = "Save", existing = FALSE)

write.FCS(ff.meta, filename = paste0(output.file, ".fcs"))

# Figures
pdf.name <- paste0(output.file, ".pdf")
pdf(file = pdf.name,
    width = 8.3, 
    height = 11.7)
par(mfrow=c(1,1))
par(cex.main=1.4)
fSOM$prettyColnames <- as.character(channel.names)
PlotStars(fSOM)
# Plot metaclusters
PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metacluster), 
          backgroundColor = grDevices::colorRampPalette(rainbow(5)),
          main = "MST with metaclusters")
# Plot metaclusters
PlotStars(fSOM, view= "grid", backgroundValues = as.factor(metacluster), 
          backgroundColor = grDevices::colorRampPalette(rainbow(5)),
          main = "MST with metaclusters")
# Query for CD34+ CD10+ population (MRD)
# query <- c("FL3 INT" = "high", "FL6 INT" = "high")
# query_res <- QueryStarPlot(UpdateNodeSize(fSOM, reset=TRUE), query, plot = FALSE)
# cellTypes <- factor(rep("Others", 49), levels=c("Others", "CD34+ CD10+"))
# cellTypes[query_res$selected] <- "CD34+ CD10+"
# PlotStars(fSOM, backgroundValues = cellTypes, backgroundColor = c("#FFFFFF00", "#0000FF22"))
dev.off()
