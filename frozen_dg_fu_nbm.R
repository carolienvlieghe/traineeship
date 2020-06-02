#######################################
###Frozen flowSOM with tagged files ###
#######################################

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(FlowSOM)
library(rstudioapi)

# Assign in- & output
merged.fcs <- choose.files(multi = FALSE, caption = "Select merged FCS file to process")

# Vector with dataset names
title.frozenfsom <- c("Dg", "FU","merged NBM")

# dispersion coefficient
dispersion.coef <- 50  #default=30 how big you want the nodes

# ## From original script
# # Define here the TIME parameter number - needed to be able to target each fcs files
# # Why? files are tagged?
# TIME.parameter <- 16
# coef.gap <- 20 #maximum time distance between 2 cells in a file
# #Node of interrest extraction settings
# #target NOI based on CD38 vs CD10
# CD45.parameter <- 5 #CD10
# SSC.parameter <- 10 #CD38
# DIAG.pct.threshold <- 1 #min DIAG cell percentage threshold
# Ratio.threshold <- 10 #min ratio DIAG/NBM threshold
# DIAG.cell.nb.threshold <- 10 #min DIAG cell number threshold

# ADC resolution of the cytometer
adc.resolution.Area <- 20
# adc.resolution.Height <- 18  
# adc.resolution.Width <- 10    
# adc.resolution.TIME <- 20    
# adc.resolution.N.parameter <- 10

# Cytometer name
cytometer.name <- "Navios BeD"

# name of the Rdata file containing the NBM flowSOM result
name.file.flowSOM.fixed.res <- choose.files(caption = "Upload Rdata file", multi = FALSE)

# Enter here the max scale for FSC and SSC of the cytometer
max.scale.FSC.SSC <- 2^20

# Read FCS file
ff <- read.FCS(filename = merged.fcs)

# get index for channels to transform (FL)
FL.index <- grep("FL", colnames(ff))
TAG.index <- grep("TAG", colnames(ff))

# Identify the datasets in the merged file through the tags
DIAG_indices <- which(exprs(ff)[,TAG.index]==250)
FU_indices <- which(exprs(ff)[,TAG.index]==500)
NBM_indices <- which(exprs(ff)[,TAG.index]==750)

# Logicle transformation of flowframe
ff.trf <- transform(ff,transformList(colnames(ff)[FL.index],logicleTransform(w = 0.9, t = 1048576)))
# Check transformation
# save (pretty: CD) channelnames of flowframe in variable
FL.names <- as.character(parameters(ff)$desc[FL.index])
windows(title = "Logicle transformation")
par(mfrow= c(4,3), pty="s")
nb.plots.per.windows <- 10
aa <- 0
for (FL in FL.index) {
  aa <- aa+1
  d1 <- density(exprs(ff.trf[,FL]))
  plot(d1, xlab = colnames(ff.trf[,FL]), xlim=c(min(exprs(ff.trf[,FL])),max(exprs(ff.trf[,FL]))), main= FL.names[aa])
}

#  normalize with mean = 0, SD = 1
normsd <- function(x) {
  return ((x - mean(x))/sd(x))
}
FL.norm <- normsd(exprs(ff.trf)[,FL.index])
# replace the expression matrix with the normalized values
exprs(ff.trf)[,FL.index] <- FL.norm

######################
### Frozen flowSOM ###
######################

load(file = name.file.flowSOM.fixed.res) # Load the Rdata of the flowSOM result of the reference NBM (= fSOM.res)

#reference MST plot
fSOM.fixed.NBM <- fSOM.res[[1]]

# Plot The MST with node numbers
# PlotNumbers(fSOM.fixed.NBM, view = "MST")

# Built the MST tree from the flowSOM result of the reference fSOM and get the coördinates (--> can be added to flowFrame)

# vectors to save coördinates
node.fSOM <- vector()
node.fSOM.x <- NULL
node.fSOM.y <- NULL
# number of events
nb.cells <- nrow(ff.trf[NBM_indices])
# x/y coordinates per node
layout.fcs <- fSOM.fixed.NBM$MST$l
# which node does an event belong to: node number per event
mapping.fcs <- fSOM.fixed.NBM$map$mapping[,1]
# Size of every node
vertex.size.fcs <- fSOM.fixed.NBM$MST$size
# number of clusters
nb.cluster <- (fSOM.fixed.NBM$map$xdim) * (fSOM.fixed.NBM$map$ydim)
nb.nodes <- nb.cluster
nb.points.per.line <- nb.cells /nb.nodes
for (i in 1:nb.cluster) {
  tempo <- mapping.fcs[mapping.fcs==i] # list with node numbers repeated the number of times it is present in the mapping
  tempo_indices <- which(mapping.fcs==i) # indices of the events that belong to a certain node number
  
  random.nb <- runif(n=tempo,min = 0, max = 1 ) # creates random numbers between 0 and 1 for every event in with certain node number
  
  random.radius <-  random.nb * vertex.size.fcs[i]/2
  random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
  
  node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
  node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
}

# Start PDF file to save plots
windows()
par(mfrow = c(1,1))
plot(node.fSOM.x, node.fSOM.y, main = "Frozen MST layout")
points(layout.fcs, col="red")

# Determine plot boundaries, in units of the data
xmin.frozen <- par("usr")[1]
xmax.frozen <- par("usr")[2]
ymin.frozen <- par("usr")[3]
ymax.frozen <- par("usr")[4]

### scaling node.fSOM.x, node.fSOM.y, node.fSOM.density and node.fSOM.cluster on 2^20 bit #####
# Is this necessary???
node.fSOM.x.min <- xmin.frozen
node.fSOM.x.max <- xmax.frozen
node.fSOM.y.min <- ymin.frozen
node.fSOM.y.max <- ymax.frozen

node.fSOM.x.coef <- ((2^adc.resolution.Area)-1)/(node.fSOM.x.max-node.fSOM.x.min)
node.fSOM.y.coef <- ((2^adc.resolution.Area)-1)/(node.fSOM.y.max-node.fSOM.y.min)

node.fSOM.x.ind.term <- 1- (node.fSOM.x.coef * node.fSOM.x.min)
node.fSOM.y.ind.term <- 1- (node.fSOM.y.coef * node.fSOM.y.min)

node.fSOM.x.scaled <- (node.fSOM.x.coef * node.fSOM.x) + node.fSOM.x.ind.term
node.fSOM.y.scaled <- (node.fSOM.y.coef * node.fSOM.y) + node.fSOM.y.ind.term

# windows()
# par(mfrow=c(1,1))
# plot(node.fSOM.x.scaled,node.fSOM.y.scaled, main="Rescaled MST layout + Vertex size")

# Make separate flowframes of each data set without TAG (why???)
ff.dg.tag <- ff.trf[DIAG_indices,]
ff.dg <- ff.dg.tag[,1:(ncol(ff.dg.tag)-1)]
ff.FU.tag <- ff.trf[FU_indices,]
ff.FU <- ff.FU.tag[,1:(ncol(ff.FU.tag)-1)]
ff.NBM.tag <- ff.trf[NBM_indices,]
ff.NBM <- ff.NBM.tag[,1:(ncol(ff.NBM.tag)-1)]

# Map new data to an existing flowSOM object
fSOM.dg.frozen <- NewData(fSOM.fixed.NBM, ff.dg)
fSOM.FU.frozen <- NewData(fSOM.fixed.NBM, ff.FU)
fSOM.NBM.frozen <- NewData(fSOM.fixed.NBM, ff.NBM)

ff.nodes <- list()
ff.meta <- list()
ff.list <- list(ff.dg, ff.FU, ff.NBM)
fSOM.list <- list(fSOM.dg.frozen, fSOM.FU.frozen, fSOM.NBM.frozen)
indices.list <- list(DIAG_indices, FU_indices, NBM_indices)

# windows()
# n2mfrow(3)
# par(mfrow = c(3,1), pty="s")
# par(mar=c(5.1,4.1,4.1,2.1))
# par(mar=c(1,1,1,1))

for (m in 1:3) {
  nb.cells <- nrow(ff.list[[m]]) # number of events in dataset
  # Built MST with coordinates
  node.fSOM <- vector()
  node.fSOM.x <- NULL
  node.fSOM.y <- NULL
  layout.fcs <- fSOM.list[[m]]$MST$l
  # not necessary
  # plot(layout.fcs, cex.main=0.8, 
       # main= paste(title.frozenfsom[m],"-FROZEN on NBM MST layout", sep = "")) # creates 3 identical MST's
  
  xmin.frozen.NBM <- par("usr")[1]
  xmax.frozen.NBM <- par("usr")[2]
  ymin.frozen.NBM <- par("usr")[3]
  ymax.frozen.NBM <- par("usr")[4]
  
  mapping.fcs <- fSOM.list[[m]]$map$mapping[,1]
  vertex.size.fcs <- fSOM.list[[m]]$MST$size
  nb.cluster <- (fSOM.list[[m]]$map$xdim) * (fSOM.list[[m]]$map$ydim)
  nb.nodes <- nb.cluster
  nb.points.per.line <- nb.cells /nb.nodes
  
  for (i in 1:nb.cluster) {
    tempo <- mapping.fcs[mapping.fcs==i]
    tempo_indices <- which(mapping.fcs==i)
    random.nb <- runif(n=tempo,min = 0, max = 1 )
    random.radius <-  random.nb * vertex.size.fcs[i]/2
    random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
    
    node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
    node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
  }
  
  # Plot the nodes of the dataset onto the original MST
  # plot(node.fSOM.x,node.fSOM.y, cex.main = 0.8,
  #      main=paste(title.frozenfsom[m],"-FROZEN on NBM MST layout + Vertex size",sep = ""))
  # points(layout.fcs, col="red") 
  
  node.fSOM.x.min <- xmin.frozen
  node.fSOM.x.max <- xmax.frozen
  node.fSOM.y.min <- ymin.frozen
  node.fSOM.y.max <- ymax.frozen
  
  node.fSOM.x.coef <- ((2^adc.resolution.Area)-1)/(node.fSOM.x.max-node.fSOM.x.min)
  node.fSOM.y.coef <- ((2^adc.resolution.Area)-1)/(node.fSOM.y.max-node.fSOM.y.min)
  
  node.fSOM.x.ind.term <- 1- (node.fSOM.x.coef * node.fSOM.x.min)
  node.fSOM.y.ind.term <- 1- (node.fSOM.y.coef * node.fSOM.y.min)
  
  node.fSOM.x.scaled <- (node.fSOM.x.coef * node.fSOM.x) + node.fSOM.x.ind.term
  node.fSOM.y.scaled <- (node.fSOM.y.coef * node.fSOM.y) + node.fSOM.y.ind.term
  
  # plot(node.fSOM.x.scaled,node.fSOM.y.scaled, 
  #      cex.main = 0.8,
  #      main=paste(title.frozenfsom[m],"nodes on frozen MST layout",sep = " "))
  # points(layout.fcs, col="red") 
  
  # metacluster
  metacluster <- metaClustering_consensus(fSOM.list[[m]]$map$codes, k = 5)
  data.metacluster<- metacluster[fSOM.list[[m]]$map$mapping[,1]]
  
  ff.nodes[[m]] <- cbind(node.fSOM.x.scaled, node.fSOM.y.scaled, data.metacluster)
  colnames(ff.nodes[[m]]) <- c("FROZEN flowSOM X", "Frozen flowSOM Y", "Metaclustering Consensus")
  
  ff.new <- ff[indices.list[[m]],]
  ff.new <- fr_append_cols(ff.new, ff.nodes[[m]])
  mtrx <- exprs(ff.new)
  if (m==1){
    merged.mtrx <- mtrx
  } else {
    merged.mtrx <- rbind(merged.mtrx, mtrx)
  }
}

exprs(ff.new) <- merged.mtrx
ff.frozen <- ff.new

# Add range keyword to description slot and ranges to annotated datafram
for (i in c((ncol(ff)+1):ncol(ff.frozen))){
  # adjust rowname of adf
  r.name <- paste0("$P",i)
  rownames(ff.frozen@parameters@data)[i] <- r.name
  # adjust desc of adf
  desc <- ff.frozen@parameters@data$name[i]
  ff.frozen@parameters@data$desc[i] <- desc
  # Add $PnN to the desc slot
  name.kw <- paste0("$P",i,"N")
  ff.frozen@description[name.kw] <- ff.frozen@parameters@data$name[i]
  desc.kw <- paste0("$P",i,"S")
  ff.frozen@description[desc.kw] <- ff.frozen@parameters@data$desc[i]
  # Adjust range of adf and add $PnR keyword to the desc slot
  range.kw <- paste0("$P",i,"R")
  if (ff.frozen@parameters@data$name[i]=="Metaclustering Consensus"){
    ff.frozen@parameters@data$range[i] <- ff.frozen@parameters@data$maxRange[i] + 1
    ff.frozen@parameters@data$maxRange[i] <- ff.frozen@parameters@data$maxRange[i] + 1
    ff.frozen@parameters@data$minRange[i] <- 0
    ff.frozen@description[range.kw] <- ff.frozen@parameters@data$range[i]
  } else {
    ff.frozen@parameters@data$range[i] <- 1100
    ff.frozen@parameters@data$maxRange[i] <- 1100
    ff.frozen@parameters@data$minRange[i] <- 0
    ff.frozen@description[range.kw] <- ff.frozen@parameters@data$range[i]
  }
}


output.file <- selectFile(caption = "Save file as ...", label = "Save", existing = FALSE)
write.FCS(ff.frozen, filename = paste0(output.file, ".fcs"))

# Figures
channel.names <- ff.trf@parameters@data$desc # for legend plots
pdf.name <- paste0(output.file, ".pdf")
pdf(file = pdf.name,
    width = 8.3, 
    height = 11.7)
par(mfrow=c(1,1))
par(cex.main=1.4)
fSOM.res$prettyColnames <- as.character(channel.names)
for (i in 1:3) {
  PlotStars(fSOM.list[[i]], view = "MST", backgroundValues = as.factor(metacluster), 
            backgroundColor = grDevices::colorRampPalette(rainbow(5)),
            main = title.frozenfsom[[i]])
}
dev.off()

# # Plot metaclusters
# PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metacluster), 
#           backgroundColor = grDevices::colorRampPalette(rainbow(5)),
#           main = "MST with metaclusters")
# # Plot metaclusters
# PlotStars(fSOM, view= "grid", backgroundValues = as.factor(metacluster), 
#           backgroundColor = grDevices::colorRampPalette(rainbow(5)),
#           main = "MST with metaclusters")
# # Query for CD34+ CD10+ population (MRD)
# query <- c("FL3 INT" = "high", "FL6 INT" = "high")
# query_res <- QueryStarPlot(UpdateNodeSize(fSOM, reset=TRUE), query, plot = FALSE)
# cellTypes <- factor(rep("Others", 49), levels=c("Others", "CD34+ CD10+"))
# cellTypes[query_res$selected] <- "CD34+ CD10+"
# PlotStars(fSOM, backgroundValues = cellTypes, backgroundColor = c("#FFFFFF00", "#0000FF22"))
# dev.off()


# PlotGroups()
# remove variable vector list matrix dataFrame from RAM memory
# remove(name.file.flowSOM.fixed.res, fSOM.res, fSOM.fixed.NBM, ff.nodes, ff.dg.tag, ff.dg, ff.FU.tag, ff.FU, ff.NBM.tag , ff.NBM,
# fSOM1.FROZEN, fSOM2.FROZEN, fSOM.fixed.NBM,ff.list, fSOM.list)

#### END FIXED FLOWSOM ###




####################
### Free flowSOM ###
####################

# # Assign the FL columns to be used for the flowSOM analysis
# parameters.to.use <- c(7:11,14:15)
# # Assign parameters for flowSOM
# x.dim <- 11
# y.dim <- 11
# nb.clusters.fsom <- 5
# fl.names <- ff.trf@parameters@data$desc # for legend plots
# 
# set.seed(123)
# 
# fSOM.res.free <- FlowSOM(ff.trf, 
#                     compensate = FALSE,
#                     transform = FALSE,
#                     scale = FALSE,
#                     scaled.center = FALSE,
#                     scaled.scale = FALSE,
#                     colsToUse = parameters.to.use,
#                     xdim = x.dim,
#                     ydim = y.dim,
#                     nClus = 100,
#                     rlen = 10)
# fSOM <- fSOM.res.free[[1]]
# # Save coördinates of nodes and add it to flowframe
# nodes.mapping <- fSOM$map$mapping[,1] # value between 1 - 100 (--> which node does the event belong to)
# nb.cells <- nrow(ff.trf)
# node.fSOM <- vector()
# node.fSOM.x <- NULL
# node.fSOM.y <- NULL
# layout.fcs <- fSOM$MST$l # x y coordinates per node??
# # plot(layout.fcs, main= "Frozen flowSOM")
# mapping.fcs <- fSOM$map$mapping[,1]
# vertex.size.fcs <- fSOM$MST$size # size of every node
# nb.clusters <- (fSOM$map$xdim) * (fSOM$map$ydim)
# nb.nodes <- nb.clusters
# nb.points.per.line <- nb.cells/nb.nodes
# for (i in 1:nb.clusters) {
#   tempo <- mapping.fcs[mapping.fcs==i]
#   tempo_indices <- which(mapping.fcs==i)
#   
#   random.nb <- runif(n=tempo,min = 0, max = 1 )
#   
#   random.radius <-  random.nb * vertex.size.fcs[i]/2
#   random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
#   
#   node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
#   node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
# }
# # plot(node.fSOM.x, node.fSOM.y, main= "Frozen flowSOM + size")
# 
# matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
# matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=10, xmax=1014, ymin=10, ymax=1014)
# colnames(matrix.node.fSOM.scaled) <- c("Free FlowSOM X","Free FlowSOM Y")
# 
# ff.fsom <- cbind2 (ff, matrix.node.fSOM.scaled)
# 
# # get metaclustering data for fcs file
# metacluster <- fSOM.res [[2]]
# metacluster <- metaClustering_consensus(fSOM$map$codes, k = 5)
# data.metacluster <- metacluster[fSOM$map$mapping[,1]] # = metaClustering_perCell <- metaClustering[fSOM$map$mapping[,1]]
# data.metacluster <- as.matrix(data.metacluster)
# colnames(data.metacluster) <- "Metaclustering Consensus FREE"
# 
# ff.meta <-  cbind2(ff.fsom, data.metacluster)
# 
# for (i in 1:3) {
#   mtrx <- exprs(ff.meta[indices.list[[i]],18:20])
#   if (i==1) {
#     free.mtrx <- mtrx
#   } else {
#     free.mtrx <- rbind(free.mtrx, mtrx)
#     }
# }
# 
# ff.total <- fr_append_cols(ff.frozen, free.mtrx)
# 
# ff.new <- ff[indices.list[[m]],]
# ff.new <- fr_append_cols(ff.new, ff.nodes[[m]])
# mtrx <- exprs(ff.new)
# if (m==1){
#   merged.mtrx <- mtrx
# } else {
#   merged.mtrx <- rbind(merged.mtrx, mtrx)
# }
# 
# output.file <- selectFile(caption = "Save file as ...", label = "Save", existing = FALSE)
# 
# write.FCS(ff.total, filename = paste0(output.file, ".fcs"))
# 
# # Figures
# pdf.name <- paste0(output.file, ".pdf")
# pdf(file = pdf.name,
#     width = 8.3, 
#     height = 11.7)
# par(mfrow=c(1,1))
# par(cex.main=1.4)
# fSOM$prettyColnames <- as.character(fl.names)
# PlotStars(fSOM)
# # Plot metaclusters
# PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metacluster), 
#           backgroundColor = grDevices::colorRampPalette(rainbow(5)),
#           main = "MST with metaclusters")
# # Plot metaclusters
# PlotStars(fSOM, view= "grid", backgroundValues = as.factor(metacluster), 
#           backgroundColor = grDevices::colorRampPalette(rainbow(5)),
#           main = "MST with metaclusters")
# # Query for CD34+ CD10+ population (MRD)
# query <- c("FL3 INT" = "high", "FL6 INT" = "high")
# query_res <- QueryStarPlot(UpdateNodeSize(fSOM, reset=TRUE), query, plot = FALSE)
# cellTypes <- factor(rep("Others", 49), levels=c("Others", "CD34+ CD10+"))
# cellTypes[query_res$selected] <- "CD34+ CD10+"
# PlotStars(fSOM, backgroundValues = cellTypes, backgroundColor = c("#FFFFFF00", "#0000FF22"))
# dev.off()
# 
