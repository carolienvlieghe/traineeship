##################### Script 1011 tSNE & FlowSOM for Rstudio #####################
##################################################################################

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

##################################################################################
# Loading packages
library("Rtsne")
library("GEOmap")
library("pca3d")
library(rgl)
library(MASS)
library(flowCore)
library(FlowSOM)
library(umap)

##################################################################################
### Assign in & output ###
input.folder <- "D:/school/Stage officieel/csv_out/"
output.folder <- "D:/school/Stage officieel/tsne_out/"
dir.create(path=output.folder)
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S")
file.name.fcs <- list.files(path = input.folder, pattern = "\\.fcs$", full.names = TRUE)
for (fcs in file.name.fcs) {
  file.name.base <- basename(fcs)
  base <- substr(file.name.base, start = 19, stop = nchar(file.name.base)-4)
  output.file.path <- paste0(output.folder, date.format, " - flowSOM - ", base, ".fcs")
  pdf.name <- paste0(output.folder,date.format," - plots - ", base, ".pdf")
  
  ##################################################################################
  ### Read fcs file to analyze ###
  ff <- read.FCS(fcs) # opens flowFrame
  
  ##################################################################################
  ### Transformation & normalization ###
  # Assign channels that need to be transformed
  # channels <- colnames(ff[,6:15])
  # Logicle transformation on FL channels, look at range, m can be default unless range exceeds 4.5 decades try loop?
  # lgcl <- estimateLogicle(ff, channels)
  # transformed <- transform(ff, lgcl)
  
  # Better transformation with??:
  transformed <- transform(ff,transformList(colnames(ff)[6:15],logicleTransform()))
  
  # normalize function
  # normalize with median of best negative population: the fcs file is a subpopulation were all cells were CD5 negative!
  # norm_median <- function(x) {
  #   return ( x - median(exprs(transformed[,13])))
  # }
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  # normalize with mean = 0, SD = 1
  # normsd <- function(x) {
  #   return ((x - mean(x))/sd(x))
  # }
  # apply normalization function
  inputmatrix <- exprs(transformed[,1:15]) # converts flowframe into matrix
  # FL's
  ffnormalized <- normalize(inputmatrix[,6:15])
  # scatter
  scatternorm <- normalize(inputmatrix[,1:5])
  # Combine scatter and FL
  ffmatrix <- cbind(scatternorm, ffnormalized)
  # scale to [0,1024]
  linear <- ffmatrix*1024
  # make flowframe of scaled [0,1024] data
  param <- parameters(ff)
  ff.scaled <- flowFrame(linear, parameters = param)
  
  # assign names and description names according to kaluza: script csv2fcs adjusted
  fcs.description.name <- colnames(linear)
  fcs.description.desc <- as.character(parameters(ff)$desc[1:15])
  
  
  ##################################################################################
  ################################## Start #########################################
  ##################################################################################
  
  # open pdf file to save plots
  # pdf.output.path <- paste0(output.folder, date.format, " - flowSOM - ", base, ".pdf")
  # pdf(pdf.output.path) # opens new pdf file
  
  # Assign parameters to use for tSNE, flowSOM and umap
  # Which parameters to use: don't use parameters used in preprocessing
  parameters.to.use <- c(8:11, 14)
  
  ### Check transformation with individual density plots
  
  # Save density plots in pdf
  pdf(file = pdf.name,
      width = 8.3, 
      height = 11.7)
  
  nb.plot.per.row <- 4
  nb.plot.per.column <- 3
  par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s")
  nb.plot.per.window <- nb.plot.per.row * nb.plot.per.column
  aa <- 0
  for (i in c(6:15)) {
    if(aa==nb.plot.per.window){
      windows()
      par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s")
      aa <- 0
    }
    aa <- aa +1
    x1 <- linear[,i]
    d1 <- density(x1)
    plot(d1, xlab = fcs.description.desc[i], main = "Check logicle transformation", cew.main=0.7, xlim=c(0,1024))
  }
  
  ### Paired density plots on 1 page for FL's used in tsne
  # This slows the script and could be potentially cause trouble when saving in PDF, so if we want this: save as jpeg
  # panel.col <- function(x,y,...){
  #   points(x, y, pch=".", cex=0.51, col= densCols(x,y, colramp = colorRampPalette(c("blue2","green2","red2","yellow"))))
  # }
  # windows()
  # pairs(linear [,parameters.to.use], lower.panel=NULL,
  #       panel=panel.col,
  #       labels=fcs.description.desc[parameters.to.use], gap=0.2,
  #       xlim=c(0,1024),ylim=c(0,1024),
  #       xaxt = "n", yaxt = "n",
  #       main= "Add all plot - cells gated" )
  
  ###############
  ### FlowSOM ###
  ###############
  ### The columns used for the flowSOM have to be displayed as logicle 
  flowsom.res <- FlowSOM(ff.scaled,
                         scale=FALSE,
                         colsToUse = parameters.to.use,
                         xdim = round(sqrt(120)),
                         ydim = round(sqrt(120)),
                         scaled.center = FALSE,
                         scaled.scale = FALSE,
                         nClus = 10,
                         rlen = 10,
                         compensate = FALSE,
                         transform = FALSE)
  
  fSOM <- flowsom.res[[1]]
  # change names for legend
  fSOM$prettyColnames <- fcs.description.desc
  
  fSOM <- UpdateNodeSize(fSOM, maxNodeSize = 7, reset=TRUE)
  
  plot.new() # start new page for plot
  par(mfrow=c(1,1))
  par(cex.main=1.4)
  PlotStars(fSOM)
  summary(fSOM)
  
  # dev.off()
  # append fsom in fcs
  layout.fcs <- fSOM$MST$l
  mapping.fcs <- fSOM$map$mapping[,1]
  
  node.fSOM.x  <-  layout.fcs[,1]
  node.fSOM.y  <- layout.fcs[,2]
  
  vertex.size.fcs <- fSOM$MST$size
  
  nb.cluster <- (fSOM$map$xdim) * (fSOM$map$ydim)
  
  nb.nodes <- nb.cluster
  nb.cells <- nrow(ff.scaled)
  nb.points.per.line <- nb.cells /nb.nodes
  dispersion.coef <- 60
  
  for (i in 1:nb.cluster) {
    tempo <- mapping.fcs[mapping.fcs==i]
    tempo_indices <- which(mapping.fcs==i)
    
    random.nb <- runif(n=tempo,min = 0, max = 1 )
    
    random.radius <-  random.nb * vertex.size.fcs[i]/2
    random.angle <- runif(n=length(tempo),min = 0, max = 2*pi)
    
    node.fSOM.x [tempo_indices] <- (random.radius*cos(random.angle))/dispersion.coef + layout.fcs[i,1]
    node.fSOM.y [tempo_indices] <- (random.radius*sin(random.angle))/dispersion.coef + layout.fcs[i,2]
    
  }
  
  matrix.node.fSOM <- cbind(node.fSOM.x, node.fSOM.y)
  matrix.node.fSOM.scaled <- igraph::norm_coords(matrix.node.fSOM, xmin=10, xmax=1014, ymin=10, ymax=1014)
  colnames(matrix.node.fSOM.scaled) <- c("FlowSOM.1","FlowSOM.2")
  
  ff.fsom <- cbind2 (ff, matrix.node.fSOM.scaled)
  
  #Emdedd the FlowSOM Metaclusterization result in the fcs file
  metacluster <- flowsom.res [[2]]
  metacluster <- metaClustering_consensus(fSOM$map$codes, k = 6)
  #metacluster <- MetaClustering(fSOM$map$codes, method = "metaClustering_Consensus", max = 10) #== flowSOM.res[[2]]
  
  ###
  par(mfrow=c(1,1))
  par(cex.main=1.4)
  PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metacluster), main = "PlotStars + AutoClusters") 
  dev.off()
  ###
  
  data.metacluster <- metacluster[fSOM$map$mapping[,1]]
  data.metacluster <- as.matrix(data.metacluster)
  colnames(data.metacluster) <- "Metaclustering Consensus"
  
  # windows()
  # t-SNE Dot Plot
  # plot(umap.1, umap.2, main="tSNE Plot",pch=20, cex=0.3, col= data.metacluster)
  # 
  # colorPalette <-  colorRampPalette(c("blue2", "green2","red2", "yellow") )
  # color.ben <- colorPalette(max(data.metacluster)) [data.metacluster]
  
  meta.cluster <- as.matrix(data.metacluster)
  colnames(meta.cluster) <- "FlowSOM Metaclusterisation"
  
  ff.meta <-  cbind2(ff.fsom, data.metacluster)
  
  
  ############
  ### umap ###
  ############
  # plot umap maybe with same colors as flowSOM metaclusters?
  # Pre-defined settings
  umap.custom.settings <- umap.defaults
  umap.custom.settings$n_neighbors <- 15 
  umap.custom.settings$n_nepochs <- 2000 
  
  # make matrix of parameters to use
  umap.mtrx <- as.matrix(linear[,parameters.to.use])
  # calculate umap
  umap_out <- umap(umap.mtrx, umap.custom.settings)
  umap.1 <- umap_out$layout[,1]
  umap.2 <- umap_out$layout[,2]
  
  #scale the t-sne results on 1024 channel
  umap.1.scaled <- scales::rescale(umap.1, to = c(0, 1024))
  umap.2.scaled <- scales::rescale(umap.2, to = c(0, 1024))
  
  matrix.umap.scaled <- cbind(umap.1.scaled,umap.2.scaled)
  ff.umap <- cbind2(ff.meta, matrix.umap.scaled)
  
  
  ############
  ### tSNE ###
  ############
  # takes a long time
  #tSNE calculation on selected parameters
  rtsne_out <- Rtsne(linear[,parameters.to.use])
  tsne.1<-rtsne_out$Y[,1]
  tsne.2<-rtsne_out$Y[,2]
  
  # scale the t-sne results on 1024 channel
  matrix.tsne <- cbind(tsne.1, tsne.2)
  matrix.tsne.scaled <- igraph::norm_coords(matrix.tsne, xmin=10, xmax=1014, ymin=10, ymax=1014)
  colnames(matrix.tsne.scaled) <- c("tsne.1","tsne.2")
  
  tsne.1.scaled <- matrix.tsne.scaled[,1]
  tsne.2.scaled <- matrix.tsne.scaled[,2]
  ff.total <- cbind2(ff.umap, matrix.tsne.scaled)
  
  
  write.FCS(ff.total, filename = output.file.path)
}

