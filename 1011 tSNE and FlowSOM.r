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
# Read FCS file
input.folder <- "D:/school/Stage officieel/csv_out/"
fcs.path <- list.files(path= input.folder, pattern = "\\.fcs", full.names=TRUE)
ff <- read.FCS(fcs.path, column.pattern ="TIME", invert.pattern = TRUE) 
# ff = FlowFrame, nrows = events, ncol = parameters, class of ff = FlowFrame, type = s4 --> for plots should be numeric, use exprs() to make matrix
# column.pattern and invert.pattern make sure TIME is not included.

#### Logicle transformation
summary(ff)
channels <- colnames(ff[,1:15])
print(channels)


lgcl <- estimateLogicle(ff, channels)
transformed <- transform(ff, lgcl)
summary(transformed)

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
inputmatrix <- exprs(transformed[,1:15]) # converts flowframe into matrix
summary(inputmatrix)
ffnormalized <- normalize(inputmatrix)
summary(ffnormalized)
logicle <- ffnormalized
linear <- logicle*1024


##################################################################################
# Predefined settings
##################################################################################
# Enter here the column index/numbering to be used for tsne and FlowSOM calculation
tsne.parameter <- c(6:11, 13:14)
tsne.perplexity <- 30 # default =30
tsne.theta <- 0.1 # 0 for exact t-sne, default = 0.5
tsne.max_iter <- 1000 #default = 1000
tsne.scaled.center <- FALSE
tsne.scaled.scale <- FALSE


FlowSOM.parameter <- c(6:11, 13:14)
FlowSOM.autocluster.nb <- 20 #number of clusters to find
FlowSOM.nodes.total.nb <- 120
set_seed <- 1234
FlowSOM.scale <- FALSE
FlowSOM.scaled.center <- FALSE
FlowSOM.scaled.scale <- FALSE


umap.custom.settings <- umap.defaults
umap.custom.settings$n_neighbors <- 15 #15
umap.custom.settings$n_nepochs <- 2000 #200


##################################################################################
# Enter here the ADC resolution from your flow cytometer 
# ie: BD Diva= 18, BCI Navios=20, BCI CytoFLex=24, BCI MoFLo XDP and Astrios=32, BCI CyAn=16, 
#<N> parameter=10, parameters coming from Kaluza logicle Matrix =10 
ADC.resolution <- 10
# Removed scaling to linear --> already done by extraction with plugin


###################################################################################
# Enter here the path of the output folder             
output.folder <- "D:/school/Stage officieel/tsne_out/"
dir.create(path=output.folder)
# Enter here the ouput fcs file name             
output.fcs.file.name <- "Carolien.fcs"


##################################################################################
# Enter here the number of cells to be used for computation : downsampling
nb.cells.to.be.taken <- 100000000  #max 100 000 cells in R 32 bit


###################################################################################
# Script start here                                                                 #
###################################################################################

###################################################################################
# Downsampling
if(nb.cells.to.be.taken < nrow(ff)) {
                                  indices <- sample(nrow(linear), nb.cells.to.be.taken)
                                  table.to.be.computed <- linear[linear_indices,]
		        table.to.be.used.as.fcs.file <- linear[linear_indices,]
}else{
  table.to.be.computed <- linear
  table.to.be.used.as.fcs.file <- linear
  }

#Check Compatibility test between table.to.be.computed and table.to.be.used.as.fcs.file
if((colnames(table.to.be.computed) != colnames(table.to.be.used.as.fcs.file)) & 
(nrow(table.to.be.computed)!= nrow(table.to.be.used.as.fcs.file)) ) {
 stop("The tsne_lgcl parameters do not fit with the tsne_linear parameters")} 

#loading packages
library("Rtsne")
library("GEOmap")
library(MASS)
library(flowCore)

#Extraction of the parameter names
fcs.description.name <- colnames(table.to.be.computed)
fcs.description.desc <-  fcs.description.name

#check logicle coef in 1D plot
windows()
nb.plot.per.row <- 4
nb.plot.per.column <- 4
par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) )
nb.plot.per.windows <- nb.plot.per.row * nb.plot.per.column
aa <- 0

for(i in tsne.parameter) {
if(aa==nb.plot.per.windows){
windows()
par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) )
aa <- 0
}
aa <- aa+1
x1 <- table.to.be.computed[,i]
d1 <- density(x1)
plot(d1, xlab= fcs.description.name[i], main="Check logicle transformation", cex.main=0.7, xlim = c(0,1024))
}

#pairs density plot
panel.col <- function(x,y,...){
 points(x, y, pch=".", cex=0.51, 
         col= densCols(x,y,
                       colramp = colorRampPalette(c("blue2","green2","red2","yellow"))))
}

windows()
pairs(table.to.be.computed [,tsne.parameter], lower.panel=NULL,
      panel=panel.col,
      labels=fcs.description.name[tsne.parameter], gap=0.2,
      xlim=c(0,1024),ylim=c(0,1024),
      xaxt = "n", yaxt = "n",
      main= "Add all plot - cells gated" )

#windows()
#tSNE calculation on selected parameters
#rtsne_out <- Rtsne(table.to.be.computed[,tsne.parameter])
#tsne.1<-rtsne_out$Y[,1]
#tsne.2<-rtsne_out$Y[,2]

#t-SNE Dot Plot
#plot(tsne.1, tsne.2, main="tSNE Plot",pch=20, cex=0.3)

#t-SNE Density plot
#zz2 <- kde2d(tsne.1, tsne.2, n = 200)
#windows()
#image(zz2, col=c('white','blue','cyan','green','yellow','green','red','darkred'), nlevel=11, legend.only=TRUE,
#      main="R Image plot - t-sne - File #1", xlab="t-sne 1",ylab="t-sne 2") #, xlim=x.range, ylim=y.range)
#contour(zz2, drawlabels = FALSE, nlevels = 11, col= "black", add=TRUE) 
#axis(2)
#axis(1)
#box()


#scale the t-sne results on 1024 channel
#matrix.tsne <- cbind(tsne.1, tsne.2)
#matrix.tsne.scaled <- igraph::norm_coords(matrix.tsne, xmin=10, xmax=1014, ymin=10, ymax=1014)
#colnames(matrix.tsne.scaled) <- c("tsne.1","tsne.2")

#tsne.1.scaled <- matrix.tsne.scaled[,1]
#tsne.2.scaled <- matrix.tsne.scaled[,2]

#scale the flowFrame on 1024 channel
ff.scaled.csv <- table.to.be.used.as.fcs.file 
ff.scaled <- flowFrame(ff.scaled.csv)


#Add the t-sne matrix to the fcs file 
#ff.new <- cbind2(ff.scaled,
#                 matrix.tsne.scaled)

ff.new <- ff.scaled

#FlowSOM calculation
library(FlowSOM)
csv.FlowSOM <- as.matrix(table.to.be.computed)
#csv.FLowSOM.scaled <-  igraph::norm_coords(csv.FlowSOM, xmin=10, xmax=1014, ymin=10, ymax=1014)
csv.FLowSOM.scaled <- csv.FlowSOM


ff.FlowSOM.scaled <- flowFrame(csv.FLowSOM.scaled)

flowSOM.res <- FlowSOM(ff.FlowSOM.scaled, 
			scale=FlowSOM.scale,  
			colsToUse=FlowSOM.parameter,  
			#maxMeta=5, 
			xdim = round(sqrt(FlowSOM.nodes.total.nb)), 
			ydim = round(sqrt(FlowSOM.nodes.total.nb)),
			scaled.center = FlowSOM.scaled.center,
			scaled.scale = FlowSOM.scaled.scale,
			nClus = FlowSOM.autocluster.nb,
			rlen=10,
			compensate = FALSE,
			transform = FALSE)

fSOM <- flowSOM.res [[1]]

#FlowSOM.parameter <- tsne.parameter
colsToCluster <- FlowSOM.parameter
fcs.description.name <- colnames(table.to.be.used.as.fcs.file)
fcs.description.desc <-  colnames(table.to.be.used.as.fcs.file)



###
windows()
par(mfrow=c(1,1))
par(cex.main=1.4)
PlotStars(fSOM)
###

###
windows()
target.channel <- FlowSOM.parameter 
marker.name.desc <- fcs.description.desc[target.channel]
marker.name.name <- fcs.description.name[target.channel]
nb.plot.per.row <- 2
nb.plot.per.column <- 2
par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5), cex.main=0.8 )
nb.plot.per.windows <- nb.plot.per.row * nb.plot.per.column
aa <- 0

for(i in 1: length(colsToCluster)) {
  if(aa==nb.plot.per.windows){
    windows()
    par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) , cex.main=0.8)
    aa <- 0
  }
  aa <- aa+1
  
  marker.name <- colsToCluster [i]
  PlotMarker(fSOM,marker.name, 
             main = paste("FlowSOM plot Marker", " <",marker.name.name[i],"> ", 
                          #marker.name.desc[i],
                          sep=""))
}


#plot all nodes in a dot plot: check if there are nodes in CD19+CD45RO+
windows()
PlotCenters(fSOM, 7, 8)

###

    #ff_o <- ff.new 
    #ff <- ff.FlowSOM.scaled
      #s <- seq_len(nrow(ff_o))
      #fsom_f <- NewData(fsom, ff)
    #m <- matrix(0, nrow = nrow(ff_o), ncol = 1)
    #m[s, ] <- fsom_f$map$mapping[, 1]
    #colnames(m) <- "FlowSOM"
    #ff_o <- flowCore::cbind2(ff_o, m)

    #flowCore::write.FCS(ff_o, filename = gsub("\\.fcs", "_FlowSOM.fcs", original_files[i]))
 
layout.fcs <- fSOM$MST$l
mapping.fcs <- fSOM$map$mapping[,1]

  node.fSOM.x  <-  layout.fcs[,1]
  node.fSOM.y  <- layout.fcs[,2]

vertex.size.fcs <- fSOM$MST$size

nb.cluster <- (fSOM$map$xdim) * (fSOM$map$ydim)

nb.nodes <- nb.cluster
nb.cells <- nrow(csv.FLowSOM.scaled)
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

ff.new.1 <- cbind2 (ff.new, matrix.node.fSOM.scaled)

#Emdedd the FlowSOM Metaclusterization result in the fcs file
metacluster <- flowSOM.res [[2]]
metacluster <- metaClustering_consensus(fSOM$map$codes, k = FlowSOM.autocluster.nb)
#metacluster <- MetaClustering(fSOM$map$codes, method = "metaClustering_Consensus", max = 10) #== flowSOM.res[[2]]

###
windows()
par(mfrow=c(1,1))
par(cex.main=1.4)
PlotStars(fSOM, view= "MST", backgroundValues = as.factor(metacluster), main = "PlotStars + AutoClusters") 
###

data.metacluster <- metacluster[fSOM$map$mapping[,1]]
data.metacluster <- as.matrix(data.metacluster)
colnames(data.metacluster) <- "Metaclustering Consensus"

#windows()
#t-SNE Dot Plot
#plot(tsne.1, tsne.2, main="tSNE Plot",pch=20, cex=0.3, col= data.metacluster)

#colorPalette <-  colorRampPalette(c("blue2", "green2","red2", "yellow") )
#color.ben <- colorPalette(max(data.metacluster)) [data.metacluster]

#meta.cluster <- as.matrix(data.metacluster)
#colnames(meta.cluster) <- "FlowSOM Metaclusterisation"

ff.new.2 <-  cbind2(ff.new.1, data.metacluster)

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*
#ff t-sne calculation
# Note: this script should be sourced as: source('<path to file>', chdir=T)

FAST_TSNE_SCRIPT_DIR <- "C:"
#FAST_TSNE_SCRIPT_DIR <<- getwd()

#t-sne matrix

X <- as.matrix(table.to.be.computed[,tsne.parameter])
X <- scale(X, center= tsne.scaled.center, scale = tsne.scaled.scale)


message("FIt-SNE R wrapper loading.")
message("FIt-SNE root directory was set to ",  FAST_TSNE_SCRIPT_DIR)

# Compute FIt-SNE of a dataset
#       dims - dimensionality of the embedding. Default 2.
#       perplexity - perplexity is used to determine the
#           bandwidth of the Gaussian kernel in the input
#           space.  Default 30.
#       theta - Set to 0 for exact.  If non-zero, then will use either
#           Barnes Hut or FIt-SNE based on nbody_algo.  If Barnes Hut, then
#           this determins the accuracy of BH approximation.
#           Default 0.5.
#       max_iter - Number of iterations of t-SNE to run.
#           Default 1000.
#       fft_not_bh - if theta is nonzero, this determins whether to
#            use FIt-SNE or Barnes Hut approximation. Default is FIt-SNE.
#            set to be True for FIt-SNE
#       ann_not_vptree - use vp-trees (as in bhtsne) or approximate nearest neighbors (default).
#            set to be True for approximate nearest neighbors
#       exaggeration_factor - coefficient for early exaggeration
#           (>1). Default 12.
#       no_momentum_during_exag - Set to 0 to use momentum
#           and other optimization tricks. 1 to do plain,vanilla
#           gradient descent (useful for testing large exaggeration
#           coefficients)
#       stop_early_exag_iter - When to switch off early exaggeration.
#           Default 250.
#       start_late_exag_iter - When to start late
#           exaggeration. set to -1 to not use late exaggeration
#           Default -1.
#       late_exag_coeff - Late exaggeration coefficient.
#          Set to -1 to not use late exaggeration.
#           Default -1
#       nterms - If using FIt-SNE, this is the number of
#                      interpolation points per sub-interval
#       intervals_per_integer - See min_num_intervals              
#       min_num_intervals - Let maxloc = ceil(max(max(X)))
#           and minloc = floor(min(min(X))). i.e. the points are in
#           a [minloc]^no_dims by [maxloc]^no_dims interval/square.
#           The number of intervals in each dimension is either
#           min_num_intervals or ceil((maxloc -
#           minloc)/intervals_per_integer), whichever is
#           larger. min_num_intervals must be an integer >0,
#           and intervals_per_integer must be >0. Default:
#           min_num_intervals=50, intervals_per_integer =
#           1
#
#       sigma - Fixed sigma value to use when perplexity==-1
#            Default -1 (None)
#       K - Number of nearest neighbours to get when using fixed sigma
#            Default -30 (None)
#
#       initialization - N x no_dims array to intialize the solution
#            Default: None
#
#       load_affinities - 
#            If 1, input similarities are loaded from a file and not computed
#            If 2, input similarities are saved into a file.
#            If 0, affinities are neither saved nor loaded
#
#       perplexity_list - if perplexity==0 then perplexity combination will
#            be used with values taken from perplexity_list. Default: NULL
#       df - Degree of freedom of t-distribution, must be greater than 0.
#       Values smaller than 1 correspond to heavier tails, which can often 
#       resolve substructure in the embedding. See Kobak et al. (2019) for
#       details. Default is 1.0
#
fftRtsne <- function(X, 
                     dims = 2, perplexity = tsne.perplexity, theta = tsne.theta,
                     max_iter = tsne.max_iter,
                     fft_not_bh = TRUE,
                     ann_not_vptree = TRUE,
                     stop_early_exag_iter = 250,
                     exaggeration_factor = 12.0, no_momentum_during_exag = FALSE,
                     start_late_exag_iter = -1.0, late_exag_coeff = 1.0,
                     mom_switch_iter = 250, momentum = 0.5, final_momentum = 0.8, learning_rate = 200,
                     n_trees = 50, search_k = -1, rand_seed = -1,
                     nterms = 3, intervals_per_integer = 1, min_num_intervals = 50, 
                     K = -1, sigma = -30, initialization = NULL,
                     data_path = NULL, result_path = NULL,
                     load_affinities = NULL,
                     fast_tsne_path = NULL, nthreads = 0, perplexity_list = NULL, 
                     get_costs = FALSE, df = 1.0) {
  
  version_number <- '1.1.0'
  
  if (is.null(fast_tsne_path)) {
    if (.Platform$OS.type == "unix") {
      fast_tsne_path <- file.path(FAST_TSNE_SCRIPT_DIR, "bin", "fast_tsne")
    } else {
      fast_tsne_path <- file.path(FAST_TSNE_SCRIPT_DIR, "bin", "FItSNE.exe")
    }
  }
  
  if (is.null(data_path)) {
    data_path <- tempfile(pattern = 'fftRtsne_data_', fileext = '.dat')
  }
  if (is.null(result_path)) {
    result_path <- tempfile(pattern = 'fftRtsne_result_', fileext = '.dat')
  }
  if (is.null(fast_tsne_path)) {
    fast_tsne_path <- system2('which', 'fast_tsne', stdout = TRUE)
  }
  fast_tsne_path <- normalizePath(fast_tsne_path)
  if (!file_test('-x', fast_tsne_path)) {
    stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
  }
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  if (!is.numeric(theta) || (theta < 0.0) || (theta > 1.0) ) { stop("Incorrect theta.")}
  if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
  if (!is.matrix(X)) { stop("Input X is not a matrix")}
  if (!(max_iter > 0)) { stop("Incorrect number of iterations.")}
  if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter < 0) { stop("stop_early_exag_iter should be a positive integer")}
  if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
  if (!is.numeric(df)) { stop("df should be numeric")}
  if (!is.wholenumber(dims) || dims <= 0) { stop("Incorrect dimensionality.")}
  if (search_k == -1) {
    if (perplexity > 0) {
      search_k <- n_trees * perplexity * 3
    } else if (perplexity == 0) {
      search_k <- n_trees * max(perplexity_list) * 3
    } else { 
      search_k <- n_trees * K
    }
  }
  
  if (fft_not_bh) {
    nbody_algo <- 2
  } else {
    nbody_algo <- 1
  }
  
  if (is.null(load_affinities)) {
    load_affinities <- 0
  } else {
    if (load_affinities == 'load') {
      load_affinities <- 1
    } else if (load_affinities == 'save') {
      load_affinities <- 2
    } else {
      load_affinities <- 0
    }
  }
  
  if (ann_not_vptree) {
    knn_algo <- 1
  } else {
    knn_algo <- 2
  }
  tX <- as.numeric(t(X))
  
  f <- file(data_path, "wb")
  n <- nrow(X)
  D <- ncol(X)
  writeBin(as.integer(n), f, size = 4)
  writeBin(as.integer(D), f, size = 4)
  writeBin(as.numeric(theta), f, size = 8) #theta
  writeBin(as.numeric(perplexity), f, size = 8)
  
  if (perplexity == 0) {
    writeBin(as.integer(length(perplexity_list)), f, size = 4)
    writeBin(perplexity_list, f) 
  }
  
  writeBin(as.integer(dims), f, size = 4)
  writeBin(as.integer(max_iter), f, size = 4)
  writeBin(as.integer(stop_early_exag_iter), f, size = 4)
  writeBin(as.integer(mom_switch_iter), f, size = 4)
  writeBin(as.numeric(momentum), f, size = 8)
  writeBin(as.numeric(final_momentum), f, size = 8)
  writeBin(as.numeric(learning_rate), f, size = 8)
  writeBin(as.integer(K), f, size = 4) #K
  writeBin(as.numeric(sigma), f, size = 8) #sigma
  writeBin(as.integer(nbody_algo), f, size = 4)  #not barnes hut
  writeBin(as.integer(knn_algo), f, size = 4) 
  writeBin(as.numeric(exaggeration_factor), f, size = 8) #compexag
  writeBin(as.integer(no_momentum_during_exag), f, size = 4) 
  writeBin(as.integer(n_trees), f, size = 4) 
  writeBin(as.integer(search_k), f, size = 4) 
  writeBin(as.integer(start_late_exag_iter), f, size = 4) 
  writeBin(as.numeric(late_exag_coeff), f, size = 8) 
  
  writeBin(as.integer(nterms), f, size = 4) 
  writeBin(as.numeric(intervals_per_integer), f, size = 8) 
  writeBin(as.integer(min_num_intervals), f, size = 4) 
  writeBin(tX, f) 
  writeBin(as.integer(rand_seed), f, size = 4) 
  writeBin(as.numeric(df), f, size = 8)
  writeBin(as.integer(load_affinities), f, size = 4) 
  if (!is.null(initialization)) { writeBin( c(t(initialization)), f) }		
  close(f) 
  
  flag <- system2(command = fast_tsne_path, 
                  args = c(version_number, data_path, result_path, nthreads))
  if (flag != 0) {
    stop('tsne call failed')
  }
  f <- file(result_path, "rb")
  n <- readBin(f, integer(), n = 1, size = 4)
  d <- readBin(f, integer(), n = 1, size = 4)
  Y <- readBin(f, numeric(), n = n * d)
  Y <- t(matrix(Y, nrow = d))
  if (get_costs) {
    readBin(f, integer(), n = 1, size = 4)
    costs <- readBin(f, numeric(), n = max_iter, size = 8)
    Yout <- list(Y = Y, costs = costs)
  } else {
    Yout <- Y
  }
  close(f)
  file.remove(data_path)
  file.remove(result_path)
  Yout
}


fftRtsne_Res <- fftRtsne(X)
#plot(fftRtsne_Res)

fftsne.1 <- fftRtsne_Res[,1]
fftsne.2 <- fftRtsne_Res[,2]

windows()
#t-SNE Dot Plot
plot(fftsne.1, fftsne.2, main="ff tSNE Plot",pch=20, cex=0.3)

#t-SNE Density plot
zz2 <- kde2d(fftsne.1,fftsne.2, n = 200)
windows()
image(zz2, col=c('white','blue','cyan','green','yellow','green','red','darkred'), nlevel=11, legend.only=TRUE,
      main="R Image plot - fft-sne - File #1", xlab="t-sne 1",ylab="t-sne 2") #, xlim=x.range, ylim=y.range)
contour(zz2, drawlabels = FALSE, nlevels = 11, col= "black", add=TRUE) 
axis(2)
axis(1)
box()


#scale the t-sne results on 1024 channel
matrix.fftsne <- cbind(fftsne.1, fftsne.2)
matrix.fftsne.scaled <- igraph::norm_coords(matrix.fftsne, xmin=10, xmax=1014, ymin=10, ymax=1014)
colnames(matrix.fftsne.scaled) <- c("ff.tsne.1","ff.tsne.2")

fftsne.1.scaled <- matrix.fftsne.scaled[,1]
fftsne.2.scaled <- matrix.fftsne.scaled[,2]

windows()
#t-SNE Dot Plot
plot(fftsne.1.scaled, fftsne.2.scaled, main="ff tSNE Plot scaled",pch=20, cex=0.3)


ff.new.3 <-  cbind2(ff.new.2, matrix.fftsne.scaled)


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*
#umap calculation
library(umap)
umap_out <- umap(X, umap.custom.settings)

umap.1 <- umap_out$layout[,1]
umap.2 <- umap_out$layout[,2]
#plot the tSNE

windows()
plot(umap.1,umap.2, main="umap Plot",pch=20, cex=0.3)
#scale the t-sne results on 1024 channel
xmin.graph.umap <-  par("usr")[1]
xmax.graph.umap <- par("usr")[2]

ymin.graph.umap <- par("usr")[3]
ymax.graph.umap <- par("usr")[4]

umap.1.scaled <- scales::rescale(umap.1, to = c(1, 2^10), from = c( xmin.graph.umap, xmax.graph.umap))
umap.2.scaled <- scales::rescale(umap.2, to = c(1, 2^10), from = c( ymin.graph.umap, ymax.graph.umap))

windows()
plot(umap.1.scaled,umap.2.scaled, main="umap Plot - rescale",pch=20, cex=0.3)



zz2 <- kde2d(umap.1.scaled, umap.2.scaled, n = 200)
windows()
image(zz2, col=c('white','blue','cyan','green','yellow','green','red','darkred'), nlevel=11, legend.only=TRUE,
      main="R Image plot - umap - File #1", xlab="t-sne 1",ylab="t-sne 2") #, xlim=x.range, ylim=y.range)
contour(zz2, drawlabels = FALSE, nlevels = 11, col= "black", add=TRUE) 
axis(2)
axis(1)
box()

matrix.umap.scaled <- cbind(umap.1.scaled,umap.2.scaled)

ff.new.4 <-  cbind2(ff.new.3, matrix.umap.scaled)


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*

#Define the complete output fcs file name
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S-")

sFilename <- paste(output.folder, date.format, output.fcs.file.name, sep="")
Filename.tempo <- sFilename
print(sFilename)

write.FCS(ff.new.4, filename = sFilename)


#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*


 #windows()
#color.ben <- c("blue2", "cyan2",  "darkslateblue",  "green4", "lightblue3"   ,"red2", "lightgoldenrod3", "indianred4",  "ivory4")
#plot(fftsne.1, fftsne.2, main="ff tSNE Plot",pch=20, cex=1, col= color.ben[data.metacluster])

windows()
plot(fftsne.1, fftsne.2, main="ff tSNE Plot - FlowSOM clustering",pch=20, cex=1, 
     col= rainbow (length(unique(metacluster)), alpha=0.3) [data.metacluster])

mylims <- par("usr")

legend.position <- (mylims[2] - mylims[1])/40
library(plotrix)
color.legend(xl = mylims[2], xr=mylims[2] +legend.position ,
             yt = mylims[4], yb = mylims[3],
             legend = 1:length(unique(metacluster)),
             rect.col = rainbow (length(unique(metacluster)), alpha=0.3),
             align = "lt", gradient = "y", cex=0.6)





#t-sne plot colorisation per marker intensity
library(grDevices)

colorPalette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", 
                                   "#7FFF7F", "yellow", "#FF7F00","red", "#7F0000"))



windows()

target.channel <-1:(ncol(table.to.be.computed) )
marker.name.desc <- fcs.description.desc[target.channel]
marker.name.name <- fcs.description.name[target.channel]
nb.plot.per.row <- 2
nb.plot.per.column <- 2
par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5), cex.main=0.8 )
nb.plot.per.windows <- nb.plot.per.row * nb.plot.per.column
aa <- 0

for(i in 1: length(target.channel)) {
  if(aa==nb.plot.per.windows){
    windows()
    par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) , cex.main=0.8)
    aa <- 0
  }
  aa <- aa+1
  marker.name <- marker.name.desc [i]
  color.Palette.tsne <- colorPalette(100)[as.numeric(cut(table.to.be.computed[,i],breaks =100))]
  plot(fftsne.1,fftsne.2, pch=20, cex = 0.5, col= color.Palette.tsne, 
     main = paste("ff t-sne: marker expression for ", colnames(table.to.be.computed)[i]),
     cex.main=0.8, axes=TRUE, xlab="ff t-sne 1", ylab="ff t-sne 2")
  
mylims <- par("usr")
legend.position <- (mylims[2] - mylims[1])/60
library(plotrix)
color.legend(xl = mylims[1], xr= mylims[1]+ (mylims[2]-mylims[1])/2 ,
             yt = mylims[3]+legend.position, yb = mylims[3],
             legend = c(round(min(table.to.be.computed[,i]),1), 
                        round(max(table.to.be.computed[,i]), 1)) ,
             rect.col = c("blue", "#007FFF", "cyan", 
                                   "#7FFF7F", "yellow", "#FF7F00","red", "#7F0000"),
             align = "lt", gradient = "x", cex=0.5) 

}

#-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*
#Plot the heatmap
cluster.DF <- cbind(logicle, data.metacluster)
nb.clusters <- length(unique(data.metacluster)) 
cluster.mean <- list()
ncol.cluster.DF <- ncol(cluster.DF)


for(i in 1:nb.clusters){
  cluster.DF.indices <- which(cluster.DF[,ncol.cluster.DF]==i)
  cluster.mean [[i]] <- apply(cluster.DF[cluster.DF.indices, ], 2, median)
}
cluster.mean.DF <- do.call(cbind, cluster.mean)
colnames(cluster.mean.DF) <- paste("cluster ", 1:nb.clusters, sep="")

lab.Row <- colnames(logicle)

lab.Row <- strsplit(lab.Row, " ")
fcs.desc <- c()
for(j in 1 : length(lab.Row) ) {
  fcs.desc [j] <- lab.Row [[j]] [1]
}

lab.Row <- fcs.desc

library(gplots)

windows()

heatmap.2(as.matrix(cluster.mean.DF[1:nrow(cluster.mean.DF)-1, ] ), 
		col= bluered(15), keysize=1, trace= "none", 
		symbreaks = FALSE, symkey = FALSE, main = "gplots::heatmaps clustering results",
		cex.main= 0.6, margins = c(10,6),
		key.par = list(mgp= c(1.5, 0.5, 0), mar = c(4.5, 2.5, 3, 1.5)),
		cexRow=1, cexCol = 1,
		labRow = lab.Row, dendogram = "column", Rowv = FALSE, Colv= TRUE) 

# > meta_mfi <- FlowSOM::MetaclusterMFIs(flowSOM.res)
# > str(meta_mfi)
# num [1:20, 1:15] 963 965 950 953 963 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:20] "1" "2" "3" "4" ...
# ..$ : chr [1:15] "FS.PEAK" "FS.INT" "FS.TOF" "SS.INT" ...
# > heatmap.2(meta_mfi)


# AggregateFlowFrames()
# of
# flowSet(ff, ff2)
# FlowSOMSubset()
# fSOM$metaData

# table(GetClusters(fSOM)) / 
# PlotVariable() 
