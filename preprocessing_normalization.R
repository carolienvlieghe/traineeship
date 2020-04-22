############################
### Normalization script ###
############################

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(ggcyto)

# Assign in- & output
output.folder <- "D:/school/Stage officieel/norm_out/"
input.folder <- "D:/school/Stage officieel/csv_out/"

# read set of fcs files to process
set <- read.flowSet(files = NULL, path = input.folder)

# Assign number of flowframes (nb.ff) in the flowset, channels to use & markernames
nb.ff <- length(set)
channels <- 6:15 # all fluorochromes
markers <- c("CD58", "CD81", "CD34", "CD22", "CD38", "CD10", "CD19", "CD5", "CD20", "CD45") # B-ALL tube
channel.names <- grep('FL', colnames(set), value = TRUE)

# Transform data to Logicle (biexponential transformation): check density for good transformation (change parameters of function if not good)
# Only the fluorochrome channels are transformed
trf.set <- transform(set,transformList(colnames(set)[channels],logicleTransform(w = 1, t = 1048576)))

# Density plots of all samples to check the logicle transformation, 1 sample per window
for(i in c(1:nb.ff)){
  windows(title = paste("Logicle transformation of sample", i, sep = " "))
  par(mfrow= c(4,3), pty="s")
  nb.plots.per.windows <- 10
  aa <- 0
  for (j in channels) {
    if(aa==nb.plots.per.windows) {
      windows(title = paste("Logicle transformation of sample", i, sep = " "))
      par(mfrow= c(4,3), pty="s")
      aa <- 0
    }
    aa <- aa+1
    d1 <- density(exprs(trf.set[[i]][,j]))
    plot(d1, xlab = colnames(trf.set[[i]][,j]), xlim=c(-1,5), main= markers[aa])
  }
}

# E.g. Scatterplot of two markers:
# plot(x = exprs(trf.set[[1]][,13]), y = exprs(trf.set[[1]][,15]))

# assign values for the landmark for the peak of the negative population (for markers that are only positive, chose a landmark for the peak of the neg. pop)
reference.neg.pop.mode <- c(1, #FL1 CD58
                            3.5, # FL2 CD81 only pos
                            1, # FL3 CD34
                            1.5, # FL4 CD22 dim to pos
                            1.5, # FL5 CD38 dim to pos
                            1, # FL6 CD10
                            2.5, # FL7 CD19 only pos
                            1, # FL8 CD5
                            1, # FL9 CD20
                            1.5) # FL10 CD45

# select peaks
# ff = flowframe, FL = fluorochrome
a <- 1
for (ff in c(1:nb.ff)){
  aa <- 1
  negative.population.mode <- c()
  # Select the peaks with the cursor
  for(FL in channels) {
    windows(title = paste("selection of negative peaks of sample", j, sep = " "))
    x1 <- exprs(trf.set[[ff]][,FL]) # all fluorescence intensities of parameter i
    d1 <- density(x1) # compute density of these values
    plot(d1, xlab= colnames(exprs(trf.set[[ff]][,FL])), main= markers[aa], sub = "Click on negative peak (for CD19 & CD81: chose most positive peak)", cex.main=0.7, xlim=c(-1,5)) # plot these densities
    coord.1 <- locator(1, type = "l", col="red") # locator reads the position of the graphics cursor, max nb of points to locate, type l = line
    # Possibility to turn back when you clicked wrong?
    negative.population.mode[FL] <- coord.1$x # x-coördinate of the line is added for every parameter
    abline(v = negative.population.mode[FL],col="red",lwd=1,lty=5) #adds lines through the cuurrent plot --> v = xvalue vertical line
    aa <- aa+1
    dev.off()
  }
  negative.population.mode <- negative.population.mode[!is.na(negative.population.mode)]
  aa <- 1
  # Compute the new values
  for (FL in channels){
    new.value <- each_row(trf.set[[ff]][,FL], function(x) {x + (reference.neg.pop.mode[a] - negative.population.mode[a])})
    new.value <- as.matrix(new.value)
    if ( FL == 6){
      new <- new.value
    } else {
      new <- cbind(new, new.value)
    }
    a <- a+1
  }
  a <- 1
  # Assign the new values to the flowframe
  colnames(new) <- channel.names
  new <- cbind(exprs(trf.set[[ff]][,1:5]), new)
  new <- cbind(new, exprs(trf.set[[ff]][,16]))
  exprs(trf.set[[ff]]) <- new
}
# Samples are now normalized 


# Inverse transform (necessary? better for output)
logicle <- logicleTransform(w = 1, t = 1048576)
inv <- inverseLogicleTransform(trans = logicle)
inv.set <- transform(trf.set,transformList(colnames(trf.set)[channels],inv))


# Save flowset as separate fcs files in assigned output folder
write.flowSet(inv.set, outdir = output.folder)
