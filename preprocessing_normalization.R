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
fcs.files <- choose.files(multi = TRUE, caption = "Select fcs files to normalize")
output.folder <- choose.dir(caption = "Select the folder to store the normalized fcs files")

# Number of fcs files
nb.fcs <- length(fcs.files)

# Channels to normalize
channels <- 6:15 # FL's
marker.names <- c("CD58", "CD81", "CD34", "CD22", "CD38", "CD10", "CD19", "CD5", "CD20", "CD45") # B-ALL tube

### start ###

# assign values for the landmark for the peak of the negative population (for markers that are only positive, chose a landmark for the peak of the neg. pop)
reference.neg.pop.mode <- c(1.2, #FL1 CD58
                            3.5, # FL2 CD81 only pos
                            1.2, # FL3 CD34
                            1.8, # FL4 CD22 dim to pos
                            3, # FL5 CD38 dim to pos
                            1.2, # FL6 CD10
                            2.5, # FL7 CD19 only pos
                            1.2, # FL8 CD5
                            1.2, # FL9 CD20
                            2.3) # FL10 CD45

# Parameters for logicle and inverse transformation
logicle <- logicleTransform(w = 1, t = 1048576)
inv <- inverseLogicleTransform(trans = logicle)

index <- 1
for (fcs in 1:nb.fcs){
  # open fcs file --> store it as flowframe ff
  ff <- read.FCS(filename = fcs.files[fcs])
  # save channelnames of flowframe in variable
  channel.names <- grep('FL', colnames(ff), value = TRUE)
  # tranform flowframe using the logicle transformation (biexponential)
  ff.trf <- transform(ff, transformList(colnames(ff)[channels],logicleTransform(w = 1, t = 1048576)))
  # create density plots to check the logicle transformation
  windows(title = paste("Logicle transformation of sample", identifier(ff), sep = " "))
  par(mfrow= c(4,3), pty="s")
  nb.plots.per.windows <- 10
  name <- 0 # points out which markername to use as title for the plot
  for (FL in channels) {
    name <- name+1
    d1 <- density(exprs(ff.trf[,FL]))
    plot(d1, xlab = colnames(ff.trf[,FL]), xlim=c(-0.5,4.5), main= marker.names[name])
  }
  # select peaks
  # ff = flowframe, FL = fluorochrome
  name2 <- 1
  negative.population.mode <- c()
  # Select the peaks with the cursor
  for(FL in channels) {
    windows(title = paste("selection of negative peaks of sample", identifier(ff), sep = " "))
    x1 <- exprs(ff.trf[,FL]) # all fluorescence intensities of parameter i
    d1 <- density(x1) # compute density of these values
    plot(d1, xlab= colnames(exprs(ff.trf[,FL])), main= marker.names[name2], sub = "Select the negative peaks. For CD19, CD81 and CD38 select the most positive population, for CD45 select the middle population", cex.main=0.7, xlim=c(-0.5,4.5)) # plot these densities
    coord.1 <- locator(1, type = "l", col="red") # locator reads the position of the graphics cursor, max nb of points to locate, type l = line
    # Possibility to turn back when you clicked wrong?
    negative.population.mode[FL] <- coord.1$x # x-co�rdinate of the line is added for every parameter
    abline(v = negative.population.mode[FL],col="red",lwd=1,lty=5) #adds lines through the cuurrent plot --> v = xvalue vertical line
    name2 <- name2+1
    dev.off()
  }
  negative.population.mode <- negative.population.mode[!is.na(negative.population.mode)]
  name2 <- 1
  # Compute the new values
  for (FL in channels){
    new.value <- each_row(ff.trf[,FL], function(x) {x + (reference.neg.pop.mode[index] - negative.population.mode[index])})
    new.value <- as.matrix(new.value)
    if ( FL == 6){
      new <- new.value
    } else {
      new <- cbind(new, new.value)
    }
    index <- index+1
  }
  index <- 1
  # Assign the new values to the flowframe
  colnames(new) <- channel.names
  exprs(ff.trf)[,6:15] <- new
  
  # Inverse transform (necessary? better for output )
  ff.inv <- transform(ff.trf,transformList(colnames(ff.trf[channels]),inv))
  
  # Save flowset as separate fcs files in assigned output folder
  file.name <- paste0(output.folder,"\\normalized_", identifier(ff),".fcs")
  write.FCS(ff.inv, filename = file.name)
}


# Samples are now normalized 
# Check normalized samples
#check csv corrected
# for(i in c(1:nb.ff)){
#   windows(title = paste("Logicle transformation of sample", i, sep = " "))
#   par(mfrow= c(4,3), pty="s")
#   nb.plots.per.windows <- 10
#   aa <- 0
#   for (j in channels) {
#     if(aa==nb.plots.per.windows) {
#       windows(title = paste("Logicle transformation of sample", i, sep = " "))
#       par(mfrow= c(4,3), pty="s")
#       aa <- 0
#     }
#     aa <- aa+1
#     d1 <- density(exprs(trf.set[[i]][,j]))
#     plot(d1, xlab = colnames(trf.set[[i]][,j]), xlim=c(-1,5), main= markers[aa])
#     abline(v = negative.population.mode[i],col="red",lwd=1,lty=5)
#     abline(v = reference.neg.pop.mode[i],col="blue",lwd=1,lty=5)
#   }
# }