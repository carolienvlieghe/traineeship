############################
### Normalization script ###
############################

# Normalize the normal bone marrow samples before merging

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)

# Assign in- & output
fcs.files <- choose.files(multi = TRUE, caption = "Select fcs files to normalize")
output.folder <- choose.dir(caption = "Select the folder to store the normalized fcs files")

# Number of fcs files to process
nb.fcs <- length(fcs.files)

### start ###

# Assign values as landmark for the peak of the negative population (for markers that are only positive, choose a landmark for the positive peak)
reference.neg.pop.mode <- c(1.05, #FL1 CD58
                            3.4, # FL2 CD81 only pos
                            1.06, # FL3 CD34
                            1.75, # FL4 CD22 dim to pos
                            3.02, # FL5 CD38 dim to pos
                            1.02, # FL6 CD10
                            2.62, # FL7 CD19 only pos
                            1.05, # FL8 CD5
                            1.25, # FL9 CD20
                            2.48) # FL10 CD45

# Define the logicle and inverse logicle function in a variable
logicle <- logicleTransform(w = 0.9, t = 1048576)
inv <- inverseLogicleTransform(trans = logicle)

# Popup window with explanation
winDialog(type = "ok", "Select the  negative peaks for each channel with the cursor.\n For CD81, CD38 & CD19 select the positive peak.\n For CD45 select the middle peak.\n The blue line refers to the reference peak")

# Loop over every fcs file to perform normalization
for (fcs in 1:nb.fcs){
  index <- 1
  # open fcs file --> store it as flowframe ff
  ff <- read.FCS(filename = fcs.files[fcs])
  # get index for channels to transform (FL)
  FL.index <- grep("FL", colnames(ff))
  # save (pretty: CD) channelnames of flowframe in variable
  FL.names <- as.character(parameters(ff)$desc[FL.index])
  # tranform flowframe using the logicle transformation (biexponential)
  ff.trf <- transform(ff, transformList((colnames(ff)[FL.index]), logicle))
  # create density plots to check the logicle transformation
  windows(title = paste("Logicle transformation of sample", identifier(ff), sep = " "))
  par(mfrow= c(4,3), pty="s")
  nb.plots.per.windows <- 10
  name <- 0 # points out which markername to use as title for the plot
  for (FL in FL.index) {
    name <- name+1
    d1 <- density(exprs(ff.trf[,FL]))
    plot(d1, xlab = colnames(ff.trf[,FL]), xlim=c(-0.5,4.5), main= FL.names[name])
  }
  # select peaks
  # ff = flowframe, FL = fluorochrome
  name <- 1
  negative.population.mode <- c()
  # Select the peaks with the cursor
  for(FL in FL.index) {
    windows(title = paste("Sample:", identifier(ff), sep = " "))
    x1 <- exprs(ff.trf[,FL]) # all values of parameter FL
    d1 <- density(x1) # compute density of these values
    plot(d1, xlab= colnames(exprs(ff.trf[,FL])), main= FL.names[name], sub = "Select the negative/positive peaks", cex.main=0.7, xlim=c(-0.5,4.5)) # plot these densities
    abline(v = reference.neg.pop.mode[name],col="blue",lwd=1,lty=5) # add blue line corresponding with reference population
    coord.1 <- locator(1, type = "l", col="red") # locator reads the position of the graphics cursor, max nb of points to locate, type l = line
    # Possibility to turn back when you clicked wrong?
    negative.population.mode[FL] <- coord.1$x # x-coordinate of the line is added for every parameter
    abline(v = negative.population.mode[FL],col="red",lwd=1,lty=5) #adds lines through the cuurrent plot --> v = xvalue vertical line
    name <- name+1
    dev.off()
  }
  negative.population.mode <- negative.population.mode[!is.na(negative.population.mode)] #NA's are removed
  # Compute the new values
  for (FL in FL.index){
    new.value <- each_row(ff.trf[,FL], function(x) {x + (reference.neg.pop.mode[index] - negative.population.mode[index])})
    new.value <- as.matrix(new.value)
    if ( FL == 6){
      new <- new.value
    } else {
      new <- cbind(new, new.value)
    }
    index <- index+1
  }
  # Replace the old values with the new values in the expression matrix of the flowframe
  colnames(new) <- colnames(ff.trf)[FL.index]
  exprs(ff.trf)[,FL.index] <- new
  
  # Inverse transform (necessary? better for output)
  ff.inv <- transform(ff.trf,transformList(colnames(ff.trf[,FL.index]),inv))
  
  # Save flowset as separate fcs files in assigned output folder
  file.name <- paste0(output.folder,"\\normalized_", identifier(ff),".fcs")
  write.FCS(ff.inv, filename = file.name)
}
