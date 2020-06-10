####################################
######## TAG & merge script ########
####################################

# Assign a TAG to each dataset to identify the sample type in Kaluza
# 250 = Diagnosis
# 500 = Follow up
# 750 = Normal Bone Marrow

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(rstudioapi)

# Select fcs files according to the filetype
sample.type <- c("NBM","FU","Dg")
fcs.file <- list()
fcs.list <- list()
# Loop over the sample types to select the corresponding fcs files through a dialog window
for (i in sample.type) {
  fcs.file <- list(choose.files(multi = FALSE, caption = paste("Select", i, "file", sep = " ")))
  names(fcs.file) <- c(i)
  if (fcs.file=="character(0)") {
    # when no file is uploaded for a certain sample type (cancel)
    # nothing has to happen, proceed to next item in loop
  } else {
    fcs.list <- append(fcs.list, fcs.file)
  }
}

# Number of fcs files to process
nb.fcs <- length(fcs.list)

# Start (loop over every selected file)
for (fcs in 1:nb.fcs){
  # Read fcs file, store it in flowframe ff
  ff <- read.FCS(filename = as.character(fcs.list[fcs]))
  # Assign the TAG for the sample types
  if (names(fcs.list[fcs])=="NBM") {
    TAG.fcs <- as.matrix(rep(750, nrow(ff)))
  } else {
    if (names(fcs.list[fcs]) == "FU") {
      TAG.fcs <- as.matrix(rep(500, nrow(ff)))
    } else {
      TAG.fcs <- as.matrix(rep(250, nrow(ff)))
    }
  }
  colnames(TAG.fcs) <- "TAG"
  # Make new column in the flowframe containing the tag to identify the sample type
  ff.new <- fr_append_cols(ff, TAG.fcs)
  # Add range keyword & add range info to Annotated dataframe --> important for visualization in Kaluza
  range.kw <- paste0("$P",ncol(ff.new),"R")
  ff.new@description[range.kw] <- '1024'
  ff.new@parameters@data$range[ncol(ff.new)] <- '1024'
  ff.new@parameters@data$maxRange[ncol(ff.new)] <- '1024'
  ff.new@parameters@data$minRange[ncol(ff.new)] <- '1024'
  
  # Merge the data of every file
  mtrx <- exprs(ff.new)
  if ( fcs == 1) {
    column.names <- colnames(ff.new)
    merged.ff <- mtrx
  } else {
    colnames(mtrx) <- column.names
    merged.ff <- rbind(merged.ff, mtrx)
  }
}

# Replace expression matrix of the flowframe with updated matrix containing the tags
# Easy way to maintain phenodata of original files:
colnames(ff.new) <- column.names
exprs(ff.new) <- merged.ff

# Save the file in an interactive way --> dialog window
file.name <- file.name <- paste0(selectFile(caption = "Save file as ...", 
                                            label = "Save", existing = FALSE), ".fcs")
write.FCS(ff.new, file.name)
