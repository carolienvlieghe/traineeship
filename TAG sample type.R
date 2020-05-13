############################
######## TAG script ########
############################

## Script to assign a specific TAG to a file
## Purpose: ability to merge files in Kaluza and still identify sample


### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(ggcyto)
library(flowWorkspace)

# Assign in- & output
fcs.files <- choose.files(multi = TRUE, caption = "Select fcs files to tag")
output.folder <- choose.dir(caption = "Select the folder to store the tagged fcs files")

# Number of fcs files to process
nb.fcs <- length(fcs.files)

# Start
for (fcs in 1:nb.fcs){
  # open fcs file --> store it as flowframe ff
  ff <- read.FCS(filename = fcs.files[fcs])
  # Select the sample type
  # NBM = normal bone marrow
  # FU = follow up
  # Dg = diagnosis
  sample.type <- select.list(choices = c("NBM", "FU", "Dg"), 
                             preselect = NULL, multiple = FALSE, 
                             title = paste0("Select Sample type of ", identifier(ff)), 
                             graphics = TRUE)
  # Assign a preset value to the sample type
  if (sample.type == "NBM") {
    TAG.fcs <- as.matrix(rep(750, nrow(ff)))
  } else {
    if (sample.type == "FU") {
      TAG.fcs <- as.matrix(rep(500, nrow(ff)))
    } else {
      TAG.fcs <- as.matrix(rep(250, nrow(ff)))
    }
  }
  colnames(TAG.fcs) <- "TAG"
  # Make new column in the flowframe containing the tag to identify the sample type
  ff.new <- fr_append_cols(ff, TAG.fcs)
  description(ff.new)$`$P17R` <- 1024
  file.name <- paste0(output.folder, "\\sample_", fcs, ".fcs")
  write.FCS(ff.new, file.name)
}