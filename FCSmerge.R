##################################
### Simple merge of flowframes ###
##################################

# The script concatenates the expression matrices without tags

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(rstudioapi)

# Select files to merge
input.folder <- choose.files(multi = TRUE, caption = "Select fcs files to merge")
answer <- winDialog(type = "yesno", "Do you want to select another file?")
# while loop makes it possible to select files from different directories
while (answer=="YES") {
  input.folder <- append(input.folder, choose.files(multi = TRUE, caption = "Select fcs files to merge"))
  answer <- winDialog(type = "yesno", "Do you want to select another file?")
}

# read set of fcs files to process
# the columnnames of the different expression matrices need to be the same
# this means only files run with the same panel can be merged 
set <- read.flowSet(files = input.folder)
nb.ff <- length(set)

for (ff in 1:nb.ff) {
  mtrx <- exprs(set[[ff]])
  if ( ff == 1){
    new.ff <- mtrx
  } else {
    new.ff <- rbind(new.ff, mtrx)
  }
} 

# replace random ff of set with new merged matrix
# This is an easy way to keep Annotated Dataframe
exprs(set[[1]]) <- new.ff
# Save the file in interactive way --> dialog window
file.name <- paste0(selectFile(caption = "Save file as ...", 
                               label = "Save", existing = FALSE)
                    , ".fcs")
# Write & save fcs file
write.FCS(x = set[[1]], filename = file.name)