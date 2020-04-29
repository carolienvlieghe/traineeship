############################
### Simple merge of flowframes ###
############################
### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)

# Assign in- & output
input.folder <- "D:/path_to_inputfolder/"
output.folder <- "D:/path_to_outputfolder/"
dir.create(path = output.folder)

# read set of fcs files to process
set <- read.flowSet(files = NULL, path = input.folder)
# channel.names <- grep('FL', colnames(set), value = TRUE)
nb.ff <- length(set)

for (ff in 1:nb.ff) {
  mtrx <- exprs(set[[ff]])
  if ( ff == 1){
    new.ff <- mtrx
  } else {
    new.ff <- rbind(new.ff, mtrx)
  }
} 

# save merged matrix in random ff of set and write FCS --> to keep phenodata
exprs(set[[1]]) <- new.ff
write.FCS(x = set[[1]], filename = paste0(output.folder, "15_NBM.fcs"))
