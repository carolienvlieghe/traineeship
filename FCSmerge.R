##################################
### Simple merge of flowframes ###
##################################

# The script concatenates the expression matrices

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

# Load packages
library(flowCore)
library(rstudioapi)

# Select files to merge, use while loop to be able to select files from different directories
input.folder <- choose.files(multi = TRUE, caption = "Select fcs files to merge")
answer <- winDialog(type = "yesno", "Do you want to select another file?")
while (answer=="YES") {
  input.folder <- append(input.folder, choose.files(multi = TRUE, caption = "Select fcs files to merge"))
  answer <- winDialog(type = "yesno", "Do you want to select another file?")
}

# read set of fcs files to process
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

# save merged matrix in random ff off set (easy way to keep Annotated Dataframe) and write FCS
exprs(set[[1]]) <- new.ff
# ask for filename
file.name <- paste0(selectFile(caption = "Save file as ...", 
                               label = "Save", existing = FALSE)
                    , ".fcs")
write.FCS(x = set[[1]], filename = file.name)
