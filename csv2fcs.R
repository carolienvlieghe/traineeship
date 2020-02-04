################ Conversion script: CSV to FCS ################
###############################################################
# Save events from population of interest as CSV file in Kaluza

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

### Load libraries ###
library(XML)
library(methods)
library(flowCore)
library(xml2)
library(purrr)
library(rvest)
library(Biobase)

### Assign in & output ###
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S")
output.folder <- "D:/school/Stage officieel/Output/"
dir.create(path = output.folder)
input.folder <- "D:/school/Stage officieel/Input/"

# Zip and unzip the .analysis container to extract xml and fcs
new.file <- gsub("analysis","zip", fcs.path)
file.rename(from = fcs.path, to = new.file) # file is now a zipped folder
zip.path <- list.files(path= input.folder, pattern = "\\.zip$", full.names=TRUE)
zip.path.name <- basename(zip.path)
zip.path.name.1 <- substr(zip.path.name,start=1,stop = nchar(zip.path.name)-4)
file.name.zip <- substr(zip.path, start = 1, stop = nchar(zip.path)-4) 
dir.create(path = file.name.zip) # make folder to extract files from zipped folder
unzip(zip.path, exdir = file.name.zip) # unzip the folder

# List the .analysis, csv and fcs files in the input folder
fcs.path <- list.files(path= input.folder, pattern = "\\.analysis$", full.names=TRUE)
csv.path <- list.files(path= input.folder, pattern = "\\.csv$", full.names=TRUE)
file.name.fcs <- list.files(path = file.name.zip, pattern = "\\.fcs$", full.names = TRUE)

# assign the output for the new fcs file 
csv.base <- basename(csv.path)
base <- substr(csv.base, start = 1, stop = nchar(csv.base)-4) # keep only the name of the file without .csv
output.file <- paste(date.format, base, sep = " - ") # assign a new name to the new fcs file
output.file.path <- paste0(output.folder, output.file, ".fcs") # assign the whole path to the new fcs file


### Read the files ###
# Read the csv file
csv.file <- read.csv(csv.path) # class = dataframe

# Read fcs file dataset 1 in FCS 2.0 format (LMD from navios consists of 2 datasets, 1: FCS2.0 and 2: FCS3.0)
fcs2 <- read.FCS(file.name.fcs,dataset=1)


### Make flowFrame from csv ###
csvmatrix <- as.matrix(csv.file) # turn csv dataframe into a matrix 
ff <- new("flowFrame",exprs=csvmatrix) # make new flowFrame from csvmatrix

# Change rownames of annotated dataframe (look at it with: pData(parameters(flowFrame))) to fcs keyword $Pn
rows <- c()
for (i in c(1:16)) {
  rows <- append(rows, paste0("$P",i))
}
rownames(pData(parameters(ff))) <- rows


### Write to flowFrame to FCS###
write.FCS(x = ff, output.file.path)
