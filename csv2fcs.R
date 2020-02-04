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

# List the .analysis and csv files in the input folder
fcs.path <- list.files(path= input.folder, pattern = "\\.analysis$", full.names=TRUE)
csv.path <- list.files(path= input.folder, pattern = "\\.csv$", full.names=TRUE)

# assign the output for the new fcs file 
csv.base <- basename(csv.path)
base <- substr(csv.base, start = 1, stop = nchar(csv.base)-4) # keep only the name of the file without .csv
output.file <- paste(date.format, base, sep = " - ") # assign a new name to the new fcs file
output.file.path <- paste0(output.folder, output.file, ".fcs") # assign the whole path to the new fcs file

### Read the files ###
# Read the csv file
csv.file <- read.csv(csv.path) # class = dataframe

### Make flowFrame ###
csvmatrix <- as.matrix(csv) # turn csv dataframe into a matrix 
ff <- new("flowFrame",exprs=csvmatrix) # make new flowFrame from csvmatrix

### Write to flowFrame to FCS###
write.FCS(x = ff, output.file.path)