################ Conversion script: CSV to FCS ################
###############################################################
# Save events from population of interest as CSV file in Kaluza
# Copy csv to inputfolder

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

### Load libraries ###
library(methods)
library(flowCore)
library(Biobase)

### Assign in & output ###
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S")
output.folder <- "D:/school/Stage officieel/csv_out/"
dir.create(path = output.folder)
input.folder <- "D:/school/Stage officieel/csv_in/"

# List the and csv files in the input folder (make loop to convert multiple csv files)
csv.path <- list.files(path= input.folder, pattern = "\\.csv$", full.names=TRUE)

# assign the output for the new fcs file 
csv.base <- basename(csv.path)
base <- substr(csv.base, start = 1, stop = nchar(csv.base)-4) # keep only the name of the file without .csv
output.file <- paste(date.format, base, sep = " - ") # assign a new name to the new fcs file
output.file.path <- paste0(output.folder, output.file, ".fcs") # assign the whole path to the new fcs file


### Read the files ###
# Read the csv file
csv.file <- read.csv(csv.path) # class = dataframe


### Make flowFrame from csv ###
# Replace . with whitespace in colnames of scatter parameters
col <- c()
for (i in c(1:5)) {
  new <- sub("\\.", " ", colnames(csv.file[i]))
  col <- append(col, new)
}
col <- append(col, colnames(csv.file[6:length(csv.file)]))
colnames(csv.file) <- col

# turn csv dataframe into a matrix
csvmatrix <- as.matrix(csv.file)  
# turn matrix into flowFrame
ff <- new("flowFrame",exprs=csvmatrix)

# number of variables
no_var <- c(1:length(parameters(ff)$name))
no_fl <- c(no_var[7]:length(no_var)-1)

# Change rownames of annotated dataframe (look at it with: pData(parameters(flowFrame))) to fcs keyword $Pn
rows <- c()
for (i in no_var) {
  rows <- append(rows, paste0("$P",i))
}
rownames(pData(parameters(ff))) <- rows

# Change fluorochrome names to match the names in kaluza
new <- as.character(parameters(ff)$name[1:5])

FL <- c()
a <- 1
for (i in no_fl){
  new <- append(new, paste0("FL",a," INT"))
  a <- a + 1
} 
new <- append(new, as.character(parameters(ff)$name[length(no_var)]))
new <- as.factor(new)
colnames(ff) <- new


### Write to flowFrame to FCS###
write.FCS(x = ff, output.file.path)
