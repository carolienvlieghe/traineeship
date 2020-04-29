################ Conversion script: CSV to FCS ################
###############################################################
# Save events from population of interest as CSV file in Kaluza
# Copy csv to inputfolder

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

### Load libraries ###
library(flowCore)

### Assign in & output folders ###
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S")
# Select the csv files you want to convert in the dialog window and list all files in csv.path variable
csv.path <- choose.files(multi = TRUE, caption = "Select csv files to convert")
# Open dialog window to chose output folder
output.folder <- choose.dir(caption = "Select folder to store new fcs files")
dir.create(path = output.folder)

# Start loop to convert each csv file to fcs format
for (csv in csv.path) {
  # assign the output for the new fcs file 
  csv.base <- basename(csv)
  base <- substr(csv.base, start = 1, stop = nchar(csv.base)-4) # keep only the name of the file without .csv
  output.file <- paste(date.format, base, sep = " - ") # assign a new name to the new fcs file
  output.file.path <- paste0(output.folder, output.file, ".fcs") # assign the whole path to the new fcs file
  
  
  ### Read the files ###
  # Read the csv file
  csv.file <- read.csv(csv) # class = dataframe
  # Kaluza scales the values to 1024 before displaying
  # MFI/1024 corresponds with original MFI in pre-analysis

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
  
  # Adjust ranges
  adc.resolution.Area <- 20
  adc.resolution.Height <- 18
  adc.resolution.TIME <- 20
  adc.resolution.Width <- 10
  
  for (i in no_var) {
    if ((ff@parameters@data$name[i]=="FS H") | (ff@parameters@data$name[i]=="FS-H") |
        (ff@parameters@data$name[i]=="FSC-H") | (ff@parameters@data$name[i]=="FSC H") |
        (ff@parameters@data$name[i]=="FS-PEAK") | (ff@parameters@data$name[i]=="FS PEAK")) {
      ff@parameters@data$range[i] <- as.character(2^adc.resolution.Height)
      ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Height)
      ff@parameters@data$minRange[i] <- 0
    } else {
      if (ff@parameters@data$name[i]=="TIME") {
        ff@parameters@data$range[i] <- as.character(2^adc.resolution.TIME)
        ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.TIME)
        ff@parameters@data$minRange[i] <- 0
      } else {
        if ((ff@parameters@data$name[i]=="FS W") | (ff@parameters@data$name[i]=="FS-W") |
            (ff@parameters@data$name[i]=="FSC-Width") | (ff@parameters@data$name[i]=="FSC Width") |
            (ff@parameters@data$name[i]=="FS-TOF") | (ff@parameters@data$name[i]=="FS TOF")) {
          ff@parameters@data$range[i] <- as.character(2^adc.resolution.Width)
          ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Width)
          ff@parameters@data$minRange[i] <- 0
        } else {
          if ((ff@parameters@data$name[i]=="SS H") | (ff@parameters@data$name[i]=="SS-H") |
              (ff@parameters@data$name[i]=="SSC-H") | (ff@parameters@data$name[i]=="SSC H") |
              (ff@parameters@data$name[i]=="SS-PEAK") | (ff@parameters@data$name[i]=="SS PEAK")) {
            ff@parameters@data$range[i] <- as.character(2^adc.resolution.Height)
            ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Height)
            ff@parameters@data$minRange[i] <- 0
          } else {
            if ((ff@parameters@data$name[i]=="SS W") | (ff@parameters@data$name[i]=="SS-W") |
                (ff@parameters@data$name[i]=="SSC-Width") | (ff@parameters@data$name[i]=="SSC Width") |
                (ff@parameters@data$name[i]=="SS-TOF") | (ff@parameters@data$name[i]=="SS TOF")) {
              ff@parameters@data$range[i] <- as.character(2^adc.resolution.Width)
              ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Width)
              ff@parameters@data$minRange[i] <- 0
            } else {
              ff@parameters@data$range[i] <- as.character(2^adc.resolution.Area)
              ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Area)
              ff@parameters@data$minRange[i] <- 0
            }
          }
        }
      }
    }
  }
  # try to extract this from range info and make it into loop
  description(ff)$`$P1R` <- "262144"
  description(ff)$`$P2R` <- "1048576"
  description(ff)$`$P3R` <- "1024"
  description(ff)$`$P4R` <- "1048576"
  description(ff)$`$P5R` <- "1024"
  description(ff)$`$P6R` <- "1048576"
  description(ff)$`$P7R` <- "1048576"
  description(ff)$`$P8R` <- "1048576"
  description(ff)$`$P9R` <- "1048576"
  description(ff)$`$P10R` <- "1048576"
  description(ff)$`$P11R` <- "1048576"
  description(ff)$`$P12R` <- "1048576"
  description(ff)$`$P13R` <- "1048576"
  description(ff)$`$P14R` <- "1048576"
  description(ff)$`$P15R` <- "1048576"
  description(ff)$`$P16R` <- "1048576"
  ### Write to flowFrame to FCS###
  write.FCS(x = ff, output.file.path)
}