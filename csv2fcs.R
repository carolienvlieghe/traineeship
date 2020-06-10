#####################################
### Conversion script: CSV to FCS ###
#####################################
# Save events from population of interest as CSV file in Kaluza
# Copy csv to inputfolder

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window

### Load libraries ###
library(flowCore)

### Assign in- & output folders ###
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M")
# Select the csv files you want to convert in the dialog window and list all files in csv.path variable
csv.path <- choose.files(multi = TRUE, caption = "Select csv files to convert")
# Open dialog window to choose output folder
output.folder <- choose.dir(caption = "Select folder to store new fcs files")

# Start loop to convert each csv file to fcs format
for (csv in csv.path) {
  # assign the output for the new fcs file 
  csv.base <- basename(csv)
  base <- substr(csv.base, start = 1, stop = nchar(csv.base)-4) # remove path & extension from filename
  output.file <- paste(date.format, base, sep = " - ") # assign a new name to the new fcs file
  output.file.path <- paste0(output.folder,"\\", output.file, ".fcs") # assign the whole path to the new fcs file
  
  
  ### Read the files ###
  # Read the csv file
  csv.file <- read.csv(csv) # class = dataframe

  ### Make flowFrame from csv ###
  # Replace . with whitespace in colnames of scatter parameters to correspond with names in Kaluza
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
  ff <- flowFrame(csvmatrix)
  
  # number of variables & fluorochromes
  nb.var <- c(1:length(parameters(ff)$name))
  nb.fl <- grep("CD", colnames(ff))
  
  # Look at annotated dataframe (look at it with: pData(parameters(flowFrame))) to fcs keyword $Pn
  # Change fluorochrome names to match the names in kaluza
  a <- 1
  for (i in nb.fl){
    colnames(ff)[i] <- paste0("FL", a, " INT")
    a <- a + 1
  } 

  # Resolution maximum Area, height and width specific for Navios
  adc.resolution.Area <- 20
  adc.resolution.Height <- 18
  adc.resolution.TIME <- 20
  adc.resolution.Width <- 10
  
  # Adjust range in Annotated Dataframe --> imported for visualization in Kaluza
  for (i in nb.var) {
    if ((ff@parameters@data$name[i]=="FS H") | (ff@parameters@data$name[i]=="FS-H") |
        (ff@parameters@data$name[i]=="FSC-H") | (ff@parameters@data$name[i]=="FSC H") |
        (ff@parameters@data$name[i]=="FS-PEAK") | (ff@parameters@data$name[i]=="FS PEAK")) {
      ff@parameters@data$range[i] <- as.character(2^adc.resolution.Height)
      ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Height)
      ff@parameters@data$minRange[i] <- 0
      range.kw <- paste0("$P",i,"R") # range keyword
      ff@description[range.kw] <- as.character(2^adc.resolution.Height) # add range
    } else {
      if (ff@parameters@data$name[i]=="TIME") {        ff@parameters@data$range[i] <- as.character(2^adc.resolution.TIME)
        ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.TIME)
        ff@parameters@data$minRange[i] <- 0
        range.kw <- paste0("$P",i,"R") # range keyword
        ff@description[range.kw] <- as.character(2^adc.resolution.TIME) # add range
      } else {
        if ((ff@parameters@data$name[i]=="FS W") | (ff@parameters@data$name[i]=="FS-W") |
            (ff@parameters@data$name[i]=="FSC-Width") | (ff@parameters@data$name[i]=="FSC Width") |
            (ff@parameters@data$name[i]=="FS-TOF") | (ff@parameters@data$name[i]=="FS TOF")) {
          ff@parameters@data$range[i] <- as.character(2^adc.resolution.Width)
          ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Width)
          ff@parameters@data$minRange[i] <- 0
          range.kw <- paste0("$P",i,"R") # range keyword
          ff@description[range.kw] <- as.character(2^adc.resolution.Width) # add range
        } else {
          if ((ff@parameters@data$name[i]=="SS H") | (ff@parameters@data$name[i]=="SS-H") |
              (ff@parameters@data$name[i]=="SSC-H") | (ff@parameters@data$name[i]=="SSC H") |
              (ff@parameters@data$name[i]=="SS-PEAK") | (ff@parameters@data$name[i]=="SS PEAK")) {
            ff@parameters@data$range[i] <- as.character(2^adc.resolution.Height)
            ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Height)
            ff@parameters@data$minRange[i] <- 0
            range.kw <- paste0("$P",i,"R") # range keyword
            ff@description[range.kw] <- as.character(2^adc.resolution.Height) # add range
          } else {
            if ((ff@parameters@data$name[i]=="SS W") | (ff@parameters@data$name[i]=="SS-W") |
                (ff@parameters@data$name[i]=="SSC-Width") | (ff@parameters@data$name[i]=="SSC Width") |
                (ff@parameters@data$name[i]=="SS-TOF") | (ff@parameters@data$name[i]=="SS TOF")) {
              ff@parameters@data$range[i] <- as.character(2^adc.resolution.Width)
              ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Width)
              ff@parameters@data$minRange[i] <- 0
              range.kw <- paste0("$P",i,"R") # range keyword
              ff@description[range.kw] <- as.character(2^adc.resolution.Width) # add range
            } else {
              ff@parameters@data$range[i] <- as.character(2^adc.resolution.Area)
              ff@parameters@data$maxRange[i] <- as.character(2^adc.resolution.Area)
              ff@parameters@data$minRange[i] <- -346
              range.kw <- paste0("$P",i,"R") # range keyword
              ff@description[range.kw] <- as.character(2^adc.resolution.Area) # add range
            }
          }
        }
      }
    }
  }
  
  # Save flowframes as fcs file 
  write.FCS(x = ff, output.file.path)
}