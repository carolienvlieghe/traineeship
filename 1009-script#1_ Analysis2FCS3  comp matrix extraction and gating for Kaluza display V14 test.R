###############################################################################################
# This script is able to create a fcs3.0 from a Kaluza analysis file. 
# The modified compensation matrix from the Klauza analysis is paste into the new fcs3.0 file
# This script works only with analysis file create with Kaluza v2.0. It will not work with 
# analysis file created with version Kaluza 1.3 or Kaluza 1.5
###############################################################################################


rm(list=ls())
cat("\014")
if(!is.null(dev.list())) dev.off()


library(XML)
library(methods)
library(flowCore)
library(xml2)
library(purrr)
library(rvest)


#create an output folder where all new fcs files will be sent

date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S-")
output.folder <- "C:/Users/carol/OneDrive/Documenten/school/Stage officieel/Extraction_output/"
dir.create(path = output.folder)


#C:\Users\BEDUPONT\Documents\Demo Kaluza 1.2\304-NL\analysis files
# rename .analysis  in .zip
input.folder <- "C:/Users/carol/OneDrive/Documenten/school/Stage officieel/Extraction_input/"



#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-



fcs.path <- list.files(path= input.folder, pattern = "\\.analysis$", full.names=TRUE)

nb.files <- length(fcs.path)

fcs.path.new <- list()
for (i in 1:nb.files){
  old.file <- fcs.path[[i]]
  new.file <- gsub("analysis","zip", old.file)
  file.rename(from = old.file, to = new.file)
  
}
#@END rename analysis file in zip


###############################################################################################
#loop to extract the fcs 3.0 and xml slot
#  + create the fcs file
###############################################################################################

zip.path <- list.files(path= input.folder, pattern = "\\.zip$", full.names=TRUE)

zip.path.name <- basename(zip.path)
zip.path.name.1 <- substr(zip.path.name,start=1,stop = nchar(zip.path.name)-4)

nb.files.zip <- length(zip.path)


###############################################################################################

for (a in 1:nb.files.zip){
  
  file.name.zip.1 <- substr(zip.path[[a]], start = 1, stop = nchar(zip.path[[a]])-4)
  file.name.zip <- file.name.zip.1
  dir.create(path = file.name.zip)
  unzip(zip.path[[a]], exdir = file.name.zip)
  
  #Extract the xml solt to get the compensation matrix
  file.name.xml <- list.files(path = file.name.zip, pattern = "\\.xml$", full.names = TRUE)
  
  doc <- read_xml(file.name.xml)
  result <- xmlParse(file = file.name.xml, useInternalNodes = FALSE, 
                                ignoreBlanks=FALSE,
                                replaceEntities = TRUE,
                                asTree = TRUE,
                                isSchema = FALSE)
  
  ######
  doc.berlin <- xmlInternalTreeParse(file.name.xml)
  poc <-xmlRoot(doc.berlin)
  p <- poc[[1]]
  PanelList <- xmlApply(p[[1]], function(x) xmlGetAttr(x[["DataSet"]], "N"))
  ######
  
  
  
  #count the number of fcs files in the zip file
  files.name.fcs <- list.files(path= file.name.zip, pattern = "\\.fcs$", full.names=TRUE)
  
  nb.files.fcs.zip <- length(files.name.fcs)
  
  meta <- list()
  
  for (aa in 1:nb.files.fcs.zip){
  
  fcs.test <- read.FCS(files.name.fcs[[aa]],dataset=1)
  
  
  #@START comp matrix extraction for any instrument
  #detector name extraction with xml2 package
  doc.1 <- xml_find_all(doc, ".//Compensation/CompensatedDetectors" )
  doc.2 <- xml_find_all(doc.1[aa], ".//D" )
  doc.3 <- html_text(doc.2) #comp.matrix.parameter.name
  nb.parameter.comp.matrix <- length(doc.3)
  
  number.coef <- (nb.parameter.comp.matrix^2) - nb.parameter.comp.matrix
  
  #compensation matrix extraction with xml2 package
  
  doc.10 <- xml_find_all(doc, ".//Compensation/S" )
  doc.11 <- xml_find_all(doc.10[aa], ".//SV" )
  doc.13 <- html_attrs(doc.11)
  SV.table <- do.call(rbind, doc.13)
  SV.table.1 <- SV.table
  
  SV.table.1.row <- SV.table.1[,1]
  indices.matrix.row <- match(SV.table.1.row, doc.3)
  #is.na(indices.matrix.row) <- 0
  
  SV.table.1.col <- SV.table.1[,2]
  indices.matrix.col <- match(SV.table.1.col, doc.3)
  #is.NA.indices.matrix.col <- is.na(indices.matrix.col)
  
  #indices.numeric.coef <- !is.na(indices.matrix.row) & !is.na(indices.matrix.col)
  numeric.coef <- as.numeric(SV.table.1[,3])
  #numeric.1 <- numeric.coef[indices.numeric.coef]
  
  #SV.table.3 <- cbind(indices.matrix.row[indices.numeric.coef],
  #                    indices.matrix.col[indices.numeric.coef], 
  #                    numeric.1)
  
  SV.table.3 <- cbind(indices.matrix.row,
                      indices.matrix.col, 
                      numeric.coef)
  
  
  #Building the comp matrix
  #build a A x A zero diag matrix
  comp.matrix <- diag(x=1, nrow= nb.parameter.comp.matrix, ncol = nb.parameter.comp.matrix)
  
  #matrix Kaluza format
  for (z in 1: nrow(SV.table.3)){
    row <- SV.table.3[z,1]
    column <- SV.table.3[z,2]
    comp.matrix[row, column] <- SV.table.3[z,3]
  }
  
  #matrix R + fcs file format
  matrix.coef <- comp.matrix.fcs <- t(comp.matrix)
  
  comp.matrix.parameter.name <- doc.3
  
  
  ###Method 2 - comp matrix extraction
  comp.berlin <- getNodeSet(p[[1]][[aa]], "Protocol/Compensation", addFinalizer=TRUE)
  comp.berlin.1 <- getNodeSet(p[[1]][[aa]], "Protocol/Compensation/S/SV", addFinalizer=TRUE)
  
  
###############################################################################################
#                                      Berlin                                                 #
###############################################################################################  
  
    #Extration of the main gate : must be on linear parameter FSC vs SSC by example.
  #@STRART the good one !!!
  gates.berlin <- getNodeSet(p[[1]][[aa]], "Protocol/Gates", addFinalizer=TRUE)
  
  # M="FL1 INT" : parameter; S="C" is Logicle;  S="O" is Log; S="I" is Linear      i<-1
  # B = Boolean; R = rectangle ; P= polygon ; L = histogram
  gate.mat <- lapply(1:xmlSize(gates.berlin[[1]]), function(i){
    if (xmlName(gates.berlin[[1]][[i]])=="B") {xmlValue(gates.berlin[[1]][[i]][["ExpressionText"]])} 
    else {
      if (xmlName(gates.berlin[[1]][[i]])=="R") {
        matrix(as.numeric(unlist(strsplit((matrix(xmlSApply(gates.berlin[[1]][[i]][["P"]], 
          function(x) xmlGetAttr(x,"O")),nrow=2,ncol=1)),","))),nrow=2,ncol=2,byrow=TRUE,
          dimnames=list(c("min","max"),as.list(xmlSApply(gates.berlin[[1]][[i]][["A"]], 
          function(x) xmlGetAttr(x,"M"))))) }
      else {
      if (xmlName(gates.berlin[[1]][[i]])=="P") {
        matrix(as.numeric(unlist(strsplit((matrix(xmlSApply(gates.berlin[[1]][[i]][["P"]], 
         function(x) xmlGetAttr(x,"O")),nrow=xmlSize(gates.berlin[[1]][[i]][["P"]]),ncol=1)),","))),
         nrow=xmlSize(gates.berlin[[1]][[i]][["P"]]),ncol=2,byrow=TRUE,
         dimnames=list(as.list(xmlSApply(gates.berlin[[1]][[i]][["P"]],xmlName)),as.list(xmlSApply(gates.berlin[[1]][[i]][["A"]], 
         function(x) xmlGetAttr(x,"M"))))) }
        else {
        if (xmlName(gates.berlin[[1]][[i]])=="L") {
          extraction.1 <-  matrix(xmlSApply(gates.berlin[[1]][[i]][["P"]], 
                                            function(x) xmlGetAttr(x,"O")),nrow=2,ncol=1) 
          
          extraction.2 <- strsplit(extraction.1, ",")
          extraction.3 <- lapply(extraction.2, "[[",1)
          extraction.4 <- unlist(as.numeric(extraction.3))
          extraction.5 <- matrix(extraction.4, nrow = 2, ncol = 1,
                                 dimnames=list(c("min","max"),
                                               as.list(xmlSApply(gates.berlin[[1]][[i]][["A"]], 
                                                      function(x) xmlGetAttr(x,"M")))))
          #parameter name extraction
          #extraction.6 <- xmlSApply(gates.berlin[[1]][[i]][["A"]], function(x) xmlGetAttr(x,"M"))
          return(extraction.5)
        }
          
          #matrix(as.numeric(unlist(strsplit((matrix(xmlSApply(gates.berlin[[1]][[i]][["P"]], 
           # function(x) xmlGetAttr(x,"O")),nrow=2,ncol=1)),","))),
            #nrow=2,ncol=2,byrow=TRUE,
            #dimnames=list(c("min","max"),as.list(xmlSApply(gates.berlin[[1]][[i]][["A"]], 
            #function(x) xmlGetAttr(x,"M"))))) }
        else {
          "NA"; warning("Ellipsoid gate not yet supported")
          if (xmlName(gates.berlin[[1]][[i]])=="E") {
            
            
            #Ellipse coordinates extraction - Center Points
            extraction.1 <-  matrix(xmlSApply(gates.berlin[[1]][[i]][["P"]], 
                                              function(x) xmlGetAttr(x,"O")),nrow=1,ncol=1) 
            
            extraction.center.points <- as.numeric(unlist(strsplit(extraction.1, ","))) #Xc Yc
            
            
            # a & b Coefficient extraction
            extraction.coef.b <-  as.numeric((xmlSApply(gates.berlin[[1]][[i]][["Coefficients"]], 
                                              function(x) xmlValue(x))[1])) #a
            
            extraction.coef.a <-  as.numeric((xmlSApply(gates.berlin[[1]][[i]][["Coefficients"]], 
                                                        function(x) xmlValue(x))[2])) #c
            
            a.b.vector <- c(extraction.coef.a, extraction.coef.b)
            
            #angle extraction
            extraction.angle <- -as.numeric(xmlGetAttr(gates.berlin[[1]][[i]],"L")) * (pi/180)
            
            if(length(extraction.angle)==0) { extraction.angle <- 0}
            
            angle.ellipse <- c(extraction.angle, 0)
            
            #parameter names extraction
            extraction.parameter.1 <- xmlSApply(gates.berlin[[1]][[i]][["A"]], 
                                                function(x) xmlGetAttr(x,"M"))
            
            #matrix with points
            extraction.final.matrix <- rbind(extraction.center.points, a.b.vector, angle.ellipse)
            colnames(extraction.final.matrix) <- extraction.parameter.1
            
            return(extraction.final.matrix)
            
          }
         }
        }
       }
      }
    })
  
#####################################################
  
  #Extraction of region name and region matrix from all Rectangular and Polygonal region
  region.name <- xmlSApply(gates.berlin[[1]], xmlGetAttr, "N")
  
  region.name.index.R.P <- which(names(region.name) != "B")
  
  region.name.R.P <- region.name[region.name.index.R.P]
  region.matrix.R.P <- gate.mat[region.name.index.R.P]
  
  tempo.table <- rep(0,length(region.name.R.P))
  region.table <- cbind(names(region.name.R.P), unname(region.name.R.P), as.numeric(tempo.table))
  colnames(region.table) <- c("Region type","Regions", "Selection 0 or 1")
  
  #this code to be able to copy the gating strategu to all fcs files from the analasis 
  if(aa==1) {region.table <- edit(region.table, title= "Enter here your gating strategy") ####!!!!!!!!!!!
             region.table.file <- paste(output.folder, "region.table.rds", sep="")
             saveRDS(region.table, file= paste(output.folder, "region.table.rds", sep=""))}
  
  region.table <- readRDS(region.table.file)
  
  #this code to enable a gating strategy for each files from the analysis
  #region.table <- edit(region.table, title= "Enter here your gating strategy")
  
  warning("!! Only use Rectangle, Polygone or Histogram regions to define your populations in Kaluza !! ")
  warning("!! Do not draw any region in Radar Plot - it will crash the script !!")
  #colnames of each region.matrix.R.P have to be rename with the FCS3.0 parameter name !
  
###############################################################################################  

  tempo.Navios <- (fcs.test@description$`$CYT` == "Gallios") |
                  (fcs.test@description$`$CYT` == "Navios") |
                  (fcs.test@description$`$CYT` == "Navios EX") |
                  (fcs.test@description$`$CYT` == "Cytomics FC 500")
  
  tempo.GalliosK <- fcs.test@description$`$CYT` == "Gallios (Kaluza)"
  
  tempo.BD <- substr(fcs.test@description$`CREATOR`, start = 1, stop = 11) #== "BD FACSDiva"
  
  tempo.CytoFLEX <- (fcs.test@description$`$CYT` == "CytoFLEX") |
             (fcs.test@description$`$CYT` == "CytoFLEX LX") | 
             (fcs.test@description$`$CYT` == "CytoFLEX S")
  
  
  
  #@START script for Gallios : Navios / FC500
  if (any(tempo.Navios) == TRUE){
    
    #Open the fcs3.0 slot and save the fcs with it's new comp matrix
  fcs <- read.FCS(files.name.fcs[[aa]],dataset=2) #open slot fcs3.0
  fcs2 <- read.FCS(files.name.fcs[[aa]],dataset=1) #open slot fcs2.0
  
  
  #########
  fcs.parameters.data.new.name <- c()
  fcs.description.name<-unname(fcs@parameters@data$name)
  #fcs.description.desc<-unname(fcs2@parameters@data$name)
  
  fcs2.description.name<-unname(fcs2@parameters@data$name)
  
  #number of parameters
  np <- length(fcs.description.name)
  np2 <- length(fcs2.description.name)
  
  #rename parameter name as <parameter name>-parameter description because write.FCS writes only the parameter name !
  for (k in 1:(np)){
    for (l in 1:np2){
      
      if(fcs2@description$`$CYT` == "Cytomics FC 500") {
        name.fcs3 <- substr(fcs@parameters@data$name[k],start=1,stop = nchar(fcs@parameters@data$name[k])-0)
        name.fcs2 <-substr(fcs2@parameters@data$name[l],start=1,stop = nchar(fcs2@parameters@data$name[l])-4)
      }
      else{
      name.fcs3 <- substr(fcs@parameters@data$name[k],start=1,stop = nchar(fcs@parameters@data$name[k])-2)
      name.fcs2 <-substr(fcs2@parameters@data$name[l],start=1,stop = nchar(fcs2@parameters@data$name[l])-8)
      }
      
      if(fcs@parameters@data$name[k]=="FS-H"){
        fcs@parameters@data$desc[k]<-"FS PEAK"
        fcs.parameters.data.new.name[k]<-"FS PEAK"
      }
      else{
        if((fcs@parameters@data$name[k]=="FS-A") | (fcs@parameters@data$name[k]=="FS")){
          fcs@parameters@data$desc[k]<-"FS INT"
          ifelse(fcs2@description$`$CYT` == "Cytomics FC 500" , 
                 fcs.parameters.data.new.name[k]<-"FS", 
                 fcs.parameters.data.new.name[k]<-"FS INT")
        }
        else{
          if(fcs@parameters@data$name[k]=="FS-W"){
            fcs@parameters@data$desc[k]<-"FS-W"
            fcs.parameters.data.new.name[k]<-"FS TOF"
          }
          else{
            if(fcs@parameters@data$name[k]=="SS-H"){
              fcs@parameters@data$desc[k]<-"SS PEAK"
              fcs.parameters.data.new.name[k]<-"SS PEAK"
            }
            else{
              if((fcs@parameters@data$name[k]=="SS-A") | (fcs@parameters@data$name[k]=="SS")){
                fcs@parameters@data$desc[k]<-"SS INT"
                ifelse(fcs2@description$`$CYT` == "Cytomics FC 500" , 
                       fcs.parameters.data.new.name[k]<-"SS", 
                       fcs.parameters.data.new.name[k]<-"SS INT")
              }
              else{
                if(fcs@parameters@data$name[k]=="SS-W"){
                  fcs@parameters@data$desc[k]<-"SS-W"
                  fcs.parameters.data.new.name[k] <- "SS TOF"
                }
                else{
                  if(fcs@parameters@data$name[k]=="TIME"){
                    fcs@parameters@data$desc[k]<-"TIME"
                    fcs.parameters.data.new.name[k] <- "TIME"
                  }
                  else{
                    if (name.fcs3==name.fcs2){
                      fcs@parameters@data$desc[k]<- fcs2@parameters@data$desc[l]
                      fcs.parameters.data.new.name[k]<- substr(fcs2@parameters@data$name[l],
                                                               start=1,
                                                               stop = nchar(fcs2@parameters@data$name[l])-4)
                    }
                    else{
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  old.fcs.parameter.name <- unname(fcs@parameters@data$name)
  
  #fcs@parameters@data$name <- fcs.parameters.data.new.name
  
  fcs.name <- unname(fcs.parameters.data.new.name)
  fcs.desc <- unname(fcs@parameters@data$desc)
  
  fcs.description.name <- unname(fcs.parameters.data.new.name)
  fcs.description.desc <- unname(fcs@parameters@data$desc)
  
  
  #############################################################################################
  #Compensate the fcs file with the new matrix
  if(fcs2@description$`$CYT` == "Cytomics FC 500") {
  tempo.comp <- grep("^FL[0-9]",old.fcs.parameter.name)
  tempo2.comp <- which(old.fcs.parameter.name%in%comp.matrix.parameter.name)
  comp.matrix.fcs.name <- old.fcs.parameter.name[tempo2.comp]
  }
  else{
  comp.matrix.fcs.name <- paste(comp.matrix.parameter.name, "-A",sep="")  
  }
  colnames(comp.matrix.fcs) <- comp.matrix.fcs.name
  fcs.comp <- compensate(fcs, comp.matrix.fcs)
  
  
  
  #############################################################################################
  #create biexponentiel ff to be able to gate on
  #target.column <- c(5:14) #
  ##############################################################
  
  ###############################################################
  # Enter here the individual optimized m coefficient for each ##
  #                          parameter                          #
  #lgcl.coef<-c(4.9, #FITC
  #             4.2, # PE
   #            4.2, # ECD
    #           4.0, # PC55
     #          4.1, # PC7
      #         4.4, # APC
       #        4.2, # A700 4.2
        #       4.9, # AA750
         #      4.5, # PB 4.6
          #     4.4, # krO
           #    4.4,
            #   4.4)
  ###############################################################
  
  #y<-0
  #ff.logicle.gated <- fcs.comp
  
  #for(z in target.column){
   # y<-y+1
    #lgcl <- logicleTransform( w =0.9, t= 2^20, m =lgcl.coef[y])  
    ## Remark:  low W=large peak // low m=thin peak  // t= data resolution
    #tl<-transformList(fcs@parameters@data$name[z],lgcl)
    #ff.logicle.gated <-transform(ff.logicle.gated,tl)
  #}
  
  #par(mfrow= c(3,3))
  #for(i in target.column) {
   # x1 <- ff.logicle.gated@exprs[,i]
    #d1 <- density(x1)
    #plot(d1, xlab= ff.logicle.gated@parameters@data$desc[i], 
     #    xlim=c(0, ff.logicle.gated@parameters@data$maxRange[i]),
      #   main="Check logicle display", cex.main=0.7)
  #}
  #check comp and logicle
  #plotDens(ff.logicle.gated, c(10,11))
  
  
  #############################################################################################
  
  #gating
  
  #rename the parameter to fit with the fcs3.0 original parameter name
  for (index in 1:length(region.name.R.P)) {
    
    indices.gate.matrix <- match(colnames(region.matrix.R.P[[index]]),
                                 fcs.parameters.data.new.name)
    
    new.name.gate.matrix <- old.fcs.parameter.name[indices.gate.matrix]
    colnames(region.matrix.R.P[[index]]) <- (new.name.gate.matrix)
  }
  
  
  #based on the table region.table apply the gating stategy
  
  #If no gating fcs == fcs
  if (sum(as.numeric(region.table[,3])) !=0){
    indices.region <- which(region.table[,3]=="1")
    sFilename.gating.strategy <- paste(region.name.R.P[indices.region], collapse = "_")
    
    for (loop.region in indices.region){
      gate.matrix <- as.matrix(region.matrix.R.P[[loop.region]])
      
      if(region.table[loop.region,1] == "R"){
        #transform rectangular gate.matrix into a polygon 
        #x.matrix <- c(rep(gate.matrix[1,1],2), rep(gate.matrix[2,1],2))
        #y.matrix <- c(gate.matrix[1:2,2], rev(gate.matrix[1:2,2]))
        #gate.matrix.R <- cbind(x.matrix, y.matrix)
        #colnames(gate.matrix.R) <- colnames(gate.matrix)
        #R.gate <- rectangleGate(filterId = " ", (gate.matrix))
        #R.gate <- polygonGate(filterId = " ", gate.matrix.R)
        #fcs <- Subset(fcs, R.gate)
        
        x.R.indice <- match(colnames(gate.matrix)[1],old.fcs.parameter.name)
        y.R.indice <- match(colnames(gate.matrix)[2],old.fcs.parameter.name)
        
        x.R <- (fcs.comp@exprs[,x.R.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.R.indice] < gate.matrix[2,1])
        y.R <- (fcs.comp@exprs[,y.R.indice] > gate.matrix[1,2])  & (fcs.comp@exprs[,y.R.indice] < gate.matrix[2,2])
        
        x.y.R <- x.R & y.R
        
        fcs <- fcs[x.y.R]
        fcs.comp <- fcs.comp[x.y.R]
        
        }else{
          if(region.table[loop.region,1] == "P"){  
          #P.gate <- polygonGate(filterId = " ", gate.matrix)
          #fcs <- Subset(fcs.comp, P.gate)
            x.P.indice <- match(colnames(gate.matrix)[1],old.fcs.parameter.name)
            y.P.indice <- match(colnames(gate.matrix)[2],old.fcs.parameter.name)
            
            P.gate <- sp::point.in.polygon (point.x = fcs.comp@exprs[,x.P.indice], 
                                            point.y = fcs.comp@exprs[,y.P.indice],
                                            pol.x = gate.matrix[,1], 
                                            pol.y = gate.matrix[,2], 
                                            mode.checked = FALSE)
            
            fcs <- fcs [P.gate==1]
            fcs.comp <- fcs.comp [P.gate==1]
          }else{
            if(region.table[loop.region,1] == "L"){  
              x.L.indice <- match(colnames(gate.matrix)[1],old.fcs.parameter.name) 
              x.L <- (fcs.comp@exprs[,x.L.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.L.indice] < gate.matrix[2,1])
              fcs <- fcs[x.L]
              fcs.comp <- fcs.comp[x.L]
            }
        }
    }
   }
  }    
  
  if (sum(as.numeric(region.table[,3])) ==0) {sFilename.gating.strategy <- "ungated"}
  
  #############################################################################################
  # Set the compensation matrix to a R/Diva/fcs3.1 format #
  #############################################################################################
  fcs@description$'$SPILLOVER'<- matrix.coef
  #############################################################################################
  
  nb.FLA <- nrow(matrix.coef)
  #name.FLA <- paste("FL", 1:10, "-A",sep="")
  
  tempo.Navios.comp <- grep("^FL[0-9]",fcs.description.name)
  if(length (tempo.Navios.comp) != 0) {
  name.FLA.1 <- grep("^FL[0-9]",fcs.description.name)
  name.FLA <- fcs.description.name[tempo.Navios.comp]
  } else {
    name.FLA <- comp.matrix.parameter.name
  }
    
  FLA.descr <- name.FLA[1]
  for (z in 2:length(name.FLA)) {
    FLA.descr <- paste(FLA.descr,name.FLA[z],sep = ",")
  }
  

  matrix.coef2.vector<-c(t(matrix.coef))
  matrix.coef<-matrix.coef2.vector[1]
  for (z in 2:length(matrix.coef2.vector)) {
    matrix.coef <- paste(matrix.coef,matrix.coef2.vector[z],sep = ",")
  }
  
  
  #############################################################################################
  # START building the fcs file with new parameters #
  #############################################################################################
  
  #calculate number of events and number of parameters
  
  cytometer.name <- "Navios BeD"
  adc.resolution.Area <- 20
  adc.resolution.Height <- 18
  adc.resolution.TIME <- 20
  adc.resolution.Width <- 10
  
  
  #csv.File2be.Converted <- fcs@exprs
  csv.File2be.Converted <- fcs.comp@exprs
  
  #FC500 fc3.0 tous les paramètres sont cochés alors que dans 
  #le fc2 seuls les paramètres d'intérêt sont présents
  #On retire donc les paramètres du fc3.0 qui ne sont pas présents dans le slot fc2
  #On renome les noms et descriptions des paramètres. 
  
  target.parameter.name.not.NA <- !is.na(fcs.description.name)
  csv.File2be.Converted <- csv.File2be.Converted [,target.parameter.name.not.NA]
  
  fcs.name <- fcs.name[target.parameter.name.not.NA]
  fcs.desc <- fcs.desc[target.parameter.name.not.NA]
  
  fcs.description.name <- fcs.description.name[target.parameter.name.not.NA]
  fcs.description.desc <- fcs.description.desc[target.parameter.name.not.NA]
  
  
  
  cEvents <- nrow(csv.File2be.Converted)
  cParams <- ncol(csv.File2be.Converted)
  
  #file.name <- zip.path[[a]]
  #file.name.1 <- unlist(strsplit(file.name, "/"))
  #file.name.2 <- file.name.1[length(file.name.1)]
  #file.name.3 <- substr(file.name.2, start = 1, stop = nchar(file.name.2)-4)
  
  #sFilename <- paste(output.folder,date.format,file.name.3,"-new.fcs", sep="")
  
  
  
  sFilename <- paste(output.folder, date.format, "_",
                     fcs@description$`@SAMPLEID1`,"_",
                     fcs@description$`@SAMPLEID2`,"_",
                     fcs@description$`@SAMPLEID3`,"_",
                     fcs@description$`@SAMPLEID4`,
                     "_gating_", sFilename.gating.strategy,
                     "-new.fcs", sep="")
  
  
  #sFilename.2 <- paste(date.format,file.name.3,"-new.fcs", sep="")  
  sFilename.2 <- paste(date.format, "_",
                     fcs@description$`@SAMPLEID1`,"_",
                     fcs@description$`@SAMPLEID2`,"_",
                     fcs@description$`@SAMPLEID3`,"_",
                     fcs@description$`@SAMPLEID4`,
                     "_gating_", sFilename.gating.strategy,
                     "-new.fcs", sep="")
  #END Cytomet is Gallios or Navios or FC500
  } else {
    
  ################### Gallios Kaluza 
    #START for Gallios Kaluza 
    if(any(tempo.GalliosK) == TRUE ) {
      
      fcs <- fcs.test
      
      fcs.name <- unname(fcs.test@parameters@data$name)
      fcs.desc <- unname(fcs.test@parameters@data$desc)
      
      fcs.description.name <- unname(fcs.test@parameters@data$name)
      fcs.description.desc <- unname(fcs.test@parameters@data$desc)
      
      #############################################################################################
      # Set the compensation matrix to a R/Diva/fcs3.1 format #
      #############################################################################################
      
      #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      #extraction of the right parameter columns of the comp matrix
      parameters.comp.matrix.1 <- grep("^Time",fcs.description.name)
      parameters.comp.matrix.2 <- grep("^FSC",fcs.description.name)
      parameters.comp.matrix.3 <- grep("^SSC",fcs.description.name)
      parameters.comp.matrix.4 <- grep("^FS",fcs.description.name)
      parameters.comp.matrix.5 <- grep("^SS",fcs.description.name)
      
      parameters.comp.matrix.1 <- c(parameters.comp.matrix.1, 
                                    parameters.comp.matrix.2,
                                    parameters.comp.matrix.3,
                                    parameters.comp.matrix.4,
                                    parameters.comp.matrix.5)
      
      #name.FLA <- fcs.description.name[-parameters.comp.matrix.1]
      
      #matrix R + fcs file format
      #tempo.BD.comp <- grep("^[A-Z]",fcs.description.name)
      #if(length (tempo.GalliosK.comp) != 0) {
      #  name.FLA.1 <- grep("^FL[0-9]",fcs.description.name)
      #  name.FLA <- fcs.description.name[tempo.GalliosK.comp]
      #} else {
      #  name.FLA <- comp.matrix.parameter.name
      #}
      #nb.FLA <- nrow(matrix.coef)
      matrix.coef <- comp.matrix.fcs <- t(comp.matrix)
      
      #Compensate the fcs file with the new matrix
      #comp.matrix.fcs.name <- paste(comp.matrix.parameter.name, "-A",sep="")
      #colnames(comp.matrix.fcs) <- comp.matrix.fcs.name
      
      #colnames(comp.matrix.fcs) <- name.FLA
      #fcs.comp <- compensate(fcs.test, comp.matrix.fcs)
      
      #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      
      nb.FLA <- nrow(comp.matrix)
      #name.FLA <- paste("FL", 1:10, "-A",sep="")
      
      nb.FLA <- nrow(matrix.coef)
      #name.FLA <- paste("FL", 1:10, "-A",sep="")
      
      tempo.GalliosK.comp <- grep("^FL[0-9]",fcs.description.name)
      if(length (tempo.GalliosK.comp) != 0) {
        name.FLA.1 <- grep("^FL[0-9]",fcs.description.name)
        name.FLA <- fcs.description.name[tempo.GalliosK.comp]
      } else {
        name.FLA <- comp.matrix.parameter.name
      }
      
      colnames(comp.matrix.fcs) <- name.FLA
      fcs.comp <- compensate(fcs.test, comp.matrix.fcs)
      
      FLA.descr <- name.FLA[1]
      for (z in 2:length(name.FLA)) {
        FLA.descr <- paste(FLA.descr,name.FLA[z],sep = ",")
      }
      
      matrix.coef2.vector <- c((comp.matrix))
      matrix.coef<-matrix.coef2.vector[1]
      for (z in 2:length(matrix.coef2.vector)) {
        matrix.coef <- paste(matrix.coef,matrix.coef2.vector[z],sep = ",")
      }
      
      #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      #based on the table region.table apply the gating stategy
      
      #If no gating fcs == fcs
      if (sum(as.numeric(region.table[,3])) !=0){
        indices.region <- which(region.table[,3]=="1")
        sFilename.gating.strategy <- paste(region.name.R.P[indices.region], collapse = "_")
        
        for (loop.region in indices.region){
          gate.matrix <- as.matrix(region.matrix.R.P[[loop.region]])
          
          if(region.table[loop.region,1] == "R"){
            #transform rectangular gate.matrix into a polygon 
            #x.matrix <- c(rep(gate.matrix[1,1],2), rep(gate.matrix[2,1],2))
            #y.matrix <- c(gate.matrix[1:2,2], rev(gate.matrix[1:2,2]))
            #gate.matrix.R <- cbind(x.matrix, y.matrix)
            #colnames(gate.matrix.R) <- colnames(gate.matrix)
            #R.gate <- rectangleGate(filterId = " ", (gate.matrix))
            #R.gate <- polygonGate(filterId = " ", gate.matrix.R)
            #fcs <- Subset(fcs, R.gate)
            
            x.R.indice <- match(colnames(gate.matrix)[1],fcs.description.name)
            y.R.indice <- match(colnames(gate.matrix)[2],fcs.description.name)
            
            x.R <- (fcs.comp@exprs[,x.R.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.R.indice] < gate.matrix[2,1])
            y.R <- (fcs.comp@exprs[,y.R.indice] > gate.matrix[1,2])  & (fcs.comp@exprs[,y.R.indice] < gate.matrix[2,2])
            
            x.y.R <- x.R & y.R
            
            fcs <- fcs[x.y.R]
            fcs.comp <- fcs.comp[x.y.R]
            
          }else{
            if(region.table[loop.region,1] == "P"){  
              #P.gate <- polygonGate(filterId = " ", gate.matrix)
              #fcs <- Subset(fcs.comp, P.gate)
              x.P.indice <- match(colnames(gate.matrix)[1],fcs.description.name)
              y.P.indice <- match(colnames(gate.matrix)[2],fcs.description.name)
              
              P.gate <- sp::point.in.polygon (point.x = fcs.comp@exprs[,x.P.indice], 
                                              point.y = fcs.comp@exprs[,y.P.indice],
                                              pol.x = gate.matrix[,1], 
                                              pol.y = gate.matrix[,2], 
                                              mode.checked = FALSE)
              
              fcs <- fcs [P.gate==1]
              fcs.comp <- fcs.comp [P.gate==1]
            }else{
              if(region.table[loop.region,1] == "L"){  
                x.L.indice <- match(colnames(gate.matrix)[1],fcs.description.name) 
                x.L <- (fcs.comp@exprs[,x.L.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.L.indice] < gate.matrix[2,1])
                fcs <- fcs[x.L]
                fcs.comp <- fcs.comp[x.L]
              }
            }
          }
        }
      }    
      
      if (sum(as.numeric(region.table[,3])) ==0) {sFilename.gating.strategy <- "ungated"}
      
      #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
      
    cytometer.name <- "Gallios (Kaluza) BeD"
    adc.resolution.Area <- 20
    adc.resolution.Height <- 18
    adc.resolution.TIME <- 20
    adc.resolution.Width <- 10
    
    sFilename <- paste(output.folder, date.format, #"_",
                       fcs.test@description$`@SAMPLEID1`,"_",
                       fcs.test@description$`@SAMPLEID2`,"_",
                       fcs.test@description$`@SAMPLEID3`,"_",
                       fcs.test@description$`@SAMPLEID4`,
                       "-new.fcs", sep="")
    
    
    #sFilename.2 <- paste(date.format,file.name.3,"-new.fcs", sep="")  
    sFilename.2 <- paste(date.format, #"_",
                         fcs.test@description$`@SAMPLEID1`,"_",
                         fcs.test@description$`@SAMPLEID2`,"_",
                         fcs.test@description$`@SAMPLEID3`,"_",
                         fcs.test@description$`@SAMPLEID4`,
                         "-new.fcs", sep="")
    
    csv.File2be.Converted <- fcs@exprs
    
    #csv.File2be.Converted <- fcs@exprs #
    csv.File2be.Converted <- fcs.comp@exprs #
    
    
    
    fcs.name <- fcs.name
    fcs.desc <- fcs.desc
    
    fcs.description.name <- fcs.description.name
    fcs.description.desc <- fcs.description.desc
    
    
    
    }else{
    #BD FACSDiva
    if(length(tempo.BD) != 0) {
      #############################################################################################
      fcs <- fcs.test
      
      fcs.name <- unname(fcs.test@parameters@data$name)
      fcs.desc <- unname(fcs.test@parameters@data$desc)
      
      fcs.description.name <- unname(fcs.test@parameters@data$name)
      fcs.description.desc <- unname(fcs.test@parameters@data$desc)
      
      #extraction of the right parameter columns of the comp matrix
      parameters.comp.matrix.1 <- grep("^Time",fcs.description.name)
      parameters.comp.matrix.2 <- grep("^FSC",fcs.description.name)
      parameters.comp.matrix.3 <- grep("^SSC",fcs.description.name)
      
      
      parameters.comp.matrix.1 <- c(parameters.comp.matrix.1, 
                                  parameters.comp.matrix.2,
                                  parameters.comp.matrix.3)
      
      
      
      name.FLA <- fcs.description.name[-parameters.comp.matrix.1]
      
      #matrix R + fcs file format
      #tempo.BD.comp <- grep("^[A-Z]",fcs.description.name)
      #if(length (tempo.GalliosK.comp) != 0) {
      #  name.FLA.1 <- grep("^FL[0-9]",fcs.description.name)
      #  name.FLA <- fcs.description.name[tempo.GalliosK.comp]
      #} else {
      #  name.FLA <- comp.matrix.parameter.name
      #}
      
      
      nb.FLA <- nrow(matrix.coef)
      matrix.coef <- comp.matrix.fcs <- t(comp.matrix)
      
      #Compensate the fcs file with the new matrix
      #comp.matrix.fcs.name <- paste(comp.matrix.parameter.name, "-A",sep="")
      #colnames(comp.matrix.fcs) <- comp.matrix.fcs.name
      colnames(comp.matrix.fcs) <- name.FLA
      fcs.comp <- compensate(fcs.test, comp.matrix.fcs)
      
      # Set the compensation matrix to a R/Diva/fcs3.1 format #
            #name.FLA <- paste("FL", 1:10, "-A",sep="")
            #name.FLA.1 <- fcs.description.name[-parameters.comp.matrix.1]
      #name.FLA <- name.FLA.1
      #name.FLA <- comp.matrix.parameter.name
      
      FLA.descr <- name.FLA[1]
      for (z in 2:length(name.FLA)) {
        FLA.descr <- paste(FLA.descr,name.FLA[z],sep = ",")
      }
      
      
      matrix.coef2.vector <- c(t(matrix.coef))
      matrix.coef <- matrix.coef2.vector[1]
      for (z in 2:length(matrix.coef2.vector)) {
        matrix.coef <- paste(matrix.coef,matrix.coef2.vector[z],sep = ",")
      }
      
      
      
      
      
      #based on the table region.table apply the gating stategy
      
      #If no gating fcs == fcs
      if (sum(as.numeric(region.table[,3])) !=0){
        indices.region <- which(region.table[,3]=="1")
        sFilename.gating.strategy <- paste(region.name.R.P[indices.region], collapse = "_")
        
        for (loop.region in indices.region){
          gate.matrix <- as.matrix(region.matrix.R.P[[loop.region]])
          
          if(region.table[loop.region,1] == "R"){
            #transform rectangular gate.matrix into a polygon 
            #x.matrix <- c(rep(gate.matrix[1,1],2), rep(gate.matrix[2,1],2))
            #y.matrix <- c(gate.matrix[1:2,2], rev(gate.matrix[1:2,2]))
            #gate.matrix.R <- cbind(x.matrix, y.matrix)
            #colnames(gate.matrix.R) <- colnames(gate.matrix)
            #R.gate <- rectangleGate(filterId = " ", (gate.matrix))
            #R.gate <- polygonGate(filterId = " ", gate.matrix.R)
            #fcs <- Subset(fcs, R.gate)
            
            x.R.indice <- match(colnames(gate.matrix)[1],fcs.description.name)
            y.R.indice <- match(colnames(gate.matrix)[2],fcs.description.name)
            
            x.R <- (fcs.comp@exprs[,x.R.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.R.indice] < gate.matrix[2,1])
            y.R <- (fcs.comp@exprs[,y.R.indice] > gate.matrix[1,2])  & (fcs.comp@exprs[,y.R.indice] < gate.matrix[2,2])
            
            x.y.R <- x.R & y.R
            
            fcs <- fcs[x.y.R]
            fcs.comp <- fcs.comp[x.y.R]
            
          }else{
            if(region.table[loop.region,1] == "P"){  
              #P.gate <- polygonGate(filterId = " ", gate.matrix)
              #fcs <- Subset(fcs.comp, P.gate)
              x.P.indice <- match(colnames(gate.matrix)[1],fcs.description.name)
              y.P.indice <- match(colnames(gate.matrix)[2],fcs.description.name)
              
              P.gate <- sp::point.in.polygon (point.x = fcs.comp@exprs[,x.P.indice], 
                                              point.y = fcs.comp@exprs[,y.P.indice],
                                              pol.x = gate.matrix[,1], 
                                              pol.y = gate.matrix[,2], 
                                              mode.checked = FALSE)
              
              fcs <- fcs [P.gate==1]
              fcs.comp <- fcs.comp [P.gate==1]
            }else{
              if(region.table[loop.region,1] == "L"){  
                x.L.indice <- match(colnames(gate.matrix)[1],fcs.description.name) 
                x.L <- (fcs.comp@exprs[,x.L.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.L.indice] < gate.matrix[2,1])
                fcs <- fcs[x.L]
                fcs.comp <- fcs.comp[x.L]
              }
            }
          }
        }
      }    
      
      if (sum(as.numeric(region.table[,3])) ==0) {sFilename.gating.strategy <- "ungated"}

      
      cytometer.name <- "BD BeD"
      adc.resolution.Area <- 18
      adc.resolution.Height <- 18
      adc.resolution.TIME <- 18
      adc.resolution.Width <- 18
      
      
      
      
      #sFilename <- paste(output.folder, date.format, "_",
       #                  substr(fcs.test@description$`$FIL`, start = 1, 
        #                        stop = nchar(fcs.test@description$`$FIL`)-4),
         #                "-new.fcs", sep="")
      
      sFilename <- paste(output.folder, 
                         #date.format, "_",
                         zip.path.name.1 [a],"_", aa,
                         "-new.fcs", sep="")
      
      
      #sFilename.2 <- paste(date.format,file.name.3,"-new.fcs", sep="")  
      #sFilename.2 <- paste(date.format, "_",
       #                    substr(fcs.test@description$`$FIL`, start = 1, 
        #                          stop = nchar(fcs.test@description$`$FIL`)-4),
         #                  "-new.fcs", sep="")
      
      sFilename.2 <- paste(date.format, "_",
                         zip.path.name.1 [a],"_", aa,
                         "-new.fcs", sep="")
      
      
      #csv.File2be.Converted <- fcs@exprs #
      csv.File2be.Converted <- fcs.comp@exprs #
      
      
      
      fcs.name <- fcs.name
      fcs.desc <- fcs.desc
      
      fcs.description.name <- fcs.description.name
      fcs.description.desc <- fcs.description.desc
      
    }else{
    ##############CytoFLEX####################
    if(any(tempo.CytoFLEX) == TRUE) {
      
      fcs <- fcs.test
      
      fcs.name <- unname(fcs.test@parameters@data$name)
      fcs.desc <- unname(fcs.test@parameters@data$desc)
      
      fcs.description.name <- unname(fcs.test@parameters@data$name)
      fcs.description.desc <- unname(fcs.test@parameters@data$desc)
      
      #########################################################################################
      #extraction of the right parameter columns of the comp matrix
      #name.FLA <- comp.matrix.parameter.name
      FL.indices <- grep("^FL[0-9]",fcs.description.name)
      FL.H.indices <- grep("-H",fcs.description.name[FL.indices])
      FL.A.indices <- grep("-A",fcs.description.name[FL.indices])
      
      if((any(FL.H.indices)==TRUE) & (any(FL.A.indices)==TRUE)){  #if there are -H parameter AND -A parameterthen big matrix
      
      name.FLA.H <- paste0(comp.matrix.parameter.name,"-H")
      name.FLA.A <- paste0(comp.matrix.parameter.name,"-A")
      name.FL.H.A <- c(name.FLA.H, name.FLA.A)

      FLA.descr <- paste(name.FL.H.A, collapse=",")
      name.FLA <- FLA.descr
      
      matrix.coef.1 <- comp.matrix.fcs <- (comp.matrix) #####ne pas prendre la transposée
      
      nb.column <- ncol(comp.matrix)
      nb.FLA <- NCOL(comp.matrix) *2

      comp.matrix.1 <- t(comp.matrix) # !!!!!!!!!
      comp.matrix.2 <- matrix(data= rep(0, nb.column^2), ncol = nb.column)
      
      comp.matrix.3 <- cbind(comp.matrix.1, comp.matrix.2) #matrix has to be square
      comp.matrix.4 <- cbind(comp.matrix.2, comp.matrix.1)
      
      comp.matrix.5 <- rbind(comp.matrix.3, comp.matrix.4)
      colnames(comp.matrix.5) <- name.FL.H.A
      fcs.comp <- compensate(fcs.test,(comp.matrix.5))
      
      }else{
        if((any(FL.H.indices)==TRUE) & (any(FL.A.indices)==FALSE)){
          name.FLA.H <- paste0(comp.matrix.parameter.name,"-H")
          #name.FLA.A <- paste0(comp.matrix.parameter.name,"-A")
          name.FL.H.A <- c(name.FLA.H)
          
          FLA.descr <- paste(name.FL.H.A, collapse=",")
          name.FLA <- FLA.descr
          
          matrix.coef.1 <- comp.matrix.fcs <- t(comp.matrix) #####ne pas prendre la transposée
          
          nb.column <- ncol(comp.matrix)
          nb.FLA <- NCOL(comp.matrix)
          
          #comp.matrix.1 <- comp.matrix
          #comp.matrix.2 <- matrix(data= rep(0, nb.column^2), ncol = nb.column)
          
          #comp.matrix.3 <- cbind(comp.matrix.1, comp.matrix.2) #matrix has to be square
          #comp.matrix.4 <- cbind(comp.matrix.2, comp.matrix.1)
          
          comp.matrix.5 <- t(comp.matrix) #!!!!!!!!!!!!!!!!!
          colnames(comp.matrix.5) <- name.FL.H.A
          fcs.comp <- compensate(fcs.test,(comp.matrix.5))
          
        }else{
          if((any(FL.H.indices)==FALSE) & (any(FL.A.indices)==TRUE)){
            #name.FLA.H <- paste0(comp.matrix.parameter.name,"-H")
            name.FLA.A <- paste0(comp.matrix.parameter.name,"-A")
            name.FL.H.A <- c(name.FLA.A)
            
            FLA.descr <- paste(name.FL.H.A, collapse=",")
            name.FLA <- FLA.descr
            
            matrix.coef.1 <- comp.matrix.fcs <- t(comp.matrix) #####ne pas prendre la transposée
            
            nb.column <- ncol(comp.matrix)
            nb.FLA <- NCOL(comp.matrix)
            
            #comp.matrix.1 <- comp.matrix
            #comp.matrix.2 <- matrix(data= rep(0, nb.column^2), ncol = nb.column)
            
            #comp.matrix.3 <- cbind(comp.matrix.1, comp.matrix.2) #matrix has to be square
            #comp.matrix.4 <- cbind(comp.matrix.2, comp.matrix.1)
            
            comp.matrix.5 <- t(comp.matrix)
            colnames(comp.matrix.5) <- name.FL.H.A
            fcs.comp <- compensate(fcs.test,(comp.matrix.5))
          }
        }
      }
           
      
      
      
      
      
      # Set the compensation matrix to a R/Diva/fcs3.1 format #
      matrix.coef2.vector <- c((comp.matrix.5)) #####ne pas prendre la transposée
      matrix.coef <- paste(matrix.coef2.vector, collapse = ",")

      #based on the table region.table apply the gating stategy
      
      #If no gating fcs == fcs
      if (sum(as.numeric(region.table[,3])) !=0){
        indices.region <- which(region.table[,3]=="1")
        sFilename.gating.strategy <- paste(region.name.R.P[indices.region], collapse = "_")
        
        for (loop.region in indices.region){
          gate.matrix <- as.matrix(region.matrix.R.P[[loop.region]])
          
          if(region.table[loop.region,1] == "R"){
            #transform rectangular gate.matrix into a polygon 
            #x.matrix <- c(rep(gate.matrix[1,1],2), rep(gate.matrix[2,1],2))
            #y.matrix <- c(gate.matrix[1:2,2], rev(gate.matrix[1:2,2]))
            #gate.matrix.R <- cbind(x.matrix, y.matrix)
            #colnames(gate.matrix.R) <- colnames(gate.matrix)
            #R.gate <- rectangleGate(filterId = " ", (gate.matrix))
            #R.gate <- polygonGate(filterId = " ", gate.matrix.R)
            #fcs <- Subset(fcs, R.gate)
            
            x.R.indice <- match(colnames(gate.matrix)[1],fcs.description.name)
            y.R.indice <- match(colnames(gate.matrix)[2],fcs.description.name)
            
            x.R <- (fcs.comp@exprs[,x.R.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.R.indice] < gate.matrix[2,1])
            y.R <- (fcs.comp@exprs[,y.R.indice] > gate.matrix[1,2])  & (fcs.comp@exprs[,y.R.indice] < gate.matrix[2,2])
            
            x.y.R <- x.R & y.R
            
            fcs <- fcs[x.y.R]
            fcs.comp <- fcs.comp[x.y.R]
            
          }else{
            if(region.table[loop.region,1] == "P"){  
              #P.gate <- polygonGate(filterId = " ", gate.matrix)
              #fcs <- Subset(fcs.comp, P.gate)
              x.P.indice <- match(colnames(gate.matrix)[1],fcs.description.name)
              y.P.indice <- match(colnames(gate.matrix)[2],fcs.description.name)
              
              P.gate <- sp::point.in.polygon (point.x = fcs.comp@exprs[,x.P.indice], 
                                              point.y = fcs.comp@exprs[,y.P.indice],
                                              pol.x = gate.matrix[,1], 
                                              pol.y = gate.matrix[,2], 
                                              mode.checked = FALSE)
              
              fcs <- fcs [P.gate==1]
              fcs.comp <- fcs.comp [P.gate==1]
            }else{
              if(region.table[loop.region,1] == "L"){  
                x.L.indice <- match(colnames(gate.matrix)[1],fcs.description.name) 
                x.L <- (fcs.comp@exprs[,x.L.indice] > gate.matrix[1,1])  & (fcs.comp@exprs[,x.L.indice] < gate.matrix[2,1])
                fcs <- fcs[x.L]
                fcs.comp <- fcs.comp[x.L]
              }
            }
          }
        }
      }    
      
      if (sum(as.numeric(region.table[,3])) ==0) {sFilename.gating.strategy <- "ungated"}
      
      
      
      
            
      cytometer.name <- "CytoFLEX BeD"
      adc.resolution.Area <- 24
      adc.resolution.Height <- 24
      adc.resolution.TIME <- 32
      adc.resolution.Width <- 13
      
      sFilename <- paste(output.folder, date.format, "_",
                         substr(fcs.test@description$`$FIL`, start = 1, 
                                stop = nchar(fcs.test@description$`$FIL`)-4),
                         "-new.fcs", sep="")
      
      
      #sFilename.2 <- paste(date.format,file.name.3,"-new.fcs", sep="")  
      sFilename.2 <- paste(date.format, "_",
                           substr(fcs.test@description$`$FIL`, start = 1, 
                                  stop = nchar(fcs.test@description$`$FIL`)-4),
                           "-new.fcs", sep="")
      
      
      #test
      #comp.matrix.0 <- fcs.test@description$`$SPILLOVER`
      #ncol(comp.matrix)
      #fcs.test.comp <- compensate(fcs.test,comp.matrix.5)
      
      #csv.File2be.Converted <- fcs@exprs #
      csv.File2be.Converted <- fcs.comp@exprs #
      
      
    }
####################### END CytoFLEX ################    
    else{
    print("This instrument is not supported yet - contact your Beckman Coulter representative")
            break
       }
      }
     }
    }
    
#### Build the fcs file    
  cEvents <- nrow(csv.File2be.Converted)
  cParams <- ncol(csv.File2be.Converted)
    
  
  #create the HEADER and TEXT sections of the FCS file as strings
  #Delimiter choice
  I<-rawToChar(as.raw(124)) #124="|" // 12="/f" BD Diva // #92="\" R
  sText.1<-I
  sText.1<- paste(sText.1,"$BEGINANALYSIS",I,0,I,sep="")
  sText.1<- paste(sText.1,"$ENDANALYSIS",I,0,I,sep="")
  sText.1<- paste(sText.1,"$BEGINSTEXT",I,0,I,sep="")
  sText.1<- paste(sText.1,"$ENDSTEXT",I,0,I,sep="")
  #$BEGINDATA and $ENDATA hereafter
  sText.3<- paste("$FIL",I,sFilename.2,I,sep="")
  #sText.3<- paste(sText.3,"$SYS",I,"Windows XP 5.1",I,sep="")
  cEvents.diva<-sprintf("%-19.f",cEvents)
  sText.3<- paste(sText.3,"$TOT",I,cEvents.diva,I,sep="")
  sText.3<- paste(sText.3,"$PAR",I,cParams,I,sep="")
  sText.3<- paste(sText.3,"$MODE",I,"L",I,sep="")
  sText.3<- paste(sText.3,"$BYTEORD",I,"4,3,2,1",I,sep="")
  sText.3<- paste(sText.3,"$DATATYPE",I,"F",I,sep="")
  sText.3<- paste(sText.3,"$NEXTDATA",I,"         0",I,sep="")
  #sText.3<- paste(sText.3,"CREATOR",I,"BD FACSDiva Software Version 6.1.2",I,sep="")
  #sText.3<- paste(sText.3,"TUBE NAME",I,fcs.name,I,sep="")
  #sText.3<- paste(sText.3,"$SRC",I,"Compensation Controls",I,sep="")
  #sText.3<- paste(sText.3,"EXPERIMENT NAME",I,"20151201",I,sep="")
  #sText.3<- paste(sText.3,"GUID",I,"55193328-3c99-46c6-97ca-1413fe274753",I,sep="")
  #sText.3<- paste(sText.3,"$DATE",I,fs@description$`$DATE`,I,sep="")
  #sText.3<- paste(sText.3,"$BTIM",I,"16:30:56",I,sep="")
  #sText.3<- paste(sText.3,"$ETIM",I,"16:35:56",I,sep="")
  #sText.3<- paste(sText.3,"WINDOW EXTENSION",I,"0.00",I,sep="")
  #sText.3<- paste(sText.3,"EXPORT USER NAME",I,"BENOIT",I,sep="")
  #sText.3<- paste(sText.3,"EXPORT TIME",I,"01-DEC-2015-16:30:56",I,sep="")
  #sText.3<- paste(sText.3,"FSC ASF",I,"0.00",I,sep="")
  #sText.3<- paste(sText.3,"AUTOBS",I,"TRUE",I,sep="")
  sText.3<- paste(sText.3,"$CYT",I, cytometer.name,I,sep="")
  sText.3<- paste(sText.3,"$TIMESTEP",I,"0.01",I,sep="")
  #sText.3<-paste(sText.3,"SPILL",I,nb.FLA,",",FLA.descr,",",matrix.coef,I,sep="")
  
  for (n in 1: cParams){
    sText.3<-paste(sText.3,"$P",n,"N",I,fcs.description.name[n],I,sep="")
    sText.3<-paste(sText.3,"$P",n,"S",I,fcs.description.desc[n],I,sep="")
    if ((fcs.description.name[n]=="FS H") | (fcs.description.name[n]=="FS-H") |
        (fcs.description.name[n]=="FSC-H") | (fcs.description.name[n]=="FSC H") |
        (fcs.description.name[n]=="FSC-H") | (fcs.description.name[n]=="FSC H") |
        (fcs.description.name[n]=="FS-PEAK") | (fcs.description.name[n]=="FS PEAK")){
      range<-as.character(2^adc.resolution.Height)
    } else {
      if (fcs.description.name[n]=="TIME") {
        range<-as.character(2^adc.resolution.TIME)
      } else {
        if ((fcs.description.name[n]=="FS W") | (fcs.description.name[n]=="FS-W") |
            (fcs.description.name[n]=="FSC-Width") | (fcs.description.name[n]=="FSC Width") |
            (fcs.description.name[n]=="FS TOF") | (fcs.description.name[n]=="FS-TOF")){
          range<-as.character(2^adc.resolution.Width)
        } else {
          if ((fcs.description.name[n]=="SS H") | (fcs.description.name[n]=="SS-H") |
              (fcs.description.name[n]=="SSC-H") | (fcs.description.name[n]=="SSC H") |
              (fcs.description.name[n]=="SSC-H") | (fcs.description.name[n]=="SSC H") |
              (fcs.description.name[n]=="SS-PEAK") | (fcs.description.name[n]=="SS PEAK")) {
            range<-as.character(2^adc.resolution.Height)
          } else {
            if ((fcs.description.name[n]=="SS W") | (fcs.description.name[n]=="SS-W") |
                (fcs.description.name[n]=="SSC-Width") | (fcs.description.name[n]=="SSC Width") |
                (fcs.description.name[n]=="SS TOF") | (fcs.description.name[n]=="SS-TOF")){
              range<-as.character(2^13)
            } else {
              range<-as.character(2^adc.resolution.Area)
            }
          }
        }
      }
    }
    sText.3<-paste(sText.3,"$P",n,"R",I,range,I,sep="")
    sText.3<-paste(sText.3,"$P",n,"B",I,"32",I,sep="")
    sText.3<-paste(sText.3,"$P",n,"E",I,"0,0",I,sep="")
    #sText.3<-paste(sText.3,"$P",n,"V",I,"357",I,sep="")
    #sText.3<-paste(sText.3,"$P",n,"G",I,"1.0",I,sep="") #$PxG=0.01 for time
    #sText.3<-paste(sText.3,"$P",n,"DISPLAY",I,"LOG",I,sep="") #no keyword $PxDISPLAY for time
    #sText.3<-paste(sText.3,"$P",n,"BS",I,"-1",I,sep="") #BS=0 for Time
    #sText.3<-paste(sText.3,"$P",n,"MS",I,"0",I,sep="") #MS=0 for Time
  }
  sText.1.3<-paste(sText.1,sText.3,sep = "")
  nTextLength<-nchar(sText.1.3) +
    nchar(paste("$BEGINDATA",I,"0000",I,"$ENDDATA",I,"0000000000000000000",I,sep = ""))
  nTextLength<-nTextLength+5 #+5 <-5 blanks between ENDTEXT and STARTDATA
  #$BEGINDATA 4 bytes
  sText.2<- paste("$BEGINDATA",I,sprintf("%4.f",(256 + nTextLength)),I,sep="")
  #$ENDDATA 19 bytes - factor 4 because 4 bytes=32bit used/data -
  sText.2<- paste(sText.2,"$ENDDATA",I,
                  sprintf("%-19.f",(256 + nTextLength -1 + (cEvents * cParams * 4))),I,sep="")
  #Add 198 blank at the end of the SHeader
  space <- " "
  sText.198 <- NULL
  for (ii in 1:198){
    sText.198<- paste(space,sText.198,sep="")
  }
  
  sText.194 <- NULL
  for (ii in 1:194){
    sText.194<- paste(space,sText.194,sep="")
  }
  
  
  
  
  sText.4blank <- NULL
  for (ii in 1:4){
    sText.4blank<- paste(space,sText.4blank,sep="")
  }
  sText.7blank <- NULL
  for (ii in 1:7){
    sText.7blank<- paste(space,sText.7blank,sep="")
  }
  
  if ((256 + nTextLength -1 + (cEvents * cParams * 4)) >= 99999999){
    sHeader<-paste("FCS3.0", #START TEXT=#256
                   sText.4blank,
                   sprintf("%08.f",256), # START TEXT
                   sprintf("%08.f",(256 + nTextLength - 6)), #END TEXT
                   sprintf("%08.f",(256 + nTextLength)), #START DATA
                   #sprintf("%08.f",(256 + nTextLength -1 + (cEvents * cParams * 4))), #END DATA
                   sprintf("%08.f",0), #END DATA
                   paste(sText.7blank,0,sep=""),#7 blank + 0
                   paste(sText.7blank,0,sep=""), #7 blank + 0
                   sText.198,sep="")#198 = 6+4+ blank
  }else{
    sHeader<-paste("FCS3.0", #START TEXT=#256
                   sText.4blank,
                   sprintf("%08.f",256), # START TEXT
                   sprintf("%08.f",(256 + nTextLength - 6)), #END TEXT
                   sprintf("%08.f",(256 + nTextLength)), #START DATA
                   sprintf("%08.f",(256 + nTextLength -1 + (cEvents * cParams * 4))), #END DATA
                   #sprintf("%08.f",0), #END DATA
                   paste(sText.7blank,0,sep=""),#7 blank + 0
                   paste(sText.7blank,0,sep=""), #7 blank + 0
                   sText.198,sep="")#198 = 6+4+ blank  
  }
  
  
  
  
  #Add 5 blanks at the end of the TEXT
  sText.5 <- NULL
  for (ii in 1:5){
    sText.5<- paste(space,sText.5,sep="")
  }#
  #sText.5 <- paste(sText.5,I,sep = "") #### new
  #BD Diva sText right order
  sText<-paste(sText.1,sText.2,sText.3,sText.5,sep ="")
  #convert the dataframe/table csv in a vector to be able to use the writeBin function
  csv.vector<-c(t(csv.File2be.Converted))
  csv.vector <- as.numeric (csv.vector)
  
  #Write out the FCS file
  con<-file(sFilename,open="wb") #open for reading in binary mode
  writeChar(sHeader,con,eos=NULL)
  writeChar(sText,con,eos=NULL)
  writeBin(csv.vector,con,size=4,endian="big") #size=4 number of byte =>32 bit // little=1234 big=4321
  writeChar("00000000",con,eos=NULL)
  close(con)
}
unlink(file.name.zip, recursive = TRUE) #delete the tempo folder containing all fcs files + xml
}

###############################################################################################



#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
#                                rename .zip  in .analysis
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
  
  fcs.path <- list.files(path= input.folder, pattern = "\\.zip$", full.names=TRUE)
  
  nb.files <- length(fcs.path)
  
  fcs.path.new <- list()
  for (i in 1:nb.files){
    old.file <- fcs.path[[i]]
    new.file <- gsub("zip","analysis", old.file)
    file.rename(from = old.file, to = new.file)
    
  }
  #@END rename analysis file in zip
  


