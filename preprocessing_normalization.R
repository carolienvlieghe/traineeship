############################
### Normalization script ###
############################

### Clear Rstudio windows ###
rm(list=ls()) # removes all object from Rstudio environment window
cat("\014") # clears Rstudio console window
if(!is.null(dev.list())) dev.off() # clears the Rstudio plot window


library(flowCore)
# assign in- & output
output.folder <- "D:/school/Stage officieel/norm_out/"
input.folder <- "D:/school/Stage officieel/csv_in/"
file.name.csv <- list.files(path = input.folder, pattern = "\\.csv$", full.names = TRUE)

# read csv file
linear <- as.matrix(read.csv(file.name.csv))
ff <- flowFrame(linear)
ff_t <- transform(ff,transformList(colnames(ff)[6:15],logicleTransform()))
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
inputmatrix <- exprs(ff_t[,1:15]) # converts flowframe into matrix
# FL's
ffnormalized <- normalize(inputmatrix[,6:15])
# scatter
scatternorm <- inputmatrix[,1:5]
# Combine scatter and FL
ffmatrix <- cbind(scatternorm, ffnormalized*1024)
# scale to [0,1024]
logicle <- ffmatrix

## define settings
target.column <- 6:15
TAG.channel <- 500
final.cell.number <- 10000000000
CD45.parameter.nb <- 15
SSINT.parameter.nb <- 4
# ADC resolution of the cytometer
adc.resolution.Area <- 20
adc.resolution.Height <- 18
adc.resolution.Width <- 10
adc.resolution.TIME <- 20
# the max scale for FSC and SSC of the cytometer
max.scale.FSC.SSC <- 2^20

tempo <- colnames(linear[,1:15])
colnames(logicle) <- tempo

###################################################################################
#Enter here the name of the exported RPlugIn Logicle Matrix i.eFlowSOM_lgcl
table.to.be.computed <- logicle
################################################################################### 

###################################################################################
#Enter here the name of the exported RPlugIn Linear Matrix i.e FlowSOM
table.to.be.used.as.fcs.file <- logicle
################################################################################### 

nb.parameter <- ncol(table.to.be.used.as.fcs.file)


library(flowDensity)

  
###############################################################
  
TAG.fcs <- rep(TAG.channel, nrow(logicle))

#############    check the biexpoentiel   #####################
#check logicle coef in 1D plot
windows()
nb.plot.per.row <- 2
nb.plot.per.column <- 2
par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) )
nb.plot.per.windows <- nb.plot.per.row * nb.plot.per.column
aa <- 0

for(i in target.column) {
  if(aa==nb.plot.per.windows){
    windows()
    par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) )
    aa <- 0
  }
  aa <- aa+1
    
  x1 <- table.to.be.computed[,i]
  d1 <- density(x1)
  plot(d1, xlab= colnames(table.to.be.computed)[i], xlim=c(0,1024),
       main="Check logicle transformation", cex.main=0.7)
}
dev.off()
  
#Selection of the negative mode of each negative lymphocyte population
#using the mouse
negative.population.mode <- c()

for(i in target.column) {
  windows()
  x1 <- logicle[,i]
  d1 <- density(x1)
  plot(d1, xlab= colnames(table.to.be.computed)[i], main="Click on negative mode peak", cex.main=0.7, xlim=c(0,1024))
  coord.1 <- locator(1, type = "l", col="red")
  negative.population.mode[i] <- coord.1$x
  abline(v = negative.population.mode[i],col="red",lwd=1,lty=5)
  dev.off()
}
negative.population.mode <- negative.population.mode[!is.na(negative.population.mode)]
cat(negative.population.mode)
reference.neg.pop.mode <- c(200, #1.300439, #FITC
                            200, #0.994293, # PE
                            200, #1.4777779, #ECD
                            200,#1.525965, #PC5.5
                            200,#1.587294, #PC7
                            200,#1.254920, #APC
                            200, #A700
                            200,#2.916238, #AA750
                            200,#1.123684, #PB
                            700,#2.944265) #KrO
                            200,#BV610
                            200,#BV660
                            200)#BV780
                            
                            
virtual.gain <- reference.neg.pop.mode - negative.population.mode

csv.logicle.corrected <- lapply(c(1:10), function(i)
  table.to.be.computed[,i+5] +
    (reference.neg.pop.mode[i] - negative.population.mode[i]))

csv.logicle.corrected <- do.call(cbind, csv.logicle.corrected)

if(nrow(csv.logicle.corrected) >final.cell.number){
  sample_indices <- sample(1:nrow(csv.logicle.corrected), final.cell.number)
  csv.logicle.corrected <- csv.logicle.corrected[sample_indices,]
  table.to.be.used.as.fcs.file <- table.to.be.used.as.fcs.file[sample_indices,]	
}

TIME <- 1:nrow(csv.logicle.corrected)
TIME <- TIME/10

csv.corrected <- cbind(csv.logicle.corrected,TIME)
csv.corrected <- cbind(logicle[,1:5], csv.corrected)
colnames(csv.corrected) <- c(colnames(logicle), "TIME")

#check csv corrected
for(i in target.column) {
  if(aa==nb.plot.per.windows){
    windows()
    par(mfrow= c(nb.plot.per.row, nb.plot.per.column), pty="s", mar= c(5,4,4,2)-c(1,1,1,1.5) )
    aa <- 0
  }
  aa <- aa+1
  
  x1 <- csv.corrected[,i]
  d1 <- density(x1)
  plot(d1, xlab= colnames(table.to.be.computed)[i], xlim= c(0,1024),
       main="Check corrected csv ", cex.main=0.7)
  abline(v = negative.population.mode[i],col="red",lwd=1,lty=5)
  abline(v = reference.neg.pop.mode[i],col="blue",lwd=1,lty=5)
}

#End check

parameter.name.linear <- colnames(logicle)
parameter.name.logicle <- paste0("<N>-", parameter.name.linear)
#parameter.name.logicle[c(FSCINT.parameter.nb ,SSINT.parameter.nb)] <- paste0(parameter.name.linear[c(FSCINT.parameter.nb ,SSINT.parameter.nb)], "-<N>")

csv.DF.numeric <- cbind(csv.corrected, TAG.fcs)
colnames(csv.DF.numeric) <- c(parameter.name.logicle, "TIME" ,"TAG.fcs")

print(colnames(csv.DF.numeric) )
print(nrow(csv.DF.numeric) )

ff <- flowFrame(csv.DF.numeric)

#Define the complete output fcs file name
date <- Sys.time()
date.format <- format(date, format= "%Y%m%d-%H%M%S-")

sFilename <- paste("D:/school/Stage officieel/", "norm","_normalized.fcs" ,sep="")
write.FCS(ff, filename = sFilename)
summary(ff)
