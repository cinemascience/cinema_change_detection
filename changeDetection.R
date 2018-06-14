#---
#A command line change detection code
#---

#--
#inputs are as follows:
# R < changeDetection.R workingDirectory
#--

# clear all previous environment variables and memory
rm(list=ls())
#gc()

#current directory should already be working directory
#print(getwd())

#install required libraries
if(!require(argparse)){install.packages("argparse")}
if(!require(biglm)){install.packages("biglm")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(jsonlite)){install.packages("jsonlite")}

# load required libraries
library("argparse") #to parse the commandline input
library("biglm") #for Myers et al. change detection code
library("ggplot2") #to draw the resulting plot
library("jsonlite") #to save the json file

# parse the arguments to understand next steps
parser <- ArgumentParser()

parser$add_argument("-f", "--csvFile", help='csv input file with data')
parser$add_argument("-p", "--path", help='path to .cdb')
parser$add_argument("-r", "--restrict", nargs='+', help='restrict rows based on list of column names and values')
parser$add_argument("-x", "--xVal", help='the x paramter column to iterate over')
parser$add_argument("-y", "--yVal", help='the y parameter column to find change over')
parser$add_argument("-cd", "--cdParams", type="double", nargs='+', help='parameters for Myers et al. change detection algorithm')
parser$add_argument("-o", "--output", help='output phrase to mark files by')

args <- parser$parse_args()

#print(args$csvFile)
#print(args$path)
#print(args$restrict)
#print(args$x)
#print(args$y)
#print(args$cdParams)


# read input csvFile
csvInputFile <- read.csv(file=args$csvFile)

rowValues <- vector(mode = "numeric", length=0)
# cull out the needed rows according to the restrictions given
for (i in seq(1,length(args$restrict),2)) #interate through the restrictions
{
    csvInputFile <- csvInputFile[csvInputFile[,args$restrict[i]]==as.integer(args$restrict[i+1]),]
    rowValues <- as.numeric(rownames(csvInputFile))
    #csvInputFile <- csvInputFile[csvInputFile[,"phi"]==-180,]
}

featureMetric <- csvInputFile[,args$y]
inputParam <- 1:dim(csvInputFile)[1]
inputValues <- csvInputFile[,args$x]

# --- change detection using Myers et al.

# source file for Myers et al. technique
source('cinema_change_detection/run_lm.R')
#source('run_lm.R')

#set parameters for Myers et al. techniques
B <- as.integer(args$cd[1])
alpha <- as.double(args$cd[2])
deltasq <- as.double(args$cd[3])

# run Myers et al. linear fits algorithm
z <- run_lm(inputParam, featureMetric, B, alpha, deltasq)

# Dataframe to store starting and ending points of the linear regression model
#for each time step.
seg <- data.frame()
i <- 1
seg[i, 1] <- i
seg[i, 2] <- z$Time[i]
seg[i, 3] <- z$beta0[i] + z$beta1[i] * 1
seg[i, 4] <- z$beta0[i] + z$beta1[i] * z$Time[i]

for(i in 2:z$Partition)
{
    seg[i, 1] <- z$Time[i - 1] + 1
    seg[i, 2] <- z$Time[i]
    seg[i, 3] <- z$beta0[i] + z$beta1[i] * 1
    seg[i, 4] <- z$beta0[i] + z$beta1[i] * (z$Time[i] - z$Time[i - 1])
}

colnames(seg) <- c("x", "xend", "y", "yend")

changePoints <- vector(mode = "numeric", length=0)
changePoints <- c(changePoints, seg[,1])
changePoints <- c(changePoints, length(inputParam))
changePointsGraph <- vector(mode = "numeric", length=0)

changePointsGraph <- inputValues[changePoints]

changePointsY <- vector(mode = "numeric", length=0)
for (i in 1:length(changePoints))
{
    changePointsY <- c(changePointsY, featureMetric[changePoints[i]])
}

changePoints <- changePoints - 1

df <- data.frame(v1=changePoints, v2=changePointsY)
colnames(df) <- c("changePoints", "changePointsY")

#changePointsToZero <- inputParam - 1
changePointsToZero <- inputValues
featureMetricToZero <- featureMetric

dfData <- data.frame(v1=changePointsToZero, v2=featureMetricToZero)
colnames(dfData) <- c("changePointsToZero", "featureMetricToZero")

seg[,1] <- inputValues[seg[,1]]
seg[,2] <- inputValues[seg[,2]]

#--- change detection calculated

#----saving results to csv, as image and as json---"

#if resulting directory doesn't exist, create it
if (!dir.exists(file.path(args$path, "CCD")))
    dir.create(file.path(args$path, "CCD"))

# Plotting the data, the partition lines, and the linear regression model lines in each partition.
changePlot <- ggplot(dfData, aes(x = changePointsToZero, y = featureMetricToZero)) + geom_point(size = 1.2, alpha=1) + geom_line(size=0.5, alpha=1) + labs(x = args$x, y = args$y) + geom_segment(data = seg, mapping = aes(x = x, y = y, xend = xend, yend = yend), size = 2, colour = rgb(199/255,27/255,0/255), alpha = 1) + geom_point(data=df, mapping = aes(x=changePointsGraph, y=changePointsY), color=rgb(57/255,131/255,235/255), alpha=1, size = 4) + theme(panel.background = element_rect(fill="grey99", color="black"), axis.text = element_text(size=18), axis.title = element_text(size=20), panel.grid.major = element_line(color = "grey80"), panel.grid.minor = element_line(color="grey90"))

ggsave(paste(args$path, "/CCD/", "CCD_",args$output, ".png", sep=""), plot = changePlot, width=512, height=128, units ="mm", dpi=100)

#write the corresponding json file
titleText <- ""
for (i in seq(1,length(args$restrict),2)) #interate through the restrictions
{
    titleText <- paste(titleText, args$restrict[i], ":", args$restrict[i+1], ", ")
}
titleText <- substr(titleText, 1, nchar(titleText)-2)
descText <- paste("There are", toString(length(changePoints)), "change points in this study.", sep=" ")
#imageText <- paste(args$path, "/CCD/", "CCD_",args$x,"_",args$y, ".png", sep="")
imageText <- paste("CCD_", args$output, ".png", sep="")
jsonData <- list (title = titleText, desc = descText, parameters = list(args$x, args$y), date = Sys.time(), image = imageText, cinema = list(changecolumn=paste("CCD_", args$output, sep="")))

json<-toJSON(jsonData, auto_unbox=TRUE, pretty = TRUE)
write(json, file= paste(args$path, "/CCD/", "CCD_",args$output, ".json", sep=""))


# order segments into row for for csv file
csvOutputFile <- read.csv(file=args$csvFile, check.names=FALSE) #need to get full db now

changePoints <- changePoints + 1

changePtsIter = 1

changePts <- rep(0,dim(csvOutputFile)[1])

for (i in 1:length(changePoints))
{
    changePts[rowValues[changePoints[i]]] = 1;
}

#print(which(csvOutputFile$phi == -180))

#changePts <- 1:dim(csvInputFile)[1]

#changePoints <- changePoints + 1

#for(i in 1:length(changePts))
#{ if (changePoints[changePtsIter] == i) {
#    changePts[i] = 1
#    if (changePtsIter < length(changePoints)) { changePtsIter =
#        changePtsIter + 1 } } else { changePts[i] = 0; }
#}

#write output csvFile

csvOutputFile <- cbind(csvOutputFile, changePts)
colnames(csvOutputFile)[which(names(csvOutputFile) == "changePts")] <- paste("CCD_",args$output, sep="")
write.csv(csvOutputFile, args$csvFile, row.names=FALSE)


