#---
#A command line multi-variate change detection code
#---

#--
# example usage command
# Rscript multiVariatechangeDetection.R -f /path/to/data.csv -p /path/to/nyx.cdb -r phi -180 theta -45 iso 1 -x time -y "entropy" "image canny count" "sample" -cd 0.05 199 2 0.05 -o iso1_Time_vs_Variables
#--

# clear all previous environment variables and memory
rm(list=ls())
#gc()

#current directory should already be working directory
#print(getwd())

#install required libraries
if(!require(argparse)){install.packages("argparse")}
if(!require(ecp)){install.packages("ecp")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(reshape)){install.packages("reshape")}
if(!require(jsonlite)){install.packages("jsonlite")}

# load required libraries
library("argparse") #to parse the commandline input
library("ecp") #to execute the multi-variate change detection
library("ggplot2") #to draw the resulting plot
library("reshape") #helper function for getting data into plot
library("jsonlite") #to save the json file

# parse the arguments to understand next steps
parser <- ArgumentParser()

parser$add_argument("-f", "--csvFile", help='csv input file with data')
parser$add_argument("-p", "--path", help='path to .cdb')
parser$add_argument("-r", "--restrict", nargs='+', help='restrict rows based on list of column names and values')
parser$add_argument("-x", "--xVal", help='the x paramter column to iterate over')
parser$add_argument("-y", "--yVal", nargs='+', help='the y parameter column(s) to find change over')
parser$add_argument("-cd", "--cdParams", type="double", nargs='+', help='parameters for ecp (divisive) change detection algorithm')
parser$add_argument("-o", "--output", help='output phrase to mark files by')

args <- parser$parse_args()

# read input csvFile
csvInputFile <- read.csv(file=args$csvFile)

rowValues <- vector(mode = "numeric", length=0)

# cull out the needed rows according to the restrictions given
if (length(args$restrict) > 0)
{
    for (i in seq(1,length(args$restrict),2)) #interate through the restrictions
    {
        csvInputFile <- csvInputFile[csvInputFile[,args$restrict[i]]==as.integer(args$restrict[i+1]),]
        rowValues <- as.numeric(rownames(csvInputFile))
    }
} else
{
    rowValues <- inputParam <- 1:dim(csvInputFile)[1]
}

#get feature data
featureData <- vector("list", length(args$y))
for (i in seq(1, length(args$y)))
{
    noSpaces <- gsub(" ",".",args$y[i])
    featureMetric <- csvInputFile[,noSpaces]
    featureData[[i]] <- featureMetric
}
inputValues <- csvInputFile[,args$x]

# --- change detection using ecp library

#set parameters for ecp division technique
sigLevel <- as.double(args$cd[1])
iterations <- as.integer(args$cd[2])
minSize <- as.integer(args$cd[3])
alpha <- as.double(args$cd[4])

normData <- vector("list", length(args$y))
for (i in seq(1, length(args$y)))
{
    if ((max(featureData[[i]]) - min(featureData[[i]])) != 0)
    {
        normD = (featureData[[i]] - min(featureData[[i]])) / (max(featureData[[i]]) - min(featureData[[i]]))
        normData[[i]] <- normD
    } else
    {
        normD <- rep(0.5,dim(csvInputFile)[1])
        normData[[i]] <- normD
    }
}

data = matrix(normData[[1]])
if (length(args$y) >= 2)
{
    for (i in seq(2, length(args$y)))
    {
        rows = dim(data)[1]
        cols = dim(data)[2]
        data = matrix(c(data, normData[[i]]), nrow=rows, ncol=cols+1)
    }
}

changeResults = e.divisive(X=data,sig.lvl=sigLevel,R=iterations,k=NULL,min.size=minSize,alpha=alpha)

#------ change detection concluded

#----saving results to csv, as image and as json---"

#if resulting directory doesn't exist, create it
if (!dir.exists(file.path(args$path, "CCD")))
dir.create(file.path(args$path, "CCD"))

#edit last change point location if needed
changePoints <- changeResults$estimates
if (changePoints[length(changePoints)] > length(inputValues))
{
    changePoints[length(changePoints)] = length(inputValues)
}

#draw graph#
changeDF <- data.frame(id=inputValues)
for (i in seq(1, length(args$y)))
{
    changeDF[args$y[i]] <- normData[[i]]
}
headings <- args$y[1]
if (length(args$y) >= 2)
{
    for (i in seq(2, length(args$y)))
    {
        headings <- c(headings,args$y[i])
    }
}
changeDF.long <- melt(changeDF, id="id", measure = headings)
changePlot <- ggplot(changeDF.long, aes(id, value, color = variable)) + geom_line() + labs(x= args$x, y = "Normalized Value")
for (i in seq(1, length(changeResults$estimates)))
{
    changePlot = changePlot + geom_vline(xintercept = inputValues[changePoints], linetype="dotted", color = "black", size = 1)
}
ggsave(paste(args$path, "/CCD/", "CCD_",args$output, ".png", sep=""), plot = changePlot, width=512, height=128, units ="mm", dpi=100)

#write the corresponding json file#
titleText <- ""
if (length(args$restrict) > 0)
{
    for (i in seq(1,length(args$restrict),2)) #interate through the restrictions
    {
        titleText <- paste(titleText, args$restrict[i], ":", args$restrict[i+1], ", ")
    }
    titleText <- substr(titleText, 1, nchar(titleText)-2)
} else
{
    titleText <- paste("Change Point Detection with no Restrictions")
}

descText <- paste("There are", toString(length(changeResults$estimates)), "change points in this study.", sep=" ")
#imageText <- paste(args$path, "/CCD/", "CCD_",args$x,"_",args$y, ".png", sep="")
imageText <- paste("CCD_", args$output, ".png", sep="")
jsonData <- list (title = titleText, desc = descText, parameters = list(args$x, args$y), date = Sys.time(), image = imageText, cinema = list(changecolumn=paste("CCD_", args$output, sep="")))

json<-toJSON(jsonData, auto_unbox=TRUE, pretty = TRUE)
write(json, file= paste(args$path, "/CCD/", "CCD_",args$output, ".json", sep=""))

#add change points to csv file#


csvOutputFile <- read.csv(file=args$csvFile, check.names=FALSE) #need to get full db now
changePtsIter = 1
changePts <- rep(0,dim(csvOutputFile)[1])

for (i in 1:length(changePoints))
{
    changePts[rowValues[changePoints[i]]] = 1;
}

#write output csvFile
csvOutputFile <- cbind(csvOutputFile, changePts)
colnames(csvOutputFile)[which(names(csvOutputFile) == "changePts")] <- paste("CCD_",args$output, sep="")
write.csv(csvOutputFile, args$csvFile, row.names=FALSE)
