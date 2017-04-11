#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file

## 06102016_NmrAlignment_wrapper.R
## Marie Tremblay-Franco
## marie.tremblay-franco@toulouse.inra.fr

runExampleL <- FALSE

##------------------------------
## Options
##------------------------------
strAsFacL <- options()$stringsAsFactors
options(stringsAsFactors = FALSE)


##------------------------------
## Libraries loading
##------------------------------
	# ParseCommandArgs function
library(batch)
	# Alignment
library(speaq)


# R script call
source_local <- function(fname)
{
	argv <- commandArgs(trailingOnly = FALSE)
	base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
	source(paste(base_dir, fname, sep="/"))
}
# Function import
source_local("NmrAlignment_script.R")


##------------------------------
## Errors ?????????????????????
##------------------------------


##------------------------------
## Constants
##------------------------------
topEnvC <- environment()
flagC <- "\n"


##------------------------------
## Script
##------------------------------
if(!runExampleL)
    argLs <- parseCommandArgs(evaluate=FALSE)


## Parameters Loading
##-------------------
	# Inputs
		## Library of spectra to align
if (!is.null(argLs[["zipfile"]])){
	fileType="zip"
	zipfile= argLs[["zipfile"]]
	directory=unzip(zipfile, list=F)
	directory=paste(getwd(),strsplit(directory[1],"/")[[1]][2],sep="/")
} else if (!is.null(argLs[["tsvfile"]])){
	fileType="tsv"
	directory <- read.table(argLs[["tsvfile"]],check.names=FALSE,header=TRUE,sep="\t")
}


		## Spectral width
leftBorder <- argLs[["left_border"]]
rightBorder <- argLs[["right_border"]]

		##Exclusion zone(s)
exclusionZones <- argLs[["zone_exclusion_choices.choice"]]
exclusionZonesBorders <- NULL
if (!is.null(argLs$zone_exclusion_left))
{
   for(i in which(names(argLs)=="zone_exclusion_left"))
   {
     exclusionZonesBorders <- c(exclusionZonesBorders,list(c(argLs[[i]],argLs[[i+1]])))
   }
}

		## Reference spectrum
reference <- argLs[["reference"]]

		## Size of a small nDivRange
nDivRange <- argLs[["nDivRange"]]

		## Intensity threshold for peak removal
baselineThresh <- argLs[["baselineThresh"]]


	# Outputs
logOut <- argLs[["logOut"]]
alignedSpectra <- argLs[["alignedSpectra"]]
graphOut <- argLs[["graphOut"]]


## Checking arguments
##-------------------
error.stock <- "\n"
if(length(error.stock) > 1)
  stop(error.stock)


## Computation
##------------
directory.alignement <- nmr.alignment(fileType=fileType,directory=directory,leftBorder=leftBorder,rightBorder=rightBorder,exclusionZones=exclusionZones,
                                  exclusionZonesBorders=exclusionZonesBorders, reference=reference, nDivRange=nDivRange,
                                  baselineThresh=baselineThresh, maxshift=50, verbose=FALSE)
directory.raw <- directory.alignement[[1]]
directory.aligned <- directory.alignement[[2]]

## Saving
##-------
	# Aligned spectra
t.directory.aligned <- t(directory.aligned)
rownames(t.directory.aligned) <- colnames(directory.aligned)
# colnames(t.directory.aligned) <- c("Bucket",colnames(t.directory.aligned))
write.table(t.directory.aligned,file=alignedSpectra,row.names=TRUE,quote=FALSE,sep="\t")


excludedZone <- NULL
for (c in 1:length(exclusionZonesBorders))
{
  excludedZone <- c(excludedZone,exclusionZonesBorders[[c]])
  excludedZone <- sort(excludedZone)
}

## Graphical output: overlay of raw and estimated spectra
pdf(graphOut,onefile=TRUE)
par(mfrow=c(2,1))

raw.spectra <- data.frame(directory.raw)
colnames(raw.spectra) <- substr(colnames(raw.spectra),2,7)

aligned.spectra <- data.frame(directory.aligned)
colnames(aligned.spectra) <- substr(colnames(aligned.spectra),2,7)

drawSpec(raw.spectra,xlab="", ylab="Raw spectra", main="")
drawSpec(aligned.spectra,xlab="", ylab="Aligned spectra", main="")

nbZones <- length(excludedZone)/2
if (nbZones != 0)
{
  n <- length(excludedZone)
  drawSpec(raw.spectra[,1:which(round(as.numeric(colnames(raw.spectra)),2) == excludedZone[n])[1]],xlab="", ylab="Raw spectra", main="")
  drawSpec(aligned.spectra[,1:which(round(as.numeric(colnames(aligned.spectra)),2) == excludedZone[n])[1]],xlab="", ylab="Aligned spectra", main="")

  n <- n - 1
  while (n >= nbZones & nbZones > 1)
  {
    drawSpec(raw.spectra[,(which(round(as.numeric(colnames(raw.spectra)),2) == excludedZone[n])[1]):(which(round(as.numeric(colnames(raw.spectra)),2) == excludedZone[n-1])[1])],xlab="", ylab="Raw spectra", main="")
    drawSpec(aligned.spectra[,(which(round(as.numeric(colnames(aligned.spectra)),2) == excludedZone[n])[1]):(which(round(as.numeric(colnames(aligned.spectra)),2) == excludedZone[n-1])[1])],xlab="", ylab="Aligned spectra", main="")
    n <- n - 2
  }

  drawSpec(raw.spectra[,(which(round(as.numeric(colnames(raw.spectra)),2) == excludedZone[1])[1]):ncol(raw.spectra)],xlab="", ylab="Raw spectra", main="")
  drawSpec(aligned.spectra[,(which(round(as.numeric(colnames(aligned.spectra)),2) == excludedZone[1])[1]):ncol(aligned.spectra)],xlab="", ylab="Aligned spectra", main="")
}
drawSpec(raw.spectra[,(which(round(as.numeric(colnames(raw.spectra)),2) == 2.4)[1]):(which(round(as.numeric(colnames(raw.spectra)),2) == 2.8)[1])],xlab="", ylab="Raw spectra", main="")
drawSpec(aligned.spectra[,(which(round(as.numeric(colnames(aligned.spectra)),2) == 2.4)[1]):(which(round(as.numeric(colnames(aligned.spectra)),2) == 2.8)[1])],xlab="", ylab="Aligned spectra", main="")

dev.off()


## Ending
##---------------------
cat("\nEnd of 'NMR alignment' Galaxy module call: ", as.character(Sys.time()), sep = "")
options(stringsAsFactors = strAsFacL)
rm(list = ls())
