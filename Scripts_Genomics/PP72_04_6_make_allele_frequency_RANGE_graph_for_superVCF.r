rm(list=ls())

palette(colors())

library(stringr)
library("ggplot2")


# get the current folder (folder from where the script is executed)
library("rstudioapi", lib.loc="~/Library/R/3.3/library")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

thousand.folders <- list.dirs(path=getwd(), full.names = TRUE)

for (folder in thousand.folders ) {
  print(folder)
  setwd(folder)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  myFiles <- Sys.glob(      "*superVCF"       )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  print(myFiles)

  Dataset <- ""
  
  if (length(myFiles) > 0) {
    for (myFile in myFiles) {
      myY_lim =""
      myShortSampleName <- strsplit(myFile,'-on-',fixed=TRUE)[[1]][1]
      mySamplename <- sub(".superVCF", "", myFile)
      myFilename <- paste(mySamplename, ".superVCF.allele_freq_cov_", sprintf("%03d", 0),".pdf", sep ="")
      pdf(myFilename, width = 4, height = 4, useDingbats=F)
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1,1,myShortSampleName, adj = c(0.5, NA), cex = 3)
      dev.off()
      
      
      
      for (myCovThreshold in seq(10,130, by = 20)) {
        Dataset <- read.delim(myFile,h=T)
      
        mySamplename <- sub(".superVCF", "", myFile)
        
        
        Dataset <- subset(Dataset, grepl(",", Dataset$Occurrence_Alleles))
        Dataset$number_of_allele <- str_count(Dataset$Occurrence_Alleles, ",")+1
        Dataset <- subset(Dataset, Dataset$number_of_allele == 2)
        Dataset <- subset(Dataset, !grepl("I|D", Dataset$Types_Alleles) )
        Dataset <- subset(Dataset, Dataset$Repeat == 0)   # no repeats
        
        Dataset$Total_occ <- sapply(Dataset$Occurrence_Alleles,    function(x) sum( as.numeric(unlist(strsplit(as.character(x), split=","))) ) )
        
        Dataset <- subset(Dataset, Dataset$Total_occ >= myCovThreshold)   # cov min
        myFilename <- paste(mySamplename, ".superVCF.allele_freq_cov_", sprintf("%03d", myCovThreshold),".pdf", sep ="")
        
        
        Dataset$random_allele_occ <- sapply(Dataset$Occurrence_Alleles,   function(x)  sample(as.numeric(unlist(strsplit(as.character(x), split=","))), 1 ) )
  
        Dataset$freq <- Dataset$random_allele_occ/Dataset$Total_occ
        if (myY_lim =="") {
          myGraph <- ggplot(Dataset, aes(x = Dataset$freq)) + geom_histogram(binwidth = 0.033, alpha = 0.5, fill="black", position = 'identity')  +
            labs(x = "Allele frequency", y = "Counts" )
          myY_lim <- ggplot_build(myGraph)$layout$panel_ranges[[1]]$y.range
        }
        pdf(myFilename, width = 4, height = 4, useDingbats=F)
          myGraph <- ggplot(Dataset, aes(x = Dataset$freq)) + geom_histogram(binwidth = 0.033, alpha = 0.5, fill="black", position = 'identity')  +
            labs(x = "", y = "" ) +
            scale_y_continuous(limits = myY_lim) +
            theme(text = element_text(size=12)) 
          print(myGraph)
        
        dev.off()
      }
    }
  }
}



