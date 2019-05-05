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
  myFiles <- Sys.glob(      "*subVCF"       )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  print(myFiles)

  Dataset <- ""
  
  if (length(myFiles) > 0) {
    for (myFile in myFiles) {
      Dataset <- read.delim(myFile,h=T)
    
      mySamplename <- strsplit(myFile,'.',fixed=TRUE)[[1]][1]
      myFilename <- paste(mySamplename, ".allele_freq.pdf", sep ="")
      
      Dataset <- subset(Dataset, grepl(",", Dataset$Occurrence_Alleles))
      Dataset$number_of_allele <- str_count(Dataset$Occurrence_Alleles, ",")+1
      Dataset <- subset(Dataset, Dataset$number_of_allele == 2)
      Dataset <- subset(Dataset, !grepl("I|D", Dataset$Types_Alleles) )
      Dataset <- subset(Dataset, Dataset$Repeat == 0)   # no repeats
      
      Dataset$Total_occ <- sapply(Dataset$Occurrence_Alleles,    function(x) sum( as.numeric(unlist(strsplit(as.character(x), split=","))) ) )
                        
      Dataset$random_allele_occ <- sapply(Dataset$Occurrence_Alleles,   function(x)  sample(as.numeric(unlist(strsplit(as.character(x), split=","))), 1 ) )

      Dataset$freq <- Dataset$random_allele_occ/Dataset$Total_occ
      
      pdf(myFilename, width = 4, height = 4, useDingbats=F)
      print(ggplot(Dataset, aes(x = Dataset$freq)) + geom_histogram(binwidth = 0.033, alpha = 0.5, fill="black", position = 'identity') + ggtitle(mySamplename) +
        labs(x = "Allele frequency", y = "Counts" ) +
        theme(text = element_text(size=12)) )
      
      dev.off()
      
    }
  }
}



