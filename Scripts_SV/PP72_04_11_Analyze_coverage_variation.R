rm(list=ls())

palette(colors())

library(stringr)
library("ggplot2")


# get the current folder (folder from where the script is executed)
library("rstudioapi", lib.loc="~/Library/R/3.3/library")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
myFiles <- Sys.glob(      "Ratio_Coverage_genome*.txt"       )
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print(myFiles)

if (length(myFiles) > 0) {
  for (myFile in sort(myFiles, decreasing = F)) {
    Dataset <- read.delim(myFile,h=T)
    
    attach(Dataset)
    hist(Dataset[,5])
  }
}

hist(Dataset[,4], breaks = 400)
hist(Dataset[,5], breaks = 400)
hist(Dataset[,6], breaks = 400)
hist(Dataset[,7], breaks = 400)
hist(Dataset[,8], breaks = 400)
hist(Dataset[,9], breaks = 400)
hist(Dataset[,10], breaks = 400)
hist(Dataset[,11], breaks = 400)
hist(Dataset[,12], breaks = 400)
hist(Dataset[,13], breaks = 400)
hist(Dataset[,16], breaks = 400)
hist(Dataset[,17], breaks = 400)
hist(Dataset[,18], breaks = 400)



  MeltedData <- NULL
  mygraph <- ggplot()
  ColorVect <- NULL
    
  if (length(myFiles) > 0) {
    for (myFile in sort(myFiles, decreasing = F)) {
      Dataset <- read.delim(myFile,h=T)
    
      mySamplename <- strsplit(myFile,'.',fixed=TRUE)[[1]][1]
      myIsolate <- strsplit(myFile,'-',fixed=TRUE)[[1]][1]
      myType <- substr(strsplit(myFile,'-',fixed=TRUE)[[1]][2],1,nchar(as.character(strsplit(myFile,'-',fixed=TRUE)[[1]][2]))-1)
      myFilename <- paste(myIsolate, ".R.density_allele_freq.pdf", sep ="")
      
      print(mySamplename)
      Dataset <- subset(Dataset, grepl(",", Dataset$Occurrence_Alleles))
      Dataset$number_of_allele <- str_count(Dataset$Occurrence_Alleles, ",")+1
      Dataset <- subset(Dataset, Dataset$number_of_allele == 2)
      Dataset <- subset(Dataset, !grepl("I|D", Dataset$Types_Alleles) )
      Dataset <- subset(Dataset, Dataset$Repeat != 0)   # no repeats
      Dataset <- subset(Dataset, Dataset$Coding == 1)   # coding
      
      myDim <- c(mySamplename, dim(Dataset))
      print(myDim)
      write.table(myDim, file=paste("Number_of_positions_for", mySamplename,"__", myFilename, ".txt", sep=""), sep="\t")
      
      Dataset$Total_occ <- sapply(Dataset$Occurrence_Alleles,    function(x) sum( as.numeric(unlist(strsplit(as.character(x), split=","))) ) )
                        
      Dataset$random_allele_occ <- sapply(Dataset$Occurrence_Alleles,   function(x)  sample(as.numeric(unlist(strsplit(as.character(x), split=","))), 1 ) )

      Dataset$freq <- Dataset$random_allele_occ/Dataset$Total_occ
      
      MeltedData <- rbind(MeltedData, data.frame(mySamplename, Dataset$freq))
      
      ColorVect <- c(ColorVect, as.vector((subset(myColorAssign, myColorAssign$Isolate == myIsolate & myColorAssign$Type == myType ))[[3]]) )
      
    }
    colnames(MeltedData) <- c("mySamplename", "freq")
    mygraph <- ggplot(MeltedData, aes(x=freq)) +
      stat_density(aes(group=mySamplename, colour=mySamplename), size = 1, position="identity",geom="line") +
      ggtitle(myIsolate) +
      labs(x = "Allele frequency", y = "Density" ) +
      scale_color_manual(values = ColorVect) +
      theme(text = element_text(size=8),
            legend.position = "none",
            panel.border = element_rect(colour = "grey40", fill=NA, size=1),
                    plot.title = element_text(size = 8),
                    panel.background = element_rect(fill = NA),
                    plot.background = element_rect(colour = NA),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    axis.title = element_text(size = 8),
                    axis.title.y = element_text(angle=90,vjust =1.5),
                    axis.title.x = element_text(vjust = -0.1),
                    axis.text = element_text(), 
                    axis.ticks = element_line()) 
      

    
    pdf(myFilename, width = 2, height = 2, useDingbats=F)
      print(mygraph)
    dev.off()  
      
  }




