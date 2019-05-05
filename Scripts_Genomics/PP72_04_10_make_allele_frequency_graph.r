rm(list=ls())

palette(colors())
library(stringr)
library("ggplot2")

# get the current folder (folder from where the script is executed)
library("rstudioapi", lib.loc="~/Library/R/3.3/library")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Assign the colors to each isolate
gg_color_hue_dark <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue_light <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 85, c = 100)[1:n]
}

myColorsDark = gg_color_hue_dark(5)
myColorsLight = gg_color_hue_light(5)

myColorAssign <- rbind( data.frame(Isolate=c("A1","A4","A5","B3","C2"),
                            Type=c("WG"),
                            myColors=myColorsDark),
                        data.frame(Isolate=c("A1","A4","A5","B3","C2"),
                            Type=c("RS"),
                            myColors=myColorsLight)     )

thousand.folders <- list.dirs(path=getwd(), full.names = TRUE)

for (folder in thousand.folders ) {
  print(folder)
  setwd(folder)
  
  #get Folder ID (GG, GL, MultiRep2_version2)
  foldername <- gsub("/", "",  gsub("PRV4", "",  tail( strsplit(as.character(folder),'_',fixed=TRUE)[[1]], 1) ) )
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  myFiles <- Sys.glob(      "*subVCF"       )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  print(myFiles)

  MeltedData <- NULL
  mygraph <- ggplot()
  ColorVect <- NULL
    
  if (length(myFiles) > 0) {
    for (myFile in sort(myFiles, decreasing = F)) {
      Dataset <- read.delim(myFile,h=T)
    
      mySamplename <- strsplit(myFile,'.',fixed=TRUE)[[1]][1]
      myIsolate <- strsplit(myFile,'-',fixed=TRUE)[[1]][1]
      myType <- substr(strsplit(myFile,'-',fixed=TRUE)[[1]][2],1,nchar(as.character(strsplit(myFile,'-',fixed=TRUE)[[1]][2]))-1)
      myFilename <- paste(myIsolate, ".nR-C.density_allele_freq.pdf", sep ="")
      
      print(mySamplename)
      Dataset <- subset(Dataset, grepl(",", Dataset$Occurrence_Alleles))
      Dataset$number_of_allele <- str_count(Dataset$Occurrence_Alleles, ",")+1
      Dataset <- subset(Dataset, Dataset$number_of_allele == 2)
      Dataset <- subset(Dataset, !grepl("I|D", Dataset$Types_Alleles) )
      Dataset <- subset(Dataset, Dataset$Repeat == 0)   # no repeats
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
}


