rm(list=ls())

palette(colors())

library("UpSetR")

# get the current folder (folder from where the script is executed)
library("rstudioapi", lib.loc="~/Library/R/3.3/library")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

thousand.folders <- list.dirs(path=getwd(), full.names = TRUE)

for (folder in thousand.folders ) {
  print(folder)
  setwd(folder)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #myFiles <- Sys.glob(      "*norepcod-noInDels_*.ForVenn.txt"       )
  
  myFiles <- Sys.glob(      "*all-but-noInDels_*.ForVenn.txt"       )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  print(myFiles)

  Dataset <- ""
  
  if (length(myFiles) > 0) {
    for (myFile in myFiles) {
      myTable <- read.delim(myFile,h=F)
    
      myfilename <- strsplit(myFile,'.',fixed=TRUE)[[1]][1]
      mySampleName <- strsplit(myfilename,'-on-',fixed=TRUE)[[1]][1]
      myIsolateName <- strsplit(mySampleName,'-',fixed=TRUE)[[1]][1]
      myGenomeName <- strsplit( strsplit(myfilename,'-on-',fixed=TRUE)[[1]][2], '_',fixed=TRUE)[[1]][1]
      myGenomeMethod <- paste(strsplit(myfilename,'_',fixed=TRUE)[[1]][3],
                              substr(strsplit(myfilename,'_',fixed=TRUE)[[1]][4], 1, 100),
                              sep="_")
    
      myList <- as.list(myTable)
      names(myList) <- mySampleName
    
      if (Dataset != "") {
        Dataset <- c(Dataset, myList)
      } else {
        Dataset <- myList
      }
    }
    pdf(paste("UpSet_Graph_", myIsolateName,"-on-", myGenomeName,"_", myGenomeMethod,".pdf", sep=""), width = 10, height = 6, onefile=FALSE, useDingbats=F)
    upset(fromList(Dataset), nsets = 7, number.angles = 0, point.size = 3.5, line.size = 0.5, 
          mainbar.y.label = "Intersection size", sets.x.label = "        Number of poly-allelic sites", 
          text.scale = c(1.8, 2.5, 1.8, 1.8, 2.5, 2.5),  order.by = "freq", nintersects= 15, mb.ratio = c(0.65, 0.35))
    dev.off()
    # text.scale = c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
    
    all_commons <- Reduce(intersect, Dataset)
    mySize <- length(all_commons)
    write.table(all_commons, file=paste("Common_elments_",mySize,"_",myGenomeMethod,".txt" ,sep=""), sep="\t")
  }
  
  
}


