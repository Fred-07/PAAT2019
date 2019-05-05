rm(list=ls())

#------------Function------------
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#------------------------------------




palette(colors())

# load libraries
library("RColorBrewer");
library(ggplot2);
library(reshape2);
library(plyr)
library(readr)


# get the current folder (folder from where the script is executed)
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Dataset <- read.delim("List_allelic_freq_for_Deletion_20170713.txt", row.names=1, header=F)

Dataset <- t(Dataset)
nbcol <- ncol(Dataset)
col.names <- as.vector(colnames(Dataset))
breaks.ncol <- as.numeric(1:nbcol)

Dataset <- melt(Dataset)

# Remove big value 
Dataset <- subset(Dataset, Dataset$value < 20000)

Dataset.summary <- summarySE(Dataset, measurevar="value", groupvars="Var2")

pdf(paste("Deletion_Allelic-freq__violin-plot.pdf", sep=""), width = 3, height = 2, useDingbats=T)
myGraph <- ggplot()+
  geom_violin(aes(Var2, value), data = Dataset, colour = I("black"), size =0.25, fill="grey") +
  geom_crossbar(data=Dataset.summary,aes(x=Var2,ymin=value, ymax=value,y=value,group=Var2), colour = I("blue"), width = 0.5, fatten=1) +
  scale_x_discrete(drop=FALSE) +
  labs(x = "Isolates", y = "Allelic frequency" ) +
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
myGraph
dev.off()


  
