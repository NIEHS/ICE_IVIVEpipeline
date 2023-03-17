EADboxplot <- function(EAD.out, invivo = NULL, label = "EAD",chemDisplay = "CASRN", EADUnit = "mg/kg/dose", modelType = "1C", species = "human", route = "iv", date = "2021", doses.per.day = 1, wth = 0.2, axis.text = 8, title.text = 10, shape.size = 1, dpi = 300)
{
  # wth -- the width of the boxplot (for plot)
  # shape.size -- the size of shape indicating in vivo data points (for plot)
  # axis.text -- the font size of axis label
  # title.text -- the font size of title label
  library(dplyr)  #example function: %>%, filter(), slice(), arrange(), select(),rename(),distinct(), mutate(), transmute(), summarise(), sample_n()
  library(tidyr)  #example function:gather()
  library(scales) #example function:trans_breaks()
  library(ggplot2)

  ## prepare EAD data for plotting
  seldata <- EAD.out[grep(label,names(EAD.out))]
  selEAD <- cbind(EAD.out$CASRN, EAD.out$ChemicalName, seldata)
  names(selEAD)[1:2] <- c("CASRN","ChemicalName")
  EAD.up <- selEAD %>% gather("Category", "EAD", 3:ncol(selEAD)) # this is to unpivot the data to make plotting easy
  maxy <- max(EAD.up$EAD, na.rm=T)

  if (is.null(invivo))
  {
    plot1 <- ggplot() + geom_boxplot(data=EAD.up, aes_string(x=chemDisplay, y="EAD"), colour = "darkgreen",outlier.colour = "NA", outlier.shape = 16, width=wth) +
      theme_bw() + #set background to white instead of gray
      ggtitle(paste(label," (",modelType,"_",species,"_",route,")", sep="","\n"))+
      theme(axis.text.x=element_text(size=axis.text, angle =65, hjust = 1, colour="black")) +
      theme(axis.text.y=element_text(size=axis.text, angle = 0, hjust = 1, colour="black")) +
      scale_x_discrete(label = function(x) abbreviate(x, minlength=12)) +
      scale_y_log10(paste(label," ","(", EADUnit,")",sep="","\n"), breaks=trans_breaks("log10", function(x) 10 ^ x)(c(1e-6, maxy))) +
      theme(plot.title = element_text(hjust=0.5,size=title.text, face="bold"),legend.text = element_text(size = title.text), axis.title.x = element_text(size=title.text,angle=0,face="bold",colour="black"),
            axis.title.y = element_text(size=title.text,angle=90,face="bold",colour="black"))  +
      labs(x='Chemical')
    plot1
    ggsave(plot1, file=paste(label, "_", modelType, "_", species, "_", route, "_dpd", doses.per.day, "_",date, ".jpeg"), width=7,height=5, dpi=dpi,units = "in", device='jpeg')
  }


  if(!is.null(invivo))
  {
    ## prepare in vivo data for plotting
    selchem <- EAD.out%>% select(CASRN, ChemicalName)
    invivo0 <- merge.data.frame(selchem, invivo, by= c("CASRN"), all = TRUE) ## merge() two datasets based on c("column Name")
    invivo1 <- invivo0[!(is.na(invivo0$ChemicalName)),]  # remove extra chemicals in the invivo dataset
    invivo.up <- invivo1  %>% gather("Category", "invivoLevel", 3:ncol(invivo1)) # this is to unpivot the data to make plotting easy
    plot1 <- ggplot() + geom_boxplot(data=EAD.up, aes_string(x=chemDisplay, y="EAD"), colour = "darkgreen",outlier.colour = "NA", outlier.shape = 16, width=wth) +
    geom_point(data=invivo.up, aes_string(x=chemDisplay, y="invivoLevel", shape="Category"), size=shape.size, alpha=0.5, colour="mediumblue") +
    theme_bw() + #set background to white instead of gray
    ggtitle(paste(label," (",modelType,"_",species,"_",route,")", " and in vivo levels ",sep="","\n"))+
    theme(axis.text.x=element_text(size=axis.text, angle = 65, hjust = 1, colour="black")) +
    theme(axis.text.y=element_text(size=axis.text, angle = 0, hjust = 1, colour="black")) +
    scale_x_discrete(label = function(x) abbreviate(x, minlength=12)) +
    scale_y_log10(paste(label," or in vivo level", " ","(", EADUnit,")",sep="","\n"), breaks=trans_breaks("log10", function(x) 10 ^ x)(c(1e-6, maxy))) +
    theme(plot.title = element_text(hjust=0.5,size=title.text, face="bold"),legend.text = element_text(size = title.text), axis.title.x = element_text(size=title.text,angle=0,face="bold",colour="black"),
          axis.title.y = element_text(size=title.text,angle=90,face="bold",colour="black"))  +
    labs(x='Chemical')
    plot1
    ggsave(plot1, file=paste(label, "_invivo_", modelType, "_", species, "_", route, "_dpd", doses.per.day, "_", date, ".jpeg"), width=7,height=5, dpi=dpi,units = "in", device='jpeg')
  }
  plot1

}


