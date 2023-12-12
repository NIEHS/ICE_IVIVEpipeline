####functions
#Name: CalcEAD
#Author: ILS (comptox@ils-inc.com)
#Date:05 April, 2019
#Version: 1.0
#License: MIT
#Summary: takes 2 required inputs: plasma concentration at steady state, in vitro activity to calculate the estimated administerd dose
#Css:dataframe containing the CASRN and the circulating chemical concentration, assumed the dosage is 1mg/kg/day along with optional values for correction factors
#Css must include "CASRN" as the first column and the second is the corresponding circulating concentration
#inVitro:dataframe of ACC or AC50 values along with "CASRN" as the first column (ID column): ex: CASRN, ACC1, ACC2, ACCn
#adj.fu= (optional, valid only for 1 compartment models) Column name corresponding to fraction unbound measures to be used in adjusting the EADin Css
#adj.armitage = (optional) Column name corresponding to css that lists the armitage adjustment factor
#Library required: none


CalcEAD <- function(Css, inVitro, adj.fu=NULL, adj.armitage=NULL, ...){
  options(stringsAsFactors = FALSE)

  if(is.null(adj.fu)){ Css$adj.fu<-1}else{Css$adj.fu<-Css[,adj.fu]} #this tests if the fu adjustment is provided/wanted, if not replaces with 1 if it is uses user provided column
  if(is.null(adj.armitage)){ Css$adj.arm<-1}else{Css$adj.arm<-Css[,adj.armitage]}
  inVitro<-inVitro[,c(colnames(inVitro) %in% sub("chem.*", "", colnames(inVitro), ignore.case = TRUE))]#this drops out any"chemical name" field that would be passed along
  inVitro.m <- merge(Css, inVitro, by= c("CASRN"), all = TRUE) ## merge() two datasets based on the chemical ID
  conc<-as.numeric(inVitro.m[,colnames(Css)[2]])#this is the circulating concentraiton
  #conc<-as.numeric(inVitro.m[,grep("%",colnames(inVitro.m))])#this is the circulating concentraiton
  adj.fu<-as.numeric(inVitro.m$adj.fu)
  adj.arm<-as.numeric(inVitro.m$adj.arm)
  i<-ncol(Css)+1 #this skips the Css object info

  if(nrow(inVitro.m) ==1){
    EAD <- t(sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) adj.arm*x/conc/adj.fu))

  } else {
    EAD <- sapply(inVitro.m[,i:ncol(inVitro.m)], function(x) adj.arm*x/conc/adj.fu)
  }
  colnames(EAD) <- paste("EAD", colnames(EAD), sep = "_")

  Output <- cbind(inVitro.m, EAD)
  return(Output)

}
