####functions
#Name: SteadyState
#Author: ILS (comptox@ils-inc.com)
#Date:26 April, 2019
#Version: 1.0
#License: MIT
#Summary: takes 3 inputs(required input "inputData") to calculate the plasma concentration at steady state of the chemicals
#inputData:dataframe containing chemicals and their clearance parameters. Must have the following structure/column names:
#CASRN, Chemical, fu, Clint, MW
#Note that Clint is the invitro (hepatocyte) clearance rate with units of ul/min/10^6 cells
#nsamples: number of times to sample the Monte Carlo distribution of the population variance
#species: one of human or rat
#ConcentrationUnit: what units should the Css be in? choice of "uM" or "mg/L"; defaults to "uM"
#Library required: none


steadyState <- function(inputData, nsamples=300, species = "human", ConcentrationUnit="uM"){
  options(stringsAsFactors = FALSE)
  oriNames<-colnames(inputData)#saving for output
  colnames(inputData)<-tolower(colnames(inputData))#to deal with any capitalization issues
  all.perc <- c("50perc", "95perc")#these are the population percentages to sample for MC
  species = tolower(species) #this will take care of any issues in formatting

  Dose         <- rep(1/24, nrow(inputData))  # daily dose of 1 mg/kg/day converted to a unit of mg/kg/h
  fu          <- inputData$fu  # unitless
  CLinvitroMean <-inputData$clint # This is the in vitro (hepatocyte) clearance with the unit "ul/min/10^6 cells"
  #calculate the clearance rates
  if(species =="rat"){
    CLrenalMean  <- 0.08*fu  # rat glomerular filtration rate is 0.08(L/h)
    CLintrinMean <- 0.06555 * CLinvitroMean # converting factor from CLinvitroMean (ul/min/10^6 cells, which is the units from OPERA) to CLintrinMean (L/h) for rat, 0.06555 = (115*9.5)*60/1e6
  }else
  {
    CLrenalMean <-  6.7*fu                 # L/h, will match so long as the front part in the column name is "renal", case sensitive. human glomerular filtration rate is 6.7 (L/h)
    CLintrinMean <- 10.2* CLinvitroMean         # converting factor from CLinvitroMean (ul/min/10^6 cells, which is what predicted from OPERA model) to CLintrinMean (L/h) for human, 10.2= (103*1650)*60/1e6
  }

  MW           <- inputData$mw  # molecular weight
  Css          <- NULL  # Steady state blood concentration
  CssAll       <- NULL  # CssAll is accumulated Css value obtained from each round of calculation
  BWAll        <- NULL  # BW is body weight, BWAll is accumulated BW value obtained from each round of calculation
  QlAll        <- NULL  # Ql is liver blood flow, QlAll is accumulated Ql value obtained from each round of calculation
  CLintrinAll  <- NULL  # CLintrin is intrinsic clearance (L/h), CLintrinAll is accumulated CLintrin value obtained from each round of calculation
  CLrenalAll   <- NULL  # CLrenal is renal clearance (L/h), CLrenalAll is accumulated CLrenal value obtained from each round of calculation

  for (j in 1:nsamples) {           # The P-PK modeling and Monte Carlo simulation
    if(species == "rat")
    { BW <- rnorm(1, 0.25, 0.25*0.2)  # assume a normal distribution, mean 0.25, std 0.25*0.2
    Ql <- rnorm(1,13750, 13750*0.2)*60/10^6  # assume normal distribution, mean 13750 (uL/min), Sd = 13750*0.2, also change the unit from ul/min to L/h
    } else
    {
      BW <- rnorm(1, 70, 70*0.2)  # assume a normal distribution, mean 70, std 705*0.2
      Ql <- rnorm(1,1.65E6, 1.65E6*0.2)*60/10^6  # assume normal distribution, mean 1.65E6 (uL/min), Sd = 1.65E6*0.2, also change the unit from ul/min to L/h
    }
    CLrenal <- rnorm(length(CLrenalMean), CLrenalMean, CLrenalMean*0.2)  # assume the coefficeient of variance (cv) is 0.2 for renal clearance, cv=sd/mean
    CLintrin <- rnorm(length(CLintrinMean), CLintrinMean, CLintrinMean*0.2)  # assume the coefficeient of variance (cv) is 0.2 for intrinsic clearance, cv=sd/mean
    Css <- Dose*BW/(CLrenal+(Ql*fu*CLintrin)/(Ql+fu*CLintrin))  # unit is mg/L

    if (ConcentrationUnit =="mg/L")
    {
      Css <- Css
    }else{
      Css <- Css*1000/MW  # convert unit from mg/L to uM; set as default
      ConcentrationUnit ="uM" #this sets it as the default if it doesnt match the alternate option
    }
    CssAll <- cbind(CssAll,Css)
    BWAll[j] <- BW
    QlAll[j] <- Ql
    CLrenalAll <- cbind(CLrenalAll,CLrenal)
    CLintrinAll <- cbind(CLintrinAll,CLintrin)
  }
  CssAll_temp <- as.data.frame(t(apply(CssAll, 1, quantile, probs= c(.5, .95)))) #get the 50th and 95th percentile for Css at 1mg/kg/day
  CssAll <- cbind(inputData, CssAll_temp,ConcentrationUnit)
  colnames(CssAll)<-c(oriNames, colnames(CssAll_temp), "Css_Unit") #this will fix the forcing to lower of the input and allow it to match what user supplied
  return(CssAll)
}