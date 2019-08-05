####functions
#Name: glu_MaxConc
#Author: ILS (comptox@ils-inc.com)
#Date:26 April, 2019
#Version: 1.0
#License: MIT
##Summary: uses a 3 compartment model (liver, kidney, gut) to calculate the max concentration of a chemical.
#Incorporates addtional gut glucoronidation for BPA-like compounds. Defaults to BPA values for parameters not supplied
#physParam: single row data fram object with the following fields: species, bw, QC.mLmin, QliverC, QkidneyC, VplasmaC, VliverC, VkidneyC, enterocytes
#chemParam: data frame with the following columns:CASRN, ChemicalName, Species,mw, vmaxliver, kmliver, kuptakeC, kGIingC, kmgutg,vmaxgutgc, fgutg, fliverg
#route: route of admistration. one of "oral" or "injection"
#interv= dosing interval (h)
#ndays = total length of exposure (days)
#dose= dose (amount) of substance given, units in mg/kg
#ConcentrationUnit= choice of "uM" or "mg/L"; the output units for plasma substance concentration. Defaults to "uM"

#Library required: deSolve package
#dose= dose given, for oral dosage mg/kg for injection dose in
#Library required: deSolve, parallel (if parallel=TRUE)

library(deSolve)

glu_MaxConc<-function(physParam, chemParam, Species="human", route="iv", ConcentrationUnit="uM", interv=24, ndays=3, dose=1, ...){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # physicology specific parameters
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 library(deSolve)
   physParam <- subset(physParam, tolower(species) == tolower(Species))
  species<-tolower(physParam$species)
  if(species == "rat"){
    GFR <- 0.08  # rat glomerular filtration rate (L/h) # Wetmore et al, 2013
    #clintC <- 0.06555 # converting factor from CLinvitro (ul/min/10^6 cells) to CLintrinsic (L/h) for rat # Wetmore et al, 2013

  } else #default to human values for any other species
  { GFR <- 6.7  # human glomerular filtration rate (L/h) # Wetmore et al, 2012
  #clintC <- 10.2 # converting factor from CLinvitro (ul/min/10^6 cells) to CLintrinsic (L/h) for hum # Wetmore et al, 2012
  # CLintrinsic= CLinvitro*(103*1650*60)/1e6
  }
  bw <- physParam$bw #(kg)|Body weight
  # Blood flow rate #
  QC<- physParam$QC.mLmin*60/1000      #(L/h)|converting from (mL/min) input; Cardiac output according to Davies et al
  # Fraction of blood flows
  QliverC <- physParam$QliverC  #(fraction of QCliver
  QkidneyC <- physParam$QkidneyC  #(fraction of QC)|kidney
  # Fraction of tissue volumes
  VplasmaC <- physParam$VplasmaC  #(fraction of bw)|plasma
  VliverC <- physParam$VliverC   #(fraction of bw)|liver
  VkidneyC <- physParam$VkidneyC  #(fraction of bw)|kidney
  VbodyC <- 1-(VliverC+VplasmaC+VkidneyC) #(fraction of bw)|rest body
  enterocytes <- physParam$enterocytes #(L)|Sum of enterocytes weights in duodenum, jujunum and ileum (Gertz 2011)
  #VbodygC <- VplasmaC                #(fraction of bw)|fractional volume of the distribution for glucuronidaed chemical, set to plasma
  VbodygC <- 1                        #190413(fraction of bw)|fractional volume of the distribution for glucuronidaed chemical, set to whole body

  #Scaled blood flows
  Qliver <- QliverC*QC               #(L/h)|Blood flow to the liver
  Qkidney <- QkidneyC*QC             #(L/h)|Blood flow to the kidney
  Qbody <- QC-Qliver-Qkidney         #(L/h)|Blood flow to the rest body
  #Scaled tissue volumes
  Vliver <- VliverC*bw               #(L)|Volume of the liver
  Vplasma <- VplasmaC*bw             #(L)|Volume of the plasma
  Vkidney <- VkidneyC*bw             #(L)|Volume of the kidney
  Vbody<- VbodyC*bw                  #(L)|Volume of the rest body
  Vbodyg<- VbodygC*bw                #(L)|Volume of the distribution for glucuronidaed chemical
  ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Chemical specific parameters
  #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cmax <- NULL
  cmaxall <- NULL
  EADall <- NULL
  colnames(chemParam)<-tolower(colnames(chemParam))#this is in the event we do a user provided option
  for (m in 1:nrow(chemParam)) {
    fu <- chemParam$fu[m]
      fliverg <- chemParam$fliverg[m]                  #correction factor of glucuronidation in the liver
    if (is.na(fliverg)) {fliverg <- 1}

    #updating to use the vaue for BPA vs Clint coming from OPERA
    vmaxliver <- chemParam$vmaxliver[m]
    kmliver <- chemParam$kmliver[m]

    #Scaled kinetic parameters
    kuptakeC <- chemParam$kuptakec[m] #(1/h/bw^0.25)|oral uptake of chemical from the gut (mainly small intestine) into the liver (optional parameter)
    kGIingC <- chemParam$kgiingc[m]                  #(1/h/bw^0.25)|transport of glucuronidated chemical from enterocytes into serum  (optional parameter)


    if(species == "rat"){
      if (is.na(kuptakeC) || is.null(kuptakeC)) {kuptakeC <- 0.38}
      if (is.na(vmaxliver) || is.null(vmaxliver)) {vmaxliver <- 660012.4089*fliverg} #units of nmol/h, using values from yang et al 2013 converted for units }
      if (is.na(kmliver) || is.null(kmliver)) {kmliver <- 16200} #units of nmol/L,  using values from yang et al 2013
      if (is.na(kGIingC) || is.null(kGIingC)) {kGIingC <- 1} #using values from yang et al 2013
    } else
    { if (is.na(kuptakeC) || is.null(kuptakeC)) {kuptakeC <- 2.1}
      if (is.na(vmaxliver) || is.null(vmaxliver)) {vmaxliver <- 15192576*fliverg}#units of nmol/h, using values from Karrer et 2018 converted for units
      if (is.na(kmliver) || is.null(kmliver)) {kmliver <- 45800} #units of nmol/L, #units of nmol/h, using values from Karrer et 2018
      if (is.na(kGIingC) || is.null(kGIingC)) {kGIingC <- 50} #using values from Karrer et 2018

    }

    kuptake <- kuptakeC/bw^0.25                      #(1/h)|Uptake of parent chemical from gut (mainly small intestine) into the liver
    kGIing <- kGIingC/bw^0.25                        #(1/h)|Uptake of glucuronidaed chemical from gut (mainly small intestine) into serum


    # Scaling of metabolic parameters
    kmgutg <- chemParam$kmgutg[m]
    if (is.na(kmgutg)) {kmgutg <- 58400}

    vmaxgutgc <- chemParam$vmaxgutgc[m]
    if (is.na(vmaxgutgc)) {vmaxgutgc <- 361}
    vmaxgutgCnew <- vmaxgutgc*bw/(bw^0.75)

    kurinechem <- fu*GFR    #(L/h)|Clearance of parent chemical via urine
    kurinechemg <- fu*GFR*6 #(L/h)|Clearance of glucuronidated chemical via urine

    fgutg <- chemParam$fgutg[m]                      #correction factor of glucuronidation in the gut
    if (is.na(fgutg)) {fgutg <- 1}

    vmaxgutg <- vmaxgutgCnew*fgutg*bw^0.75
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #  Dosing Parameters
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++----
    if(route =="oral"){
      oraldose<-dose
      ivdose<-0
    }
    if(route == "iv"){
      oraldose<-0
      ivdose<-dose
    }
    # Oral dosing
    D.o <- oraldose*1E6 #(ng/kg bw/d)|oral dose is equally distributed among the dosings
    dose.O <- D.o/chemParam$mw[m] #(nmol/kg/d)|oral dose
    uptake.O <- bw*dose.O  #(nmol/day)|total amount of uptake per day
    period.O <- 3/60 #(h)|uptake period
    koa <- uptake.O/period.O #(nmol/h)|uptake rate

    # i.v. dosing
    D.iv <- ivdose*1E6 #(ng/kg bw/d)|oral dose is equally distributed among the dosings
    dose.iv <- D.iv/chemParam$mw[m] #(nmol/kg/d)|oral dose
    uptake.iv <- bw*dose.iv  #(nmol/day)|total amount of uptake per day
    period.iv <- 0.015 #(h)|uptake period
    kiv <- uptake.iv/period.iv #(nmol/h)|uptake rate
    ###############################################

    para <- unlist(c(data.frame(QC, Qliver, Qkidney, Qbody, Vliver, Vplasma, Vkidney, Vbody, Vbodyg,  enterocytes, pliver=chemParam$pliver[m], pkidney=chemParam$pkidney[m], pbody=chemParam$pbody[m],kmliver, vmaxliver, kmgutg, kuptake, kGIing, kurinechem, kurinechemg, vmaxgutg, koa, kiv)))
    #For checking the values of parameters, writing them to disk
    # Names <- "Value"
    # write.table(rbind(Names,bw,QC,QliverC,Qliver,QkidneyC,Qkidney,Qbody,VliverC,Vliver,VplasmaC,Vplasma,VkidneyC,Vkidney,VbodyC,Vbody,VbodygC,Vbodyg,fu,chemParam$pliver[m],
    #                   chemParam$pkidney[m],chemParam$pbody[m],enterocytes,kuptakeC,kuptake,kmgutg,fgutg,kmliver, vmaxliver,fliverg,kGIingC,kGIing,vmaxgutgCnew,vmaxgutg,
    #                   kurinechem, kurinechemg,uptake.O,koa,uptake.iv,kiv), file = "test.csv", row.names = T, col.names = FALSE, sep = ",")
    #

    yini <- unlist(c(data.frame(
      Input.O = 0,  #(nmol/h)
      Input.iv= 0,  #(nmol/h)
      AGI = 0,      #(nmol)|Amount of parent chemical in gut (mainly small intestine)
      AAO = 0,      #(nmol)|Amount of parent chemical taken up from gut (mainly small intestine) into serum
      AGImet = 0,   #(nmol)|Amount of glucuronidaed chemical formed in gut (mainly small intestine)
      AGIchemg = 0, #(nmol)|Amount of glucuronidaed chemical in gut (mainly small intestine)
      AGIin = 0,    #(nmol)|Amount of glucuronidaed chemical taken up from gut (mainly small intestine) into serum
      Aplasma = 0,  #(nmol)|Amount of parent chemical in plasma
      Aliver = 0,   #(nmol)|Amount of parent chemical in liver
      Amet_liver = 0, #(nmol)|Amount of glucuronidated chemical in liver
      Akidney=0,    #(nmol)|Amount of chem in kidney
      Abody = 0,    #(nmol)|Amount of parent chemical in rest body tissue
      Aurinechem = 0, #(nmol)|Cumulative amount of parent chemical excreted into urine
      Achemg = 0,   #(nmol)|Amount of glucuronidaed chemical in the system
      Aurinechemg = 0, #(nmol)|Amount of glucuronidaed chemical in the bladder
      Cgut = 0     #(nmol/L)|Concentration of parent chemical in the gut (mainly small intestine)
      # CVliver = 0,  #(nmol/L)|Venous blood concentration of parent chemical leaving the liver
      # CV=0,         #(nmol/L)|Concentration of parent chemical in the venous plasma
      # CA=0,         #(nmol/L)|concentration of parent chemical in the arterial plasma
      # Cchemg=0      #(nmol/L)|Concentration of glucuronidaed chemical in the system
    )))

    dInput.O <- 0
    dInput.iv <- 0

    PBTKmod <- function(times, y, parms){
      with (as.list(c(y, parms)),
            {
              #Oral dosing
              for (i in 1:(ndays*(24/interv)))  {
                t1<- times-interv*(i-1)
                if(t1<=period.O && t1>=0){onoff.O=1} else{onoff.O=0}  # error corrected after introducing t1
                dInput.O <- dInput.O + koa*onoff.O #(nmol/h)
              }
              #i.v. dosing
              for (j in 1:(ndays*(24/interv)))  {
                t2<- times-interv*(j-1)
                if(t2<=period.iv && t2>=0){onoff.iv=1} else{onoff.iv=0}
                dInput.iv <- dInput.iv + kiv*onoff.iv #(nmol/h)
              }

              Cgut <- AGI/enterocytes #(nmol/L)|Concentration of parent chemical in the gut (mainly small intestine)
              RGImet <- vmaxgutg*Cgut/(kmgutg+Cgut)#(nmol/h)|Rate of chemical glucuronidation in the gut
              RAO <- kuptake*AGI#(nmol/h)|Uptake rate of parent chemical from the gut (mainly small intestine) intoliver
              RGI <- dInput.O-RGImet-RAO #(nmol/h)|Rate of parent chemical amount change in the gut(mainly small intestine)
              #Roral <-RAO#(nmol/h)|Rate of parent chemical peroral uptake
              #Amount of glucuronidaed chemical in GI tract
              RGIin <- kGIing*AGIchemg#(nmol/h)|Uptake rate of glucuronidaed chemical from gut (mainly small intestine) into serum
              RGIchemg <- RGImet - RGIin #(nmol/h)|Rate of glucuronidaed chemical amount change in the gut (mainly small intestine)
              ### C's and CV's ###
              Cliver <- Aliver/Vliver #(nmol/L)|Concentration of parent chemical in the liver
              CVliver <- Aliver/(Vliver*chemParam$pliver[m])#(nmol/L)|Venous blood concentration of parent chemical leaving the liver
              Ckidney <- Akidney/Vkidney#(nmol/L)|Concentraitoin of chem in kidney
              CVkidney <- Akidney/(Vkidney*chemParam$pkidney[m])#(nmol/L)|Venous blood concentration of chem leaving the kidney
              Cbody <- Abody/Vbody#(nmol/L)|Concentraitoin of parent chemical in the rest body
              CVbody <- Abody/(Vbody*chemParam$pbody[m])#(nmol/L)|Venous blood concentration of parent chemical leaving the rest body
              CV <- (CVliver*Qliver+CVkidney*Qkidney+CVbody*Qbody)/QC
              #(nmol/L)|Concentration of parent chemical in the venous plasma.
              CA <- Aplasma/Vplasma #(nmol/L)|concentration of parent chemical in the arterial plasma

              #Excretion of parent chemical in urine
              Rurinechem <- kurinechem*CA #(nmol/h)|Rate of parent chemical excreted into the urine

              #Amount of parent chemical in the plasma
              Rplasma <- QC*(CV-CA)+dInput.iv #(nmol/h)|Rate of parent chemical amount change in the plasma.

              #Amount of parent chemical in the liver
              Rmet_liver <- vmaxliver*CVliver/(kmliver+CVliver) #(nmol/h)|Rate of parent chemical metabolism in the liver, CLint(L/h),estimated from CLinvitro or use QPPR, when CVliver << Km

              #Amount of chem in the kidney
              Rkidney <- Qkidney*(CA-CVkidney)-Rurinechem#(nmol/h)|Rate of chem amount change in the kidney
              #Amount of parent chemical in rest body
              Rbody <- Qbody*(CA-CVbody)#(nmol/h)|Rate of parent chemical amount change in rest body
              #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              # Model for glucuronidaed chemical
              #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              #Fate of glucuronidaed chemical formed in the liver- Rmet_liver#(nmol/h)|Taken up into systemic circulation
              #Fate of glucuronidaed chemical formed in GI (mainly SI) - RGIin  #(nmol/h)|Taken up into systemic circulation
              Cchemg <- Achemg/(Vbodyg+1E-34) #(nmol/L)|Concentration of glucuronidaed chemical in the system
              #Urinary excretion of glucuronidaed chemical
              Rurinechemg <- kurinechemg*Cchemg #(nmol/h)|Rate of glucuronidaed chemical amount change in the bladder
              Rchemg <- Rmet_liver + RGIin - Rurinechemg
              #(nmol/h)|Rate of glucuronidaed chemical amount change in the system
              Rliver <- Qliver*(CA-CVliver) + RAO - Rmet_liver #(nmol/h)|Rate of parent chemical amount change in the liver
              dydt <-c(dInput.O, dInput.iv, RGI, RAO, RGImet, RGIchemg, RGIin, Rplasma, Rliver, Rmet_liver, Rkidney, Rbody, Rurinechem, Rchemg, Rurinechemg, Cgut)
              #to troubleshoot...anything that is blank/numeric(0) wont go into dydt object
              #print(c(list(dInput.O, dInput.iv, RGI, RAO, RGImet, RGIchemg, RGIin, Rplasma, Rliver, Rmet_liver, Rkidney, Rbody, Rurinechem, Rchemg, Rurinechemg, Cgut, CVliver, CV, CA, Cchemg)))
              conc <- c(CV=CV,Cchemg=Cchemg,CVliver=CVliver,CA=CA)
              res <- list(dydt, conc)
              return(res)
            })}
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Solve the system of differential equations
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    simu <- seq(0, ndays*24*60, .5)/60 #(h) time # simulation for 3 day, 3*24*60
      v <- ode(y=yini, func=PBTKmod, times=simu, parms=para, method="lsoda")
      y <- as.data.frame(v)

    #plots for checking behavior in the different compartments/components
    # par(mfrow=c(3,4))
    # plot(y$time , y$Input.O, xlab = c('Time'),  col = 1, type = 'l')
    # plot(y$time , y$Input.iv, xlab = c('Time'),  col = 2, type = 'l')
    # plot(y$time , y$AAO, xlab = c('Time'),  col = 2, type = 'l')
    # plot(y$time , y$AGImet, xlab = c('Time'),  col = 3, type = 'l')
    # plot(y$time , y$Aplasma, xlab = c('Time'),  col = 4, type = 'l')
    # plot(y$time , y$Aliver, xlab = c('Time'),  col = 5, type = 'l')
    # plot(y$time , y$Achemg, xlab = c('Time'),  col = 6, type = 'l')
    # plot(y$time , y$AGIin, xlab = c('Time'),  col = 7, type = 'l')
    # plot(y$time , y$Cgut, xlab = c('Time'),  col = 8, type = 'l')
    # plot(y$time , y$CV, xlab = c('Time'),  col = 9, type = 'l')
    # plot(y$time , y$CVliver, xlab = c('Time'),  col = 10, type = 'l')
    # plot(y$time , y$Cchemg, xlab = c('Time'),  col = 11, type = 'l')
    # #plot(y$time , y$CA, xlab = c('Time'),  col = 12, type = 'l')

    ################## generating file containing cmax for all the chemicals ################################

    if (ConcentrationUnit =="mg/L")
    {
      cmax[m] <- max(y$CV/1000)*chemParam$mw[m]*1e-03 #mg/L|maximum blood concentration during the simulation
    }else{
      cmax[m] <- max(y$CV/1000) #uM|maximum blood concentration during the simulation
      ConcentrationUnit ="uM" #this sets it as the default if it doesnt match the alternate option
    }
    cmaxout <- as.data.frame(cbind(chemParam$casrn[m],chemParam$chemicalname[m],route,dose, interv,ndays, cmax[m], ConcentrationUnit), stringsAsFactors=FALSE)
    names(cmaxout)<- c("CASRN", "ChemicalName", "route", "dose, mg/kg", "interv", "ndays", "cmax", "cmax_unit")
    cmaxall <- rbind(cmaxall, cmaxout, stringsAsFactors=FALSE)
  }
  return(cmaxall)
}
