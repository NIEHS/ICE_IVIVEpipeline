# ICE_IVIVEpipeline
The workflow allows the flexibility to select from three different rat and human PK models: a 1 compartment model that incorporates Monte Carlo simulation to simulate the population variance (1C), a 3 compartment model leveraging the  [EPA's httk package](https://github.com/USEPA/CompTox-ExpoCast-httk) (Solve_3comp), and a pbpk model that is tailored for compounds with additional glucuronidation(3CGlu). The workflow is to predict the daily equivalent administered dose (EAD, mg/kg/day) that would lead to steady state blood concentration equivalent to the bioactive concentration from in vitro assays and compared to the predicted lowest effective levels (LELs) of in vivo assays, which is user provided

# Required libraries
	library(plyr) # splitting, combining and applying data
	library(deSolve) # Solves for initial value problems of differential equations
	library(tidyr) # helps tidy data easily
	library(ggplot2) # for creating elegant complex plots
	library(scales) # scaling functions for visualizing
	library(foreach) # for copying functions and libraries to each cluster 
	library(doParallel) # for parallelization 

# Input files
	ChemicalData_rnotebook.txt: example data file of the chemical property data needed to run the workflow
	invitroData_xc.txt: example data file of the in vitro assay data needed to run the workflow. Note that the first column must be a CASRN and the following columns assay response values in units matching those that the plasma concentration of chemical are in

# Code files
	steadyState.R: the 1 compartment tk model that calculates a steady state plasma concentration of a given chemical
	glu_MaxConc.R: 3 compartment model that calculates the maximum plasma concentration of a chemical assuming BPA-like properties for glucoronidation
	CalcEAD.R: calculates the estimated administered dose based on a given plasma concentration and invitro tests
	EADboxplot.R: generates boxplots of the EAD values for each chemical based on the in vitro activity concentrations
	ICE_IVIVE.Rmd: the R notebook with the workflow code

