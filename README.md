# ICE_IVIVEpipeline
The workflow allows the flexibility to select from three different rat and human PK models: a 1 compartment model that incorporates Monte Carlo simulation to simulate the population variance (1C), 3 compartment model leveraging the  [EPA's httk package](https://github.com/USEPA/CompTox-ExpoCast-httk) (Solve_3comp) and (solve_pbtk). The workflow is to predict the daily equivalent administered dose (EAD, mg/kg/day) that would lead to steady state blood concentration equivalent to the bioactive concentration from in vitro assays and compared to the predicted lowest effective levels (LELs) of in vivo assays, which is user provided


# Required libraries
	library(tidyverse)
	library(deSolve) # Solves for initial value problems of differential equations
	library(doParallel) # for parallelization 
	library(httk) #this is needed for models: solve_3comp, solve_pbtk. The code is compatible with httk_2.0.2
	

# Input files
	ChemicalData_Rnotebook.txt: example data file of the chemical property data needed to run the workflow
	InvitroData_Rnotebook.txt: example data file of the in vitro assay data needed to run the workflow. Note that the first column must be a CASRN and the following columns assay response values in units matching those that the plasma concentration of chemical are in

# Code files
	steadyState.R: the 1 compartment tk model that calculates a steady state plasma concentration of a given chemical
	CalcEAD.R: calculates the estimated administered dose based on a given plasma concentration and invitro tests
	EADboxplot.R: generates boxplots of the EAD values for each chemical based on the in vitro activity concentrations
	ICE_IVIVE.Rmd: the R notebook with the workflow code

