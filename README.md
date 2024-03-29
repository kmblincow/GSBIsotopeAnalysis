# GSBIsotopeAnalysis
This repository houses the code to perform the analyses in the Blincow et al. (2022) manuscript which details the trophic ecology of Giant Sea Bass using stable isotope and gut content analyses.

Blincow KM, Swalethorp R, Ramírez-Valdez A, Semmens BX (2022) Giant appetites: exploring the trophic ecology of California’s largest kelp forest predator, the giant sea bass Stereolepis gigas. Mar Ecol Prog Ser 695:157-171. https://doi.org/10.3354/meps14123

#Sample Map   
##(Fig. 1)  
Files:  
SampleMap.R

#Linear Models for Bulk C/N and Total Length  
##(Table 1, Fig. 2)  
Files:   
JAGS_linearmodel_intonly.R  
LooTableData.csv  
FinalGSBBulk.csv  
  
#SIBER Analysis  
##(Fig. 3)  
Files:  
SIBERAnalysis.R  
FinalGSBBulk.csv  
  
#Mixing Model  
##(Fig.4)  
Files:  
MixingModel.R  
FinalGSBBulk.csv  
PPSources_means.csv (derived from Lit_PP_Isotopes.csv - shown in Table A2)  
gsb_TEF_stCrnge.csv  

#CSIA Analysis  
##(Fig. 5, Fig. A2)  
CSIAonly_exploration.R  
CSIA_length.csv  
  
#Additional Files  
##(Fig. A1)  
SummaryPlots_Mar2021.R (summary plots of bulk isotope data)  
GutContentData.R (generates summary statistics for stomach content data, relies on gutcontents_split.csv)  
