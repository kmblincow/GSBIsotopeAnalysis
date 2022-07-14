# GSBIsotopeAnalysis
This repository houses the code to perform the analyses in the Blincow et al. (2022) manuscript which details the trophic ecology of Giant Sea Bass using stable isotope and gut content analyses.

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
PPSources.csv (derived from Lit_PP_Isotopes.csv - shown in Table A2)
gsb_TEF.csv

#CSIA Analysis
##(Fig. 5, Fig. A2)
CSIAonly_exploration.R
CSIA_length.csv

#Additional Files
##(Fig. A1)
SummaryPlots_Mar2021.R (summary plots of bulk isotope data)
GutContentData.R (generates summary statistics for stomach content data, relies on gutcontents_split.csv)
