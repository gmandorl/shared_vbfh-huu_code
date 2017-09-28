#!/bin/bash


# QCDcorrection=nom
# QCDcorrection=up
# QCDcorrection=down


# JEScorrection=nom
# JEScorrection=up
# JEScorrection=down


QCDcorrectionARRAY=(nom up down nom nom up down)
JEScorrectionARRAY=(nom nom nom up down up down)

for i in $(seq 0 1 6); do QCDcorrection=${QCDcorrectionARRAY[$i]}; JEScorrection=${JEScorrectionARRAY[$i]};

echo QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection




 hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root
 
 
# hadd -f DYJetstoLL_amc_Pt_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetsToLL_Pt-0To50_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetsToLL_Pt-50To100_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetsToLL_Pt-100To250_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetsToLL_Pt-250To400_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DYJetsToLL_Pt-400To650_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetsToLL_Pt-650ToInf_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root 

# hadd -f DYJetstoLL_amc_Pt_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetsToLL_Pt-0To50_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetsToLL_Pt-50To100_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetsToLL_Pt-100To250_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetsToLL_Pt-250To400_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DYJetsToLL_Pt-400To650_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root 

# hadd -f DYInclusivetoLL_M_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DY0JetsToLL_M_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DY1JetsToLL_M_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DY2JetsToLL_M_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DY3JetsToLL_M_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root DY4JetsToLL_M_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root




 hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   ST_tW_top_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root     ST_tW_antitop_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root     ST_s-channel_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root


done



# hadd -f DYJetstoLL_M_mu_QCDScalenom_JESnom_v25_reskim.root   DY0JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root   DY1JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root   DY2JetstoLL_amc_mu_QCDScalenom_JESnom_v25_reskim.root  DY3JetstoLL_amc_2J_mu_QCDScalenom_JESnom_v25_reskim.root DY4JetstoLL_amc_2J_mu_QCDScalenom_JESnom_v25_reskim.root

#hadd -f DYJetstoLL_M_mu_QCDScalenom_JESnom_v25_reskim.root DY0JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root DY1JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root DY2JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root DY3JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root DY4JetsToLL_M_mu_QCDScalenom_JESnom_v25_reskim.root

# hadd -f DYJetstoLL_amc_mu_QCDScalenom_JESnom_v25_reskim.root   DYJetstoLL_HT*reskim.root   DYJetstoLL_amc_0J_mu_QCDScalenom_JESnom_v25_reskim.root   DYJetstoLL_amc_1J_mu_QCDScalenom_JESnom_v25_reskim.root  DYJetstoLL_amc_2J_mu_QCDScalenom_JESnom_v25_reskim.root


# hadd -f ST_mu_QCDScalenom_JESnom_v25_reskim.root   ST_tW_top_mu_QCDScalenom_JESnom_v25_reskim.root     ST_tW_antitop_mu_QCDScalenom_JESnom_v25_reskim.root     ST_s-channel_mu_QCDScalenom_JESnom_v25_reskim.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScalenom_JESnom_v25_reskim.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScalenom_JESnom_v25_reskim.root










