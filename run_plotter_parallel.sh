#!/bin/bash


QCDcorrection=$1
JEScorrection=$2

echo QCDcorrection $QCDcorrection "     " JEScorrection $JEScorrection;


echo "DYJetstoLL_amc_2J"
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetstoLL_amc_0J_v25_reskim.root DYJetstoLL_amc_0J mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetstoLL_amc_1J_v25_reskim.root DYJetstoLL_amc_1J mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetstoLL_amc_2J_v25_reskim.root DYJetstoLL_amc_2J mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
echo "End DYJetstoLL_amc_2J"
 
 

 



 echo "Data"
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/SingleMuon_reminiaod_v25.root SingleMuon mu  1 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
 echo "End Data"


  echo "Single Top"
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/ST_tW_top_v25_reskim.root ST_tW_top mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/ST_tW_antitop_v25_reskim.root ST_tW_antitop mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/ST_s_v25_reskim.root ST_s-channel mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/ST_t_top_v25_reskim.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/ST_t_antitop_v25_reskim.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
  echo "End Single Top"


  echo "TT and Diboson"
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/TT_v25_reskim.root TT mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/WW_v25_reskim.root WW mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/WZ_v25_reskim.root WZ mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/ZZ_v25_reskim.root ZZ mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/WJetsToLnu_madgraph_v25_reskim.root WJetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
 echo "End TT and Diboson"


 echo "Signal"
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/VBF_HToMuMu_v25_reskim.root VBF_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/GluGlu_HToMuMu_v25_reskim.root GluGlu_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
 echo "End Signal"


#######done





#  DY SAMPLE FOR bdt TRAINING

# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetsToLL_M-50_v25_reskim.root DYJetsToLL_M mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DY0JetsToLL_M-50_v25_reskim.root DY0JetsToLL_M mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DY1JetsToLL_M-50_v25_reskim.root DY1JetsToLL_M mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DY2JetsToLL_M-50_v25_reskim.root DY2JetsToLL_M mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DY3JetsToLL_M-50_v25_reskim.root DY3JetsToLL_M mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DY4JetsToLL_M-50_v25_reskim.root DY4JetsToLL_M mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;


# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetsToLL_Pt-0To50_v25_reskim.root    DYJetsToLL_Pt-0To50_amc mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetsToLL_Pt-50To100_v25_reskim.root  DYJetsToLL_Pt-50To100_amc mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetsToLL_Pt-100To250_v25_reskim.root DYJetsToLL_Pt-100To250_amc mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetsToLL_Pt-250To400_v25_reskim.root DYJetsToLL_Pt-250To400_amc mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /afs/cern.ch/work/g/gimandor/public/Hmumu_nTuples/DYJetsToLL_Pt-400To650_v25_reskim.root DYJetsToLL_Pt-400To650_amc mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_Pt-650ToInf_v25_reskim.root DYJetsToLL_Pt-650ToInf_amc mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;




#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT100to200_v25_reskim.root DYJetstoLL_HT100_200 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
 #./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT200to400_v25_reskim.root DYJetstoLL_HT200_400 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
 #./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT400to600_v25_reskim.root DYJetstoLL_HT400_600 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
 #./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT600to800_v25_reskim.root DYJetstoLL_HT600_800 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
 #./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT800to1200_v25_reskim.root DYJetstoLL_HT800_1200 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
 #./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT1200to2500_v25_reskim.root DYJetstoLL_HT1200_2500 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
 #./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT2500toInf_v25_reskim.root DYJetstoLL_HT2500_Inf mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir






