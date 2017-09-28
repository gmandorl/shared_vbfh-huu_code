#!/bin/bash
# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`  -lMLP -lXMLIO -lTMVA
#g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/include -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/include  -lboost_thread -llwtnn -lboost_system -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2-ikhhed/include/eigen3 

g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  
# g++ ./plotter_vbfzll.C -g -o plot `root-config --cflags --glibs`   -lMLP -lXMLIO -lTMVA  -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/include -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/boost/1.57.0-ikhhed2/include  -lboost_thread -llwtnn -lboost_system -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lwtnn/1.0-ikhhed/lib -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2-ikhhed/include/eigen3
# QCDcorrection=nom
# QCDcorrection=up
# QCDcorrection=down


# JEScorrection=nom
# JEScorrection=up
# JEScorrection=down

QCDcorrectionARRAY=(nom up down nom nom up down)
JEScorrectionARRAY=(nom nom nom up down up down)

for i in $(seq 0 1 4 ); do QCDcorrection=${QCDcorrectionARRAY[$i]}; JEScorrection=${JEScorrectionARRAY[$i]};

echo run_plotter_parallel.sh $QCDcorrection $JEScorrection;
source run_plotter_parallel.sh $QCDcorrection $JEScorrection &


done


#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_madgraph_v25_reskim.root DYJetstoLL_madgraph mu  0 0 nom 0 nom  v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/VBF_HToMuMu_v25_reskim.root VBF_HToMuMu mu  0 0 nom 0 nom v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/TT_v25_reskim.root TT mu  0 0 nom 0 nom  v25 reskim histoFileDir

## ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetsToLL_M-50_v25_reskim.root DYJetsToLL_M mu  0 0 nom 0 nom  v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DY0JetsToLL_M-50_v25_reskim.root DY0JetsToLL_M mu  0 0 nom 0 nom  v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DY1JetsToLL_M-50_v25_reskim.root DY1JetsToLL_M mu  0 0 nom 0 nom  v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DY2JetsToLL_M-50_v25_reskim.root DY2JetsToLL_M mu  0 0 nom 0 nom  v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DY3JetsToLL_M-50_v25_reskim.root DY3JetsToLL_M mu  0 0 nom 0 nom  v25 reskim histoFileDir
#./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DY4JetsToLL_M-50_v25_reskim.root DY4JetsToLL_M mu  0 0 nom 0 nom  v25 reskim histoFileDir






# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_amc_0J_v25_reskim.root DYJetstoLL_amc_0J mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_amc_1J_v25_reskim.root DYJetstoLL_amc_1J mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_amc_2J_v25_reskim.root DYJetstoLL_amc_2J mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# 
# 
# 
# 
# 
# 
# 
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/SingleMuon_reminiaod_v25.root SingleMuon mu  1 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# 
# 
# 
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_tW_top_v25_reskim.root ST_tW_top mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_tW_antitop_v25_reskim.root ST_tW_antitop mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_s_v25_reskim.root ST_s-channel mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_t_top_v25_reskim.root ST_t-channel_top_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ST_t_antitop_v25_reskim.root ST_t-channel_antitop_4f_inclusiveDecays mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# 
# 
# 
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/TT_v25_reskim.root TT mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WW_v25_reskim.root WW mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WZ_v25_reskim.root WZ mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/ZZ_v25_reskim.root ZZ mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/WJetsToLnu_madgraph_v25_reskim.root WJetsToLNu mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;


# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/VBF_HToMuMu_v25_reskim.root VBF_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;
# ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/GluGlu_HToMuMu_v25_reskim.root GluGlu_HToMuMu mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir;




# done




####### ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/EWK_LL_JJ_v25_reskim.root EWK_LL_JJ mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir



# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT100to200_v25_reskim.root DYJetstoLL_HT100_200 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT200to400_v25_reskim.root DYJetstoLL_HT200_400 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT400to600_v25_reskim.root DYJetstoLL_HT400_600 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT600to800_v25_reskim.root DYJetstoLL_HT600_800 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT800to1200_v25_reskim.root DYJetstoLL_HT800_1200 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT1200to2500_v25_reskim.root DYJetstoLL_HT1200_2500 mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir
# # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_HT2500toInf_v25_reskim.root DYJetstoLL_HT2500_Inf mu  0 0 $QCDcorrection 0 $JEScorrection  v25 reskim histoFileDir

# # # # # # # # ./plot /gpfs/ddn/cms/user/mandorli/CMSSW_8_0_24/src/fileSkim/DYJetstoLL_madgraph_v25_reskim.root DYJetstoLL_madgraph mu  0 0 nom 0 nom  v25 reskim histoFileDir






#hadd -f DYJetstoLL_amc_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetstoLL_HT*reskim.root   DYJetstoLL_amc_0J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   DYJetstoLL_amc_1J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root  DYJetstoLL_amc_2J_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root

#hadd -f ST_mu_QCDScale${QCDcorrection}_JES${JEScorrection}_v25_reskim.root   ST_tW_top_mu_${QCDcorrection}_JES${JEScorrection}_v25_reskim.root     ST_tW_antitop_mu_${QCDcorrection}_JES${JEScorrection}_v25_reskim.root     ST_s-channel_mu_${QCDcorrection}_JES${JEScorrection}_v25_reskim.root      ST_t-channel_top_4f_inclusiveDecays_mu_${QCDcorrection}_JES${JEScorrection}_v25_reskim.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_${QCDcorrection}_JES${JEScorrection}_v25_reskim.root


# hadd -f DYJetstoLL_amc_mu_QCDScalenom_JESnom_v25_reskim.root   DYJetstoLL_HT*reskim.root   DYJetstoLL_amc_0J_mu_QCDScalenom_JESnom_v25_reskim.root   DYJetstoLL_amc_1J_mu_QCDScalenom_JESnom_v25_reskim.root  DYJetstoLL_amc_2J_mu_QCDScalenom_JESnom_v25_reskim.root



# hadd -f ST_mu_QCDScalenom_JESnom_v25_reskim.root   ST_tW_top_mu_QCDScalenom_JESnom_v25_reskim.root     ST_tW_antitop_mu_QCDScalenom_JESnom_v25_reskim.root     ST_s-channel_mu_QCDScalenom_JESnom_v25_reskim.root      ST_t-channel_top_4f_inclusiveDecays_mu_QCDScalenom_JESnom_v25_reskim.root     ST_t-channel_antitop_4f_inclusiveDecays_mu_QCDScalenom_JESnom_v25_reskim.root










