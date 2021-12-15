python setup.py 0 0

cd output/sm_cards/LIMITS/

combine -M FitDiagnostics workspace.root --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0 --expectSignal 1 --rMin 0 --rMax 2 --saveShapes --saveOverallShapes --saveWithUncertainties
#combine -M FitDiagnostics workspace.root --rMin 0 --rMax 2

combine -M Significance ggtautau_3prong_1_2015_120.txt --expectSignal 1 --rMin 0 --rMax 2

echo 'copying fitDiagnostics.root ...'
mv fitDiagnostics.root fitDiagnostics/fitDiagnostics.root
scp fitDiagnostics/fitDiagnostics.root ajofrehe@lxplus.cern.ch:/afs/cern.ch/work/a/ajofrehe/PhD/gTau/analyzer/skim_ntuples/fitDiagnostics/

source do_parabola.sh

cd -

#python getDataset.py -c 2
