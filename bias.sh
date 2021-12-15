python setup.py

cd output/sm_cards/LIMITS/

#combine -M FitDiagnostics workspace.root --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0 --saveShapes  --expectSignal 1 --rMin 0 --rMax 2 --saveOverallShapes --saveWithUncertainties
#combine -M FitDiagnostics workspace.root --rMin 0 --rMax 2

combine -M FitDiagnostics workspace.root --expectSignal 1 -t 10000 -s 12345 --rMin 0 --rMax 2 --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0


cd -

#python getDataset.py -c 2
