for HF in {1..2}
do
  python setup.py $HF 0
  cd output/sm_cards/LIMITS/
  combine -M FitDiagnostics workspace.root --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0 --expectSignal 1 --rMin 0 --rMax 2 --saveShapes --saveOverallShapes --saveWithUncertainties
  mv fitDiagnostics.root fitDiagnostics/fitDiagnostics$HF.root
  cd -
done

#HF nominal:
python setup.py -1 0
cd output/sm_cards/LIMITS/
combine -M FitDiagnostics workspace.root --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0 --expectSignal 1 --rMin 0 --rMax 2 --saveShapes --saveOverallShapes --saveWithUncertainties
mv fitDiagnostics.root fitDiagnostics/fitDiagnosticsHFnominal.root
cd -

for highNch in {5..8}
do
  python setup.py 0 $highNch
  cd output/sm_cards/LIMITS/
  combine -M FitDiagnostics workspace.root --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0 --expectSignal 1 --rMin 0 --rMax 2 --saveShapes --saveOverallShapes --saveWithUncertainties
  mv fitDiagnostics.root fitDiagnostics/fitDiagnostics$highNch.root
  cd -
done

#Nch nominal:
python setup.py 0 -1
cd output/sm_cards/LIMITS/
combine -M FitDiagnostics workspace.root --robustFit=1 --setRobustFitAlgo=Minuit2 --setRobustFitStrategy=0 --setRobustFitTolerance=0.2 --X-rtd MINIMIZER_analytic --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND  --cminPreScan --cminDefaultMinimizerStrategy 0 --expectSignal 1 --rMin 0 --rMax 2 --saveShapes --saveOverallShapes --saveWithUncertainties
mv fitDiagnostics.root fitDiagnostics/fitDiagnosticsNCHnominal.root
cd -

scp output/sm_cards/LIMITS/fitDiagnostics/fitDiagnostics* ajofrehe@lxplus.cern.ch:/afs/cern.ch/work/a/ajofrehe/PhD/gTau/analyzer/skim_ntuples/fitDiagnostics/
