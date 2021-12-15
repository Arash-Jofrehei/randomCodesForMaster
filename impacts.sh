cd output/sm_cards/LIMITS/

combineTool.py -M Impacts -d workspace.root -m 125 --doInitialFit --robustFit 1 --expectSignal 1 --rMin 0 --rMax 2
combineTool.py -M Impacts -d workspace.root -m 125 --doFits --robustFit 1 --expectSignal 1 --rMin 0 --rMax 2
combineTool.py -M Impacts -d workspace.root -m 125 -o impacts.json
plotImpacts.py -i impacts.json -o impacts

echo 'copying impacs ...'
scp impacts.pdf ajofrehe@lxplus.cern.ch:/eos/user/a/ajofrehe/www/gtau/3prong533/singlePlot/
scp impacts.pdf ajofrehe@lxplus.cern.ch:/eos/user/a/ajofrehe/www/gtau/3prong533-postfit/singlePlot/
cd -

#python getDataset.py -c 2
