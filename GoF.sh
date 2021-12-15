python setup.py

cd output/sm_cards/LIMITS/

combine -M GoodnessOfFit ggtautau_3prong_1_2015_120.txt --algo=saturated --expectSignal 1 --rMin 0 --rMax 2

combine -M GoodnessOfFit ggtautau_3prong_1_2015_120.txt --algo=saturated -t 1000 -s 12345 --expectSignal 1 --rMin 0 --rMax 2


cd -

#python getDataset.py -c 2
