#include <iostream>
using namespace std;
const int nsamples = 150;

float WF_time[nsamples*25];
float WF_val[nsamples*25];
float hodoX[2];
float hodoY[2];
float WCX[1];
float WCY[1];
float hodo_X[2];
float hodo_Y[2];
float WC_X;
float WC_Y;
float temp_amp[25];
float temp_time[25];
float max_time[25];
float EA_X;
float EA_Y;
float digiMax[25];
unsigned int run;
unsigned int spill;
unsigned int event;
unsigned int Run;
unsigned int Spill;
unsigned int Event;

ROOT::Math::Interpolator inter(nsamples+1, ROOT::Math::Interpolation::kCSPLINE);
double fit_function(double *v,double *par)
{
  return par[0]*inter.Eval(par[2]*(v[0]-263.1)+263.1-par[1]-6)/3308.7;
}

void template_fit(){
  TCanvas *canvas = new TCanvas("template fit","template fit");
  TFile *hodoFile = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_55_hodo_B3_100.root");
  TTree *hodoTree = (TTree*) hodoFile->Get("h4");
  hodoTree->SetBranchAddress("X",hodoX);
  hodoTree->SetBranchAddress("Y",hodoY);
  hodoTree->SetBranchAddress("WF_time",WF_time);
  hodoTree->SetBranchAddress("WF_val",WF_val);
  hodoTree->SetBranchAddress("run",&run);
  hodoTree->SetBranchAddress("spill",&spill);
  hodoTree->SetBranchAddress("event",&event);
  TFile *WCfile = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_55_WC_B3_100.root");
  TTree *WCtree = (TTree*) WCfile->Get("h4");
  WCtree->SetBranchAddress("X",WCX);
  WCtree->SetBranchAddress("Y",WCY);
  TFile *template_recos = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/template_recos_B3_100.root","recreate");
  template_recos->cd();
  TTree *template_tree = new TTree("template_tree","template_tree");
  template_tree->Branch("hodo_X",&hodo_X,"hodo_X[2]/F");
  template_tree->Branch("hodo_Y",&hodo_Y,"hodo_Y[2]/F");
  template_tree->Branch("WC_X",&WC_X,"WC_X/F");
  template_tree->Branch("WC_Y",&WC_Y,"WC_Y/F");
  template_tree->Branch("temp_amp",&temp_amp,"temp_amp[25]/F");
  template_tree->Branch("temp_time",&temp_time,"temp_time[25]/F");
  template_tree->Branch("max_time",&max_time,"max_time[25]/F");
  template_tree->Branch("digiMax",&digiMax,"digiMax[25]/F");
  template_tree->Branch("EA_X",&EA_X,"EA_X/F");
  template_tree->Branch("EA_Y",&EA_Y,"EA_Y/F");
  template_tree->Branch("run",&Run,"run/i");
  template_tree->Branch("spill",&Spill,"spill/i");
  template_tree->Branch("event",&Event,"event/i");
  TH1F *amp = new TH1F("amplitude obtained by templates","amplitude obtained by templates",100,0,4600);
  TH1F *template_time = new TH1F("time obtained by templates","time obtained by templates",120,-30,30);
  TH1F *waveform = new TH1F("waveform - B3 - 100 GeV - 1*1 cm^2 at center","waveform - B3 - 100 GeV - 1*1 cm^2 at center;time(ns)",nsamples,-0.125,937.375);
  TFile *template_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/template_July2017_B3_100.root");
  TProfile *mean_waveform = (TProfile*) template_file->Get("mean waveform - B3 - 100 GeV - 1*1 cm^2 at center");
  TProfile *amp_x = new TProfile("template amplitude vs X","template amplitude vs X;X(mm)",30,-15,15);
  TProfile *amp_y = new TProfile("template amplitude vs Y","template amplitude vs Y;Y(mm)",30,-15,15);
  const int nTempBins = 3750; //25*150
  double x[nsamples+1],y[nsamples+1];
  for (int i = 0;i < nsamples+1;i++){
    x[i] = 6.25 * i;
    if (i == nsamples) y[i] = 0;
    else y[i] = mean_waveform->GetBinContent(i+1);
  }
  inter.SetData(nsamples+1,x,y);
  TF1 *func = new TF1("fit",fit_function,200,380,3);
  func->SetParNames("amplitude","time","shrink ratio");
  func->SetLineColor(2);
  func->SetParameters(3500,0,1);
  func->SetParLimits(0,30,20000);
  func->SetParLimits(1,-30,30);
  func->SetParLimits(2,0.95,1.1);
  
  //float channel_peak[25] = {2790,3640,3480,3970,3770,3270,3280,3930,4410,4270,3390,3810,3400,3780,2850,4260,4260,3700,3730,3720,4060,3860,3600,4260,4290};
  float channel_peak[25] = {2850,3740,3620,4044,3940,3320,3400,4010,4590,4530,3510,3970,3520,3800,2980,4310,4260,3740,3790,3750,4160,4020,3650,4330,4400};
  float channel_x[25] = {0,0,0,0,0,2,2,2,2,2,1,1,1,1,1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2};
  float channel_y[25] = {2,1,0,-1,-2,2,1,0,-1,-2,-2,-1,0,1,2,-2,-1,0,1,2,2,1,0,-1,-2};
  float energy_sum;
  float position_weight[25];
  float position_weight_sum;
  
  const Long64_t nentries = hodoTree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0,ientry;
  float WC_dummy_X[nentries];
  float WC_dummy_Y[nentries];
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = WCtree->LoadTree(jentry);
    nb = WCtree->GetEntry(jentry);
    WC_dummy_X[jentry] = WCX[0];
    WC_dummy_Y[jentry] = WCY[0];
  }
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = hodoTree->LoadTree(jentry);
    nb = hodoTree->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << jentry << endl;
    hodo_X[0] = hodoX[0];
    hodo_X[1] = hodoX[1];
    hodo_Y[0] = hodoY[0];
    hodo_Y[1] = hodoY[1];
    EA_X = -10000;
    EA_Y = -10000;
    WC_X = WC_dummy_X[jentry];
    WC_Y = WC_dummy_Y[jentry];
    Run = run;
    Spill = spill;
    Event = event;
    for(int channel = 0;channel < 25;channel++){
      for (int i = 0;i < nsamples;i++) waveform->SetBinContent(i+1,WF_val[i+channel*nsamples]);
      max_time[channel] = ( 6.25 * waveform->GetMaximumBin() ) - 3.25;
      digiMax[channel] = waveform->GetMaximum();
      if (waveform->GetMaximum() < 30){
        temp_amp[channel] = 1;
        temp_time[channel] = -10000;
      }else{
        func->SetParameters(1.2*waveform->GetMaximum(),0,1);
        waveform->Fit("fit","Q","",200,380);
        //waveform->Fit("fit","Q","",269.1+func->GetParameter(1)-20,269.1+func->GetParameter(1)+30);
        //if (jentry == 44 && channel == 17) break;
        if (TMath::Abs(hodo_X[0]) < 3 && TMath::Abs(hodo_Y[0]) < 3 && channel == 17) amp->Fill(func->GetParameter(0));
        //template_time->Fill(func->GetParameter(1));
        temp_amp[channel] = func->GetParameter(0);
        temp_time[channel] = func->GetParameter(1);
      }
    }
    energy_sum = 0;
    position_weight_sum = 0;
    for (int i = 0;i < 25;i++) energy_sum += temp_amp[i]/channel_peak[i];
    for(int channel = 0;channel < 25;channel++){
      position_weight[channel] = 3 + TMath::Log10(temp_amp[channel]/(channel_peak[channel]*energy_sum));
      if (position_weight[channel] < 0) position_weight[channel] = 0;
      position_weight_sum += position_weight[channel];
    }
    EA_X = 0;
    for(int channel = 0;channel < 25;channel++) EA_X += 22.0 * channel_x[channel] * position_weight[channel]/position_weight_sum;
    EA_Y = 0;
    for(int channel = 0;channel < 25;channel++) EA_Y += 22.0 * channel_y[channel] * position_weight[channel]/position_weight_sum;
    template_tree->Fill();
  }
  template_recos->cd();
  template_tree->Write();
  //template_recos->Close();
  amp->Draw();
  /*TCanvas *canvas2 = new TCanvas("template fit2","template fit2");
  template_time->Draw();
  TCanvas *canvas3 = new TCanvas("template fit3","template fit3");
  amp_x->Draw();
  TCanvas *canvas4 = new TCanvas("template fit4","template fit4");
  amp_y->Draw();*/
}