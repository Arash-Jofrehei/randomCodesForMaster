#include <iostream>
#include "Math/Interpolator.h"

using namespace std;
const int nsamples = 150;
const int nboards = 3;
const int nchannels_per_board = 5;
const int nchannels = nboards*nchannels_per_board;
string Nchannels = to_string(nchannels);

float WF_time[nsamples*nchannels];
float WF_val[nsamples*nchannels];
float hodoX[2];
float hodoY[2];
float WCX[1];
float WCY[1];
float hodo_X[2];
float hodo_Y[2];
float WC_X;
float WC_Y;
float temp_amp[nchannels];
float temp_time[nchannels];
float max_time[nchannels];
float EA_X;
float EA_Y;
float digiMax[nchannels];
unsigned int run;
unsigned int spill;
unsigned int event;
unsigned int Run;
unsigned int Spill;
unsigned int Event;

double x_of_maximum = 267.219;
double maximum_of_template = 3925.89;

ROOT::Math::Interpolator inter(nsamples+2, ROOT::Math::Interpolation::kAKIMA);
double fit_function(double *v,double *par)
{
  //return par[0]*inter.Eval(par[2]*(v[0]-263.1)+263.1-par[1]-6)/3308.7;
  //return par[0]*inter.Eval(par[2]*(v[0]-281.8)+281.8-par[1])/3089;
  return par[0]*inter.Eval(par[2]*(v[0]-x_of_maximum-par[1])+x_of_maximum)/maximum_of_template;
}

void template_fit(){
  TCanvas *canvas = new TCanvas("template fit","template fit");
  TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_35_C3_100.root");
  TTree *h4tree = (TTree*) file->Get("h4");
  h4tree->SetBranchAddress("WF_time",WF_time);
  h4tree->SetBranchAddress("WF_val",WF_val);
  h4tree->SetBranchAddress("run",&run);
  h4tree->SetBranchAddress("spill",&spill);
  h4tree->SetBranchAddress("event",&event);
  TTree *hodoTree = (TTree*) file->Get("hodo");
  hodoTree->SetBranchAddress("X",hodoX);
  hodoTree->SetBranchAddress("Y",hodoY);
  TTree *wireTree = (TTree*) file->Get("wire");
  wireTree->SetBranchAddress("X",WCX);
  wireTree->SetBranchAddress("Y",WCY);
  TFile *template_recos = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/template_recos_C3_100.root","recreate");
  template_recos->cd();
  TTree *template_tree = new TTree("template_tree","template_tree");
  template_tree->Branch("hodo_X",&hodo_X,"hodo_X[2]/F");
  template_tree->Branch("hodo_Y",&hodo_Y,"hodo_Y[2]/F");
  template_tree->Branch("WC_X",&WC_X,"WC_X/F");
  template_tree->Branch("WC_Y",&WC_Y,"WC_Y/F");
  template_tree->Branch("temp_amp",&temp_amp,("temp_amp["+Nchannels+"]/F").c_str());
  template_tree->Branch("temp_time",&temp_time,("temp_time["+Nchannels+"]/F").c_str());
  template_tree->Branch("max_time",&max_time,("max_time["+Nchannels+"]/F").c_str());
  template_tree->Branch("digiMax",&digiMax,("digiMax["+Nchannels+"]/F").c_str());
  template_tree->Branch("EA_X",&EA_X,"EA_X/F");
  template_tree->Branch("EA_Y",&EA_Y,"EA_Y/F");
  template_tree->Branch("run",&Run,"run/i");
  template_tree->Branch("spill",&Spill,"spill/i");
  template_tree->Branch("event",&Event,"event/i");
  TH1F *amp = new TH1F("amplitude obtained by templates","amplitude obtained by templates",100,0,4600);
  TH1F *template_time = new TH1F("time obtained by templates","time obtained by templates",120,-30,30);
  TH1F *template_shrink = new TH1F("shrink of templates","shrink of templates",60,0.95,1.22);
  TH1F *waveform = new TH1F("waveform - C3 - 100 GeV - 6*6 mm^2 at center","waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);
  TFile *template_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/template_July2017_C3_100.root");
  TProfile *mean_waveform = (TProfile*) template_file->Get("mean waveform - C3 - 100 GeV - 6*6 mm^2 at center");
  TProfile *amp_x = new TProfile("template amplitude vs X","template amplitude vs X;X(mm)",30,-15,15);
  TProfile *amp_y = new TProfile("template amplitude vs Y","template amplitude vs Y;Y(mm)",30,-15,15);
  TAxis *xaxis = mean_waveform->GetXaxis();
  double x[nsamples+2],y[nsamples+2];
  x[0] = xaxis->GetBinCenter(1) - 6.25;
  y[0] = 0;
  for (int i = 0;i < nsamples;i++){
    //x[i] = 6.25 * i;
    x[i+1] = xaxis->GetBinCenter(i+1);
    y[i+1] = mean_waveform->GetBinContent(i+1);
  }
  x[nsamples+1] = x[nsamples]+6.25;
  y[nsamples+1] = 0;
  inter.SetData(nsamples+1,x,y);
  TF1 *func = new TF1("fit",fit_function,240,340,3);
  func->SetParNames("amplitude","time","shrink ratio");
  func->SetLineColor(2);
  func->SetParameters(3500,0,1);
  func->SetParLimits(0,30,20000);
  func->SetParLimits(1,-30,30);
  func->SetParLimits(2,0.95,1.2);
  
  //float channel_peak[nchannels] = {2850,3740,3620,4044,3940,3320,3400,4010,4590,4530,3510,3970,3520,3800,2980,4310,4260,3740,3790,3750,4160,4020,3650,4330,4400};
  //float channel_x[nchannels] = {0,0,0,0,0,2,2,2,2,2,1,1,1,1,1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2};
  //float channel_y[nchannels] = {2,1,0,-1,-2,2,1,0,-1,-2,-2,-1,0,1,2,-2,-1,0,1,2,2,1,0,-1,-2};
  float channel_peak[nchannels] = {4550,4500,3950,3950,3840,3400,4180,3750,4010,3180,2800,3940,3830,4300,4150};
  float channel_x[nchannels] = {-1,-1,-1,-1,-1,1,1,1,1,1,0,0,0,0,0};
  float channel_y[nchannels] = {-2,-1,0,1,2,-2,-1,0,1,2,2,1,0,-1,-2};
  float energy_sum;
  float position_weight[nchannels];
  float position_weight_sum;
  
  const Long64_t nentries = h4tree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0,ientry;
  float WC_dummy_X[nentries];
  float WC_dummy_Y[nentries];
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = wireTree->LoadTree(jentry);
    nb = wireTree->GetEntry(jentry);
    WC_dummy_X[jentry] = WCX[0];
    WC_dummy_Y[jentry] = WCY[0];
  }
  
  float hodo_dummy_X0[nentries];
  float hodo_dummy_Y0[nentries];
  float hodo_dummy_X1[nentries];
  float hodo_dummy_Y1[nentries];
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = hodoTree->LoadTree(jentry);
    nb = hodoTree->GetEntry(jentry);
    hodo_dummy_X0[jentry] = hodoX[0];
    hodo_dummy_Y0[jentry] = hodoY[0];
    hodo_dummy_X1[jentry] = hodoX[1];
    hodo_dummy_Y1[jentry] = hodoY[1];
  }
  
  const int nTempBins = 100*nsamples;
  TH1F *interpolated_mean_waveform = new TH1F("interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center","interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nTempBins,-0.125,937.375);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = h4tree->LoadTree(jentry);
    nb = h4tree->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << jentry << " / " << nentries << endl;
    EA_X = -10000;
    EA_Y = -10000;
    WC_X = WC_dummy_X[jentry];
    WC_Y = WC_dummy_Y[jentry];
    hodo_X[0] = hodo_dummy_X0[jentry];
    hodo_Y[0] = hodo_dummy_Y0[jentry];
    hodo_X[1] = hodo_dummy_X1[jentry];
    hodo_Y[1] = hodo_dummy_Y1[jentry];
    Run = run;
    Spill = spill;
    Event = event;
    for(int channel = 0;channel < nchannels;channel++){
      for (int i = 0;i < nsamples;i++) waveform->SetBinContent(i+1,WF_val[i+channel*nsamples]);
      max_time[channel] = ( 6.25 * waveform->GetMaximumBin() ) - 3.25;
      digiMax[channel] = waveform->GetMaximum();
      waveform->Draw();
      //cout << jentry << " / " << nentries << "    " << channel << endl;
      if (waveform->GetMaximum() < 30){
        temp_amp[channel] = 0.001;
        temp_time[channel] = -10000;
      }else{
        func->SetParameters(1.2*waveform->GetMaximum(),0,1);
        //cout << "sth" << endl;
        waveform->Fit("fit","Q","",240,340);
        //for (int i = 4000;i < 7000;i++) interpolated_mean_waveform->SetBinContent(i+1,func->GetParameter(0)*inter.Eval(func->GetParameter(2)*((i*937.5/nTempBins-0.125)-x_of_maximum-func->GetParameter(1))+x_of_maximum)/maximum_of_template);
        //interpolated_mean_waveform->SetLineColor(4);
        //interpolated_mean_waveform->Draw("same");
        //TAxis *xaxis = interpolated_mean_waveform->GetXaxis();
        //cout << xaxis->GetBinCenter(interpolated_mean_waveform->GetMaximumBin()) << endl;
        //waveform->Fit("fit","Q","",269.1+func->GetParameter(1)-20,269.1+func->GetParameter(1)+30);
        if (TMath::Abs(hodo_X[0]+6)<3&&TMath::Abs(hodo_Y[0]-7)<3 && channel == 12) amp->Fill(func->GetParameter(0));
        if (TMath::Abs(hodo_X[0]+6)<3&&TMath::Abs(hodo_Y[0]-7)<3 && channel == 12) template_time->Fill(func->GetParameter(1));
        if (TMath::Abs(hodo_X[0]+6)<3&&TMath::Abs(hodo_Y[0]-7)<3 && channel == 12) template_shrink->Fill(func->GetParameter(2));
        //template_time->Fill(func->GetParameter(1));
        temp_amp[channel] = func->GetParameter(0);
        temp_time[channel] = func->GetParameter(1);
      }
    }
    energy_sum = 0;
    position_weight_sum = 0;
    for (int i = 0;i < nsamples;i++) energy_sum += temp_amp[i]/channel_peak[i];
    for(int channel = 0;channel < nchannels;channel++){
      position_weight[channel] = 3 + TMath::Log10(temp_amp[channel]/(channel_peak[channel]*energy_sum));
      if (position_weight[channel] < 0) position_weight[channel] = 0;
      position_weight_sum += position_weight[channel];
    }
    EA_X = 0;
    for(int channel = 0;channel < nchannels;channel++) EA_X += 22.0 * channel_x[channel] * position_weight[channel]/position_weight_sum;
    EA_Y = 0;
    for(int channel = 0;channel < nchannels;channel++) EA_Y += 22.0 * channel_y[channel] * position_weight[channel]/position_weight_sum;
    template_tree->Fill();
  }
  template_recos->cd();
  template_tree->Write();
  //template_recos->Close();
  amp->Draw();
  TCanvas *canvas2 = new TCanvas("template fit2","template fit2");
  template_time->Draw();
  TCanvas *canvas3 = new TCanvas("template fit3","template fit3");
  amp_x->Draw();
  template_shrink->Draw();
  /*TCanvas *canvas4 = new TCanvas("template fit4","template fit4");
  amp_y->Draw();*/
}