#include <iostream>
#include "Math/Interpolator.h"
using namespace std;
const int nsamples = 150;
const int nboards = 5;
const int nchannels_per_board = 5;
const int nchannels = nboards*nchannels_per_board;
const int impinged_channel = 2;
const float beam_energy = 100.0;

float temp_amp[nchannels];
float temp_time[nchannels];
float EA_X;
float EA_Y;

float WF_time[nsamples*nchannels];
float WF_val[nsamples*nchannels];
float X[2];
float Y[2];

void generating_template(){
  TCanvas *canvas = new TCanvas("generating template","generating template");
  //TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_35_C3_100.root");
  TFile *file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/vfe_55_C3_100.root");
  TTree *mtree = (TTree*) file->Get("h4");
  mtree->SetBranchAddress("X",X);
  mtree->SetBranchAddress("Y",Y);
  mtree->SetBranchAddress("WF_time",WF_time);
  mtree->SetBranchAddress("WF_val",WF_val);
  TFile *template_recos = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/July2017/H4Analysis/ntuples/template_recos_new_runs_C3_100.root");
  TTree *template_tree = (TTree*) template_recos->Get("template_tree");
  template_tree->SetBranchAddress("temp_amp",temp_amp);
  template_tree->SetBranchAddress("temp_time",temp_time);
  template_tree->SetBranchAddress("EA_X",&EA_X);
  template_tree->SetBranchAddress("EA_Y",&EA_Y);
  TF1 *g = new TF1("g","gaus",0,1000);
  TH1F *waveform = new TH1F("waveform - C3 - 100 GeV - 6*6 mm^2 at center","waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);
  TProfile *mean_waveform = new TProfile("mean waveform - C3 - 100 GeV - 6*6 mm^2 at center","mean waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nsamples,-0.125,937.375);
  
  const Long64_t nentries = mtree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0,ientry;
  float dummy_amp[nentries];
  float dummy_time[nentries];
  bool dummy_constraint[nentries];
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = template_tree->LoadTree(jentry);
    nb = template_tree->GetEntry(jentry);
    dummy_amp[jentry] = temp_amp[impinged_channel];
    dummy_time[jentry] = temp_time[impinged_channel];
    if (TMath::Abs(EA_X)<2&&TMath::Abs(EA_Y)<2&&temp_amp[1]>35&&temp_amp[3]>35&&temp_amp[12]>35&&temp_amp[17]>35&&dummy_amp[jentry] > 3000.0*beam_energy/100.0) dummy_constraint[jentry] = true;
    else dummy_constraint[jentry] = false;
  }
  
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);   nbytes += nb;
    //if (TMath::Abs(X[0]+6)<3&&TMath::Abs(Y[0]-7)<3){                                                 //To be tuned
    if (dummy_constraint[jentry] == true){                                                 //To be tuned
      waveform->Reset();
      for (int i = 0;i < nsamples;i++) waveform->SetBinContent(i+1,WF_val[i+impinged_channel*nsamples]);
      for (int i = 0;i < nsamples;i++) mean_waveform->Fill(WF_time[i+impinged_channel*nsamples]-dummy_time[jentry],WF_val[i+impinged_channel*nsamples]/dummy_amp[jentry]); //comment this line and uncomment the rest of lines in the loop if it's the first time you are generating templates and they haven't been used for fits (basically there's no temp_reco file)
      /*if (waveform->GetMaximum() < 2000.0*beam_energy/100.0) continue;                          //vice versa if otherwise
      g->SetParameters(1.2*waveform->GetMaximum(),6.25*(waveform->GetMaximumBin()-1),6);
      g->SetParLimits(0,1500.0*beam_energy/100.0,10000);
      g->SetParLimits(1,240,320);
      waveform->Fit("g","Q");
      waveform->Fit("g","Q","",-18+g->GetParameter(1),18+g->GetParameter(1));
      if (g->GetParameter(1) < 260 || g->GetParameter(1) > 280) continue;
      if (g->GetParameter(0) < 2300*beam_energy/100.0) continue;
      for (int i = 0;i < nsamples;i++) mean_waveform->Fill(270+WF_time[i+impinged_channel*nsamples]-g->GetParameter(1),WF_val[i+impinged_channel*nsamples]);*/
    }
  }
  
  TAxis *xaxis = mean_waveform->GetXaxis();
  const int nTempBins = 100*nsamples;
  TH1F *interpolated_mean_waveform = new TH1F("interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center","interpolated mean waveform - C3 - 100 GeV - 6*6 mm^2 at center;time(ns)",nTempBins,-0.125,937.375);
  TAxis *xaxis2 = interpolated_mean_waveform->GetXaxis();
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
  //ROOT::Math::Interpolator inter(nsamples+1, ROOT::Math::Interpolation::kCSPLINE);
  ROOT::Math::Interpolator inter(nsamples+2, ROOT::Math::Interpolation::kAKIMA);
  inter.SetData(nsamples+2,x,y);
  for (int i = 0;i < nTempBins;i++) interpolated_mean_waveform->SetBinContent(i+1,inter.Eval(xaxis2->GetBinCenter(i+1)));
  cout << "time of maximum of template:    " << xaxis2->GetBinCenter(interpolated_mean_waveform->GetMaximumBin()) << endl;
  cout << "maximum of template:    " << interpolated_mean_waveform->GetBinContent(interpolated_mean_waveform->GetMaximumBin()) << endl;
  int bin1 = interpolated_mean_waveform->FindFirstBinAbove(interpolated_mean_waveform->GetMaximum()/2);
  int bin2 = interpolated_mean_waveform->FindLastBinAbove(interpolated_mean_waveform->GetMaximum()/2);
  cout << "FWHM of template:    " << interpolated_mean_waveform->GetBinCenter(bin2) - interpolated_mean_waveform->GetBinCenter(bin1) << endl;
  mean_waveform->Draw();
  interpolated_mean_waveform->SetLineColor(4);
  interpolated_mean_waveform->Draw("same");
  /*TFile *template_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/template_July2017_C3_100.root","recreate");
  template_file->cd();
  mean_waveform->Write();
  template_file->Close();*/
}