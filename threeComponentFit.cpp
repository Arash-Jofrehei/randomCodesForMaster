#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "string.h"
using namespace std;
int nbins = 20;
float x_start = 197;
float x_end = 217;
float y_start = 284.5;
float y_end = 304.5;
TH1F* apd_plus_electronics;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int shift = int(par[3]);
  return par[0]*apd_plus_electronics->GetBinContent(n+1+shift)+par[1]*convoluted_pulse_fscint->GetBinContent(n+1+shift)+par[2]*convoluted_pulse_WLS->GetBinContent(n+1+shift)+par[4]; //I haven't used par[4] (vertical shift) here. It unstables the fits. For marginal crystals it might be useful since I noticed that due to the noise, the baseline is incorrect for many events.
}
void threeComponentFit(){
  TCanvas *canvas = new TCanvas("threeComponentFit","threeComponentFit");
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TF1 *func = new TF1("fit",fit_function,-100,200,4);
  func->SetParNames("spike on APD","fiber scintillation","WLS","shift","vertical offset");
  func->SetParLimits(0,0,100000000);
  func->SetParLimits(1,0,100000000);
  func->SetParLimits(2,0,100000000);
  func->SetParLimits(3,20,80);
  double spike_ratio = 0;
  double fscint_ratio = 0;
  double WLS_ratio = 0;
  double integrated_charge = 1;
  TH2F *spike_hist2D_apd[4];
  spike_hist2D_apd[0] = new TH2F("spike_hist2D_apd1","nuclear counter effect contribution for APD1",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[1] = new TH2F("spike_hist2D_apd2","nuclear counter effect contribution for APD2",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[2] = new TH2F("spike_hist2D_apd3","nuclear counter effect contribution for APD3",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[3] = new TH2F("spike_hist2D_apd4","nuclear counter effect contribution for APD4",nbins,x_start,x_end,nbins,y_start,y_end);
  TH2F *fscint_hist2D_apd[4];
  fscint_hist2D_apd[0] = new TH2F("fscint_hist2D_apd1","fiber scintillation contribution for APD1",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[1] = new TH2F("fscint_hist2D_apd2","fiber scintillation contribution for APD2",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[2] = new TH2F("fscint_hist2D_apd3","fiber scintillation contribution for APD3",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[3] = new TH2F("fscint_hist2D_apd4","fiber scintillation contribution for APD4",nbins,x_start,x_end,nbins,y_start,y_end);
  TH2F *WLS_hist2D_apd[4];
  WLS_hist2D_apd[0] = new TH2F("WLS_hist2D_apd1","WLS contribution for APD1",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[1] = new TH2F("WLS_hist2D_apd2","WLS contribution for APD2",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[2] = new TH2F("WLS_hist2D_apd3","WLS contribution for APD3",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[3] = new TH2F("WLS_hist2D_apd4","WLS contribution for APD4",nbins,x_start,x_end,nbins,y_start,y_end);
  TProfile *apd[4][nbins][nbins];
  string name;
  for (int a = 0;a < 4;a++){           //reading the waveform profiles for each bin
    for(int i = 0;i < nbins;i++){
      for(int j = 0;j < nbins;j++){
        name = "APD";
        name.append(to_string(a+1));
        name.append("_");
        name.append(to_string(i));
        name.append("_");
        name.append(to_string(j));
        const char *APD_name = name.c_str();
        apd[a][i][j] = (TProfile*) apd_profile_waveform_bins->Get(APD_name);
      }
    }
  }
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      cout << endl << i+1 << "    " << j+1 << endl;
      for (int a = 0;a < 4;a++){
        func->SetParameters(50000,1,1,30,0,0);
        //fitting process depends a lot on the initial shift and more importantly the limits on the shift. Current ones work fine
        for (int k = 0;k < 100;k++){              //works with a lot less than 100 times. Sometimes fails to fit a few bins for less than ~20
          apd[a][i][j]->Fit("fit","Q","",25,60);
        }
        integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
        spike_ratio = func->GetParameter(0)/integrated_charge;
        fscint_ratio = func->GetParameter(1)/integrated_charge;
        WLS_ratio = func->GetParameter(2)/integrated_charge;
        spike_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(0));
        fscint_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(1));
        WLS_hist2D_apd[a]->SetBinContent(i+1,j+1,func->GetParameter(2));
      }
    }
  }
  func->SetParameters(50000,1,1,30,0,0);
  for (int i = 0;i < 100; i++){ // to check the fit for a specific bin
    apd[0][1][1]->SetLineWidth(0);
    apd[0][1][1]->Fit("fit","","",25,60);
  }
  integrated_charge = func->GetParameter(0)+func->GetParameter(1)+func->GetParameter(2);
  spike_ratio = func->GetParameter(0)/integrated_charge;
  fscint_ratio = func->GetParameter(1)/integrated_charge;
  WLS_ratio = func->GetParameter(2)/integrated_charge;
  cout << "spike:      " << spike_ratio*100 << endl;
  cout << "fscint:      " << fscint_ratio*100 << endl;
  cout << "WLS:      " << WLS_ratio*100 << endl;
  
  spike_hist2D_apd[0]->Draw("colz");
  //apd_plus_electronics->Draw();
  //convoluted_pulse_fscint->Draw("same");
  //convoluted_pulse_WLS->Draw("same");
}