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
float WF_time6[1024];
float WF_val6[1024];
float WF_time7[1024];
float WF_val7[1024];
float WF_time8[1024];
float WF_val8[1024];
float WF_time9[1024];
float WF_val9[1024];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TBranch        *b_nFibresOnX;
TBranch        *b_nFibresOnY;
int nbins = 20;
float x_start = 197;
float x_end = 217;
float y_start = 284.5;
float y_end = 304.5;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int WLS_shift = int(par[2]);
  int fscint_shift = int(par[3]);
  return par[0]*(par[1]*convoluted_pulse_fscint->GetBinContent(n+1+fscint_shift)+(1-par[1])*convoluted_pulse_WLS->GetBinContent(n+1+WLS_shift))+par[4];
}
void fscint_contribution(){
  TCanvas *canvas = new TCanvas("fscint_contribution","fscint_contribution");
  //TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins.root");
  //TTree *ftree = (TTree*) final->Get("h4");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  TF1 *func = new TF1("fit",fit_function,0,132,4);
  func->SetParNames("integrated charge","fiber scintillation ratio","WLS shift","fscint shift","added constant");
  //func->SetParLimits(0,0,1000000000);
  func->SetParLimits(1,0,1);
  func->SetParLimits(2,-25,25);
  func->SetParLimits(3,-25,25);
  double fscint_ratio = 0.5;
  double integrated_charge = 1;
  TH2F *ratio_hist2D_apd[4];
  ratio_hist2D_apd[0] = new TH2F("ratio_hist2D_apd1","fiber scintillation contribution for APD1 (percent)",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[1] = new TH2F("ratio_hist2D_apd2","fiber scintillation contribution for APD2 (percent)",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[2] = new TH2F("ratio_hist2D_apd3","fiber scintillation contribution for APD3 (percent)",nbins,x_start,x_end,nbins,y_start,y_end);
  ratio_hist2D_apd[3] = new TH2F("ratio_hist2D_apd4","fiber scintillation contribution for APD4 (percent)",nbins,x_start,x_end,nbins,y_start,y_end);
  bool save = false;
  TProfile *apd[4][nbins][nbins];
  string name;
  for (int a = 0;a < 4;a++){
    for(int i = 0;i < nbins;i++){
      for(int j = 0;j < nbins;j++){
        name = "APD";
        name.append(to_string(a+1));
        name.append("_");
        name.append(to_string(i));
        name.append("_");
        name.append(to_string(j));
        const char *APD_name = name.c_str();
        if (save) apd[a][i][j] = new TProfile(APD_name,APD_name,1024,-0.1,204.7);
        else apd[a][i][j] = (TProfile*) apd_profile_waveform_bins->Get(APD_name);
      }
    }
  }/*
  ftree->SetBranchAddress("WF_time6",WF_time6);
  ftree->SetBranchAddress("WF_val6",WF_val6);
  ftree->SetBranchAddress("WF_time7",WF_time7);
  ftree->SetBranchAddress("WF_val7",WF_val7);
  ftree->SetBranchAddress("WF_time8",WF_time8);
  ftree->SetBranchAddress("WF_val8",WF_val8);
  ftree->SetBranchAddress("WF_time9",WF_time9);
  ftree->SetBranchAddress("WF_val9",WF_val9);
  ftree->SetBranchAddress("x", x);
  ftree->SetBranchAddress("y", y);
  ftree->SetBranchAddress("nFibresOnX", nFibresOnX, &b_nFibresOnX);
  ftree->SetBranchAddress("nFibresOnY", nFibresOnY, &b_nFibresOnY);
  Long64_t nentries = ftree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (save){
    for (Long64_t jentry=0; jentry<nentries;jentry++){
      Long64_t ientry = ftree->LoadTree(jentry);
      nb = ftree->GetEntry(jentry);   nbytes += nb;
      if (x[1]<x_start||x[1]>=x_end||y[1]<y_start||y[1]>=y_end) continue;
      int x_index = int((nbins/(x_end-x_start))*(x[1]-x_start));
      int y_index = int((nbins/(y_end-y_start))*(y[1]-y_start));
      for (int sample = 0;sample < 1009;sample++){
        apd[0][x_index][y_index]->Fill(WF_time6[sample]-28.6,WF_val6[sample]);
        apd[1][x_index][y_index]->Fill(WF_time7[sample]-28.6,WF_val7[sample]);
        apd[2][x_index][y_index]->Fill(WF_time8[sample]-28.6,WF_val8[sample]);
        apd[3][x_index][y_index]->Fill(WF_time9[sample]-28.6,WF_val9[sample]);
      }
      for (int sample = 866;sample < 1024;sample++){
        for (int a = 0;a < 4;a++){
          apd[a][x_index][y_index]->Fill(sample/5.0,0);
        }
      }
    }
  }*/
  //TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins.root","recreate");
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      cout << endl << i+1 << "    " << j+1 << endl;
      for (int a = 0;a < 4;a++){
        //apd_profile_waveform_bins->cd();
        //apd[a][i][j]->Write();
        func->SetParameters(1,0.5,0,0,0);
        apd[a][i][j]->Fit("fit","","",0,90);
        integrated_charge = func->GetParameter(0);
        fscint_ratio = func->GetParameter(1);
        ratio_hist2D_apd[a]->SetBinContent(i+1,j+1,100*fscint_ratio);
      }
    }
  }
  func->SetParameters(1,0.5,0,0,0);
  apd[1][18][1]->Fit("fit","","",0,90);
  //ratio_hist2D_apd[1]->Draw("colz");
  //convoluted_pulse_fscint->Draw("same");
  //convoluted_pulse_WLS->Draw("same");
}