#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
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
int nbins = 40;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
double fit_function(double *v,double *par)
{
  int n = int(v[0]*5);
  int WLS_shift = int(par[3]);
  int fscint_shift = int(par[4]);
  return par[0]*convoluted_pulse_fscint->GetBinContent(n+1+fscint_shift) + par[1]*convoluted_pulse_WLS->GetBinContent(n+1+WLS_shift) + par[2];
}
void fscint_contribution(){
  TCanvas *canvas = new TCanvas("fscint_contribution","fscint_contribution");
  TFile *final = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root");
  TTree *ftree = (TTree*) final->Get("h4");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("convoluted pulse fiber scintillation");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("convoluted pulse WLS");
  TH2F *apd1_time_val = new TH2F("apd1_time_val","apd1_time_val",1024,-0.1,204.7,124000,-44000,80000);
  TH2F *apd2_time_val = new TH2F("apd2_time_val","apd2_time_val",1024,-0.1,204.7,124000,-44000,80000);
  TH2F *apd3_time_val = new TH2F("apd3_time_val","apd3_time_val",1024,-0.1,204.7,124000,-44000,80000);
  TH2F *apd4_time_val = new TH2F("apd4_time_val","apd4_time_val",1024,-0.1,204.7,124000,-44000,80000);
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
  for (Long64_t jentry=0; jentry<nentries;jentry++){
    Long64_t ientry = ftree->LoadTree(jentry);
    nb = ftree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(x[0]-207.0)<0.5&&TMath::Abs(y[0]-294.5)<0.5)
    {
      for (int sample = 0;sample < 1009;sample++){
        apd1_time_val->Fill(WF_time6[sample]-28.6,20*WF_val6[sample]);
        apd2_time_val->Fill(WF_time7[sample]-28.6,20*WF_val7[sample]);
        apd3_time_val->Fill(WF_time8[sample]-28.6,20*WF_val8[sample]);
        apd4_time_val->Fill(WF_time9[sample]-28.6,20*WF_val9[sample]);
      }
      for (int sample = 866;sample < 1024;sample++){
        apd1_time_val->Fill(sample/5.0,0);//-1000);
        apd2_time_val->Fill(sample/5.0,0);//-1000);
        apd3_time_val->Fill(sample/5.0,0);//-1000);
        apd4_time_val->Fill(sample/5.0,0);//-1000);
      }
    }
  }
  TProfile *apd1_time_val_prof = apd1_time_val->ProfileX("apd1_time_val_prof");
  TProfile *apd2_time_val_prof = apd2_time_val->ProfileX("apd2_time_val_prof");
  TProfile *apd3_time_val_prof = apd3_time_val->ProfileX("apd3_time_val_prof");
  TProfile *apd4_time_val_prof = apd4_time_val->ProfileX("apd4_time_val_prof");
  TH1D *apd1_time_val_1D = apd1_time_val_prof->ProjectionX();
  TH1D *apd2_time_val_1D = apd2_time_val_prof->ProjectionX();
  TH1D *apd3_time_val_1D = apd3_time_val_prof->ProjectionX();
  TH1D *apd4_time_val_1D = apd4_time_val_prof->ProjectionX();
  double WLS_ratio = 0.5; //    WLS/(WLS+fscint)
  double c_fscint = 1;
  double c_WLS = 1;
  TF1 *func = new TF1("fit",fit_function,0,132,2);
  func->SetParameters(1,1,0,0,0);
  func->SetParNames("fscint constant","WLS constant","adding constant","WLS shift","fscint shift");
  func->SetParLimits(0,0,1000);
  func->SetParLimits(1,0,1000);
  //func->Draw();
  apd1_time_val_1D->Fit("fit","","",0,50);
  c_fscint = func->GetParameter(0);
  c_WLS = func->GetParameter(1);
  double adding_constant = func->GetParameter(2);
  double WLS_shift = func->GetParameter(3);
  double fscint_shift = func->GetParameter(4);
  WLS_ratio = c_WLS/(c_WLS+c_fscint);
  cout << c_fscint << "      " << c_WLS << "      " << WLS_ratio << "      " << adding_constant << "      " << WLS_shift << "      " << fscint_shift << endl;
}