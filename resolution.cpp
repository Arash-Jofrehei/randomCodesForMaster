#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TMath.h"
using namespace std;

Float_t         charge_sig[18];
Float_t         x[2];
Float_t         y[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
TH1F *APD[4];
float E = 1;
int APD_start = E*10000;
int APD_end = E*50000;
int APD_bin = 201;
TF1 *crys = new TF1("crys","crystalball",APD_start,4*APD_end);
TH1F *weighted_sum_hist = new TH1F("weighted sum","weighted sum",APD_bin,4*APD_start,4*APD_end);
float mean_APDs;
float APD_weights[4] = {1,1,1,1};
float mean_APD[4];
float sigma_APD[4];

void find_APD_weights(){
  mean_APDs = 0;
  for (int i = 0;i < 4;i++){
    mean_APD[i] = APD[i]->GetMean();
    sigma_APD[i] = 0.8*(APD[i]->GetStdDev());
    crys->SetParameters(1000,mean_APD[i],sigma_APD[i],1.1,2.1);
    crys->FixParameter(3,1.1);
    crys->FixParameter(4,2.1);
    APD[i]->Draw("same");
    APD[i]->Fit("crys");
    mean_APD[i] = crys->GetParameter(1);
    sigma_APD[i] = crys->GetParameter(2);
    mean_APDs += mean_APD[i] / 4.0;
  }
  cout << "----------------------------------------------------------------" << endl << "APD weights:    ";
  for (int i = 0;i < 4;i++){
    APD_weights[i] = mean_APDs/mean_APD[i];
    cout << APD_weights[i] << "    ";
  }
  cout << endl << "----------------------------------------------------------------" << endl;
}

void resolution(){
  TCanvas *canvas = new TCanvas("resolution","resolution");
  TFile *center_run = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/analysis_4684.root");
  TTree *mtree = (TTree*) center_run->Get("h4");
  mtree->SetBranchAddress("charge_sig",charge_sig);
  mtree->SetBranchAddress("x", x);
  mtree->SetBranchAddress("y", y);
  mtree->SetBranchAddress("nFibresOnX", &nFibresOnX);
  mtree->SetBranchAddress("nFibresOnY", &nFibresOnY);
  APD[0] = new TH1F("APD1","APD1",APD_bin,APD_start,APD_end);
  APD[1] = new TH1F("APD2","APD2",APD_bin,APD_start,APD_end);
  APD[2] = new TH1F("APD3","APD3",APD_bin,APD_start,APD_end);
  APD[3] = new TH1F("APD4","APD4",APD_bin,APD_start,APD_end);
  Long64_t nentries = mtree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);   nbytes += nb;
    if (TMath::Abs(x[1]-207)<1.5&&TMath::Abs(y[1]-294.5)<1.5&&nFibresOnX[1]<3&&nFibresOnY[1]<3){
      for (int i = 0;i < 4;i++){
        APD[i]->Fill(charge_sig[6+i]);
      }
    }
  }
  find_APD_weights();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);   nbytes += nb;
    float weighted_sum = 0;
    if (TMath::Abs(x[1]-207)<2&&TMath::Abs(y[1]-294.5)<2&&nFibresOnX[1]<3&&nFibresOnY[1]<3){
      for (int i = 0;i < 4;i++){
        weighted_sum += APD_weights[i] * charge_sig[6+i];
      }
      weighted_sum_hist->Fill(weighted_sum);
    }
  }
  crys->SetParameters(1000,4*mean_APDs,4*sigma_APD[1],1.1,2.1);
  //crys->ReleaseParameter(3);
  //crys->ReleaseParameter(4);
  crys->FixParameter(3,1.1);
  crys->FixParameter(4,2.1);
  weighted_sum_hist->Draw();
  weighted_sum_hist->Fit("crys","","",215000,240000);
  double weighted_mean = crys->GetParameter(1);
  double weighted_mean_error = crys->GetParError(1);
  double weighted_sigma = crys->GetParameter(2);
  double weighted_sigma_error = crys->GetParError(2);
  double weighted_resolution = 100.0 * weighted_sigma / weighted_mean;
  double weighted_resolution_error = weighted_resolution * TMath::Sqrt(TMath::Power(weighted_mean_error/weighted_mean,2.0)+TMath::Power(weighted_sigma_error/weighted_sigma,2.0));
  cout << "weighted resolution:  " << weighted_resolution << " +- " << weighted_resolution_error << endl;
}