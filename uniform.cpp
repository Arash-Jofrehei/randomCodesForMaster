#include <TROOT.h>
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "Riostream.h"
#include "TCanvas.h"
using namespace std;
Float_t         charge_sig[18];
Float_t         x[2];
Float_t         y[2];
float           originalX[2];
float           originalY[2];
Int_t           nFibresOnX[2];
Int_t           nFibresOnY[2];
double          xbar;
double          ybar;
double          Xbar;
double          Ybar;
double          uniformed_calibrated_sum;
int x_index,y_index;
const int nbins = 80;
double c[nbins][nbins];
double n[nbins][nbins];
double calibrated_sum = 0;
double weights[18] = {1197.5/574,2523.75/1512,1197.5/779,0,0,2523.75/5300,0.94649,0.935808,1.17086,0.980063,2523.75/1850,0,0,1197.5/2561,2523.75/1433,1197.5/876,0,0};

void uniform()
{
  TCanvas *canvas = new TCanvas("c","c");
  TH1F *uni_calib = new TH1F("uni_calib","uni_calib",40,80,120);
  TH3F *Ybar_Xbar_uni = new TH3F("Ybar_Xbar_uni","Ybar_Xbar_uni",nbins,197,217,nbins,283.5,303.5,130000,-10000,120000);
  TH3F *Y1_X1_uni = new TH3F("Y1_X1_uni","Ybar_Xbar_uni",nbins,197,217,nbins,283.5,303.5,130000,-100000,1200000);
  TFile *merged_file = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/merged_crystal4apd.root"/*,"UPDATE"*/);
  TTree *mtree = (TTree*) merged_file->Get("h4");
  TTree *uniformed = new TTree("uniformed","uniformed");
  uniformed->Branch("uniformed_calibrated_sum",&uniformed_calibrated_sum,"uniformed_calibrated_sum/D");
  uniformed->Branch("originalX",originalX,"originalX[2]/F");
  uniformed->Branch("originalY",originalY,"originalY[2]/F");
  uniformed->Branch("Xbar",&Xbar,"Xbar/D");
  uniformed->Branch("Ybar",&Ybar,"Ybar/D");
  mtree->SetBranchAddress("charge_sig",charge_sig);
  mtree->SetBranchAddress("x", x);
  mtree->SetBranchAddress("y", y);
  mtree->SetBranchAddress("nFibresOnX", &nFibresOnX);
  mtree->SetBranchAddress("nFibresOnY", &nFibresOnY);
  mtree->SetBranchAddress("xbar",&xbar);
  mtree->SetBranchAddress("ybar",&ybar);
  Long64_t nentries = mtree->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (int i = 0;i < nbins;i++)
    for (int j = 0;j < nbins;j++){
      c[i][j] = 0;
      n[i][j] = 0;
    }
  for (Long64_t jentry=0; jentry<nentries/10;jentry++) {
    Long64_t ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);   nbytes += nb;
    calibrated_sum = 0;
    for (int i = 0;i < 18;i++)
      calibrated_sum += charge_sig[i] * weights[i];
    x_index = int((nbins/20.0)*(xbar - 197));
    y_index = int((nbins/20.0)*(ybar - 283.5));
    if (x_index<0||y_index<0||x_index>(nbins-1)||y_index>(nbins-1))
      continue;
    n[x_index][y_index] += 1;
    c[x_index][y_index] += calibrated_sum/100000.0;
  }
  for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      if (n[i][j] != 0)
      {
        c[i][j] /= n[i][j];
        if (c[i][j] != 0)
          c[i][j] = 0.001/c[i][j];
      }
    }
  }
  nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries/10;jentry++){
    Long64_t ientry = mtree->LoadTree(jentry);
    nb = mtree->GetEntry(jentry);   nbytes += nb;
    originalX[0] = x[0];
    originalY[0] = y[0];
    originalX[1] = x[1];
    originalY[1] = y[1];
    Xbar = xbar;
    Ybar = ybar;
    calibrated_sum = 0;
    x_index = int((nbins/20.0)*(xbar - 197));
    y_index = int((nbins/20.0)*(ybar - 283.5));
    if (x_index<0||y_index<0||x_index>(nbins-1)||y_index>(nbins-1)){
      uniformed_calibrated_sum = 0;
      //uniformed->Fill();
      continue;
    }
    for (int i = 0;i < 18;i++)
      calibrated_sum += charge_sig[i] * weights[i];
    uniformed_calibrated_sum = calibrated_sum * c[x_index][y_index];
    uni_calib->Fill(uniformed_calibrated_sum);
    Ybar_Xbar_uni->Fill(Xbar,Ybar,uniformed_calibrated_sum);
    if (x[1] > -500 && y[1] > -500){
      Y1_X1_uni->Fill(x[1],y[1],uniformed_calibrated_sum);
    }
    //uniformed->Fill();
  }
  TProfile2D *Ybar_Xbar_uni_prof = Ybar_Xbar_uni->Project3DProfile("yx");
  TProfile2D *Y1_X1_uni_prof = Y1_X1_uni->Project3DProfile("yx");
  Ybar_Xbar_uni_prof->Draw("colz");
  /*for (int i = 0;i < nbins;i++){
    for (int j = 0;j < nbins;j++){
      //if (Ybar_Xbar_uni_prof->GetBinContent(i,j) != 0)
        cout << Ybar_Xbar_uni_prof->GetBinContent(i+1,j+1) << " ";
    }
    cout << endl;
  }*/
  Y1_X1_uni_prof->Draw("colz");
  //uni_calib->Draw();
  //uniformed->Write();
}