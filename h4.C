#define h4_cxx
#include "h4.h"
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

using namespace std;

void format_h(TH1F* h, int linecolor){
  h->SetLineWidth(3);
  h->SetLineColor(linecolor);
}

TFile *center = new TFile("ntuples/analysis_4684.root");
TTree * main = (TTree*) center->Get( "h4" );
TCanvas *c1 = new TCanvas("c1","c1");
TH1F *original_sum = new TH1F("original_sum", "original_sum", 1000,8000,130000);
TH1F *weighted_sum = new TH1F("weighted_sum", "weighted_sum", 1000,8000,130000);
TH1F *apd1 = new TH1F("apd1", "apd1", 1000,8000,130000);
TH1F *apd2 = new TH1F("apd2", "apd2", 1000,8000,130000);
TH1F *apd3 = new TH1F("apd3", "apd3", 1000,8000,130000);
TH1F *apd4 = new TH1F("apd4", "apd4", 1000,8000,130000);
TH1F *weighted_apd1 = new TH1F("weighted_apd1", "weighted_apd1", 1000,8000,130000);
TH1F *weighted_apd2 = new TH1F("weighted_apd2", "weighted_apd2", 1000,8000,130000);
TH1F *weighted_apd3 = new TH1F("weighted_apd3", "weighted_apd3", 1000,8000,130000);
TH1F *weighted_apd4 = new TH1F("weighted_apd4", "weighted_apd4", 1000,8000,130000);
Double_t mean_apd[4] = {0};
Double_t weights[4] = {0};

void h4::Loop()
{
//   In a ROOT session, you can do:
//      root> .L h4.C
//      root> h4 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3){
        original_sum->Fill(charge_sig[xtal4apd_1]+charge_sig[xtal4apd_2]+charge_sig[xtal4apd_3]+charge_sig[xtal4apd_4]);
        apd1->Fill(charge_sig[xtal4apd_1]);
        apd2->Fill(charge_sig[xtal4apd_2]);
        apd3->Fill(charge_sig[xtal4apd_3]);
        apd4->Fill(charge_sig[xtal4apd_4]);
        
      }
      // if (Cut(ientry) < 0) continue;
   }
   
   original_sum->Fit("gaus","","",8000,130000);
   TF1 *fit_original_sum = original_sum->GetFunction("gaus");
   Double_t constant_original_sum = fit_original_sum->GetParameter(0);
   Double_t mean_original_sum = fit_original_sum->GetParameter(1);
   Double_t sigma_original_sum = fit_original_sum->GetParameter(2);
   
   apd1->Fit("gaus","","",8000,130000);
   TF1 *fit_apd1 = apd1->GetFunction("gaus");
   Double_t constant_apd1 = fit_apd1->GetParameter(0);
   mean_apd[0] = fit_apd1->GetParameter(1);
   Double_t sigma_apd1 = fit_apd1->GetParameter(2);
   apd1->SetLineWidth(1);
   apd1->SetLineColor(2);
   
   apd2->Fit("gaus","","",8000,130000);
   TF1 *fit_apd2 = apd2->GetFunction("gaus");
   Double_t constant_apd2 = fit_apd2->GetParameter(0);
   mean_apd[1] = fit_apd2->GetParameter(1);
   Double_t sigma_apd2 = fit_apd2->GetParameter(2);
   apd2->SetLineWidth(1);
   apd2->SetLineColor(3);
   
   apd3->Fit("gaus","","",8000,130000);
   TF1 *fit_apd3 = apd3->GetFunction("gaus");
   Double_t constant_apd3 = fit_apd3->GetParameter(0);
   mean_apd[2] = fit_apd3->GetParameter(1);
   Double_t sigma_apd3 = fit_apd3->GetParameter(2);
   apd3->SetLineWidth(1);
   apd3->SetLineColor(4);
   apd3->SetFillStyle(0); //hollow
   
   apd4->Fit("gaus","","",8000,130000);
   TF1 *fit_apd4 = apd4->GetFunction("gaus");
   Double_t constant_apd4 = fit_apd4->GetParameter(0);
   mean_apd[3] = fit_apd4->GetParameter(1);
   Double_t sigma_apd4 = fit_apd4->GetParameter(2);
   apd4->SetLineWidth(1);
   apd4->SetLineColor(6);
   apd4->SetFillStyle(3004);
   apd4->SetFillColorAlpha(kPink,0.5);
   
   c1->Clear();
   
   TF1 *crys = new TF1("crys","crystalball",8000,300000);
   apd1->Draw();
   crys->SetParameters(40,mean_apd[0],0.9*sigma_apd1,1,3);
   apd1->Fit("crys","SAME");
   mean_apd[0] = crys->GetParameter(1);
   
   apd2->Draw("SAME");
   crys->SetParameters(40,mean_apd[1],0.9*sigma_apd2,1,3);
   apd2->Fit("crys","SAME");
   mean_apd[1] = crys->GetParameter(1);
   
   apd3->Draw("SAME");
   crys->SetParameters(40,mean_apd[2],0.9*sigma_apd3,1,3);
   apd3->Fit("crys","SAME");
   mean_apd[2] = crys->GetParameter(1);
   
   apd4->Draw("SAME");
   crys->SetParameters(40,mean_apd[3],0.9*sigma_apd4,1,3);
   apd4->Fit("crys","SAME");
   mean_apd[3] = crys->GetParameter(1);
   
   Double_t weights[4] = {0};
   for (int i = 0;i<4;i++){
     weights[i] = mean_apd[i] / (mean_apd[0]+mean_apd[1]+mean_apd[2]+mean_apd[3]);
   }
   
   for (int i = 0;i<4;i++){
     weights[i] = (mean_apd[0]+mean_apd[1]+mean_apd[2]+mean_apd[3])/(4.0*mean_apd[i]);
   }
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if (nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3){
       weighted_sum->Fill(weights[0]*charge_sig[xtal4apd_1]+weights[1]*charge_sig[xtal4apd_2]+weights[2]*charge_sig[xtal4apd_3]+weights[3]*charge_sig[xtal4apd_4]);
       weighted_apd1->Fill(4*weights[0]*charge_sig[xtal4apd_1]);
       weighted_apd2->Fill(4*weights[1]*charge_sig[xtal4apd_2]);
       weighted_apd3->Fill(4*weights[2]*charge_sig[xtal4apd_3]);
       weighted_apd4->Fill(4*weights[3]*charge_sig[xtal4apd_4]);
     }
   }
   /*
   weighted_sum->SetLineColor(3);
   original_sum->SetLineWidth(1);
   weighted_sum->SetFillStyle(3004);
   weighted_sum->SetFillColorAlpha(kGreen,0.5);
   
   original_sum->Draw();
   crys->SetParameters(40,mean_original_sum,0.9*sigma_original_sum,1,3);
   original_sum->Fit("crys","SAME");
   
   weighted_sum->Draw("SAME");
   crys->SetParameters(40,mean_original_sum,0.9*sigma_original_sum,1,3);
   crys->SetLineColor(3);
   weighted_sum->Fit("crys","SAME");*/
   
   weighted_apd1->Draw();
   weighted_apd2->Draw("SAME");
   weighted_apd3->Draw("SAME");
   weighted_apd4->Draw("SAME");
   weighted_sum->Draw("SAME");
   
}
