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
TH1F *original_sum = new TH1F("original_sum", "original_sum", 1012,10011,150000);
TH1F *weighted_sum = new TH1F("weighted_sum", "weighted_sum", 1012,10011,150000);
TH1F *apd1 = new TH1F("apd1", "apd1", 1012,10011,150000);
TH1F *apd2 = new TH1F("apd2", "apd2", 1012,10011,150000);
TH1F *apd3 = new TH1F("apd3", "apd3", 1012,10011,150000);
TH1F *apd4 = new TH1F("apd4", "apd4", 1012,10011,150000);
TH1F *weighted_apd1 = new TH1F("weighted_apd1", "weighted_apd1", 1012,10011,150000);
TH1F *weighted_apd2 = new TH1F("weighted_apd2", "weighted_apd2", 1012,10011,150000);
TH1F *weighted_apd3 = new TH1F("weighted_apd3", "weighted_apd3", 1012,10011,150000);
TH1F *weighted_apd4 = new TH1F("weighted_apd4", "weighted_apd4", 1012,10011,150000);
TH1F *X0 = new TH1F("X0", "X0", 39,-19,19);
TH1F *Y0 = new TH1F("Y0", "Y0", 39,-19,19);
TH1F *X1 = new TH1F("X1", "X1", 39,-19,19);
TH1F *Y1 = new TH1F("Y1", "Y1", 39,-19,19);
TH2F *weighted_sum_X0 = new TH2F("weighted_sum_X0","weighted_sum_X0",130,-6,6,1012,10011,150000);
TH2F *weighted_sum_Y0 = new TH2F("weighted_sum_Y0","weighted_sum_Y0",130,-6,6,1012,10011,150000);
TH2F *weighted_sum_X1 = new TH2F("weighted_sum_X1","weighted_sum_X1",130,-6,6,1012,10011,150000);
TH2F *weighted_sum_Y1 = new TH2F("weighted_sum_Y1","weighted_sum_Y1",130,-6,6,1012,10011,150000);
TF1 *crys = new TF1("crys","crystalball",10011,150000);
TF1 *g = new TF1("g","gaus",10011,150000);
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
      if (nFibresOnX[0] < 3 && nFibresOnY[0] < 3 && nFibresOnX[1] < 3 && nFibresOnY[1] < 3 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3 && X[1] > -3 && X[1] < 3 && Y[1] > -3 && Y[1] < 3){
        original_sum->Fill(charge_sig[xtal4apd_1]+charge_sig[xtal4apd_2]+charge_sig[xtal4apd_3]+charge_sig[xtal4apd_4]);
        apd1->Fill(charge_sig[xtal4apd_1]);
        apd2->Fill(charge_sig[xtal4apd_2]);
        apd3->Fill(charge_sig[xtal4apd_3]);
        apd4->Fill(charge_sig[xtal4apd_4]);
        X0->Fill(X[0]);
        Y0->Fill(Y[0]);
        X1->Fill(X[1]);
        Y1->Fill(Y[1]);
        
      }
   }
   
   original_sum->Fit("g","","",10011,150000);
   Double_t constant_original_sum = g->GetParameter(0);
   Double_t mean_original_sum = g->GetParameter(1);
   Double_t sigma_original_sum = g->GetParameter(2);
   cout << mean_original_sum << "    " << sigma_original_sum << endl;
   
   apd1->Fit("g","","",10011,150000);
   Double_t constant_apd1 = g->GetParameter(0);
   mean_apd[0] = g->GetParameter(1);
   Double_t sigma_apd1 = g->GetParameter(2);
   apd1->SetLineWidth(1);
   apd1->SetLineColor(2);
   
   apd2->Fit("g","","",10011,150000);
   Double_t constant_apd2 = g->GetParameter(0);
   mean_apd[1] = g->GetParameter(1);
   Double_t sigma_apd2 = g->GetParameter(2);
   apd2->SetLineWidth(1);
   apd2->SetLineColor(3);
   
   apd3->Fit("g","","",10011,150000);
   Double_t constant_apd3 = g->GetParameter(0);
   mean_apd[2] = g->GetParameter(1);
   Double_t sigma_apd3 = g->GetParameter(2);
   apd3->SetLineWidth(1);
   apd3->SetLineColor(4);
   apd3->SetFillStyle(0); //hollow
   
   apd4->Fit("g","","",10011,150000);
   Double_t constant_apd4 = g->GetParameter(0);
   mean_apd[3] = g->GetParameter(1);
   Double_t sigma_apd4 = g->GetParameter(2);
   apd4->SetLineWidth(1);
   apd4->SetLineColor(6);
   apd4->SetFillStyle(3004);
   apd4->SetFillColorAlpha(kPink,0.5);
   
   c1->Clear();
   
   cout << "APD 1 :   " << endl;
   apd1->Draw();
   crys->SetParameters(40,mean_apd[0],0.9*sigma_apd1,1,3);
   //crys->SetParLimits(3,1.17,2);
   apd1->Fit("crys","SAME");
   mean_apd[0] = crys->GetParameter(1);
   sigma_apd1 = crys->GetParameter(2);/*
   g->SetParameter(2,sigma_apd1);
   g->FixParameter(1,mean_apd[0]);
   apd1->Fit(g,"B");
   sigma_apd1 = g->GetParameter(2);*/
   
   cout << "APD 2 :   " << endl;
   apd2->Draw("SAME");
   crys->SetParameters(40,mean_apd[1],0.9*sigma_apd2,1,3);
   //crys->SetParLimits(3,1.17,2);
   apd2->Fit("crys","SAME");
   mean_apd[1] = crys->GetParameter(1);
   sigma_apd2 = crys->GetParameter(2);/*
   g->FixParameter(1,mean_apd[1]);
   g->SetParameter(2,sigma_apd2);
   apd2->Fit(g,"B");
   sigma_apd2 = g->GetParameter(2);*/
   
   cout << "APD 3 :   " << endl;
   apd3->Draw("SAME");
   crys->SetParameters(40,mean_apd[2],0.9*sigma_apd3,1,3);
   //crys->SetParLimits(3,1.17,2);
   apd3->Fit("crys","SAME");
   mean_apd[2] = crys->GetParameter(1);
   sigma_apd3 = crys->GetParameter(2);/*
   g->FixParameter(1,mean_apd[2]);
   g->SetParameter(2,sigma_apd3);
   apd3->Fit(g,"B");
   sigma_apd3 = g->GetParameter(2);*/
   
   cout << "APD 4 :   " << endl;
   apd4->Draw("SAME");
   crys->SetParameters(40,mean_apd[3],0.9*sigma_apd4,1,3);
   //crys->SetParLimits(3,1.17,2);
   apd4->Fit("crys","SAME");
   mean_apd[3] = crys->GetParameter(1);
   sigma_apd4 = crys->GetParameter(2);/*
   g->FixParameter(1,mean_apd[3]);
   g->SetParameter(2,sigma_apd4);
   apd4->Fit(g,"B");
   sigma_apd4 = g->GetParameter(2);*/
   
   Double_t weights[4] = {0};
   for (int i = 0;i<4;i++){
     weights[i] = mean_apd[i] / (mean_apd[0]+mean_apd[1]+mean_apd[2]+mean_apd[3]);
   }
   
   for (int i = 0;i<4;i++){
     weights[i] = (mean_apd[0]+mean_apd[1]+mean_apd[2]+mean_apd[3])/(4.0*mean_apd[i]);
   }
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if (nFibresOnX[0] < 3 && nFibresOnY[0] < 3 && nFibresOnX[1] < 3 && nFibresOnY[1] < 3 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3 && X[1] > -3 && X[1] < 3 && Y[1] > -3 && Y[1] < 3){
       weighted_sum->Fill(weights[0]*charge_sig[xtal4apd_1]+weights[1]*charge_sig[xtal4apd_2]+weights[2]*charge_sig[xtal4apd_3]+weights[3]*charge_sig[xtal4apd_4]);
       weighted_apd1->Fill(4*weights[0]*charge_sig[xtal4apd_1]);
       weighted_apd2->Fill(4*weights[1]*charge_sig[xtal4apd_2]);
       weighted_apd3->Fill(4*weights[2]*charge_sig[xtal4apd_3]);
       weighted_apd4->Fill(4*weights[3]*charge_sig[xtal4apd_4]);
       weighted_sum_X0->Fill(X[0],weights[0]*charge_sig[xtal4apd_1]+weights[1]*charge_sig[xtal4apd_2]+weights[2]*charge_sig[xtal4apd_3]+weights[3]*charge_sig[xtal4apd_4]);
       weighted_sum_Y0->Fill(Y[0],weights[0]*charge_sig[xtal4apd_1]+weights[1]*charge_sig[xtal4apd_2]+weights[2]*charge_sig[xtal4apd_3]+weights[3]*charge_sig[xtal4apd_4]);
       weighted_sum_X1->Fill(X[1],weights[0]*charge_sig[xtal4apd_1]+weights[1]*charge_sig[xtal4apd_2]+weights[2]*charge_sig[xtal4apd_3]+weights[3]*charge_sig[xtal4apd_4]);
       weighted_sum_Y1->Fill(Y[1],weights[0]*charge_sig[xtal4apd_1]+weights[1]*charge_sig[xtal4apd_2]+weights[2]*charge_sig[xtal4apd_3]+weights[3]*charge_sig[xtal4apd_4]);
     }
   }
   
   weighted_sum->SetLineColor(3);
   original_sum->SetLineWidth(1);
   weighted_sum->SetFillStyle(3004);
   weighted_sum->SetFillColorAlpha(kGreen,0.5);
   
   cout << "original_sum :   " << endl;
   original_sum->Draw("SAME");
   crys->SetParameters(40,mean_original_sum,0.9*sigma_original_sum,1,3);
   //crys->SetParLimits(3,1.17,2);
   original_sum->Fit("crys","SAME");
   mean_original_sum = crys->GetParameter(1);
   sigma_original_sum = crys->GetParameter(2);/*
   g->FixParameter(1,mean_original_sum);
   g->SetParameter(2,sigma_original_sum);
   original_sum->Fit(g,"B");
   sigma_original_sum = g->GetParameter(2);*/
   
   cout << "weighted_sum :   " << endl;
   weighted_sum->Draw("SAME");
   crys->SetParameters(40,mean_original_sum,0.9*sigma_original_sum,1,3);
   //crys->SetParLimits(3,1.17,2);
   crys->SetLineColor(3);
   weighted_sum->Fit("crys","SAME");
   Double_t mean_weighted_sum = crys->GetParameter(1);
   Double_t sigma_weighted_sum = crys->GetParameter(2);/*
   g->FixParameter(1,mean_weighted_sum);
   g->SetParameter(2,sigma_weighted_sum);
   weighted_sum->Fit(g,"B");
   sigma_weighted_sum = g->GetParameter(2);*/
   /*
   mean_apd[0] = apd1->GetMean();
   mean_apd[1] = apd2->GetMean();
   mean_apd[2] = apd3->GetMean();
   mean_apd[3] = apd4->GetMean();
   mean_original_sum = original_sum->GetMean();
   mean_weighted_sum = weighted_sum->GetMean();
   sigma_apd1 = apd1->GetRMS();
   sigma_apd2 = apd2->GetRMS();
   sigma_apd3 = apd3->GetRMS();
   sigma_apd4 = apd4->GetRMS();
   sigma_original_sum = original_sum->GetRMS();
   sigma_weighted_sum = weighted_sum->GetRMS();
   */
   Double_t resolution_original_sum = 100*sigma_original_sum / mean_original_sum;
   Double_t resolution_weighted_sum = 100*sigma_weighted_sum / mean_weighted_sum;
   Double_t resolution_apd1 = 100*sigma_apd1 / mean_apd[0];
   Double_t resolution_apd2 = 100*sigma_apd2 / mean_apd[1];
   Double_t resolution_apd3 = 100*sigma_apd3 / mean_apd[2];
   Double_t resolution_apd4 = 100*sigma_apd4 / mean_apd[3];
   
   cout << "original resolution:  " << resolution_original_sum << endl << "weighted resolution:  " << resolution_weighted_sum << endl;
   cout << "apd1 resolution:  " << resolution_apd1 << endl;
   cout << "apd2 resolution:  " << resolution_apd2 << endl;
   cout << "apd3 resolution:  " << resolution_apd3 << endl;
   cout << "apd4 resolution:  " << resolution_apd4 << endl;
   cout << "sum of sigma^2 of APDs:  " << sigma_apd1*sigma_apd1 + sigma_apd2*sigma_apd2 + sigma_apd3*sigma_apd3 + sigma_apd4*sigma_apd4 << endl;
   cout << "sigma^2 of original_sum:  " << sigma_original_sum*sigma_original_sum << endl;
   cout << "sigma^2 of weighted_sum:  " << sigma_weighted_sum*sigma_weighted_sum << endl;
   //cout << weights[0] << "   " << weights[1] << "   "<< weights[2] << "   " << weights[3] << endl;
   /*
   weighted_apd1->Draw();
   weighted_apd2->Draw("SAME");
   weighted_apd3->Draw("SAME");
   weighted_apd4->Draw("SAME");
   weighted_sum->Draw("SAME");
   
   weighted_sum_X0->Draw();
   weighted_sum_Y0->Draw();
   weighted_sum_X1->Draw();
   weighted_sum_Y1->Draw();
   
   TProfile *weighted_sum_X0_prof = weighted_sum_X0->ProfileX("weighted_sum_X0_prof");
   TProfile *weighted_sum_Y0_prof = weighted_sum_Y0->ProfileX("weighted_sum_Y0_prof");
   TProfile *weighted_sum_X1_prof = weighted_sum_X1->ProfileX("weighted_sum_X1_prof");
   TProfile *weighted_sum_Y1_prof = weighted_sum_Y1->ProfileX("weighted_sum_Y1_prof");
   
   weighted_sum_X0_prof->Draw();
   weighted_sum_X0_prof->Fit("gaus","","",-5,5);
   
   weighted_sum_Y0_prof->Draw();
   weighted_sum_Y0_prof->Fit("gaus","","",-5,5);
   
   weighted_sum_X1_prof->Draw();
   weighted_sum_X1_prof->Fit("gaus","","",-5,5);
   
   weighted_sum_Y1_prof->Draw();
   weighted_sum_Y1_prof->Fit("gaus","","",-5,5);
   */
   /*
   Double_t ENERGY[] = {20,50,100,200};
   Double_t RES[] = {5.08869,2.84973,1.98741,1.74835};
   TGraph *RES_ENERGY = new TGraph(4,ENERGY,RES);
   TF1 *res_fit = new TF1("res_fit","[0]/sqrt(x) + [1]",10,210);
   RES_ENERGY->Fit(res_fit);
   RES_ENERGY->Draw();*/
}
