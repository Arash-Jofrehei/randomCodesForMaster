#define test1_cxx
#define test1_cxx
#include "test1.h"
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
void test1(){
  cout<<"Compiled"<<endl;
  }
  
  void format_h(TH1F* h, int linecolor){
     h->SetLineWidth(3);
     h->SetLineColor(linecolor);
 }
  
class Painter{
    public:
    //cout<<"drawing:"<<endl;
    TFile *center = new TFile("ntuples/analysis_4684.root");
    TTree * h4 = (TTree*) center->Get( "h4" );
    TCanvas *c1 = new TCanvas("c1","c1");
    //TH1F *original_sum = new TH1F("original_sum", "original_sum", 100,8000,300000);
    void original(){
    h4->Draw("charge_sig[xtal4apd_1]+charge_sig[xtal4apd_2]+charge_sig[xtal4apd_3]+charge_sig[xtal4apd_4]>>original_sum(100,8000,300000)","nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3");
    original_sum.Fit("gaus","","",8000,300000);
    c1->Modified();
    c1->Update();
    };
    
    /*
    
    void apd1(){
    h4->Draw("4*charge_sig[xtal4apd_1]>>apd1(100,8000,300000)","nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3");
    apd1.Fit("gaus","","",8000,300000);
    apd1->SetLineWidth(1);
    c1->Modified();
    c1->Update();
    };
    
    void apd2(){
    h4->Draw("4*charge_sig[xtal4apd_2]>>apd2(100,8000,300000)","nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3");
    apd2.Fit("gaus","","",8000,300000);
    c1->Modified();
    c1->Update();
    };
    
    void apd3(){
    h4->Draw("4*charge_sig[xtal4apd_3]>>apd3(100,8000,300000)","nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3");
    apd3.Fit("gaus","","",8000,300000);
    apd3->SetFillStyle(0);
    c1->Modified();
    c1->Update();
    };
    
    void apd4(){
    h4->Draw("4*charge_sig[xtal4apd_4]>>apd4(100,8000,300000)","nFibresOnX[0]==2 && nFibresOnY[0]==2 && X[0] > -3 && X[0] < 3 && Y[0] > -3 && Y[0] < 3");
    apd4.Fit("gaus","","",8000,300000);
    apd4->SetLineColor(6);
    apd4->SetFillStyle(3004);
    apd4->SetFillColorAlpha(kPink,0.5);
    c1->Modified();
    c1->Update();
    };*/
    
}painter;
   

/*

#define test1_cxx
#include "test1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
using namespace std;

//TCanvas *c1 = new TCanvas("c1","Charge_sig");
//  c1->SetGrid();
  
//  TH1F *charge_sig = new TH1F("charge_sig","charge_sig",1024,-5,5);
  */
void test1::Loop()
{
//   In a ROOT session, you can do:
//      root> .L test1.C
//      root> test1 t
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
      // if (Cut(ientry) < 0) continue;
   }
   
   
   //charge_sig->Fill(h4->charge_sig[xtal_4apd_1]->GetEntry(jentry));
  
  
}
