#include <string>

using namespace std;

void xtalk(){
   
   // ---- ---- ---- ---- Retrieving files and histograms ---- ---- ---- ----
   
   TFile *f_injtype1 = new TFile("runs/pixelalive_depleted_injtype1_33.root");
   TFile *f_injtype5 = new TFile("runs/pixelalive_depleted_injtype5_45.root");
   TFile *f_injtype6 = new TFile("runs/pixelalive_depleted_injtype6_48.root");
   
   float alive_eff = 0.5;
   float coupled_eff = 0.3;
   float uncoupled_eff = 0.2;
   
   string baseDir = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_";
   string shortBaseDir = "D_B(0)_O(0)_H(0)_";
   

   
   for (int ch = 15; ch >= 12; ch--){
     
     string chip = to_string(ch);
     string plotDir = "plots/xtalk/HPK/chip"+chip+"/";
     
     TCanvas *c0_pixelalive1 = (TCanvas*) f_injtype1->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive1 = (TH2F*)c0_pixelalive1->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TCanvas *c0_pixelalive5 = (TCanvas*) f_injtype5->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive5 = (TH2F*)c0_pixelalive5->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TCanvas *c0_pixelalive6 = (TCanvas*) f_injtype6->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive6 = (TH2F*)c0_pixelalive6->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     
     int nColumns = h_pixelalive5->GetXaxis()->GetNbins();
     int nRows = h_pixelalive5->GetYaxis()->GetNbins();
     
     //cout << "nRows: " << nRows << "  nColumns: " << nColumns << endl;
     
     
     
     
     
     
     // ---- ---- ---- ---- analysis ---- ---- ---- ----
     
     vector<int> dead_row,dead_col,suspicious_row,suspicious_col,confirmed_row,confirmed_col;
     TH2F *h_suspicious2D = (TH2F*)h_pixelalive5->Clone("h_suspicious2D");
     h_suspicious2D->SetTitle(("suspicious channels of chip "+chip).c_str());
     TH2F *h_confirmed2D = (TH2F*)h_pixelalive5->Clone("h_confirmed2D");
     h_confirmed2D->SetTitle(("confirmed disconnected channels of chip "+chip).c_str());
     
     
     for (int i=0; i<nRows; i++){
       for (int j=0; j<nColumns; j++){
         bool detectable = true;
         if (h_pixelalive1->GetBinContent(j+1,i+1)<alive_eff){
           dead_row.push_back(i);
           dead_col.push_back(j);
           detectable = false;
         }
         if (detectable){
           if ((i==0 && j%2 == 0) || (i==nRows-1 && j%2 == 1)) detectable = false;
         }
         if (detectable && h_pixelalive1->GetBinContent(j+1,i+1)>alive_eff && h_pixelalive5->GetBinContent(j+1,i+1)<coupled_eff){
           suspicious_row.push_back(i);
           suspicious_col.push_back(j);
           h_suspicious2D->SetBinContent(j+1,i+1,1);
         }
         else h_suspicious2D->SetBinContent(j+1,i+1,0.0001);
         if (detectable && h_pixelalive1->GetBinContent(j+1,i+1)>alive_eff && h_pixelalive5->GetBinContent(j+1,i+1)<coupled_eff && h_pixelalive6->GetBinContent(j+1,i+1)<uncoupled_eff){
           confirmed_row.push_back(i);
           confirmed_col.push_back(j);
           h_confirmed2D->SetBinContent(j+1,i+1,1);
         }
         else h_confirmed2D->SetBinContent(j+1,i+1,0.0001);
       }
     }
     
     cout << "chip " << ch << ":" << endl;
     cout << "    dead:        " << dead_row.size() << endl;
     cout << "    suspicious:  " << suspicious_row.size() << endl;
     cout << "    confirmed:   " << confirmed_row.size() << endl;
     
     
     
     
     
     
     
     
     
     // ---- ---- ---- ---- plotting ---- ---- ---- ----
     
     gStyle->SetOptStat(0);
     
     TCanvas *c_pixelalive1 = new TCanvas(("c_pixelalive1_"+chip).c_str(),("efficiency when injecting in same pixel of chip "+chip).c_str(),800,800);
     h_pixelalive1->Draw("colz");
     //c_pixelalive1->SaveAs((plotDir+"pixelalive_"+chip+".png").c_str());
     
     TCanvas *c_pixelalive5 = new TCanvas(("c_pixelalive5_"+chip).c_str(),("efficiency when injecting in coupled pixel of chip "+chip).c_str(),800,800);
     h_pixelalive5->Draw("colz");
     //c_pixelalive5->SaveAs((plotDir+"eff_coupled_"+chip+".png").c_str());
     
     TCanvas *c_pixelalive6 = new TCanvas(("c_pixelalive6_"+chip).c_str(),("efficiency when injecting in uncoupled pixel of chip "+chip).c_str(),800,800);
     h_pixelalive6->Draw("colz");
     //c_pixelalive6->SaveAs((plotDir+"eff_uncoupled_"+chip+".png").c_str());
     
     TCanvas *c_suspicious2D = new TCanvas(("c_suspicious2D_"+chip).c_str(),("suspicious channels of chip "+chip).c_str(),800,800);
     h_suspicious2D->Draw("colz");
     //c_suspicious2D->SaveAs((plotDir+"suspicious2D_"+chip+".png").c_str());
     
     TCanvas *c_confirmed2D = new TCanvas(("c_confirmed2D_"+chip).c_str(),("confirmed disconnected channels of chip "+chip).c_str(),800,800);
     h_confirmed2D->Draw("colz");
     //c_confirmed2D->SaveAs((plotDir+"confirmed2D_"+chip+".png").c_str());
     
     
     
     
     
   } //loop on chips
  
}