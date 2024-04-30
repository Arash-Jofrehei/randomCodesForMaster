#include <string>
#define alive

using namespace std;

void thr_noise_diff(){
   TFile *f_pixelalive = new TFile("runs/Run000046_PixelAlive.root");
   TFile *f_noHV = new TFile("runs/Run000045_SCurve.root");
   TFile *f_withHV = new TFile("runs/Run000047_SCurve.root");
   
   float alive_eff = 0.5;
   
   string baseDir = "Detector/Board_0/OpticalGroup_0/Hybrid_0/Chip_";
   string shortBaseDir = "D_B(0)_O(0)_H(0)_";
   
   for (int ch = 15; ch >= 12; ch--){
     
     string chip = to_string(ch);
#ifdef alive
     string plotDir = "plots/forward_bias/HPK/chip"+chip+"/alive_";
#else
     string plotDir = "plots/forward_bias/HPK/chip"+chip+"/";
#endif
     TCanvas *c0_pixelalive = (TCanvas*) f_pixelalive->Get((baseDir+chip+"/"+shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     TH2F *h_pixelalive = (TH2F*)c0_pixelalive->GetPrimitive((shortBaseDir+"PixelAlive_Chip("+chip+")").c_str());
     
     TCanvas *c0_noHV2Dthreshold = (TCanvas*) f_noHV->Get((baseDir+chip+"/"+shortBaseDir+"Threshold2D_Chip("+chip+")").c_str());
     TCanvas *c0_withHV2Dthreshold = (TCanvas*) f_withHV->Get((baseDir+chip+"/"+shortBaseDir+"Threshold2D_Chip("+chip+")").c_str());
     TH2F *h_noHV2Dthreshold = (TH2F*)c0_noHV2Dthreshold->GetPrimitive((shortBaseDir+"Threshold2D_Chip("+chip+")").c_str());
     TH2F *h_withHV2Dthreshold = (TH2F*)c0_withHV2Dthreshold->GetPrimitive((shortBaseDir+"Threshold2D_Chip("+chip+")").c_str());
     
     int nColumns = h_noHV2Dthreshold->GetXaxis()->GetNbins();
     int nRows = h_noHV2Dthreshold->GetYaxis()->GetNbins();
     
     TCanvas *c0_noHV2Dnoise = (TCanvas*) f_noHV->Get((baseDir+chip+"/"+shortBaseDir+"Noise2D_Chip("+chip+")").c_str());
     TCanvas *c0_withHV2Dnoise = (TCanvas*) f_withHV->Get((baseDir+chip+"/"+shortBaseDir+"Noise2D_Chip("+chip+")").c_str());
     TH2F *h_noHV2Dnoise = (TH2F*)c0_noHV2Dnoise->GetPrimitive((shortBaseDir+"Noise2D_Chip("+chip+")").c_str());
     TH2F *h_withHV2Dnoise = (TH2F*)c0_withHV2Dnoise->GetPrimitive((shortBaseDir+"Noise2D_Chip("+chip+")").c_str());
     
     
     // Plotting
     gStyle->SetOptStat(0);
     
     TCanvas *c_pixelalive = new TCanvas(("c_pixelalive_"+chip).c_str(),"2D threshold distribution without HV",800,800);
     h_pixelalive->Draw("colz");
     c_pixelalive->SaveAs((plotDir+"pixelalive_"+chip+".png").c_str());
     
     
     
     TCanvas *c_noHV2Dthreshold = new TCanvas(("c_noHV2DThreshold_"+chip).c_str(),"2D threshold distribution without HV",800,800);
     h_noHV2Dthreshold->Draw("colz");
     c_noHV2Dthreshold->SaveAs((plotDir+"noHV2Dthreshold_"+chip+".png").c_str());
     
     TCanvas *c_withHV2Dthreshold = new TCanvas(("c_withHV2DThreshold_"+chip).c_str(),"2D threshold distribution without HV",800,800);
     h_withHV2Dthreshold->Draw("colz");
     c_withHV2Dthreshold->SaveAs((plotDir+"withHV2Dthreshold_"+chip+".png").c_str());
     
     TCanvas *c_diff2Dthreshold = new TCanvas(("c_diff2Dthreshold_"+chip).c_str(),"2D threshold difference with/without HV",800,800);
     TH2F *h_diff2Dthreshold = (TH2F*)h_noHV2Dthreshold->Clone("h_diff2Dthreshold");
     h_diff2Dthreshold->SetTitle(("threshold difference of chip "+chip).c_str());
     h_diff2Dthreshold->Add(h_withHV2Dthreshold,-1);
     h_diff2Dthreshold->GetZaxis()->SetRangeUser(-10,250);
     h_diff2Dthreshold->Draw("colz");
     c_diff2Dthreshold->SaveAs((plotDir+"diff2Dthreshold_"+chip+".png").c_str());
     
     TCanvas *c_diff1Dthreshold = new TCanvas(("c_diff1Dthreshold_"+chip).c_str(),"threshold without HV - threshold with HV",800,800);
     TH1F *h_diff1Dthreshold = new TH1F(("h_diff1Dthreshold_"+chip).c_str(),"threshold without HV - threshold with HV",100,-1000,1000);
     for (int i=0; i<nRows; i++){
       for (int j=0; j<nColumns; j++){
#ifdef alive
         if (h_pixelalive->GetBinContent(j+1,i+1)>alive_eff) h_diff1Dthreshold->Fill(h_diff2Dthreshold->GetBinContent(j+1,i+1));
#else
         h_diff1Dthreshold->Fill(h_diff2Dthreshold->GetBinContent(j+1,i+1));
#endif
       }
     }
     h_diff1Dthreshold->Draw("colz");
     c_diff1Dthreshold->SaveAs((plotDir+"diff1Dthreshold_"+chip+".png").c_str());
     
      
     
     
     
     
     TCanvas *c_noHV2Dnoise = new TCanvas(("c_noHV2Dnoise_"+chip).c_str(),"2D noise distribution without HV",800,800);
     h_noHV2Dnoise->Draw("colz");
     c_noHV2Dnoise->SaveAs((plotDir+"noHV2Dnoise_"+chip+".png").c_str());
     
     TCanvas *c_withHV2Dnoise = new TCanvas(("c_withHV2Dnoise_"+chip).c_str(),"2D noise distribution without HV",800,800);
     h_withHV2Dnoise->Draw("colz");
     c_withHV2Dnoise->SaveAs((plotDir+"withHV2Dnoise_"+chip+".png").c_str());
     
     TCanvas *c_diff2Dnoise = new TCanvas(("c_diff2Dnoise_"+chip).c_str(),"2D noise difference with/without HV",800,800);
     TH2F *h_diff2Dnoise = (TH2F*)h_noHV2Dnoise->Clone("h_diff2Dnoise");
     h_diff2Dnoise->SetTitle(("noise difference of chip "+chip).c_str());
     h_diff2Dnoise->Add(h_withHV2Dnoise,-1);
     h_diff2Dnoise->GetZaxis()->SetRangeUser(-40,50);
     h_diff2Dnoise->Draw("colz");
     c_diff2Dnoise->SaveAs((plotDir+"diff2Dnoise_"+chip+".png").c_str());
     
     TCanvas *c_diff1Dnoise = new TCanvas(("c_diff1Dnoise_"+chip).c_str(),"noise without HV - noise with HV, chip 15",800,800);
     TH1F *h_diff1Dnoise = new TH1F(("h_diff1Dnoise_"+chip).c_str(),"noise without HV - noise with HV",100,-200,300);
     for (int i=0; i<nRows; i++){
       for (int j=0; j<nColumns; j++){
  #ifdef alive
         if (h_pixelalive->GetBinContent(j+1,i+1)>alive_eff) h_diff1Dnoise->Fill(h_diff2Dnoise->GetBinContent(j+1,i+1));
  #else
         h_diff1Dnoise->Fill(h_diff2Dnoise->GetBinContent(j+1,i+1));
  #endif
       }
     }
     h_diff1Dnoise->Draw("colz");
     c_diff1Dnoise->SaveAs((plotDir+"diff1Dnoise_"+chip+".png").c_str());
     
     
     
     
     
     
     TCanvas *c_diffThresholdNoise = new TCanvas(("c_diffThresholdNoise_"+chip).c_str(),"threshold and noise difference",800,800);
     TH2F *h_diffThresholdNoise = new TH2F(("h_diffThresholdNoise_"+chip).c_str(),"threshold and noise difference;forward bias noise - depleted noise;forward bias threshold - depleted threshold",100,-100,100,100,-100,400);
     for (int i=0; i<nRows; i++){
       for (int j=0; j<nColumns; j++){
  #ifdef alive
         if (h_pixelalive->GetBinContent(j+1,i+1)>alive_eff) h_diffThresholdNoise->Fill(h_diff2Dnoise->GetBinContent(j+1,i+1),h_diff2Dthreshold->GetBinContent(j+1,i+1));
  #else
         h_diffThresholdNoise->Fill(h_diff2Dnoise->GetBinContent(j+1,i+1),h_diff2Dthreshold->GetBinContent(j+1,i+1));
  #endif
       }
     }
     h_diffThresholdNoise->Draw("colz");
     c_diffThresholdNoise->SaveAs((plotDir+"diffThresholdNoise_"+chip+".png").c_str());
     
     
   } //loop on chips
  
}