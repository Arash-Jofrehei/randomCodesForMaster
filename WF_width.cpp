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
float x_start = -10;
float x_end = 10;
float y_start = -10;
float y_end = 10;
float max_WLS_spike_ratio = 2;
float max_fscint_spike_ratio = 2;
float WLS_ratio_jump = 0.002;
float fscint_ratio_jump = 0.002;

TH1F* convoluted_pulse_spike;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
float Stat[6];

float *stat(TH1F *h){
  Stat[1] = h->FindFirstBinAbove(h->GetMaximum()/3);
  Stat[2] = h->GetMaximumBin();
  Stat[3] = h->FindLastBinAbove(h->GetMaximum()/4);
  Stat[0] = h->GetBinContent(Stat[2]);
  Stat[4] = 100.0*h->GetBinContent((Stat[2]+Stat[3])/2)/Stat[0];
  Stat[5] = h->FindFirstBinAbove(h->GetMaximum()/2);
  return Stat;
}

float *stat(TProfile *h){
  Stat[1] = h->FindFirstBinAbove(h->GetMaximum()/3);
  Stat[2] = h->GetMaximumBin();
  Stat[3] = h->FindLastBinAbove(h->GetMaximum()/4);
  Stat[0] = h->GetBinContent(Stat[2]);
  Stat[4] = 100.0*h->GetBinContent((Stat[2]+Stat[3])/2)/Stat[0];
  Stat[5] = h->FindFirstBinAbove(h->GetMaximum()/2);
  return Stat;
}

void WF_width(){
  TCanvas *canvas = new TCanvas("WF width","WF width");
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_bins_shift.root");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  convoluted_pulse_spike = (TH1F*) convoluted_pulses->Get("normalized nuclear counter effect convoluted with APD+electronics");
  convoluted_pulse_fscint = (TH1F*) convoluted_pulses->Get("normalized fiber scintillation convoluted with APD+electronics");
  convoluted_pulse_WLS = (TH1F*) convoluted_pulses->Get("normalized WLS+CeF3 convoluted with APD+electronics");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
  TH1F *sum = new TH1F("sum","sum",1024,-0.1,204.7);
  TH1F *sum1 = new TH1F("sum1","sum1",1024,-0.1,204.7);
  TH1F *sum2 = new TH1F("sum2","sum2",1024,-0.1,204.7);
  
  TH2F *spike_hist2D_apd[4];
  spike_hist2D_apd[0] = new TH2F("spike_hist2D_apd1","nuclear counter effect contribution for APD1;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[1] = new TH2F("spike_hist2D_apd2","nuclear counter effect contribution for APD2;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[2] = new TH2F("spike_hist2D_apd3","nuclear counter effect contribution for APD3;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  spike_hist2D_apd[3] = new TH2F("spike_hist2D_apd4","nuclear counter effect contribution for APD4;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  TH2F *fscint_hist2D_apd[4];
  fscint_hist2D_apd[0] = new TH2F("fscint_hist2D_apd1","fiber scintillation contribution for APD1;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[1] = new TH2F("fscint_hist2D_apd2","fiber scintillation contribution for APD2;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[2] = new TH2F("fscint_hist2D_apd3","fiber scintillation contribution for APD3;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  fscint_hist2D_apd[3] = new TH2F("fscint_hist2D_apd4","fiber scintillation contribution for APD4;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  TH2F *WLS_hist2D_apd[4];
  WLS_hist2D_apd[0] = new TH2F("WLS_hist2D_apd1","WLS contribution for APD1;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[1] = new TH2F("WLS_hist2D_apd2","WLS contribution for APD2;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[2] = new TH2F("WLS_hist2D_apd3","WLS contribution for APD3;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  WLS_hist2D_apd[3] = new TH2F("WLS_hist2D_apd4","WLS contribution for APD4;X (mm);Y (mm)",nbins,x_start,x_end,nbins,y_start,y_end);
  
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
  /*
  float spike_scale = stat(convoluted_pulse_WLS)[0]/stat(convoluted_pulse_spike)[0];
  float fscint_scale = stat(convoluted_pulse_WLS)[0]/stat(convoluted_pulse_fscint)[0];
  float center_scale = stat(convoluted_pulse_WLS)[0]/stat(apd[0][15][15])[0];
  float fiber_scale = stat(convoluted_pulse_WLS)[0]/stat(apd[0][2][2])[0];
  
  TH1F* scaled_convoluted_pulse_spike = new TH1F("scaled convoluted pulse spike","scaled convoluted pulse spike",1024,-0.1,204.7);
  TH1F* scaled_convoluted_pulse_fscint = new TH1F("scaled convoluted pulse fscint","scaled convoluted pulse fscint",1024,-0.1,204.7);
  TH1F* scaled_convoluted_pulse_WLS = new TH1F("scaled convoluted pulse WLS","scaled convoluted pulse WLS",1024,-0.1,204.7);
  TH1F* scaled_center = new TH1F("scaled center","scaled center",1024,-0.1,204.7);
  TH1F* scaled_fiber = new TH1F("scaled fiber","scaled fiber",1024,-0.1,204.7);
  for (int i = 0;i < 1024;i++){
    scaled_convoluted_pulse_WLS->SetBinContent(i+1,convoluted_pulse_WLS->GetBinContent(i+1));
    int fscint_bin = int(stat(convoluted_pulse_WLS)[2]+(i-stat(convoluted_pulse_fscint)[2])*(stat(convoluted_pulse_WLS)[2]-stat(convoluted_pulse_WLS)[1])*1.0/(stat(convoluted_pulse_fscint)[2]-stat(convoluted_pulse_fscint)[1]));
    if (fscint_bin < 1024 && fscint_bin > -1){
      scaled_convoluted_pulse_fscint->SetBinContent(fscint_bin+1,fscint_scale*convoluted_pulse_fscint->GetBinContent(i+1));
      if (fscint_bin != 1023) scaled_convoluted_pulse_fscint->SetBinContent(fscint_bin+2,fscint_scale*convoluted_pulse_fscint->GetBinContent(i+1));
    }
    int spike_bin = int(stat(convoluted_pulse_WLS)[2]+1+(i-stat(convoluted_pulse_spike)[2])*(stat(convoluted_pulse_WLS)[2]-stat(convoluted_pulse_WLS)[1])*1.0/(stat(convoluted_pulse_spike)[2]-stat(convoluted_pulse_spike)[1]));
    if (spike_bin < 1024 && spike_bin > -1){
      scaled_convoluted_pulse_spike->SetBinContent(spike_bin+1,spike_scale*convoluted_pulse_spike->GetBinContent(i+1));
      if (spike_bin != 1023) scaled_convoluted_pulse_spike->SetBinContent(spike_bin+2,spike_scale*convoluted_pulse_spike->GetBinContent(i+1));
    }
    int center_bin = int(stat(convoluted_pulse_WLS)[2]+1+(i-stat(apd[0][15][15])[2])*(stat(convoluted_pulse_WLS)[2]-stat(convoluted_pulse_WLS)[1])*1.0/(stat(apd[0][15][15])[2]-stat(apd[0][15][15])[1]));
    if (center_bin < 1024 && center_bin > -1){
      scaled_center->SetBinContent(center_bin+1,center_scale*apd[0][15][15]->GetBinContent(i+1));
      if (center_bin != 1023) scaled_center->SetBinContent(center_bin+2,center_scale*apd[0][15][15]->GetBinContent(i+1));
    }
    int fiber_bin = int(stat(convoluted_pulse_WLS)[2]+1+(i-stat(apd[0][2][2])[2])*(stat(convoluted_pulse_WLS)[2]-stat(convoluted_pulse_WLS)[1])*1.0/(stat(apd[0][2][2])[2]-stat(apd[0][2][2])[1]));
    if (fiber_bin < 1024 && fiber_bin > -1){
      scaled_fiber->SetBinContent(fiber_bin+1,fiber_scale*apd[0][2][2]->GetBinContent(i+1));
      if (center_bin != 1023) scaled_fiber->SetBinContent(fiber_bin+2,fiber_scale*apd[0][2][2]->GetBinContent(i+1));
    }
  }
  scaled_convoluted_pulse_spike->SetLineColor(2);
  scaled_convoluted_pulse_fscint->SetLineColor(3);
  scaled_convoluted_pulse_WLS->SetLineColor(4);
  scaled_convoluted_pulse_spike->Draw();
  //scaled_convoluted_pulse_fscint->Draw("same");
  //scaled_convoluted_pulse_WLS->Draw("same");
  scaled_center->Draw("same");
  //scaled_fiber->Draw("same");
  cout << "center:    " << stat(apd[0][15][15])[2]-stat(apd[0][15][15])[1] << "   " << stat(apd[0][15][15])[3]-stat(apd[0][15][15])[2] << "   " << stat(apd[0][15][15])[4] << endl;
  cout << " fiber:    " << stat(apd[0][2][2])[2]-stat(apd[0][2][2])[1] << "    " << stat(apd[0][2][2])[3]-stat(apd[0][2][2])[2] << "   " << stat(apd[0][2][2])[4] << endl;
  cout << " spike:    " << stat(convoluted_pulse_spike)[2]-stat(convoluted_pulse_spike)[1] << "    " << stat(convoluted_pulse_spike)[3]-stat(convoluted_pulse_spike)[2] << "   " << stat(convoluted_pulse_spike)[4] << endl;
  cout << "fscint:    " << stat(convoluted_pulse_fscint)[2]-stat(convoluted_pulse_fscint)[1] << "    " << stat(convoluted_pulse_fscint)[3]-stat(convoluted_pulse_fscint)[2] << "   " << stat(convoluted_pulse_fscint)[4] << endl;
  cout << "   WLS:    " << stat(convoluted_pulse_WLS)[2]-stat(convoluted_pulse_WLS)[1] << "   " << stat(convoluted_pulse_WLS)[3]-stat(convoluted_pulse_WLS)[2] << "   " << stat(convoluted_pulse_WLS)[4] << endl;
  
  for (int i = 0;i < int(max_WLS_spike_ratio/WLS_ratio_jump);i++){
    for (int j = 0;j < int(max_fscint_spike_ratio/fscint_ratio_jump);j++){
      sum->Reset();
      sum->Add(convoluted_pulse_spike);
      sum->Add(convoluted_pulse_fscint,j*fscint_ratio_jump);
      sum->Add(convoluted_pulse_WLS,i*WLS_ratio_jump);
      if (TMath::Abs(stat(sum)[2]-stat(sum)[1] - (stat(apd[0][15][15])[2]-stat(apd[0][15][15])[1])) < 1 && TMath::Abs(stat(sum)[3]-stat(sum)[2] - (stat(apd[0][15][15])[3]-stat(apd[0][15][15])[2])) < 1 && TMath::Abs(stat(sum)[4] - stat(apd[0][15][15])[4]) < 0.6){
        cout << j*fscint_ratio_jump << "    " << i*WLS_ratio_jump << endl;
        cout << stat(sum)[2]-stat(sum)[1] - (stat(apd[0][15][15])[2]-stat(apd[0][15][15])[1]) << "                        " << stat(sum)[3]-stat(sum)[2] - (stat(apd[0][15][15])[3]-stat(apd[0][15][15])[2]) << "                        " << stat(sum)[4] - stat(apd[0][15][15])[4] << endl << endl;
      }
    }
  }*/
  
  float apd_stat[4][nbins][nbins][10];
  for (int a = 0; a < 4; a++){
    for (int i = 0; i < nbins; i++){
      for (int j = 0; j < nbins; j++){
        for (int k = 0; k < 6; k++){
          apd_stat[a][i][j][k] = stat(apd[a][i][j])[k];
        }
        apd_stat[a][i][j][8] = 1000000000;
      }
    }
  }
  
  for (int jump_WLS = 0;jump_WLS < int(max_WLS_spike_ratio/WLS_ratio_jump);jump_WLS++){
    if (jump_WLS%100 == 0) cout << jump_WLS/10 << "    percent" << endl;
    for (int jump_fscint = 0;jump_fscint < int(max_fscint_spike_ratio/fscint_ratio_jump);jump_fscint++){
      sum->Reset();
      sum->Add(convoluted_pulse_spike);
      sum->Add(convoluted_pulse_fscint,jump_fscint*fscint_ratio_jump);
      sum->Add(convoluted_pulse_WLS,jump_WLS*WLS_ratio_jump);
      float stat_sum[6];
      for (int l = 0; l < 6; l++) stat_sum[l] = stat(sum)[l];
      for (int a = 0; a < 4; a++){
        for (int i = 0; i < nbins; i++){
          for (int j = 0; j < nbins; j++){
            //if (i != 3 || j != 3) continue;
            float d = TMath::Abs(stat_sum[4] - apd_stat[a][i][j][4]) + 500*TMath::Abs(stat_sum[2]-stat_sum[1] - (apd_stat[a][i][j][2] - apd_stat[a][i][j][1])) + 500*TMath::Abs(stat_sum[2]-stat_sum[5] - (apd_stat[a][i][j][2] - apd_stat[a][i][j][6])) + 100*TMath::Abs(stat_sum[3]-stat_sum[2] - (apd_stat[a][i][j][3] - apd_stat[a][i][j][2]));
            if (d < apd_stat[a][i][j][8]){
              apd_stat[a][i][j][6] = jump_fscint*fscint_ratio_jump;
              apd_stat[a][i][j][7] = jump_WLS*WLS_ratio_jump;
              apd_stat[a][i][j][8] = d;
              apd_stat[a][i][j][9] = apd_stat[a][i][j][0]/stat_sum[0];
            }
            //if (i == (nbins-1) && j == (nbins-1) && apd_stat[a][i][j][8] == 1000000000 && jump_WLS == int(max_WLS_spike_ratio/WLS_ratio_jump)-1) cout<<"APD"<<a<<"  column "<<i<<"  row "<<j<<"  was unsuccessful"<<endl;
          }
        }
      }
    }
  }
  for (int a = 0; a < 4; a++){
    for (int i = 0; i < nbins; i++){
      for (int j = 0; j < nbins; j++){
        spike_hist2D_apd[a]->SetBinContent(i+1,j+1,apd_stat[a][i][j][9]/1000.0);
        fscint_hist2D_apd[a]->SetBinContent(i+1,j+1,apd_stat[a][i][j][9]*apd_stat[a][i][j][6]/1000.0);
        WLS_hist2D_apd[a]->SetBinContent(i+1,j+1,apd_stat[a][i][j][9]*apd_stat[a][i][j][7]/1000.0);
        if (i == j && a == 0){
          cout << i+1 << "    " << apd_stat[a][i][j][9]/1000.0 << "    " << apd_stat[a][i][j][9]*apd_stat[a][i][j][6]/1000.0 << "    " << apd_stat[a][i][j][9]*apd_stat[a][i][j][7]/1000.0 << endl;
        }
        if (i == 3 && j == 2 && a == 0){
          cout << i+1 << "    " << apd_stat[a][i][j][9]/1000.0 << "    " << apd_stat[a][i][j][9]*apd_stat[a][i][j][6]/1000.0 << "    " << apd_stat[a][i][j][9]*apd_stat[a][i][j][7]/1000.0 << endl;
        }
      }
    }
  }
  spike_hist2D_apd[0]->Draw("colz");
  TCanvas *canvas2 = new TCanvas("WF width2","WF width2");
  fscint_hist2D_apd[0]->Draw("colz");
  TCanvas *canvas3 = new TCanvas("WF width3","WF width3");
  WLS_hist2D_apd[0]->Draw("colz");
  int shift;
  TCanvas *canvas4 = new TCanvas("WF width4","WF width4");
  apd[0][3][2]->Draw();
  sum->Reset();
  sum->Add(convoluted_pulse_spike,apd_stat[0][3][2][9]);
  sum->Add(convoluted_pulse_fscint,apd_stat[0][3][2][9]*apd_stat[0][3][2][6]);
  sum->Add(convoluted_pulse_WLS,apd_stat[0][3][2][9]*apd_stat[0][3][2][7]);
  shift = stat(sum)[2]-stat(apd[0][3][2])[2];
  for (int i = 0;i < 1024;i++){
    if (i < 1024-shift){
      sum->SetBinContent(i+1,sum->GetBinContent(i+1+shift));
    }
    else sum->SetBinContent(i+1,0);
  }
  sum->SetLineColor(2);
  sum->Draw("same");
  TCanvas *canvas5 = new TCanvas("WF width5","WF width5");
  apd[0][10][10]->Draw();
  sum1->Reset();
  sum1->Add(convoluted_pulse_spike,apd_stat[0][10][10][9]);
  sum1->Add(convoluted_pulse_fscint,apd_stat[0][10][10][9]*apd_stat[0][10][10][6]);
  sum1->Add(convoluted_pulse_WLS,apd_stat[0][10][10][9]*apd_stat[0][10][10][7]);
  shift = stat(sum1)[2]-stat(apd[0][10][10])[2];
  for (int i = 0;i < 1024;i++){
    if (i < 1024-shift){
      sum1->SetBinContent(i+1,sum1->GetBinContent(i+1+shift));
    }
    else sum1->SetBinContent(i+1,0);
  }
  sum1->SetLineColor(2);
  sum1->Draw("same");
  TCanvas *canvas6 = new TCanvas("WF width6","WF width6");
  apd[0][15][15]->Draw();
  sum2->Reset();
  sum2->Add(convoluted_pulse_spike,apd_stat[0][15][15][9]);
  sum2->Add(convoluted_pulse_fscint,apd_stat[0][15][15][9]*apd_stat[0][15][15][6]);
  sum2->Add(convoluted_pulse_WLS,apd_stat[0][15][15][9]*apd_stat[0][15][15][7]);
  shift = stat(sum2)[2]-stat(apd[0][15][15])[2];
  for (int i = 0;i < 1024;i++){
    if (i < 1024-shift){
      sum2->SetBinContent(i+1,sum2->GetBinContent(i+1+shift));
    }
    else sum2->SetBinContent(i+1,0);
  }
  sum2->SetLineColor(2);
  sum2->Draw("same");
}