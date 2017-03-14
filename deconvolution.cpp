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
int nbins = 19;
float x_start = -9.5;
float x_end = 9.5;
float y_start = -9.5;
float y_end = 9.5;
TH1F* apd_plus_electronics;
TH1F* convoluted_pulse_fscint;
TH1F* convoluted_pulse_WLS;
TH1F *smoothed_histogram = new TH1F("smoothed histogram","smoothed histogram",1024,-0.1,204.7);

TH1F *smoothed(TH1F* histogram,int radius = 10,int iteration = 1,bool saveMAX = true){
  float interval[2*radius+1];
  float buffer = -1000;
  for (int epoch = 0; epoch < iteration ;epoch++){
    for (int i = 0; i < 1024; i++){
      if (i < radius || i > 1024-radius || (saveMAX == true && TMath::Abs(i+1 - histogram->GetMaximumBin()) < radius)) smoothed_histogram->SetBinContent(i+1,histogram->GetBinContent(i+1));
      else{
        for (int j = 0; j < 2*radius+1 ;j++){
          interval[j] = histogram->GetBinContent(i+1+j-radius);
        }
        if (i == 511) cout << endl << endl;
        for (int j = 0;j < 2*radius+1 ;j++){
          for (int k = j;k < 2*radius+1 ;k++){
            if (interval[j] > interval[k]){
              buffer = interval[j];
              interval[j] = interval[k];
              interval[k] = buffer;
            }
          }
        }
        smoothed_histogram->SetBinContent(i+1,interval[radius]);
      }
    }
  }
  return smoothed_histogram;
}

TH1F* deconv;

TH1F *Deconvolution(TH1F *p,TH1F *h,int size = 1024){
  float main[size];
  float temp[size];
  for (int i = 0;i < size;i++){
    main[i] = p->GetBinContent(i+1);
    temp[i] = h->GetBinContent(i+1);
  }
  float dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = 0;
    if (main[sample]>0) dummy = main[sample]/temp[0];
    for (int i = sample;i < size;i++){
      if (main[i]>=0 && temp[i-sample]>0 && main[i]/temp[i-sample] < dummy) dummy = main[i]/(1.0*temp[i-sample]);
    }
    for (int i = sample;i < size;i++){
      main[i] -= dummy*temp[i-sample];
      if (main[i]<0 && i < 40 && sample < 5) cout << main[33] << "    " << dummy << endl;
    }
    deconv->SetBinContent(sample+1,dummy);
    //cout << sample+1 << "    " << dummy << endl;
  }
  return deconv;
}
/*
TH1F *Deconvolution(TH1F *p,TH1F *h,int size = 1024){
  int buffer_size = 10;
  float dummy[buffer_size];
  float dum = 0;
  float mean_start = 0;
  int sample = 0;
  int negative_counter = 0;
  float negative_fixer;
  float negative_dum[10];
  int bin_recurrence[1024];
  for (int i = 0;i < 1024;i++) bin_recurrence[i] = 0;
  while (sample < size){
    dum = 0;
    for (int i = 0;i < buffer_size;i++) dummy[i] = p->GetBinContent(sample+1+i);
    for (int point = 0;point < sample;point++){
      for (int i = 0;i < buffer_size;i++) dummy[i] -= deconv->GetBinContent(point+1)*h->GetBinContent(i+sample-point+1);
    }
    if (sample == 0){
      for (int i = 0;i < buffer_size;i++) mean_start += dummy[i] / (buffer_size*h->GetBinContent(1+i));
      dum = mean_start;
    }
    else{
      for (int i = 0;i < buffer_size;i++){
        dummy[i] /= p->GetBinContent(1)/mean_start;
        if (h->GetBinContent(i+sample+1) == 0){
          dum = dummy[0];
          break;
        }
        dummy[i] *= h->GetBinContent(sample+1)/h->GetBinContent(i+sample+1);
        dum += dummy[i]/buffer_size;
      }
    }
    if (sample < 90 && sample > 65){
      cout << sample << "    ";
      for (int i = 0;i < 4;i++) cout << dummy[i] << "    ";
      cout << dum << endl;
    }
    if (dum <= 0){
      //negative_dum[negative_counter] = dum;
      dum = 0;
      //negative_counter++;
    }
    deconv->SetBinContent(sample+1,3*dum/4);
    if (negative_counter == 10){
      bin_recurrence[sample]+=1;
      if (bin_recurrence[sample] > 2){
        negative_counter = 0;
       sample++;
       continue;
      }
      negative_fixer = 0;
      for (int i = 0; i < 10;i++){
        negative_fixer += h->GetBinContent(1)*negative_dum[i]/(10.0*h->GetBinContent(i+2));
      }
      cout << sample-10 << "    " << deconv->GetBinContent(sample-9)+negative_fixer << endl;
      negative_counter = 0;
      if (deconv->GetBinContent(sample-9)+negative_fixer > 0) deconv->SetBinContent(sample-9,deconv->GetBinContent(sample-9)+negative_fixer);
      else{
        negative_fixer = h->GetBinContent(1)*(deconv->GetBinContent(sample-9))/(10.0*h->GetBinContent(2));
        for (int i = 0; i < 9;i++){
          negative_fixer += h->GetBinContent(1)*negative_dum[i]/(10.0*h->GetBinContent(i+3));
        }
        cout << sample-11 << "    " << deconv->GetBinContent(sample-10)+negative_fixer << endl;
        if (deconv->GetBinContent(sample-10)+negative_fixer>0) deconv->SetBinContent(sample-10,deconv->GetBinContent(sample-10)+negative_fixer);
        else{
          negative_fixer = h->GetBinContent(1)*(deconv->GetBinContent(sample-10))/(10.0*h->GetBinContent(2)) + h->GetBinContent(1)*(deconv->GetBinContent(sample-9))/(10.0*h->GetBinContent(3));
          for (int i = 0; i < 8;i++){
            negative_fixer += h->GetBinContent(1)*negative_dum[i]/(10.0*h->GetBinContent(i+4));
          }
          cout << sample-12 << "    " << deconv->GetBinContent(sample-10)+negative_fixer << endl;
          if (deconv->GetBinContent(sample-11)+negative_fixer>0) deconv->SetBinContent(sample-11,deconv->GetBinContent(sample-11)+negative_fixer);
          else deconv->SetBinContent(sample-11,0);
          sample--;
          negative_counter = 0;
        }
      }
      sample -= 10;
    }
    sample++;
  }
  return deconv;
}*/
/*
TH1F* deconv;

TH1F *Deconvolution(TH1F *p,TH1F *h,int size = 1024){
  float dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = p->GetBinContent(sample+1);
    for (int point = 0;point < sample;point++){
      dummy -= deconv->GetBinContent(point+1)*h->GetBinContent(sample-point+1);
    }
    dummy /= h->GetBinContent(1);
    cout << sample << "    " << dummy << endl;
    if (dummy < 0) dummy = 0;
    deconv->SetBinContent(sample+1,dummy);
  }
  return deconv;
}
*/
TH1F *conv;
TH1F *convolution(TH1F *f1,TH1F *f2,int size = 1024){
  double dummy;
  for (int sample = 0;sample < size;sample++){
    dummy = 0;
    for (int point = 0;point <= sample;point++){
      dummy += f1->GetBinContent(point+1)*f2->GetBinContent(1+sample-point);
    }
    conv->SetBinContent(sample+1,dummy);
  }
  return conv;
}

void deconvolution(){
  TCanvas *canvas = new TCanvas("deconvolution","deconvolution");
  deconv = new TH1F("d","d",1024,-0.1,204.7);
  conv = new TH1F("c","c",1024,-0.1,204.7);
  TH1F *test = new TH1F("test","test",1024,-0.1,204.7);
  TFile *apd_profile_waveform_bins = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/apd_profile_waveform_19bins_shift.root");
  TFile *convoluted_pulses = new TFile("/afs/cern.ch/work/a/ajofrehe/cern-summer-2016/H4Analysis/ntuples/convoluted_pulses.root");
  apd_plus_electronics = (TH1F*) convoluted_pulses->Get("normalized APD+electronics");
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
  TH1F *smoothed_test = smoothed(test,10,1,false);
  for (int i = 0;i < 1024;i++){
    test->SetBinContent(i+1,smoothed_test->GetBinContent(i+1));
  }
  TH1F *smoothed_apd_plus_electronics = smoothed(apd_plus_electronics,10,1,true);
  smoothed_apd_plus_electronics->Draw();
  for (int i = 0;i < 1024;i++){
    apd_plus_electronics->SetBinContent(i+1,smoothed_apd_plus_electronics->GetBinContent(i+1));
  }
  int offset = 0;
  for(int i = 0;i < 1024;i++){
    if ((apd_plus_electronics->GetBinContent(i+1)/apd_plus_electronics->GetBinContent(apd_plus_electronics->GetMaximumBin())) < 0.001){
      offset += 1;
    }
    if ((apd_plus_electronics->GetBinContent(i+1)/apd_plus_electronics->GetBinContent(apd_plus_electronics->GetMaximumBin())) > 0.001) break;
  }
  for(int i = 0;i < 1024-offset;i++){
    apd_plus_electronics->SetBinContent(i+1,apd_plus_electronics->GetBinContent(i+1+offset));
  }
  for(int i = 1024-offset;i < 1024;i++){
    apd_plus_electronics->SetBinContent(i+1,0);
  }
  offset = 0;
  for(int i = 0;i < 1024;i++){
    if ((apd[0][2][2]->GetBinContent(i+1)/apd[0][2][2]->GetBinContent(apd[0][2][2]->GetMaximumBin())) < 0.001){
      offset += 1;
    }
    if ((apd[0][2][2]->GetBinContent(i+1)/apd[0][2][2]->GetBinContent(apd[0][2][2]->GetMaximumBin())) > 0.001) break;
  }
  for(int i = 0;i < 1024-offset;i++){
    test->SetBinContent(i+1,apd[0][2][2]->GetBinContent(i+1+offset));
  }
  for(int i = 1024-offset;i < 1024;i++){
    test->SetBinContent(i+1,0);
  }
  //test->Smooth(1000);
  //apd_plus_electronics->Smooth(1000);
  //TH1D *dec_hist = Deconvolution(apd_plus_electronics,apd_plus_electronics);
  TH1F *dec_hist = Deconvolution(test,apd_plus_electronics);
  test->Draw();
  conv->SetLineColor(2);
  convolution(dec_hist,apd_plus_electronics);
  conv->Draw("same");
  TCanvas *canvas2 = new TCanvas("deconvolution2","deconvolution2");
  apd_plus_electronics->Draw();
  TCanvas *canvas3 = new TCanvas("deconvolution3","deconvolution3");
  dec_hist->Draw();
}