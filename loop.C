#include <vector>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
using namespace std;
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TF1.h"
#include <TStyle.h>
#include "TROOT.h"
#include <map>
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TBranch.h"
#include "TTree.h"
#include "TChain.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TREE_DATA.C"
#include "TREE_MC.C"
#include "muon_tnp_weight.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
//#include "/Applications/root/macros/AtlasStyle.C"

#define nch_cut
#define singlePlot

string anaTag = "3prong533";
string plotFormat = "pdf";


// **************************
// cut definition
  bool tau_muon_isSoft = true;
  double tau_hadron_vertexprob = 0.1;
  float HFpLeading_low = 0;
  float HFmLeading_low = 0;
  float HFpLeading_high = 4;
  float HFmLeading_high = 4;
  float pionLeadingCut = 0.5;
  float pionSubLeadingCut = 0.3;
  float pionSubSubLeadingCut = 0.3;
  float minZDCp = 0;
  float minZDCm = 0;
  float maxZDCp = 42000;
  float maxZDCm = 60000;
  float ratioZDCpm = 4.2/6;
  int min_nch = 3;
  int max_nch = 300;
  float MET_cut = 20000;
  float deltaPhi_cut = -1;//2.6;
  int MuTauCharge = -1;
  bool hasZDCinfo = false;
// end of cut definition
// **************************

double my_pi=TMath::Pi();
double deltaphi(double phi1, double phi2)
{
 double dphi=fabs(phi1-phi2);
 return (dphi<=my_pi)? dphi : 2.*my_pi-dphi;
}

string baseDir = "/eos/user/a/ajofrehe/gtau/ntuples/";

//void run(string fileData, string fileMCS){//, string fileMCB1, string fileMCB2){
  //string files[2] = {fileData,fileMCS};
  //const int nSamples = 2;

//string inputFiles[] = {"flatTuple_AOD_data_1prong_2015_280721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_SuperChic_2018_280721.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggCCbar_2015_210521.root","flatTuple_mcggBBbar_2015_060521.root","flatTuple_mcggBBbar_5f_2018_190521.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcgPb_AOD_SuperChic_2015_020721.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_090621.root"};
  
//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_data_2015_190421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_050521.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root","flatTuple_mcggCCbar_2015_300421-1.root","flatTuple_mcggBBbar_2015_060521.root"};

//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_subdata_2018_210421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_270421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_subdata_2018_210421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_pp_2pi_MadGraph_2015_280421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_2pi_MadGraph_2018_270421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_2pi_SuperChic_2018_270421.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015_220321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_embedded.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018.root", "../ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018.root", "/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggBBbar_4f_AOD_MadGraph_2015.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_080321.root"};

//string inputFiles[] = {"../ntuples/flatTuple_data_2015.root","../ntuples/flatTuple_mcggTauTau_2015.root","../ntuples/flatTuple_mcggBBbar_2015.root","../ntuples/flatTuple_mcggCCbar_2015.root","../ntuples/flatTuple_mcggTauTau_HydjetDrumMB_2018.root"};

void run(const int nSamples, string files[] = inputFiles){

  string basePlotDir = "/eos/user/a/ajofrehe/www/gtau/"+anaTag;

  //mkdir(directory.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
//void run(vector<string> files){
  //const int nSamples = files.size();
  TFile *root_file[nSamples];
  TREE_DATA *TREE;
  TREE_MC *TREEMCs[nSamples-1];
  TH1F *cutflow[nSamples];
  //SetAtlasStyle();
  setTDRStyle();
  int colors[] = {1,2,3,4,6,7,11};
  int styles[] = {20,1,1,1,1,1,1};

  Double_t Red[2]   = { 1.00, .00};
  Double_t Green[2] = { 1.00, .00};
  Double_t Blue[2]  = { 1.00, .00};
  Double_t Stops[2] = { .00, 1.00};
                                                                               
  Int_t nb=5;
  TColor::CreateGradientColorTable(2,Stops,Red,Green,Blue,nb);

  root_file[0] = new TFile((baseDir+files[0]).c_str(),"READ");
  TREE = new TREE_DATA((TTree*)(root_file[0])->Get("tree"),root_file[0]);
  cutflow[0] = (TH1F*) root_file[0]->Get("ntuplizer/cutflow_perevt");
  cutflow[0]->SetLineColor(colors[0]); cutflow[0]->SetMarkerStyle(styles[0]);
  for (int i = 1; i < nSamples; i++){
    root_file[i] = new TFile((baseDir+files[i]).c_str(),"READ");
    TREEMCs[i-1] = new TREE_MC((TTree*)(root_file[i])->Get("tree"),root_file[i]);
    cutflow[i] = (TH1F*) root_file[i]->Get("ntuplizer/cutflow_perevt");
    cutflow[i]->SetLineColor(colors[i]); cutflow[i]->SetMarkerStyle(styles[i]);
  }
  //if (nSamples>5) cutflow[5]->SetLineColor(6);
  
  TFile *pionEffFile = new TFile("pionEff.root","READ");
  TH2F *pion_eff = (TH2F*) pionEffFile->Get("pion_eff");
  TAxis *xaxis_eff = pion_eff->GetXaxis();
  TAxis *yaxis_eff = pion_eff->GetYaxis();
  //xaxis_eff->kCanExtend = false;
  //yaxis_eff->kCanExtend = false;
  
// *******************************************
  
  
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_5TeV  = "PbPb - 0.55 nb^{-1}";
  lumi_sqrtS = " (5.02 TeV/nucleon)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  //int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
  //  example_plot( iPeriod, 11 );  // left-aligned
  //  example_plot( iPeriod, 33 );  // right-aligned

  //  writeExtraText = false;       // remove Preliminary
  
  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

  //  example_plot( iPeriod, 11 );  // default: left-aligned
  //  example_plot( iPeriod, 22 );  // centered
  //  example_plot( iPeriod, 33 );  // right-aligned  

  int iPeriod = 0;
  int iPos = 11;

  //  if( iPos==0 ) relPosX = 0.12;

  int W = 800;
  int H = 800;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  // Updated by:   Dinko Ferencek (Rutgers)
  //
  int H_ref = 800; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;
/*
  TString canvName = "FigExample_";

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);

  TH1* h = new TH1F("h","h",40,70,110);
  h->GetXaxis()->SetNdivisions(6,5,0);
  h->GetXaxis()->SetTitle("m_{e^{+}e^{-}} (GeV)");  
  h->GetYaxis()->SetNdivisions(6,5,0);
  h->GetYaxis()->SetTitleOffset(1);
  h->GetYaxis()->SetTitle("Events / 0.5 GeV");  

  h->SetMaximum( 260 );
  if( iPos==1 ) h->SetMaximum( 300 );
  h->Draw();

  int histLineColor = kOrange+7;
  int histFillColor = kOrange-2;
  float markerSize  = 1.0;
*/
  // writing the lumi information and the CMS "logo"
  //CMS_lumi( canv, iPeriod, iPos );

  //canv->Update();
  //canv->RedrawAxis();
  //canv->GetFrame()->Draw();

  //canv->Print(canvName+".pdf",".pdf");
  //canv->Print(canvName+".png",".png");

  
  
// *******************************************
// histogram declaration
  
  float dataLumi = 404.318415215; // ub - 2015 data
  //float dataLumi = 425.911750;  // ub - 2018 sub data
  double PiZeroMass = 134.98; //MeV
  //string tag[] = {"", "Signal 2015 MadGraph", "Signal 2018 MadGraph", "Signal 2018 SuperChic", "Signal 2018 SuperChic + Embedded MB"};
  //string tag[] = {"", "MadGraph 2015", "CCbar 2015", "BBbar 2015", "BBbar 2018"};
  string tag[] = {"", "MC Signal", "gammaPb SuperChic 2015", "CCbar 2015", "BBbar 2015"};
  //string tag[] = {"", "MC Signal", "MC 1prong", "MC 1prong + 0#pi^{0}", "MC 1prong + 1#pi^{0}", "MC 1prong + n#pi^{0}", "MC Signal 2018"};
  //string tag[] = {"", "PbPb reco 2015", "pp reco 2015", "pp reco 2018"};
  //string tag[] = {"", "MC Signal 2015", "Signal 2018", "SuperChic 2018"};
  float nEvents[] = {1,1184600,2000000,623000,30000000,1,1,1,1};
  for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1);
  //nEvents[4] /= 10;
  //for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1)/200;
  //float crossSectionMC[4] = {570000,1500,300000,570000};
  //float crossSectionMC[] = {570000,300000,1500,1500};
  float crossSectionMC[] = {570,570,570,570,570,570};
  float SF[] = {1,0,0,0,0,0,0};
  //float pT_SF = 1677539 / 2000000.0;
  float above3GeVSuperChic = 322461;
  float nAllEventsSuperChic = 2000000;
  //float pT_SF = above3GeVSuperChic / nAllEventsSuperChic;
  float pT_SF = 0.21;
  float tau_pT_SF[] = {1,pT_SF,1,1,pT_SF,1};
  //float tau_pT_SF[] = {1,pT_SF,pT_SF,pT_SF,pT_SF,pT_SF,1};
  bool threeProng[] = {1,1,1,1,0,0,0};
  for (int s = 1; s < nSamples; s++) if(nEvents[s] != 0) SF[s] = dataLumi * crossSectionMC[s-1] / nEvents[s];
  //float SF2[] = {1,0,0,0,0,0,0};
  for (int s = 1; s < nSamples; s++) if(cutflow[s]->GetBinContent(1) != 0) SF[s] = tau_pT_SF[s] * dataLumi * crossSectionMC[s-1] / cutflow[s]->GetBinContent(1);
  //SF[1] = nEvents[0]/nEvents[1];
  //SF[2] = nEvents[0]/nEvents[2];

  int tau_mu_pt_bins = 40;
  int tau_mu_eta_bins = 13;
  int tau_mu_phi_bins = 13;
  int tau_hadron_pt_bins = 40;
  int tau_hadron_eta_bins = 13;
  int tau_hadron_phi_bins = 13;
  int tau_hadron_rhomass_bins = 12;
  int tau_hadron_nch_bins = 20;
  int tau_hadron_pvz_bins = 50;
  int calo_energy_bins = 125;
  int deltaphi_bins = 20*4;
  
  TH1F *h_cutflow[nSamples];
  TH1F *h_tau_mu_p[nSamples];
  TH1F *h_tau_mu_pz[nSamples];
  THStack *hs_tau_mu_pt = new THStack("hs_tau_mu_pt","#tau_{#mu} p_{T} [GeV]");
  TH1F *h_tau_mu_pt[nSamples];
  THStack *hs_tau_mu_eta = new THStack("hs_tau_mu_eta","#tau_{#mu} #eta");
  TH1F *h_tau_mu_eta[nSamples];
  THStack *hs_tau_mu_phi = new THStack("hs_tau_mu_phi","#tau_{#mu} #phi");
  TH1F *h_tau_mu_phi[nSamples];
  TH1F *h_tau_hadron_p[nSamples];
  TH1F *h_tau_hadron_pz[nSamples];
  THStack *hs_tau_hadron_pt = new THStack("hs_tau_hadron_pt","#tau_{hadron} p_{T} [GeV]");
  TH1F *h_tau_hadron_pt[nSamples];
  TH1F *h_gen_tau_hadron_visible_pt[nSamples];
  THStack *hs_tau_hadron_eta = new THStack("hs_tau_hadron_eta","#tau_{hadron} #eta");
  TH1F *h_tau_hadron_eta[nSamples];
  THStack *hs_tau_hadron_phi = new THStack("hs_tau_hadron_phi","#tau_{hadron} #phi");
  TH1F *h_tau_hadron_phi[nSamples];
  THStack *hs_tau_hadron_rhomass[2];
  for (int j=0; j < 2;j++) hs_tau_hadron_rhomass[j] = new THStack(("hs_tau_hadron_rhomass" + to_string(j)).c_str(),("#tau_{hadron} #rho mass" + to_string(j) + "[GeV]").c_str());
  TH1F *h_tau_hadron_rhomass[nSamples][2];
  TH2F *h_tau_hadron_rhomass2D[nSamples];
  THStack *hs_tau_hadron_nch = new THStack("hs_tau_hadron_nch","#tau_{hadron} nch");
  TH1F *h_tau_hadron_nch[nSamples];
  TH1F *h_tau_hadron_nPions[nSamples];
  THStack *hs_tau_hadron_ncand_final = new THStack("hs_tau_hadron_ncand_final","#tau_{hadron} cands final");
  TH1F *h_tau_hadron_ncand_final[nSamples];
  THStack *hs_tau_hadron_vprob = new THStack("hs_tau_hadron_vprob","#tau_{hadron} vprob (%)");
  TH1F *h_tau_hadron_vprob[nSamples];
  TH1F *h_tau_hadron_matched_pt_index[nSamples];
  THStack *hs_tau_hadron_mass = new THStack("hs_tau_hadron_mass","#tau_{hadron} mass [GeV]");
  TH1F *h_tau_hadron_mass[nSamples];
  TH1F *h_ditau_p[nSamples];
  TH1F *h_ditau_pz[nSamples];
  THStack *hs_ditau_mass = new THStack("hs_tau_hadron_mass","#tau#tau mass [GeV]");
  TH1F *h_ditau_mass[nSamples];
  THStack *hs_resVisTauPt = new THStack("hs_resTauVis","visible #tau p_{T} resolution [GeV]");
  TH1F *h_resVisTauPt[nSamples];
  TH1F *h_resVisTauEta[nSamples];
  TH1F *h_resVisTauPhi[nSamples];
  THStack *hs_tau_hadron_track_pvz[3];
  for (int j=0; j < 3;j++) hs_tau_hadron_track_pvz[j] = new THStack(("h_tau_hadron_track" + to_string(j) + "_pvz").c_str(),("#tau_{hadron} track" + to_string(j) + " pvz").c_str());
  TH1F *h_tau_hadron_track_pvz[nSamples][3];
  THStack *hs_deltaphi_tau_mu_tau_hadron = new THStack("hs_deltaphi_tau_mu_tau_hadron","#Delta#phi(#tau_{#mu}, #tau_{hadron})");
  TH1F *h_deltaphi_tau_mu_tau_hadron[nSamples];
  TH1F *h_deltaphi_tau_mu_full_tau_hadron[nSamples];
  
  
  THStack *hs_PV_N = new THStack("hs_PV_N","number of PV");
  TH1F *h_PV_N[nSamples];
  TH1F *h_sumZDCplus[nSamples];
  TH1F *h_sumZDCminus[nSamples];
  TH2F *h_sumZDC_pm[nSamples];
  TH2F *h_muEta_averageZDCside[nSamples];
  TH2F *h_averageZDCside_averageHFeta[nSamples];
  TH2F *h_tauEta_averageZDCside[nSamples];
  THStack *hs_calo_energyHFp = new THStack("hs_calo_energyHFp","energy HF+ [GeV]");
  TH1F *h_calo_energyHFp[nSamples];
  THStack *hs_calo_leadingHFp = new THStack("hs_calo_leadingHFp","energy leading tower HF+ [GeV]");
  TH1F *h_calo_leadingHFp[nSamples];
  TH1F *h_calo_energyHFp_sum[nSamples];
  TH1F *h_calo_energyHFp_size[nSamples];
  THStack *hs_calo_energyHFm = new THStack("hs_calo_energyHFm","energy HF- [GeV]");
  TH1F *h_calo_energyHFm[nSamples];
  THStack *hs_calo_leadingHFm = new THStack("hs_calo_leadingHFm","energy leading tower HF- [GeV]");
  TH1F *h_calo_leadingHFm[nSamples];
  TH1F *h_calo_energyHFm_sum[nSamples];
  TH1F *h_calo_energyHFm_size[nSamples];
  TH2F *h_calo_energyHF_pm[nSamples];
  TH2F *h_calo_leadingHF_pm[nSamples];
  TH2F *h_calo_muEta_leadingHFeta[nSamples];
  TH2F *h_calo_tauEta_leadingHFeta[nSamples];
  TH2F *h_calo_muEta_averageHFeta[nSamples];
  TH2F *h_calo_tauEta_averageHFeta[nSamples];
  THStack *hs_calo_HF_eta = new THStack("hs_calo_HF_eta","HF #eta");
  TH1F *h_calo_HF_eta[nSamples];
  TH2F *h_calo_HF_energy_eta[nSamples];
  TH2F *h_deltaphi_tau_mu_tau_hadron_mueta[nSamples];
  TH2F *h_deltaphi_tau_mu_tau_hadron_deltaeta[nSamples];
  TH2F *h_mueta_taueta[nSamples];
  TH2F *h_AP[nSamples];
  
  TH1F *h_N_genTauHadpt[nSamples];
  TH1F *h_eff_3prongReco_genTauHadpt[nSamples];
  TH1F *h_eff_tauReco_genTauHadpt[nSamples];
  TH1F *h_N_genPionPt[nSamples];
  TH1F *h_N_genCentralPionPt[nSamples];
  TH1F *h_eff_matchedPionReco_genPionPt[nSamples];
  TH1F *h_eff_matchedCentralPionReco_genPionPt[nSamples];
  TH1F *h_N_genPionEta[nSamples];
  TH1F *h_N_genHighPtPionEta[nSamples];
  TH1F *h_eff_matchedPionReco_genPionEta[nSamples];
  TH1F *h_eff_matchedHighPtPionReco_genPionEta[nSamples];
  TH1F *h_resPt_matchedPionReco[nSamples];
  TH1F *h_resR_matchedPionReco[nSamples];
  TH2F *h_resPhi_matchedPionReco[nSamples];
  TH2F *h_N_genPionPtEta[nSamples];
  TH2F *h_eff_matchedPionReco_genPionPtEta[nSamples];
  TH2F *h_eff_matchedPionReco_genPionPtEtaZoomed[nSamples];
  
  TH1F *h_NgenGamma[nSamples];
  TH1F *h_genGamma_pt[nSamples];
  TH1F *h_genGamma_eta[nSamples];
  TH1F *h_genGamma_phi[nSamples];
  TH1F *h_genGamma_deltaphi_muon[nSamples];
  TH1F *h_genGamma_deltaphi_pion[nSamples];
  TH1F *h_genGammasMinDeltaR[nSamples];
  TH1F *h_genPiZeroMinDeltaM[nSamples];
  TH1F *h_genPiZeroDeltaM[nSamples];
  TH1F *h_NgenPiZero[nSamples];
  TH1F *h_genPiZero_pt[nSamples];
  TH1F *h_genPiZero_eta[nSamples];
  TH1F *h_genPiZero_phi[nSamples];
  TH1F *h_genPiZero_deltaphi_muon[nSamples];
  TH1F *h_genPiZero_deltaphi_pion[nSamples];
  
  TH1F *h_pion_leading_pt[nSamples];
  TH1F *h_pion_subleading_pt[nSamples];
  TH1F *h_pion_subsubleading_pt[nSamples];
  
  TH1F *h_pion_leading_eta[nSamples];
  TH1F *h_pion_subleading_eta[nSamples];
  TH1F *h_pion_subsubleading_eta[nSamples];
  
  TH1F *h_pion_leading_phi[nSamples];
  TH1F *h_pion_subleading_phi[nSamples];
  TH1F *h_pion_subsubleading_phi[nSamples];
  
  TH1F *h_NrecoGamma[nSamples];
  TH1F *h_recoGamma_pt[nSamples];
  TH1F *h_recoGamma_eta[nSamples];
  TH1F *h_recoGamma_phi[nSamples];
  TH1F *h_recoGamma_deltaphi_muon[nSamples];
  TH1F *h_recoGamma_deltaphi_pion[nSamples];
  TH1F *h_recoGammasMinDeltaR[nSamples];
  TH1F *h_recoPiZeroMinDeltaM[nSamples];
  TH1F *h_recoPiZeroDeltaM[nSamples];
  TH1F *h_NrecoPiZero[nSamples];
  TH1F *h_recoPiZero_pt[nSamples];
  TH1F *h_recoPiZero_eta[nSamples];
  TH1F *h_recoPiZero_phi[nSamples];
  TH1F *h_recoPiZero_deltaphi_muon[nSamples];
  TH1F *h_recoPiZero_deltaphi_pion[nSamples];
  TH2F *h_reco_pion_energy_HCAL_ECAL[nSamples];
  
  TH1F *A_highNch_highHF[nSamples];
  TH1F *B_lowNch_highHF[nSamples];
  TH1F *C_highNch_lowHF[nSamples];
  TH1F *D_lowNch_lowHF[nSamples];
  
  
  // MET
  TH1F *h_MET[nSamples];
  
  for (int i = 0; i < nSamples; i++){
    
    h_cutflow[i] = new TH1F(("h_cutflow_" + tag[i]).c_str(),("Analysis cutflow - " + tag[i]).c_str(),8, 0, 8);
    std::string cutflow_bins_string[] = {"input from Ntuplizer", "#tau #mu", "nCh", "#tau_{hadron}", "HF", "HF & #tau_{hadron}", "HF & #tau_{hadron} & nch", "..."};
    for(size_t j=0; j< 8; j++){
      h_cutflow[0]->GetXaxis()->SetBinLabel(j+1, (cutflow_bins_string[j]).c_str());
    }
    h_cutflow[i]->Sumw2();
    /*if (i != 0) {*/h_cutflow[i]->SetLineColor(colors[i]); h_cutflow[i]->SetMarkerStyle(styles[i]);
  
    A_highNch_highHF[i] = new TH1F(("A_highNch_highHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{hadron}) high Nch - high HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    A_highNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron}) high Nch - high HF"); A_highNch_highHF[i]->Sumw2();
    A_highNch_highHF[i]->SetLineColor(colors[i]); A_highNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    B_lowNch_highHF[i] = new TH1F(("B_lowNch_highHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{hadron}) low Nch - high HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    B_lowNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron}) low Nch - high HF"); B_lowNch_highHF[i]->Sumw2();
    B_lowNch_highHF[i]->SetLineColor(colors[i]); B_lowNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    C_highNch_lowHF[i] = new TH1F(("C_highNch_lowHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{hadron}) high Nch - low HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    C_highNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron}) high Nch - low HF"); C_highNch_lowHF[i]->Sumw2();
    C_highNch_lowHF[i]->SetLineColor(colors[i]); C_highNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
    D_lowNch_lowHF[i] = new TH1F(("D_lowNch_lowHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{hadron}) low Nch - low HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    D_lowNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron}) low Nch - low HF"); D_lowNch_lowHF[i]->Sumw2();
    D_lowNch_lowHF[i]->SetLineColor(colors[i]); D_lowNch_lowHF[i]->SetMarkerStyle(styles[i]);
    
    h_pion_leading_pt[i] = new TH1F(("h_pion_leading_pt_" + tag[i]).c_str(),("leading pion pt " + tag[i]).c_str(),10, -0.5, 9.5);
    h_pion_leading_pt[i]->SetXTitle("leading pion pt [GeV]"); h_pion_leading_pt[i]->Sumw2();
    /*if (i != 0) {*/h_pion_leading_pt[i]->SetLineColor(colors[i]); h_pion_leading_pt[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_pt[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subleading_pt[i] = new TH1F(("h_pion_subleading_pt_" + tag[i]).c_str(),("subleading pion pt " + tag[i]).c_str(),10, -0.7, 9.3);
    h_pion_subleading_pt[i]->SetXTitle("subleading pion pt [GeV]"); h_pion_subleading_pt[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subleading_pt[i]->SetLineColor(colors[i]); h_pion_subleading_pt[i]->SetMarkerStyle(styles[i]);
    h_pion_subleading_pt[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subsubleading_pt[i] = new TH1F(("h_pion_subsubleading_pt_" + tag[i]).c_str(),("subsubleading pion pt " + tag[i]).c_str(),10, -0.7, 9.3);
    h_pion_subsubleading_pt[i]->SetXTitle("subsubleading pion pt [GeV]"); h_pion_subsubleading_pt[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subsubleading_pt[i]->SetLineColor(colors[i]); h_pion_subsubleading_pt[i]->SetMarkerStyle(styles[i]);
    h_pion_subsubleading_pt[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_leading_eta[i] = new TH1F(("h_pion_leading_eta_" + tag[i]).c_str(),("leading pion eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_pion_leading_eta[i]->SetXTitle("leading pion eta"); h_pion_leading_eta[i]->Sumw2();
    /*if (i != 0) {*/h_pion_leading_eta[i]->SetLineColor(colors[i]); h_pion_leading_eta[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_eta[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subleading_eta[i] = new TH1F(("h_pion_subleading_eta_" + tag[i]).c_str(),("subleading pion eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_pion_subleading_eta[i]->SetXTitle("subleading pion eta"); h_pion_subleading_eta[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subleading_eta[i]->SetLineColor(colors[i]); h_pion_subleading_eta[i]->SetMarkerStyle(styles[i]);
    h_pion_subleading_eta[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subsubleading_eta[i] = new TH1F(("h_pion_subsubleading_eta_" + tag[i]).c_str(),("subsubleading pion eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_pion_subsubleading_eta[i]->SetXTitle("subsubleading pion eta"); h_pion_subsubleading_eta[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subsubleading_eta[i]->SetLineColor(colors[i]); h_pion_subsubleading_eta[i]->SetMarkerStyle(styles[i]);
    h_pion_subsubleading_eta[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_leading_phi[i] = new TH1F(("h_pion_leading_phi_" + tag[i]).c_str(),("leading pion phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_pion_leading_phi[i]->SetXTitle("leading pion phi"); h_pion_leading_phi[i]->Sumw2();
    /*if (i != 0) {*/h_pion_leading_phi[i]->SetLineColor(colors[i]); h_pion_leading_phi[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_phi[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subleading_phi[i] = new TH1F(("h_pion_subleading_phi_" + tag[i]).c_str(),("subleading pion phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_pion_subleading_phi[i]->SetXTitle("subleading pion phi"); h_pion_subleading_phi[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subleading_phi[i]->SetLineColor(colors[i]); h_pion_subleading_phi[i]->SetMarkerStyle(styles[i]);
    h_pion_subleading_phi[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subsubleading_phi[i] = new TH1F(("h_pion_subsubleading_phi_" + tag[i]).c_str(),("subsubleading pion phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_pion_subsubleading_phi[i]->SetXTitle("subsubleading pion phi"); h_pion_subsubleading_phi[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subsubleading_phi[i]->SetLineColor(colors[i]); h_pion_subsubleading_phi[i]->SetMarkerStyle(styles[i]);
    h_pion_subsubleading_phi[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_NrecoGamma[i] = new TH1F(("h_NrecoGamma_" + tag[i]).c_str(),("total number of reco gammas " + tag[i]).c_str(),10, -0.5, 9);
    h_NrecoGamma[i]->SetXTitle("number of reco gammas"); h_NrecoGamma[i]->Sumw2();
    /*if (i != 0) {*/h_NrecoGamma[i]->SetLineColor(colors[i]); h_NrecoGamma[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_pt[i] = new TH1F(("h_recoGamma_pt_" + tag[i]).c_str(),("gamma reco p_{T} " + tag[i]).c_str(),40, 0, 8);
    h_recoGamma_pt[i]->SetXTitle("gamma reco p_{T} [GeV]"); h_recoGamma_pt[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_pt[i]->SetLineColor(colors[i]); h_recoGamma_pt[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_eta[i] = new TH1F(("h_recoGamma_eta_" + tag[i]).c_str(),("gamma reco eta " + tag[i]).c_str(),29, -3, 3);
    h_recoGamma_eta[i]->SetXTitle("gamma reco eta"); h_recoGamma_eta[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_eta[i]->SetLineColor(colors[i]); h_recoGamma_eta[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_phi[i] = new TH1F(("h_recoGamma_phi_" + tag[i]).c_str(),("gamma reco phi " + tag[i]).c_str(),20, -TMath::Pi(), TMath::Pi());
    h_recoGamma_phi[i]->SetXTitle("gamma reco phi"); h_recoGamma_phi[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_phi[i]->SetLineColor(colors[i]); h_recoGamma_phi[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_deltaphi_muon[i] = new TH1F(("h_recoGamma_deltaphi_muon_" + tag[i]).c_str(),("#Delta#phi(reco #mu, reco gamma) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_recoGamma_deltaphi_muon[i]->SetXTitle("#Delta#phi(reco #mu, reco gamma)"); h_recoGamma_deltaphi_muon[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_deltaphi_muon[i]->SetLineColor(colors[i]); h_recoGamma_deltaphi_muon[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_deltaphi_pion[i] = new TH1F(("h_recoGamma_deltaphi_pion_" + tag[i]).c_str(),("#Delta#phi(reco charged pion, reco gamma) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_recoGamma_deltaphi_pion[i]->SetXTitle("#Delta#phi(reco charged pion, reco gamma)"); h_recoGamma_deltaphi_pion[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_deltaphi_pion[i]->SetLineColor(colors[i]); h_recoGamma_deltaphi_pion[i]->SetMarkerStyle(styles[i]);
    
    h_recoGammasMinDeltaR[i] = new TH1F(("h_recoGammasMinDeltaR_" + tag[i]).c_str(),("min #DeltaR of reco gammas " + tag[i]).c_str(),100, 0, 1.75);
    h_recoGammasMinDeltaR[i]->SetXTitle("min #DeltaR of reco gammas"); h_recoGammasMinDeltaR[i]->Sumw2();
    /*if (i != 0) {*/h_recoGammasMinDeltaR[i]->SetLineColor(colors[i]); h_recoGammasMinDeltaR[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZeroMinDeltaM[i] = new TH1F(("h_recoPiZeroMinDeltaM_" + tag[i]).c_str(),("min mass difference of reco #pi^{0} wrt 134.98 MeV " + tag[i]).c_str(),50, -200, 1200);
    h_recoPiZeroMinDeltaM[i]->SetXTitle("reco #pi^{0} mass - 134.98 (MeV)"); h_recoPiZeroMinDeltaM[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZeroMinDeltaM[i]->SetLineColor(colors[i]); h_recoPiZeroMinDeltaM[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZeroDeltaM[i] = new TH1F(("h_recoPiZeroDeltaM_" + tag[i]).c_str(),("mass difference of reco #pi^{0} wrt 134.98 MeV " + tag[i]).c_str(),50, -200, 1200);
    h_recoPiZeroDeltaM[i]->SetXTitle("reco #pi^{0} mass - 134.98 (MeV)"); h_recoPiZeroDeltaM[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZeroDeltaM[i]->SetLineColor(colors[i]); h_recoPiZeroDeltaM[i]->SetMarkerStyle(styles[i]);
    
    h_NrecoPiZero[i] = new TH1F(("h_NrecoPiZero_" + tag[i]).c_str(),("number of reco PiZero candidates " + tag[i]).c_str(),10, -0.5, 9);
    h_NrecoPiZero[i]->SetXTitle("number of reco PiZeros"); h_NrecoPiZero[i]->Sumw2();
    /*if (i != 0) {*/h_NrecoPiZero[i]->SetLineColor(colors[i]); h_NrecoPiZero[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZero_pt[i] = new TH1F(("h_recoPiZero_pt_" + tag[i]).c_str(),("PiZero reco p_{T} " + tag[i]).c_str(),40, 0, 8);
    h_recoPiZero_pt[i]->SetXTitle("PiZero reco p_{T} [GeV]"); h_recoPiZero_pt[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZero_pt[i]->SetLineColor(colors[i]); h_recoPiZero_pt[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZero_eta[i] = new TH1F(("h_recoPiZero_eta_" + tag[i]).c_str(),("PiZero reco eta " + tag[i]).c_str(),29, -3, 3);
    h_recoPiZero_eta[i]->SetXTitle("PiZero reco eta"); h_recoPiZero_eta[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZero_eta[i]->SetLineColor(colors[i]); h_recoPiZero_eta[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZero_phi[i] = new TH1F(("h_recoPiZero_phi_" + tag[i]).c_str(),("PiZero reco phi " + tag[i]).c_str(),20, -TMath::Pi(), TMath::Pi());
    h_recoPiZero_phi[i]->SetXTitle("PiZero reco phi"); h_recoPiZero_phi[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZero_phi[i]->SetLineColor(colors[i]); h_recoPiZero_phi[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZero_deltaphi_muon[i] = new TH1F(("h_recoPiZero_deltaphi_muon_" + tag[i]).c_str(),("#Delta#phi(reco #mu, reco PiZero) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_recoPiZero_deltaphi_muon[i]->SetXTitle("#Delta#phi(reco #mu, reco PiZero)"); h_recoPiZero_deltaphi_muon[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZero_deltaphi_muon[i]->SetLineColor(colors[i]); h_recoPiZero_deltaphi_muon[i]->SetMarkerStyle(styles[i]);
    
    h_recoPiZero_deltaphi_pion[i] = new TH1F(("h_recoPiZero_deltaphi_pion_" + tag[i]).c_str(),("#Delta#phi(reco charged pion, reco PiZero) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_recoPiZero_deltaphi_pion[i]->SetXTitle("#Delta#phi(reco charged pion, reco PiZero)"); h_recoPiZero_deltaphi_pion[i]->Sumw2();
    /*if (i != 0) {*/h_recoPiZero_deltaphi_pion[i]->SetLineColor(colors[i]); h_recoPiZero_deltaphi_pion[i]->SetMarkerStyle(styles[i]);
    
    h_reco_pion_energy_HCAL_ECAL[i] = new TH2F(("h_reco_pion_energy_HCAL_ECAL_" + tag[i]).c_str(),("reco charged pion HCAL vs ECAL energy" + tag[i]).c_str(),20, 0, 10,20, 0, 10);
    h_reco_pion_energy_HCAL_ECAL[i]->SetXTitle("ECAL energy"); h_reco_pion_energy_HCAL_ECAL[i]->SetYTitle("HCAL energy");
    
    h_NgenGamma[i] = new TH1F(("h_NgenGamma_" + tag[i]).c_str(),("total number of gen gammas " + tag[i]).c_str(),10, -0.5, 9);
    h_NgenGamma[i]->SetXTitle("number of gen gammas"); h_NgenGamma[i]->Sumw2();
    /*if (i != 0) {*/h_NgenGamma[i]->SetLineColor(colors[i]); h_NgenGamma[i]->SetMarkerStyle(styles[i]);
    
    h_genGamma_pt[i] = new TH1F(("h_genGamma_pt_" + tag[i]).c_str(),("gamma gen p_{T} " + tag[i]).c_str(),40, 0, 8);
    h_genGamma_pt[i]->SetXTitle("gamma gen p_{T} [GeV]"); h_genGamma_pt[i]->Sumw2();
    /*if (i != 0) {*/h_genGamma_pt[i]->SetLineColor(colors[i]); h_genGamma_pt[i]->SetMarkerStyle(styles[i]);
    
    h_genGamma_eta[i] = new TH1F(("h_genGamma_eta_" + tag[i]).c_str(),("gamma gen eta " + tag[i]).c_str(),29, -3, 3);
    h_genGamma_eta[i]->SetXTitle("gamma gen eta"); h_genGamma_eta[i]->Sumw2();
    /*if (i != 0) {*/h_genGamma_eta[i]->SetLineColor(colors[i]); h_genGamma_eta[i]->SetMarkerStyle(styles[i]);
    
    h_genGamma_phi[i] = new TH1F(("h_genGamma_phi_" + tag[i]).c_str(),("gamma gen phi " + tag[i]).c_str(),20, -TMath::Pi(), TMath::Pi());
    h_genGamma_phi[i]->SetXTitle("gamma gen phi"); h_genGamma_phi[i]->Sumw2();
    /*if (i != 0) {*/h_genGamma_phi[i]->SetLineColor(colors[i]); h_genGamma_phi[i]->SetMarkerStyle(styles[i]);
    
    h_genGamma_deltaphi_muon[i] = new TH1F(("h_genGamma_deltaphi_muon_" + tag[i]).c_str(),("#Delta#phi(reco #mu, gen gamma) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_genGamma_deltaphi_muon[i]->SetXTitle("#Delta#phi(gen #mu, gen gamma)"); h_genGamma_deltaphi_muon[i]->Sumw2();
    /*if (i != 0) {*/h_genGamma_deltaphi_muon[i]->SetLineColor(colors[i]); h_genGamma_deltaphi_muon[i]->SetMarkerStyle(styles[i]);
    
    h_genGamma_deltaphi_pion[i] = new TH1F(("h_genGamma_deltaphi_pion_" + tag[i]).c_str(),("#Delta#phi(reco charged pion, gen gamma) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_genGamma_deltaphi_pion[i]->SetXTitle("#Delta#phi(reco charged pion, gen gamma)"); h_genGamma_deltaphi_pion[i]->Sumw2();
    /*if (i != 0) {*/h_genGamma_deltaphi_pion[i]->SetLineColor(colors[i]); h_genGamma_deltaphi_pion[i]->SetMarkerStyle(styles[i]);
    
    h_genGammasMinDeltaR[i] = new TH1F(("h_genGammasMinDeltaR_" + tag[i]).c_str(),("min #DeltaR of gen gammas " + tag[i]).c_str(),100, 0, 1.75);
    h_genGammasMinDeltaR[i]->SetXTitle("min #DeltaR of gen gammas"); h_genGammasMinDeltaR[i]->Sumw2();
    /*if (i != 0) {*/h_genGammasMinDeltaR[i]->SetLineColor(colors[i]); h_genGammasMinDeltaR[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZeroMinDeltaM[i] = new TH1F(("h_genPiZeroMinDeltaM_" + tag[i]).c_str(),("min mass difference of gen #pi^{0} wrt 134.98 MeV " + tag[i]).c_str(),51, -1, 1);
    h_genPiZeroMinDeltaM[i]->SetXTitle("gen #pi^{0} mass - 134.98 (MeV)"); h_genPiZeroMinDeltaM[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZeroMinDeltaM[i]->SetLineColor(colors[i]); h_genPiZeroMinDeltaM[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZeroDeltaM[i] = new TH1F(("h_genPiZeroDeltaM_" + tag[i]).c_str(),("mass difference of gen #pi^{0} wrt 134.98 MeV " + tag[i]).c_str(),51, -1, 1);
    h_genPiZeroDeltaM[i]->SetXTitle("gen #pi^{0} mass - 134.98 (MeV)"); h_genPiZeroDeltaM[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZeroDeltaM[i]->SetLineColor(colors[i]); h_genPiZeroDeltaM[i]->SetMarkerStyle(styles[i]);
    
    h_NgenPiZero[i] = new TH1F(("h_NgenPiZero_" + tag[i]).c_str(),("number of gen PiZero candidates " + tag[i]).c_str(),10, -0.5, 9);
    h_NgenPiZero[i]->SetXTitle("number of gen PiZeros"); h_NgenPiZero[i]->Sumw2();
    /*if (i != 0) {*/h_NgenPiZero[i]->SetLineColor(colors[i]); h_NgenPiZero[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZero_pt[i] = new TH1F(("h_genPiZero_pt_" + tag[i]).c_str(),("PiZero gen p_{T} " + tag[i]).c_str(),40, 0, 8);
    h_genPiZero_pt[i]->SetXTitle("PiZero gen p_{T} [GeV]"); h_genPiZero_pt[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZero_pt[i]->SetLineColor(colors[i]); h_genPiZero_pt[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZero_eta[i] = new TH1F(("h_genPiZero_eta_" + tag[i]).c_str(),("PiZero gen eta " + tag[i]).c_str(),29, -3, 3);
    h_genPiZero_eta[i]->SetXTitle("PiZero gen eta"); h_genPiZero_eta[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZero_eta[i]->SetLineColor(colors[i]); h_genPiZero_eta[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZero_phi[i] = new TH1F(("h_genPiZero_phi_" + tag[i]).c_str(),("PiZero gen phi " + tag[i]).c_str(),20, -TMath::Pi(), TMath::Pi());
    h_genPiZero_phi[i]->SetXTitle("PiZero gen phi"); h_genPiZero_phi[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZero_phi[i]->SetLineColor(colors[i]); h_genPiZero_phi[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZero_deltaphi_muon[i] = new TH1F(("h_genPiZero_deltaphi_muon_" + tag[i]).c_str(),("#Delta#phi(reco #mu, gen PiZero) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_genPiZero_deltaphi_muon[i]->SetXTitle("#Delta#phi(gen #mu, gen PiZero)"); h_genPiZero_deltaphi_muon[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZero_deltaphi_muon[i]->SetLineColor(colors[i]); h_genPiZero_deltaphi_muon[i]->SetMarkerStyle(styles[i]);
    
    h_genPiZero_deltaphi_pion[i] = new TH1F(("h_genPiZero_deltaphi_pion_" + tag[i]).c_str(),("#Delta#phi(reco charged pion, gen PiZero) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_genPiZero_deltaphi_pion[i]->SetXTitle("#Delta#phi(reco charged pion, gen PiZero)"); h_genPiZero_deltaphi_pion[i]->Sumw2();
    /*if (i != 0) {*/h_genPiZero_deltaphi_pion[i]->SetLineColor(colors[i]); h_genPiZero_deltaphi_pion[i]->SetMarkerStyle(styles[i]);
    
    h_N_genPionPt[i] = new TH1F(("h_N_genPionPt_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),40, 0, 8);
    h_N_genPionPt[i]->SetXTitle("pion gen p_{T} [GeV]"); h_N_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genPionPt[i]->SetLineColor(colors[i]); h_N_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_genCentralPionPt[i] = new TH1F(("h_N_genCentralPionPt_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),40, 0, 8);
    h_N_genCentralPionPt[i]->SetXTitle("pion gen p_{T} [GeV]"); h_N_genCentralPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genCentralPionPt[i]->SetLineColor(colors[i]); h_N_genCentralPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedPionReco_genPionPt[i] = new TH1F(("h_eff_matchedPionReco_genPionPt_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),40, 0, 8);
    h_eff_matchedPionReco_genPionPt[i]->SetXTitle("pion gen p_{T} [GeV]"); h_eff_matchedPionReco_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedPionReco_genPionPt[i]->SetLineColor(colors[i]); h_eff_matchedPionReco_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedCentralPionReco_genPionPt[i] = new TH1F(("h_eff_matchedCentralPionReco_genPionPt_" + tag[i]).c_str(),("Efficiency of reconstructing central matched pions - " + tag[i]).c_str(),40, 0, 8);
    h_eff_matchedCentralPionReco_genPionPt[i]->SetXTitle("pion gen p_{T} [GeV]"); h_eff_matchedCentralPionReco_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedCentralPionReco_genPionPt[i]->SetLineColor(colors[i]); h_eff_matchedCentralPionReco_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_genPionEta[i] = new TH1F(("h_N_genPionEta_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),30, -3, 3);
    h_N_genPionEta[i]->SetXTitle("pion gen #eta"); h_N_genPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_N_genPionEta[i]->SetLineColor(colors[i]); h_N_genPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_N_genHighPtPionEta[i] = new TH1F(("h_N_genHighPtPionEta_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),30, -3, 3);
    h_N_genHighPtPionEta[i]->SetXTitle("pion gen #eta"); h_N_genHighPtPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_N_genHighPtPionEta[i]->SetLineColor(colors[i]); h_N_genHighPtPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedPionReco_genPionEta[i] = new TH1F(("h_eff_matchedPionReco_genPionEta_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),30, -3, 3);
    h_eff_matchedPionReco_genPionEta[i]->SetXTitle("pion gen #eta"); h_eff_matchedPionReco_genPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedPionReco_genPionEta[i]->SetLineColor(colors[i]); h_eff_matchedPionReco_genPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedHighPtPionReco_genPionEta[i] = new TH1F(("h_eff_matchedHighPtPionReco_genPionEta_" + tag[i]).c_str(),("Efficiency of reconstructing high P_{T} matched pions - " + tag[i]).c_str(),30, -3, 3);
    h_eff_matchedHighPtPionReco_genPionEta[i]->SetXTitle("pion gen #eta"); h_eff_matchedHighPtPionReco_genPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedHighPtPionReco_genPionEta[i]->SetLineColor(colors[i]); h_eff_matchedHighPtPionReco_genPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_N_genPionPtEta[i] = new TH2F(("h_N_genPionPtEta_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),30, -3, 3, 40, 0, 8);
    h_N_genPionPtEta[i]->SetXTitle("pion gen #eta"); h_N_genPionPtEta[i]->SetYTitle("pion gen p_{T} [GeV]");
    
    h_eff_matchedPionReco_genPionPtEta[i] = new TH2F(("h_eff_matchedPionReco_genPionPtEta_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),30, -3, 3, 40, 0, 8);
    h_eff_matchedPionReco_genPionPtEta[i]->SetXTitle("pion gen #eta"); h_eff_matchedPionReco_genPionPtEta[i]->SetYTitle("pion gen p_{T} [GeV]");
    
    h_eff_matchedPionReco_genPionPtEtaZoomed[i] = new TH2F(("h_eff_matchedPionReco_genPionPtEtaZoomed_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),30, -3, 3, 40, 0, 8);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->SetXTitle("pion gen #eta"); h_eff_matchedPionReco_genPionPtEtaZoomed[i]->SetYTitle("pion gen p_{T} [GeV]");
    
    h_resPt_matchedPionReco[i] = new TH1F(("h_resPt_matchedPionReco_" + tag[i]).c_str(),("relative p_{T} resolution of reconstructed matched pions - " + tag[i]).c_str(),31, 0.8, 1.3);
    h_resPt_matchedPionReco[i]->SetXTitle("reco pion p_{T} / gen pion p_{T} [GeV]"); h_resPt_matchedPionReco[i]->Sumw2();
    /*if (i != 0) {*/h_resPt_matchedPionReco[i]->SetLineColor(colors[i]); h_resPt_matchedPionReco[i]->SetMarkerStyle(styles[i]);
    
    h_resR_matchedPionReco[i] = new TH1F(("h_resR_matchedPionReco_" + tag[i]).c_str(),("#DeltaR resolution of reconstructed matched pions - " + tag[i]).c_str(),20, 0, 0.02);
    h_resR_matchedPionReco[i]->SetXTitle("#Delta R matched gen & reco pions [GeV]"); h_resR_matchedPionReco[i]->Sumw2();
    /*if (i != 0) {*/h_resR_matchedPionReco[i]->SetLineColor(colors[i]); h_resR_matchedPionReco[i]->SetMarkerStyle(styles[i]);
    
    h_resPhi_matchedPionReco[i] = new TH2F(("h_resPhi_matchedPionReco_" + tag[i]).c_str(),("#phi resolution of reconstructed matched pions - " + tag[i]).c_str(),20, -0.02, 0.02, 20, -0.02, 0.02);
    h_resPhi_matchedPionReco[i]->SetXTitle("#phi matched gen pions [GeV]");
    h_resPhi_matchedPionReco[i]->SetYTitle("#phi matched reco pions [GeV]");
    
    h_N_genTauHadpt[i] = new TH1F(("h_N_genTauHadpt_" + tag[i]).c_str(),("Number of events with a reconstructed muon - " + tag[i]).c_str(),20, 0, 20);
    h_N_genTauHadpt[i]->SetXTitle("#tau_{hadron} gen p_{T} [GeV]"); h_N_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genTauHadpt[i]->SetLineColor(colors[i]); h_N_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_3prongReco_genTauHadpt[i] = new TH1F(("h_eff_3prongReco_genTauHadpt_" + tag[i]).c_str(),("Efficiency of reconstructing 3 pions - " + tag[i]).c_str(),20, 0, 20);
    h_eff_3prongReco_genTauHadpt[i]->SetXTitle("#tau_{hadron} gen p_{T} [GeV]"); h_eff_3prongReco_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_3prongReco_genTauHadpt[i]->SetLineColor(colors[i]); h_eff_3prongReco_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_tauReco_genTauHadpt[i] = new TH1F(("h_eff_tauReco_genTauHadpt_" + tag[i]).c_str(),("Efficiency of reconstructing a hadronic tau - " + tag[i]).c_str(),20, 0, 20);
    h_eff_tauReco_genTauHadpt[i]->SetXTitle("#tau_{hadron} gen p_{T} [GeV]"); h_eff_tauReco_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_tauReco_genTauHadpt[i]->SetLineColor(colors[i]); h_eff_tauReco_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_tau_mu_p[i] = new TH1F(("h_tau_mu_p_" + tag[i]).c_str(),("#tau_{#mu} total p " + tag[i]).c_str(),30, 0, 30);
    h_tau_mu_p[i]->SetXTitle("#tau_{#mu} total p [GeV]"); h_tau_mu_p[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_p[i]->SetLineColor(colors[i]); h_tau_mu_p[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_p[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_pz[i] = new TH1F(("h_tau_mu_pz_" + tag[i]).c_str(),("#tau_{#mu} p_{z} " + tag[i]).c_str(),30, 0, 30);
    h_tau_mu_pz[i]->SetXTitle("#tau_{#mu} p_{z} [GeV]"); h_tau_mu_pz[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_pz[i]->SetLineColor(colors[i]); h_tau_mu_pz[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_pz[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_pt[i] = new TH1F(("h_tau_mu_pt_" + tag[i]).c_str(),("#tau_{#mu} p_{T} " + tag[i]).c_str(),20, 0, 20);
    h_tau_mu_pt[i]->SetXTitle("#tau_{#mu} p_{T} [GeV]"); h_tau_mu_pt[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_pt[i]->SetLineColor(colors[i]); h_tau_mu_pt[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_pt[i]->SetBinErrorOption(TH1::kPoisson);
        
    h_tau_mu_eta[i] = new TH1F(("h_tau_mu_eta_" + tag[i]).c_str(),("#tau_{#mu} #eta " + tag[i]).c_str(),tau_mu_eta_bins, -2.5, 2.5);
    h_tau_mu_eta[i]->SetXTitle("#tau_{#mu} #eta"); h_tau_mu_eta[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_eta[i]->SetLineColor(colors[i]); h_tau_mu_eta[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_eta[i]->SetBinErrorOption(TH1::kPoisson);
        
    h_tau_mu_phi[i] = new TH1F(("h_tau_mu_phi_" + tag[i]).c_str(),("#tau_{#mu} #phi " + tag[i]).c_str(),tau_mu_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_mu_phi[i]->SetXTitle("#tau_{#mu} #phi"); h_tau_mu_phi[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_phi[i]->SetLineColor(colors[i]); h_tau_mu_phi[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_phi[i]->SetBinErrorOption(TH1::kPoisson);
    
    // hadron related
    
    h_tau_hadron_p[i] = new TH1F(("h_tau_hadron_p_" + tag[i]).c_str(),("#tau_{hadron} total p " + tag[i]).c_str(),10, 0, 20);
    h_tau_hadron_p[i]->SetXTitle("#tau_{hadron} total p [GeV]"); h_tau_hadron_p[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_p[i]->SetLineColor(colors[i]); h_tau_hadron_p[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_pz[i] = new TH1F(("h_tau_hadron_pz_" + tag[i]).c_str(),("#tau_{hadron} p_{z} " + tag[i]).c_str(),10, 0, 20);
    h_tau_hadron_pz[i]->SetXTitle("#tau_{hadron} p_{z} [GeV]"); h_tau_hadron_pz[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_pz[i]->SetLineColor(colors[i]); h_tau_hadron_pz[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_pt[i] = new TH1F(("h_tau_hadron_pt_" + tag[i]).c_str(),("#tau_{hadron} p_{T} " + tag[i]).c_str(),20, 0, 20);
    h_tau_hadron_pt[i]->SetXTitle("#tau_{hadron} p_{T} [GeV]"); h_tau_hadron_pt[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_pt[i]->SetLineColor(colors[i]); h_tau_hadron_pt[i]->SetMarkerStyle(styles[i]);
    
    h_gen_tau_hadron_visible_pt[i] = new TH1F(("h_gen_tau_hadron_visible_pt_" + tag[i]).c_str(),("gen #tau_{hadron} visible p_{T} " + tag[i]).c_str(),20, 0, 20);
    h_gen_tau_hadron_visible_pt[i]->SetXTitle("#tau_{hadron} p_{T} [GeV]"); h_gen_tau_hadron_visible_pt[i]->Sumw2();
    /*if (i != 0) {*/h_gen_tau_hadron_visible_pt[i]->SetLineColor(colors[i]); h_gen_tau_hadron_visible_pt[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_eta[i] = new TH1F(("h_tau_hadron_eta_" + tag[i]).c_str(),("#tau_{hadron} #eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_tau_hadron_eta[i]->SetXTitle("#tau_{hadron} #eta"); h_tau_hadron_eta[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_eta[i]->SetLineColor(colors[i]); h_tau_hadron_eta[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_phi[i] = new TH1F(("h_tau_hadron_phi_" + tag[i]).c_str(),("#tau_{hadron} #phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_hadron_phi[i]->SetXTitle("#tau_{hadron} #phi"); h_tau_hadron_phi[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_phi[i]->SetLineColor(colors[i]); h_tau_hadron_phi[i]->SetMarkerStyle(styles[i]);
        
    for (int j=0; j < 2;j++){
      h_tau_hadron_rhomass[i][j] = new TH1F(("h_tau_hadron_rhomass" + to_string(j) + "_" + tag[i]).c_str(),("#tau_{hadron} #rho_{"+to_string(j+1)+"} mass [GeV] " + tag[i]).c_str(),tau_hadron_rhomass_bins, 0.2, 1.4);
      h_tau_hadron_rhomass[i][j]->SetXTitle(("#rho_{"+to_string(j+1)+"} mass [GeV]").c_str()); h_tau_hadron_rhomass[i][j]->Sumw2();
      /*if (i != 0) {*/h_tau_hadron_rhomass[i][j]->SetLineColor(colors[i]); h_tau_hadron_rhomass[i][j]->SetMarkerStyle(styles[i]);
    }
    
    h_tau_hadron_rhomass2D[i] = new TH2F(("h_tau_hadron_rhomass2D_" + tag[i]).c_str(),("#tau_{hadron} #rho mass [GeV] " + tag[i]).c_str(),2*tau_hadron_rhomass_bins, 0.2, 1.4, 2*tau_hadron_rhomass_bins, 0.2, 1.4);
    h_tau_hadron_rhomass2D[i]->SetXTitle("#rho_{1} mass [GeV]"); h_tau_hadron_rhomass2D[i]->SetYTitle("#rho_{2} mass [GeV]");
    
    h_tau_hadron_nch[i] = new TH1F(("h_tau_hadron_nch_" + tag[i]).c_str(),("h_tau_hadron_nch_" + tag[i]).c_str(),16, -0.5, 15.5);
    h_tau_hadron_nch[i]->SetXTitle("nch"); h_tau_hadron_nch[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_nch[i]->SetLineColor(colors[i]); h_tau_hadron_nch[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_nPions[i] = new TH1F(("h_tau_hadron_nPions_" + tag[i]).c_str(),("h_tau_hadron_nPions_" + tag[i]).c_str(),7, -0.5, 6.5);
    h_tau_hadron_nPions[i]->SetXTitle("nPions"); h_tau_hadron_nPions[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_nPions[i]->SetLineColor(colors[i]); h_tau_hadron_nPions[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_ncand_final[i] = new TH1F(("h_tau_hadron_ncand_final_" + tag[i]).c_str(),("h_tau_hadron_ncand_final_" + tag[i]).c_str(),10, 0.5, 10.5);
    h_tau_hadron_ncand_final[i]->SetXTitle("ncand final"); h_tau_hadron_ncand_final[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_ncand_final[i]->SetLineColor(colors[i]); h_tau_hadron_ncand_final[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_vprob[i] = new TH1F(("h_tau_hadron_vprob_" + tag[i]).c_str(),("h_tau_hadron_vprob_" + tag[i]).c_str(),20, 0, 100);
    h_tau_hadron_vprob[i]->SetXTitle("vprob (%)"); h_tau_hadron_vprob[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_vprob[i]->SetLineColor(colors[i]); h_tau_hadron_vprob[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_matched_pt_index[i] = new TH1F(("h_tau_hadron_matched_pt_index_" + tag[i]).c_str(),("h_tau_hadron_matched_pt_index_" + tag[i]).c_str(),10, 0.5, 10.5);
    h_tau_hadron_matched_pt_index[i]->SetXTitle("p_T rank of the matched tau"); h_tau_hadron_matched_pt_index[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_matched_pt_index[i]->SetLineColor(colors[i]); h_tau_hadron_matched_pt_index[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_mass[i] = new TH1F(("h_tau_hadron_mass_" + tag[i]).c_str(),("#tau_{hadron} visible mass " + tag[i]).c_str(), 16, 0.1, 1.7);
    h_tau_hadron_mass[i]->SetXTitle("#tau_{hadron} mass [GeV]"); h_tau_hadron_mass[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_mass[i]->SetLineColor(colors[i]); h_tau_hadron_mass[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_mass[i] = new TH1F(("h_ditau_mass_" + tag[i]).c_str(),("#tau#tau visible invariant mass " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_mass[i]->SetXTitle("#tau#tau invariant mass [GeV]"); h_ditau_mass[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_mass[i]->SetLineColor(colors[i]); h_ditau_mass[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_p[i] = new TH1F(("h_ditau_p_" + tag[i]).c_str(),("#tau#tau visible total p " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_p[i]->SetXTitle("#tau#tau total p [GeV]"); h_ditau_p[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_p[i]->SetLineColor(colors[i]); h_ditau_p[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_pz[i] = new TH1F(("h_ditau_pz_" + tag[i]).c_str(),("#tau#tau visible p_{z} " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_pz[i]->SetXTitle("#tau#tau p_{z} [GeV]"); h_ditau_pz[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_pz[i]->SetLineColor(colors[i]); h_ditau_pz[i]->SetMarkerStyle(styles[i]);
    
    h_resVisTauPt[i] = new TH1F(("h_resVisTauPt_" + tag[i]).c_str(),("visible #tau p_{T} resolution " + tag[i]).c_str(), 21, -1, 1);
    h_resVisTauPt[i]->SetXTitle("visible reco #tau p_{T} - gen #tau p_{T} [GeV]"); h_resVisTauPt[i]->Sumw2();
    /*if (i != 0) {*/h_resVisTauPt[i]->SetLineColor(colors[i]); h_resVisTauPt[i]->SetMarkerStyle(styles[i]);
    
    h_resVisTauEta[i] = new TH1F(("h_resVisTauEta_" + tag[i]).c_str(),("visible #tau eta resolution " + tag[i]).c_str(), 21, -0.03, 0.03);
    h_resVisTauEta[i]->SetXTitle("visible reco #tau eta - gen #tau eta"); h_resVisTauEta[i]->Sumw2();
    /*if (i != 0) {*/h_resVisTauEta[i]->SetLineColor(colors[i]); h_resVisTauEta[i]->SetMarkerStyle(styles[i]);
    
    h_resVisTauPhi[i] = new TH1F(("h_resVisTauPhi_" + tag[i]).c_str(),("visible #tau phi resolution " + tag[i]).c_str(), 21, -0.03, 0.03);
    h_resVisTauPhi[i]->SetXTitle("visible reco #tau phi - gen #tau phi"); h_resVisTauPhi[i]->Sumw2();
    /*if (i != 0) {*/h_resVisTauPhi[i]->SetLineColor(colors[i]); h_resVisTauPhi[i]->SetMarkerStyle(styles[i]);
        
    for (int j=0; j < 3;j++){
      h_tau_hadron_track_pvz[i][j] = new TH1F(("h_tau_hadron_track" + to_string(j) + "_pvz_" + tag[i]).c_str(),("h_tau_hadron_track" + to_string(j) + "_pvz_" + tag[i]).c_str(),tau_hadron_pvz_bins, -0.5, 0.5);
      h_tau_hadron_track_pvz[i][j]->SetXTitle(("PV_{z} - trk"+to_string(j)+"_{z}").c_str()); h_tau_hadron_track_pvz[i][j]->Sumw2();
      /*if (i != 0) {*/h_tau_hadron_track_pvz[i][j]->SetLineColor(colors[i]); h_tau_hadron_track_pvz[i][j]->SetMarkerStyle(styles[i]);
    }
    
    // other
        
    h_deltaphi_tau_mu_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{hadron}) " + tag[i]).c_str(),deltaphi_bins/4, 0.75*TMath::Pi(), TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron[i]->Sumw2();
    h_deltaphi_tau_mu_tau_hadron[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_tau_hadron[i]->SetMarkerStyle(styles[i]);
    h_deltaphi_tau_mu_tau_hadron[i]->SetBinErrorOption(TH1::kPoisson);
    /*h_deltaphi_tau_mu_tau_hadron[i]->SetFillColor(i);*/
        
    h_deltaphi_tau_mu_full_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_full_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{hadron}+#pi^{0}(s)) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron}+#pi^{0}(s))"); h_deltaphi_tau_mu_full_tau_hadron[i]->Sumw2();
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_full_tau_hadron[i]->SetMarkerStyle(styles[i]);
    /*h_deltaphi_tau_mu_full_tau_hadron[i]->SetFillColor(i);*/
    
    h_deltaphi_tau_mu_tau_hadron_mueta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_mueta_" + tag[i]).c_str(),("#eta_{#mu} vs #Delta#phi(#tau_{#mu}, #tau_{hadron}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetYTitle("#eta_{#mu}");
    
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_deltaeta_" + tag[i]).c_str(),("#Delta(abs(#eta))_{#mu_#tau} vs #Delta#phi(#tau_{#mu}, #tau_{hadron}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetYTitle("#Delta(abs(#eta))_{#mu_#tau}");
    
    h_mueta_taueta[i] = new TH2F(("h_mueta_taueta_" + tag[i]).c_str(),("#tau_{hadron} #eta vs #tau_{#mu} #eta " + tag[i]).c_str(),18,-3,3, 18,-3,3);
    h_mueta_taueta[i]->SetXTitle("#tau_{hadron} #eta"); h_mueta_taueta[i]->SetYTitle("#tau_{#mu} #eta");
    
    h_PV_N[i] = new TH1F(("h_PV_N_" + tag[i]).c_str(),("h_PV_N_" + tag[i]).c_str(),5, -.5, 4.5);
    h_PV_N[i]->SetXTitle("N_{PV}"); h_PV_N[i]->Sumw2();
    /*if (i != 0) {*/h_PV_N[i]->SetLineColor(colors[i]); h_PV_N[i]->SetMarkerStyle(styles[i]);
    
    h_AP[i] = new TH2F(("h_AP_" + tag[i]).c_str(),("AP " + tag[i]).c_str(),40, -1, 1, 45,0,15);
    h_AP[i]->SetXTitle("p_{z}^{+} - p_{z}^{-} / p_{z}^{+} + p_{z}^{-}"); h_AP[i]->SetYTitle("#tau_{#mu} p_{T} [GeV]");
    
    // ZDC
    
    h_sumZDCplus[i] = new TH1F(("h_sumZDCplus_" + tag[i]).c_str(),("sum ZDC+ " + tag[i]).c_str(), 80, 0, 40000);
    h_sumZDCplus[i]->SetXTitle("sum ZDC+"); h_sumZDCplus[i]->Sumw2();
    /*if (i != 0) {*/h_sumZDCplus[i]->SetLineColor(colors[i]); //h_sumZDCplus[i]->SetMarkerStyle(styles[i]);
    
    h_sumZDCminus[i] = new TH1F(("h_sumZDCminus_" + tag[i]).c_str(),("sum ZDC- " + tag[i]).c_str(), 80, 0, 40000);
    h_sumZDCminus[i]->SetXTitle("sum ZDC-"); h_sumZDCminus[i]->Sumw2();
    /*if (i != 0) {*/h_sumZDCminus[i]->SetLineColor(colors[i]); //h_sumZDCminus[i]->SetMarkerStyle(styles[i]);
    
    h_sumZDC_pm[i] = new TH2F(("h_sumZDC_pm_" + tag[i]).c_str(),("sum ZDC+ vs ZDC- " + tag[i]).c_str(),40, 0, 20000, 40, 0, 20000);
    h_sumZDC_pm[i]->SetXTitle("sum ZDC-"); h_sumZDC_pm[i]->SetYTitle("sum ZDC+");
    
    h_averageZDCside_averageHFeta[i] = new TH2F(("h_averageZDCside_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs average ZDC side " + tag[i]).c_str(),21, -1, 1, 20, -5, 5);
    h_averageZDCside_averageHFeta[i]->SetXTitle("average ZDC side"); h_averageZDCside_averageHFeta[i]->SetYTitle("average HF #eta");
    
    h_muEta_averageZDCside[i] = new TH2F(("h_muEta_averageZDCside_" + tag[i]).c_str(),("average ZDC side vs muon #eta " + tag[i]).c_str(),19, -2.5, 2.5, 21, -1, 1);
    h_muEta_averageZDCside[i]->SetXTitle("muon #eta"); h_muEta_averageZDCside[i]->SetYTitle("average ZDC side");
    
    h_tauEta_averageZDCside[i] = new TH2F(("h_tauEta_averageZDCside_" + tag[i]).c_str(),("average ZDC side vs #tau_{hadron} #eta " + tag[i]).c_str(),19, -2.5, 2.5, 21, -1, 1);
    h_tauEta_averageZDCside[i]->SetXTitle("#tau_{hadron} #eta"); h_tauEta_averageZDCside[i]->SetYTitle("average ZDC side");
    
    //HF
    
    h_calo_energyHFp[i] = new TH1F(("h_calo_energyHFp_" + tag[i]).c_str(),("energy HF+ " + tag[i]).c_str(),40, -0.5, 19.5);
    h_calo_energyHFp[i]->SetXTitle("energy HF+ [GeV]"); h_calo_energyHFp[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFp[i]->SetLineColor(colors[i]); h_calo_energyHFp[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHFm[i] = new TH1F(("h_calo_energyHFm_" + tag[i]).c_str(),("energy HF- " + tag[i]).c_str(),40, -0.5, 19.5);
    h_calo_energyHFm[i]->SetXTitle("energy HF- [GeV]"); h_calo_energyHFm[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFm[i]->SetLineColor(colors[i]); h_calo_energyHFm[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFp[i] = new TH1F(("h_calo_leadingHFp_" + tag[i]).c_str(),("energy leading tower HF+ " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHFp[i]->SetXTitle("leading tower energy HF+ [GeV]"); h_calo_leadingHFp[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFp[i]->SetLineColor(colors[i]); h_calo_leadingHFp[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFm[i] = new TH1F(("h_calo_leadingHFm_" + tag[i]).c_str(),("energy leading tower HF- " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHFm[i]->SetXTitle("leading tower energy HF- [GeV]"); h_calo_leadingHFm[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFm[i]->SetLineColor(colors[i]); h_calo_leadingHFm[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHFp_sum[i] = new TH1F(("h_calo_energyHFp_sum_" + tag[i]).c_str(),("energy sum HF+ " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFp_sum[i]->SetXTitle("energy HF+ [GeV]"); h_calo_energyHFp_sum[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFp_sum[i]->SetLineColor(colors[i]); h_calo_energyHFp_sum[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHFm_sum[i] = new TH1F(("h_calo_energyHFm_sum_" + tag[i]).c_str(),("energy sum HF- " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFm_sum[i]->SetXTitle("energy HF- [GeV]"); h_calo_energyHFm_sum[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFm_sum[i]->SetLineColor(colors[i]); h_calo_energyHFm_sum[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHFp_size[i] = new TH1F(("h_calo_energyHFp_size_" + tag[i]).c_str(),("size HF+ " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFp_size[i]->SetXTitle("size HF+"); h_calo_energyHFp_size[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFp_size[i]->SetLineColor(colors[i]); h_calo_energyHFp_size[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHFm_size[i] = new TH1F(("h_calo_energyHFm_size_" + tag[i]).c_str(),("size HF- " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFm_size[i]->SetXTitle("size HF-"); h_calo_energyHFm_size[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFm_size[i]->SetLineColor(colors[i]); h_calo_energyHFm_size[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHF_pm[i] = new TH2F(("h_calo_energyHF_pm_" + tag[i]).c_str(),("total energy HF+ vs HF- " + tag[i]).c_str(),75, -0.5, 749.5, 75, -0.5, 749.5);
    h_calo_energyHF_pm[i]->SetXTitle("total energy HF-"); h_calo_energyHF_pm[i]->SetYTitle("total energy HF+");
    
    h_calo_leadingHF_pm[i] = new TH2F(("h_calo_leadingHF_pm_" + tag[i]).c_str(),("leading tower HF+ vs HF- " + tag[i]).c_str(),90, -0.5, 89.5, 90, -0.5, 89.5);
    h_calo_leadingHF_pm[i]->SetXTitle("leading tower HF-"); h_calo_leadingHF_pm[i]->SetYTitle("leading tower HF+");
    
    h_calo_HF_eta[i] = new TH1F(("h_calo_HF_eta_" + tag[i]).c_str(),("HF eta " + tag[i]).c_str(), 20, -5, 5);
    h_calo_HF_eta[i]->SetXTitle("#eta"); h_calo_HF_eta[i]->Sumw2();
    /*if (i != 0) {*/h_calo_HF_eta[i]->SetLineColor(colors[i]); h_calo_HF_eta[i]->SetMarkerStyle(styles[i]);
    
    h_calo_HF_energy_eta[i] = new TH2F(("h_calo_HF_energy_eta_" + tag[i]).c_str(),("energy HF vs #eta " + tag[i]).c_str(),20, -5, 5, 25, -0.5, 24.5);
    h_calo_HF_energy_eta[i]->SetXTitle("#eta"); h_calo_HF_energy_eta[i]->SetYTitle("energy HF");
    
    h_calo_muEta_leadingHFeta[i] = new TH2F(("h_calo_muEta_leadingHFeta_" + tag[i]).c_str(),("leading HF #eta vs muon #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_muEta_leadingHFeta[i]->SetXTitle("muon #eta"); h_calo_muEta_leadingHFeta[i]->SetYTitle("leading HF #eta");
    
    h_calo_tauEta_leadingHFeta[i] = new TH2F(("h_calo_tauEta_leadingHFeta_" + tag[i]).c_str(),("leading HF #eta vs #tau_{hadron} #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_tauEta_leadingHFeta[i]->SetXTitle("#tau_{hadron} #eta"); h_calo_tauEta_leadingHFeta[i]->SetYTitle("leading HF #eta");
    
    h_calo_muEta_averageHFeta[i] = new TH2F(("h_calo_muEta_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs muon #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_muEta_averageHFeta[i]->SetXTitle("muon #eta"); h_calo_muEta_averageHFeta[i]->SetYTitle("average HF #eta");
    
    h_calo_tauEta_averageHFeta[i] = new TH2F(("h_calo_tauEta_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs #tau_{hadron} #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_tauEta_averageHFeta[i]->SetXTitle("#tau_{hadron} #eta"); h_calo_tauEta_averageHFeta[i]->SetYTitle("average HF #eta");
    
    // MET
    
    h_MET[i] = new TH1F(("h_MET_" + tag[i]).c_str(),("MET " + tag[i]).c_str(),35, -0.5, 34.5);
    h_MET[i]->SetXTitle("MET"); h_MET[i]->Sumw2();
    /*if (i != 0) {*/h_MET[i]->SetLineColor(colors[i]); h_MET[i]->SetMarkerStyle(styles[i]);
    
  }

  TH2F *h_calo_energyHFp_nch = new TH2F("h_calo_energyHFp_nch","leading tower HF+ vs nch",6 , 2.5, 8.5, 12, 1, 5);
  h_calo_energyHFp_nch->SetXTitle("nch"); h_calo_energyHFp_nch->SetYTitle("leading tower HF+ [GeV]");
  
  TH2F *h_calo_energyHFm_nch = new TH2F("h_calo_energyHFm_nch","leading tower HF- vs nch",6 , 2.5, 8.5, 12, 1, 5);
  h_calo_energyHFm_nch->SetXTitle("nch"); h_calo_energyHFm_nch->SetYTitle("leading tower HF- [GeV]");

    
  TH1F *h_SB_deltaphi = new TH1F("h_SB_deltaphi","S/sqrt(S+B) for #Delta#phi(#tau_{#mu}, #tau_{hadron})",deltaphi_bins, 0, TMath::Pi());
  h_SB_deltaphi->SetXTitle("S/sqrt(S+B) for #Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_SB_deltaphi->Sumw2();
  h_SB_deltaphi->SetLineColor(colors[0]); h_SB_deltaphi->SetMarkerStyle(styles[1]);
  h_SB_deltaphi->SetBinErrorOption(TH1::kPoisson);
    
  // other
  TH1F *h_track_activity_pt     = new TH1F("h_track_activity_pt",     "h_track_activity_pt",     25, 0, 10); h_track_activity_pt->Sumw2(); h_track_activity_pt->SetXTitle("track p_{T} [GeV]");
  TH2F *h_track_activity_pt_eta = new TH2F("h_track_activity_pt_eta", "h_track_activity_pt_eta", 25, 0, 10, 10, -2.5, 2.5);
  h_track_activity_pt_eta->SetXTitle("track p_{T} [GeV]"); h_track_activity_pt_eta->SetYTitle("track #eta");
  
  TH2F *h_deltaphi_tau_mu_tau_hadron_nch = new TH2F("h_deltaphi_tau_mu_tau_hadron_nch", "nch vs #Delta#phi(#tau_{#mu}, #tau_{hadron})", 8, 0, TMath::Pi(), 5 , 0.5, 5.5);
  h_deltaphi_tau_mu_tau_hadron_nch->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron_nch->SetYTitle("nch");


// end of histogram declaration
// *******************************************

//  double mass1, mass2, mass3;
//  TFile *ntuple = new TFile("ntuple.root", "RECREATE");
//  TTree *aux;
//  aux = new TTree("tree", "tree");
//  aux->Branch("mass1", &mass1);
//  aux->Branch("mass2", &mass2);
//  aux->Branch("mass3", &mass3);

  double muon_mass = 105.6583 / 1000.;
  int entries = (TREE->fChain)->GetEntries();
  
  TLorentzVector tau_muon, tau_hadron, full_tau_hadron, gamma, tempGamma;
  int charge_counter[nSamples][3];
  for (int s = 0; s < nSamples; s++){
    for (int i = 0; i < 3; i++) charge_counter[s][i] = 0;
  }
  int pionLeadingIndex = -1;
  int pionSubLeadingIndex = -1;
  int pionSubSubLeadingIndex = -1;
  float pionLeadingPt = -1;
  float pionSubLeadingPt = -1;
  float pionSubSubLeadingPt = -1;
  float muonSF = 1;
  int muon_charge = 0;
  int tauh_charge = 0;
  float tau_weight = 1;
  bool passedNch = false;
  int candCounter = 0;
  float delta_phi = 0;
  float full_delta_phi = 0;
  cout << "Running on Data ..." << endl;
  for(int iEntry=0; iEntry<entries; iEntry++) {
    (TREE->fChain)->GetEntry(iEntry);
    
    if (!(iEntry%(entries/20))) cout << "\r" << int(100*iEntry/entries) << "%" << flush;
    TLorentzVector tau_hadron_pi1, tau_hadron_pi2, tau_hadron_pi3;
    
    bool triggered = TREE->triggered->at(0);
    tau_weight = 1;
    pionLeadingIndex = -1;
    pionSubLeadingIndex = -1;
    pionSubSubLeadingIndex = -1;
    pionLeadingPt = -1;
    pionSubLeadingPt = -1;
    pionSubSubLeadingPt = -1;
    bool passedmu  = false;
    bool passedtau = false;
    bool passedcalo = true;
    bool passedZDC = true;
    bool passedMET = true;
    bool passedGamma = true;
    passedNch = false;
    double temp_tau_pt_comparison = 0.;
    muon_charge = 0;
    double temp_nch = 0.;
    double temp_npv = 0.;
    double temp_pvz = 0.;
    double temp_pv_trk_1 = 0.;
    double temp_pv_trk_2 = 0.;
    double temp_pv_trk_3 = 0.;
    for (int i=0; i<(int)TREE->BsTauTau_mu1_pt->size(); i++) {
      if (TREE->BsTauTau_mu1_isSoft->at(i) != tau_muon_isSoft) { continue; }
      if (TREE->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }

      passedmu = true;
      if (TREE->BsTauTau_mu1_pt->at(i) > temp_tau_pt_comparison){// && (TREE->BsTauTau_nch->at(0) > 3 || TREE->BsTauTau_tau_q->at(0)*TREE->BsTauTau_mu1_q->at(i) == -1)) {
        temp_tau_pt_comparison = TREE->BsTauTau_mu1_pt->at(i);
        tau_muon.SetPtEtaPhiM (TREE->BsTauTau_mu1_pt->at(i), TREE->BsTauTau_mu1_eta->at(i), TREE->BsTauTau_mu1_phi->at(i), muon_mass);
        muon_charge = TREE->BsTauTau_mu1_q->at(i);
      }
    } // loop over the size of the muons
    if(TREE->BsTauTau_mu1_pt->size()>1) passedmu = false;

    
    temp_nch = TREE->BsTauTau_nch->at(0);
    //if (temp_nch != 3) continue;
    if (temp_nch>=min_nch && temp_nch<=max_nch) passedNch = true;
    
    for (int i=0; i<(int)TREE->reco_pion_pt->size(); i++) {
      if (TREE->reco_pion_pt->at(i) > pionLeadingPt && TREE->reco_pion_pt->at(i) > pionLeadingCut){
        pionLeadingPt = TREE->reco_pion_pt->at(i);
        pionLeadingIndex = i;
      }
    }
    if (pionLeadingIndex != -1){
      for (int i=0; i<(int)TREE->reco_pion_pt->size(); i++) {
        if (i == pionLeadingIndex) continue;
        if (TREE->reco_pion_pt->at(i) > pionSubLeadingPt && TREE->reco_pion_pt->at(i) > pionSubLeadingCut){
          pionSubLeadingPt = TREE->reco_pion_pt->at(i);
          pionSubLeadingIndex = i;
        }
      }
    }
    if (pionSubLeadingIndex != -1){
      for (int i=0; i<(int)TREE->reco_pion_pt->size(); i++) {
        if (i == pionLeadingIndex || i == pionSubLeadingIndex) continue;
        if (TREE->reco_pion_pt->at(i) > pionSubSubLeadingPt && TREE->reco_pion_pt->at(i) > pionSubSubLeadingCut){
          pionSubSubLeadingPt = TREE->reco_pion_pt->at(i);
          pionSubSubLeadingIndex = i;
        }
      }
    }
    if (pionSubSubLeadingIndex == -1) passedNch = false;
    else{
      for (int i=0; i<(int)TREE->reco_pion_pt->size(); i++) {
        if (i == pionLeadingIndex || i == pionSubLeadingIndex || i == pionSubSubLeadingIndex) continue;
        if (TREE->reco_pion_pt->at(i) > pionSubSubLeadingCut) passedNch = false;
      }
      if (passedNch) temp_nch = 3;
    }
    
    
    
    /*if (temp_nch == 3){
      pionLeadingIndex = 0;
      pionSubSubLeadingIndex = 0;
      for (int i=0; i<(int)TREE->reco_pion_pt->size(); i++) {
        int binx = xaxis_eff->FindBin(TREE->reco_pion_eta->at(i));
        int biny = yaxis_eff->FindBin(TREE->reco_pion_pt->at(i));
        float eff = pion_eff->GetBinContent(binx,biny);
        //if (eff) tau_weight /= eff;
        //cout << "binx: " << binx << " biny: " << biny << " eff: " << eff << " tau weight: " << tau_weight << endl;
        if (TREE->reco_pion_pt->at(i) > TREE->reco_pion_pt->at(pionLeadingIndex)) pionLeadingIndex = i;
        if (TREE->reco_pion_pt->at(i) < TREE->reco_pion_pt->at(pionSubSubLeadingIndex)) pionSubSubLeadingIndex = i;
      }
      for (int i=0; i<(int)TREE->reco_pion_pt->size(); i++) {
        if (i != pionLeadingIndex && i != pionSubSubLeadingIndex) pionSubLeadingIndex = i;
      }
    }*/
    
    
    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    candCounter = 0;
    float vprob = -1;
    float temp_pt = 0;
    double temp_tau_rho1 = 0.; double temp_tau_rho2 = 0.;
    for (int i=0; i<(int)TREE->BsTauTau_tau_pt->size(); i++) {
      if (muon_charge*TREE->BsTauTau_tau_q->at(i) != MuTauCharge) continue;
      //if (TREE->BsTauTau_tau_pt->at(i) < 2) continue;
      //if (temp_nch != 1 && TREE->BsTauTau_tau_pt->at(i) < 3.0) continue; //temporary commented
      //if(TMath::Abs(TREE->BsTauTau_tau_eta->at(i)) > 2) continue;
      candCounter += 1;
      if (TREE->BsTauTau_tau_pt->at(i) > temp_pt) {
        vprob = 100*TREE->BsTauTau_tau_vprob->at(i);
        if (temp_nch == 1) vprob = 100*TREE->BsTauTau_B_vprob->at(i);
        temp_pt = TREE->BsTauTau_tau_pt->at(i);
      }
      if (TREE->BsTauTau_tau_vprob->at(i) < tau_hadron_vertexprob) continue;
      //if (temp_nch == 1 && TREE->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      passedtau = true;
      if (TREE->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREE->BsTauTau_tau_pt->at(i);
        tauh_charge = TREE->BsTauTau_tau_q->at(i);
        //tau_hadron.SetPtEtaPhiM (TREE->BsTauTau_tau_pt->at(i), TREE->BsTauTau_tau_eta->at(i), TREE->BsTauTau_tau_phi->at(i), muon_mass);
        tau_hadron.SetPtEtaPhiM (TREE->BsTauTau_tau_pt->at(i), TREE->BsTauTau_tau_eta->at(i), TREE->BsTauTau_tau_phi->at(i), TREE->BsTauTau_tau_mass->at(i));
        //tau_hadron_pi1.SetPtEtaPhiM (TREE->BsTauTau_tau_pi1_pt->at(i), TREE->BsTauTau_tau_pi1_eta->at(i), TREE->BsTauTau_tau_pi1_phi->at(i), .13957018);
        //tau_hadron_pi2.SetPtEtaPhiM (TREE->BsTauTau_tau_pi2_pt->at(i), TREE->BsTauTau_tau_pi2_eta->at(i), TREE->BsTauTau_tau_pi2_phi->at(i), .13957018);
        //tau_hadron_pi3.SetPtEtaPhiM (TREE->BsTauTau_tau_pi3_pt->at(i), TREE->BsTauTau_tau_pi3_eta->at(i), TREE->BsTauTau_tau_pi3_phi->at(i), .13957018);
        if (threeProng[0]){
        temp_tau_rho1 = TREE->BsTauTau_tau_rhomass1->at(i); temp_tau_rho2 = TREE->BsTauTau_tau_rhomass2->at(i);
        temp_pv_trk_1 = TREE->BsTauTau_tau_pi1_z->at(i); temp_pv_trk_2 = TREE->BsTauTau_tau_pi2_z->at(i); temp_pv_trk_3 = TREE->BsTauTau_tau_pi3_z->at(i);
        }
      }
    } // loop over the size of the tau candidates
    
    //if ((tau_muon+tau_hadron).Pt() > MET_cut) passedMET = false;

    bool found_mu_track = false;
    bool found_pi1_track = false;
    bool found_pi2_track = false;
    bool found_pi3_track = false;
    bool high_activity = false;/*
    for (int i=0; i<(int)TREE->BsTauTau_trackPFactivity_pt->size(); i++) {
      if ( TMath::Abs(TREE->BsTauTau_trackPFactivity_pt->at(i) - tau_muon.Pt()) < 0.05 && TMath::Abs( TMath::Abs(TREE->BsTauTau_trackPFactivity_eta->at(i)) - TMath::Abs(tau_muon.Eta()) ) < 0.05)
        found_mu_track = true;
      if ( TMath::Abs(TREE->BsTauTau_trackPFactivity_pt->at(i) - tau_hadron_pi1.Pt()) < 0.05 && TMath::Abs( TMath::Abs(TREE->BsTauTau_trackPFactivity_eta->at(i)) - TMath::Abs(tau_hadron_pi1.Eta()) ) < 0.05)
        found_pi1_track = true;
      if ( TMath::Abs(TREE->BsTauTau_trackPFactivity_pt->at(i) - tau_hadron_pi2.Pt()) < 0.05 && TMath::Abs( TMath::Abs(TREE->BsTauTau_trackPFactivity_eta->at(i)) - TMath::Abs(tau_hadron_pi2.Eta()) ) < 0.05)
        found_pi2_track = true;
      if ( TMath::Abs(TREE->BsTauTau_trackPFactivity_pt->at(i) - tau_hadron_pi3.Pt()) < 0.05 && TMath::Abs( TMath::Abs(TREE->BsTauTau_trackPFactivity_eta->at(i)) - TMath::Abs(tau_hadron_pi3.Eta()) ) < 0.05)
        found_pi3_track = true;
      if (!(found_mu_track || found_pi1_track || found_pi2_track || found_pi3_track) &&  TMath::Abs(TREE->BsTauTau_trackPFactivity_eta->at(i))<2.5) {
        if (passedcalo) {
          h_track_activity_pt->Fill(TREE->BsTauTau_trackPFactivity_pt->at(i));
          h_track_activity_pt_eta->Fill(TREE->BsTauTau_trackPFactivity_pt->at(i), TREE->BsTauTau_trackPFactivity_eta->at(i));
        }
        if (TREE->BsTauTau_trackPFactivity_pt->at(i) > 0.5) {
          high_activity = true;
        }
      }
    } // check activity*/

    temp_pvz = TREE->BsTauTau_bbPV_vz->at(0);
    temp_npv = TREE->PV_N;
    double sumHFp = 0;
    double maxHFp = 0;
    int sizeHFp = 0;
    double sumHFm = 0;
    double maxHFm = 0;
    int sizeHFm = 0;
    double maxHF = 0;
    double leadingHFeta = 0;
    double averageHFeta = 0;
    h_tau_hadron_nPions[0]->Fill(TREE->BsTauTau_nPions->at(0));
    float sumZDCp = -1;
    float sumZDCm = -1;
    if (hasZDCinfo){
      sumZDCp = TREE->BsTauTau_calo_zdcSumPlus->at(0);
      sumZDCm = TREE->BsTauTau_calo_zdcSumMinus->at(0);
      if (sumZDCp > maxZDCp || sumZDCm > maxZDCm || sumZDCp < minZDCp || sumZDCm < minZDCm) passedZDC = false;
    }
    
#ifdef nch_cut
    if (triggered && passedmu && passedtau)
#else
    if (triggered && passedmu && passedtau)
#endif
    {
      for (int i=0; i<(int)TREE->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREE->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREE->BsTauTau_calo_energyHFm->at(i);
        double etaCal = TREE->BsTauTau_calo_eta->at(i);
        if (eHFp != -1){
          if (eHFp > maxHFp) maxHFp = eHFp;
          if (maxHFp > maxHF) {maxHF = maxHFp; leadingHFeta = etaCal;};
          averageHFeta += eHFp * etaCal;
          if (passedNch) h_calo_energyHFp[0]->Fill(eHFp); sizeHFp++; sumHFp += eHFp;
          if (passedNch) h_calo_HF_eta[0]->Fill(etaCal);
          if (passedNch) h_calo_HF_energy_eta[0]->Fill(etaCal,eHFp);
        }
        if (eHFm != -1){
          if (eHFm > maxHFm) maxHFm = eHFm;
          if (maxHFm > maxHF) {maxHF = maxHFm; leadingHFeta = etaCal;};
          averageHFeta += eHFm * etaCal;
          if (passedNch) h_calo_energyHFm[0]->Fill(eHFm); sizeHFm++; sumHFm += eHFm;
          if (passedNch) h_calo_HF_eta[0]->Fill(etaCal);
          if (passedNch) h_calo_HF_energy_eta[0]->Fill(etaCal,eHFm);
        }
      }
      averageHFeta /= (sumHFp+sumHFm);
      if (maxHFp > HFpLeading_high || maxHFm > HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcalo = false;
      if (passedcalo && passedNch) h_calo_energyHFp_sum[0]->Fill(sumHFp,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHFm_sum[0]->Fill(sumHFm,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHF_pm[0]->Fill(sumHFm,sumHFp,tau_weight);
      h_calo_energyHFp_nch->Fill(TREE->BsTauTau_nch->at(0),maxHFp); h_calo_energyHFm_nch->Fill(TREE->BsTauTau_nch->at(0),maxHFm,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHFp_size[0]->Fill(sizeHFp,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHFm_size[0]->Fill(sizeHFm,tau_weight);
      if (passedNch) h_calo_leadingHFp[0]->Fill(maxHFp,tau_weight);
      if (passedNch) h_calo_leadingHFm[0]->Fill(maxHFm,tau_weight);
      if (passedNch) h_calo_leadingHF_pm[0]->Fill(maxHFm,maxHFp,tau_weight);
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    }
    
#ifdef nch_cut
    if (passedmu && passedNch && passedcalo && passedMET && TREE->triggered->at(0))
#else
    if (passedmu && passedcalo && passedMET && TREE->triggered->at(0) && (temp_nch<min_nch || temp_nch>max_nch))
#endif
    {
    TLorentzVector recoPiZero;
    int Ngamma = 0;
    float minDeltaR = 1000;
    float minDeltaM = 100000; //MeV
    int NrecoPiZero = 0;
    if (!threeProng[0]){
    Ngamma = TREE->BsTauTau_nGammas->at(0);
    h_NrecoGamma[0]->Fill(Ngamma);
    full_tau_hadron = tau_hadron;
    for (int i=0; i<Ngamma; i++){
      gamma.SetPtEtaPhiM (TREE->reco_gamma_pt->at(i),TREE->reco_gamma_eta->at(i),TREE->reco_gamma_phi->at(i),0);
      //if (i == 0) recoNeutralPion = gamma;
      //else recoNeutralPion += gamma;
      h_recoGamma_pt[0]->Fill(TREE->reco_gamma_pt->at(i));
      h_recoGamma_eta[0]->Fill(TREE->reco_gamma_eta->at(i));
      if (TMath::Abs(TREE->reco_gamma_eta->at(i)) > 1.7) passedGamma = false;
      h_recoGamma_phi[0]->Fill(TREE->reco_gamma_phi->at(i));
      h_recoGamma_deltaphi_muon[0]->Fill(TMath::Abs(gamma.DeltaPhi(tau_muon)));
      h_recoGamma_deltaphi_pion[0]->Fill(TMath::Abs(gamma.DeltaPhi(tau_hadron)));
      for (int j=i+1; j<Ngamma; j++){
        tempGamma.SetPtEtaPhiM (TREE->reco_gamma_pt->at(j),TREE->reco_gamma_eta->at(j),TREE->reco_gamma_phi->at(j),0);
        float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
        if (deltaR < minDeltaR) minDeltaR = deltaR;
        recoPiZero = gamma+tempGamma;
        double pionMass = 1000*recoPiZero.M(); // MeV
        if (TMath::Abs(pionMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = pionMass-PiZeroMass;
        h_recoPiZeroDeltaM[0]->Fill(pionMass-PiZeroMass);
        if (TMath::Abs(pionMass-PiZeroMass) < 50){
          h_recoPiZero_pt[0]->Fill(recoPiZero.Pt());
          h_recoPiZero_eta[0]->Fill(recoPiZero.Eta());
          h_recoPiZero_phi[0]->Fill(recoPiZero.Phi());
          h_recoPiZero_deltaphi_muon[0]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_muon)));
          h_recoPiZero_deltaphi_pion[0]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)));
          NrecoPiZero += 1;
          full_tau_hadron += recoPiZero;
        }
      }
    }
    h_NrecoPiZero[0]->Fill(NrecoPiZero);
    h_recoGammasMinDeltaR[0]->Fill(minDeltaR);
    h_recoPiZeroMinDeltaM[0]->Fill(minDeltaM);
    } // if not three prong
    } // if in the signal region
    

    delta_phi = TMath::Abs(tau_muon.DeltaPhi(tau_hadron));
    full_delta_phi = TMath::Abs(tau_muon.DeltaPhi(full_tau_hadron));
    
    
    h_cutflow[0]->Fill(0); //input from Ntuplizer
    if (passedmu) h_cutflow[0]->Fill(1); //tau mu
    if (temp_nch>=min_nch && temp_nch<=max_nch) h_cutflow[0]->Fill(2); //nCh
    if (passedtau) h_cutflow[0]->Fill(3); //tau hadron
    if (passedcalo) h_cutflow[0]->Fill(4); //HF
    if (passedcalo && passedtau) h_cutflow[0]->Fill(5); //HF and tau hadron
    if (passedcalo && passedtau && temp_nch>=min_nch && temp_nch<=max_nch) h_cutflow[0]->Fill(6); //HF and tau hadron and nch
    //if () h_cutflow[0]->Fill(4); //vertex prob
    
    bool keepEvent = false;
#ifdef nch_cut
    if (passedmu && passedNch && passedcalo && passedMET && passedGamma && TREE->triggered->at(0))
#else
    if (passedmu && passedcalo && passedMET && passedGamma && TREE->triggered->at(0) && (temp_nch<min_nch || temp_nch>max_nch))
#endif
    {
      keepEvent = true;
    }
    if (keepEvent){
      h_tau_hadron_vprob[0]->Fill(vprob,tau_weight);
      if (!passedtau) keepEvent = false;
    }
    if (keepEvent){
      charge_counter[0][muon_charge*tauh_charge + 1] += 1;
      if (muon_charge*tauh_charge == 0) keepEvent = false;
    }
    if (keepEvent){
      h_sumZDCplus[0]->Fill(sumZDCp,tau_weight);
      h_sumZDCminus[0]->Fill(sumZDCm,tau_weight);
      h_sumZDC_pm[0]->Fill(sumZDCm,sumZDCp,tau_weight);
      float averageZDC = (sumZDCp-ratioZDCpm*sumZDCm)/(sumZDCp+ratioZDCpm*sumZDCm);
      h_muEta_averageZDCside[0]->Fill(tau_muon.Eta(),averageZDC);
      h_tauEta_averageZDCside[0]->Fill(tau_hadron.Eta(),averageZDC);
      h_averageZDCside_averageHFeta[0]->Fill(averageZDC,averageHFeta);
      if(!passedZDC) keepEvent = false;
    }
    if (keepEvent){
      h_deltaphi_tau_mu_tau_hadron[0]->Fill(delta_phi,tau_weight);
      h_deltaphi_tau_mu_full_tau_hadron[0]->Fill(full_delta_phi,tau_weight);
      h_deltaphi_tau_mu_tau_hadron_nch->Fill(delta_phi,temp_nch,tau_weight);
      h_deltaphi_tau_mu_tau_hadron_mueta[0]->Fill(delta_phi,tau_muon.Eta(),tau_weight);
      h_deltaphi_tau_mu_tau_hadron_deltaeta[0]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()),tau_weight);
      if(delta_phi < deltaPhi_cut) keepEvent = false;
    }
    if (keepEvent){
      h_tau_hadron_nch[0]->Fill(temp_nch,tau_weight);
      h_tau_hadron_ncand_final[0]->Fill(candCounter,tau_weight);
      if (!threeProng[0]) for (int i = 0; i < (int)TREE->reco_pion_ecalEnergy->size(); i++) h_reco_pion_energy_HCAL_ECAL[0]->Fill(TREE->reco_pion_ecalEnergy->at(i),TREE->reco_pion_hcalEnergy->at(i),tau_weight);
      h_tau_mu_p[0] ->Fill(tau_muon.P(),tau_weight);
      h_tau_mu_pz[0] ->Fill(tau_muon.Pz(),tau_weight);
      h_tau_mu_pt[0] ->Fill(tau_muon.Pt(),tau_weight);
      h_tau_mu_eta[0]->Fill(tau_muon.Eta(),tau_weight);
      h_tau_mu_phi[0]->Fill(tau_muon.Phi(),tau_weight);
      h_tau_hadron_p[0]->Fill(tau_hadron.P(),tau_weight);
      h_tau_hadron_pz[0]->Fill(tau_hadron.Pz(),tau_weight);
      h_tau_hadron_pt[0]->Fill(tau_hadron.Pt(),tau_weight);
      h_tau_hadron_eta[0]->Fill(tau_hadron.Eta(),tau_weight);
      h_tau_hadron_phi[0]->Fill(tau_hadron.Phi(),tau_weight);
      h_tau_hadron_rhomass[0][0]->Fill(temp_tau_rho1,tau_weight);
      h_tau_hadron_rhomass[0][1]->Fill(temp_tau_rho2,tau_weight);
      h_tau_hadron_rhomass2D[0]->Fill(temp_tau_rho1,temp_tau_rho2,tau_weight);
      h_tau_hadron_mass[0]->Fill(tau_hadron.M(),tau_weight);
      h_calo_muEta_leadingHFeta[0]->Fill(tau_muon.Eta(),leadingHFeta);
      h_calo_tauEta_leadingHFeta[0]->Fill(tau_hadron.Eta(),leadingHFeta);
      h_calo_muEta_averageHFeta[0]->Fill(tau_muon.Eta(),averageHFeta);
      h_calo_tauEta_averageHFeta[0]->Fill(tau_hadron.Eta(),averageHFeta);
      h_mueta_taueta[0]->Fill(tau_muon.Eta(),tau_hadron.Eta(),tau_weight);
      h_tau_hadron_track_pvz[0][0]->Fill(temp_pvz - temp_pv_trk_1);
      h_tau_hadron_track_pvz[0][1]->Fill(temp_pvz - temp_pv_trk_2);
      h_tau_hadron_track_pvz[0][2]->Fill(temp_pvz - temp_pv_trk_3);
      h_PV_N[0]->Fill(temp_npv);
      TLorentzVector MET;
      MET.SetPtEtaPhiM((tau_muon+tau_hadron).Pt(),(tau_muon+tau_hadron).Eta(),-(tau_muon+tau_hadron).Phi(),0);
      h_MET[0]->Fill(MET.Pt(),tau_weight);
      //h_ditau_mass[0]->Fill((tau_muon+tau_hadron+MET).M());
      h_ditau_mass[0]->Fill((tau_muon+tau_hadron).M(),tau_weight);
      h_ditau_p[0]->Fill((tau_muon+tau_hadron).P(),tau_weight);
      h_ditau_pz[0]->Fill((tau_muon+tau_hadron).Pz(),tau_weight);
      if (passedNch){
        h_pion_leading_pt[0]->Fill(TREE->reco_pion_pt->at(pionLeadingIndex),tau_weight);
        h_pion_subleading_pt[0]->Fill(TREE->reco_pion_pt->at(pionSubLeadingIndex),tau_weight);
        h_pion_subsubleading_pt[0]->Fill(TREE->reco_pion_pt->at(pionSubSubLeadingIndex),tau_weight);
        h_pion_leading_eta[0]->Fill(TREE->reco_pion_eta->at(pionLeadingIndex),tau_weight);
        h_pion_subleading_eta[0]->Fill(TREE->reco_pion_eta->at(pionSubLeadingIndex),tau_weight);
        h_pion_subsubleading_eta[0]->Fill(TREE->reco_pion_eta->at(pionSubSubLeadingIndex),tau_weight);
        h_pion_leading_phi[0]->Fill(TREE->reco_pion_phi->at(pionLeadingIndex),tau_weight);
        h_pion_subleading_phi[0]->Fill(TREE->reco_pion_phi->at(pionSubLeadingIndex),tau_weight);
        h_pion_subsubleading_phi[0]->Fill(TREE->reco_pion_phi->at(pionSubSubLeadingIndex),tau_weight);
      }
      //cout << TREE->BsTauTau_calo_zdcSumPlus->at(0) << "," << TREE->BsTauTau_calo_zdcSumMinus->at(0) << endl;
      
      h_AP[0]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2,tau_weight);
      //cout << (tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()) << " ---- " << tau_muon.Pt() << endl;
      //h_MET[0]->Fill(TREE->MET_sumEt->at(0));
      
      if (delta_phi < 2.75){// && delta_phi > 2.4){
        /*cout << "nch: " << temp_nch << endl;
        cout << "Delta Phi: " << delta_phi << endl;
        cout << "        mu Phi: " << tau_muon.Phi() << " - eta: " << tau_muon.Eta() << " - pT: " << tau_muon.Pt() << endl;
        cout << "tau hadron Phi: " << tau_hadron.Phi() << " - eta: " << tau_hadron.Eta() << " - pT: " << tau_hadron.Pt() << endl;
        cout << "      pi-1 Phi: " << tau_hadron_pi1.Phi() << " - eta: " << tau_hadron_pi1.Eta() << " - pT: " << tau_hadron_pi1.Pt() << endl;
        cout << "      pi-2 Phi: " << tau_hadron_pi2.Phi() << " - eta: " << tau_hadron_pi2.Eta() << " - pT: " << tau_hadron_pi2.Pt() << endl;
        cout << "      pi-3 Phi: " << tau_hadron_pi3.Phi() << " - eta: " << tau_hadron_pi3.Eta() << " - pT: " << tau_hadron_pi3.Pt() << endl;
        cout << "Met sum Et: " << (tau_muon+tau_hadron).Pt() << endl;
        cout << "HF+ total: " << sumHFp << " -- HF- total: " << sumHFm << endl;
        cout << "HF+ leading: " << maxHFp << " -- HF- leading: " << maxHFm << endl;
        cout << "**** **** **** ****" << endl;*/
        //cout << "tau hadron mass: " << tau_hadron.M() << endl;
      }
    } //if (keepEvent)
    
    if (passedmu && TREE->triggered->at(0) && passedtau && muon_charge*tauh_charge == -1){
      if (temp_nch >= 5 && !passedcalo) A_highNch_highHF[0]->Fill(delta_phi,tau_weight);
      if (temp_nch == 3 && !passedcalo) B_lowNch_highHF[0]->Fill(delta_phi,tau_weight);
      if (temp_nch >= 5 && passedcalo) C_highNch_lowHF[0]->Fill(delta_phi,tau_weight);
      if (temp_nch == 3 && passedcalo) D_lowNch_lowHF[0]->Fill(delta_phi,tau_weight);
    }

//      aux->Fill();

  } // loop over the entries
  
  cout << "\nMuon charge times hadronic tau charge for Data:\n -1: " << charge_counter[0][0] << "\n  0: " << charge_counter[0][1] << "\n +1: " << charge_counter[0][2] << endl;

  // *********************************************************
  // *********************************************************
  // MC STARTS HERE
  // *********************************************************
  // *********************************************************
  
  int entriesMC;
  //TLorentzVector tau_muon_mc, tau_hadron_mc;
  double temp_tau_rho1 = 0.; double temp_tau_rho2 = 0.;
  bool passedmu  = false;
  bool passedtau = false;
  bool passedcalo = true;
  bool passedGamma = true;
  double temp_tau_pt_comparison = 0.;
  double temp_nch = 0.;
  double temp_npv = 0.;
  double temp_pvz = 0.;
  int matchedpT_rank = 0;
  int acceptedGen = 0;
  float acceptedReco = 0;
  float acceptedRecoRaw = 0;
  bool matchedTauHad = false;
  int nMu3prong = 0;
  
  for (int s = 0; s < nSamples-1; s++){
  cout << "Starting " << tag[s+1] << endl;
  entriesMC = (TREEMCs[s]->fChain)->GetEntries();
  //if (s == 3) entriesMC /= 10;
  nMu3prong = cutflow[1]->GetBinContent(8);
  acceptedReco = 0;
  acceptedRecoRaw = 0;
  for(int iEntry=0; iEntry<entriesMC; iEntry++) {
    (TREEMCs[s]->fChain)->GetEntry(iEntry);
    
    if (!(iEntry%(entriesMC/20))) cout << "\r" << int(100*iEntry/entriesMC) << "%" << flush;
    //cout << "number of taus: " << (int)TREEMCs[s]->gen_tau_pt->size() << endl;
    //if (TREEMCs[s]->gen_tautau_to_mu3prong->at(0) == 0) continue;
    bool gen_tautau_to_mu3prong = TREEMCs[s]->gen_tautau_to_mu3prong->at(0);
    
    if (!threeProng[s+1]){
      bool gen_tautau_to_mu1prong = TREEMCs[s]->gen_tautau_to_mu1prong->at(0);
      bool gen_tautau_to_mu1prong_0npi = TREEMCs[s]->gen_tautau_to_mu1prong_0npi->at(0);
      bool gen_tautau_to_mu1prong_1npi = TREEMCs[s]->gen_tautau_to_mu1prong_1npi->at(0);
      bool gen_tautau_to_mu1prong_npis = TREEMCs[s]->gen_tautau_to_mu1prong_npis->at(0);
      bool passGen[] = {1,gen_tautau_to_mu1prong,gen_tautau_to_mu1prong_0npi,gen_tautau_to_mu1prong_1npi,gen_tautau_to_mu1prong_npis,gen_tautau_to_mu1prong};
      if (!passGen[s]) continue;
    }
    bool triggered = TREEMCs[s]->triggered->at(0);
    //if (s == 5) triggered = true;

    passedmu  = false;
    passedtau = false;
    passedcalo = true;
    bool passedMET = true;
    passedNch = false;
    temp_tau_pt_comparison = 0.;
    muon_charge = 0;
    muonSF = 1;
    temp_nch = 0.;
    temp_npv = 0.;
    temp_pvz = 0.;
    matchedpT_rank = 0;
    matchedTauHad = false;
    pionLeadingIndex = -1;
    pionSubLeadingIndex = -1;
    pionSubSubLeadingIndex = -1;
    pionLeadingPt = -1;
    pionSubLeadingPt = -1;
    pionSubSubLeadingPt = -1;
    
    for (int i=0; i<(int)TREEMCs[s]->BsTauTau_mu1_pt->size(); i++) {
      if (TREEMCs[s]->BsTauTau_mu1_isSoft->at(i) != tau_muon_isSoft) { continue; }
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }
      passedmu = true;
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_mu1_pt->at(i);
        tau_muon.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_mu1_pt->at(i), TREEMCs[s]->BsTauTau_mu1_eta->at(i), TREEMCs[s]->BsTauTau_mu1_phi->at(i), muon_mass);
        muon_charge = TREEMCs[s]->BsTauTau_mu1_q->at(i);
      }
    }
    
    muonSF = tnp_weight_trk_pbpb(0);
    float muonSFerror = TMath::Max(TMath::Abs(tnp_weight_trk_pbpb(-1)-tnp_weight_trk_pbpb(0)),TMath::Abs(tnp_weight_trk_pbpb(0)-tnp_weight_trk_pbpb(-2)));
    
    
    if (!TREEMCs[s]->BsTauTau_nch->empty()) temp_nch = TREEMCs[s]->BsTauTau_nch->at(0);
    if (temp_nch>=min_nch && temp_nch<=max_nch) passedNch = true;
    
    if (!TREEMCs[s]->BsTauTau_bbPV_vz->empty()) temp_pvz = TREEMCs[s]->BsTauTau_bbPV_vz->at(0);
    temp_npv = TREEMCs[s]->PV_N;
    h_tau_hadron_nPions[s+1]->Fill(TREEMCs[s]->BsTauTau_nPions->at(0));
    
    
    for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
      if (TREEMCs[s]->reco_pion_pt->at(i) > pionLeadingPt && TREEMCs[s]->reco_pion_pt->at(i) > pionLeadingCut){
        pionLeadingPt = TREEMCs[s]->reco_pion_pt->at(i);
        pionLeadingIndex = i;
      }
    }
    if (pionLeadingIndex != -1){
      for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
        if (i == pionLeadingIndex) continue;
        if (TREEMCs[s]->reco_pion_pt->at(i) > pionSubLeadingPt && TREEMCs[s]->reco_pion_pt->at(i) > pionSubLeadingCut){
          pionSubLeadingPt = TREEMCs[s]->reco_pion_pt->at(i);
          pionSubLeadingIndex = i;
        }
      }
    }
    if (pionSubLeadingIndex != -1){
      for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
        if (i == pionLeadingIndex || i == pionSubLeadingIndex) continue;
        if (TREEMCs[s]->reco_pion_pt->at(i) > pionSubSubLeadingPt && TREEMCs[s]->reco_pion_pt->at(i) > pionSubSubLeadingCut){
          pionSubSubLeadingPt = TREEMCs[s]->reco_pion_pt->at(i);
          pionSubSubLeadingIndex = i;
        }
      }
    }
    if (pionSubSubLeadingIndex == -1) passedNch = false;
    else{
      for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
        if (i == pionLeadingIndex || i == pionSubLeadingIndex || i == pionSubSubLeadingIndex) continue;
        if (TREEMCs[s]->reco_pion_pt->at(i) > pionSubSubLeadingCut) passedNch = false;
      }
      if (passedNch) temp_nch = 3;
    }
    
    /*if (temp_nch == 3){
      pionLeadingIndex = 0;
      pionSubSubLeadingIndex = 0;
      for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
        if (TREEMCs[s]->reco_pion_pt->at(i) > TREEMCs[s]->reco_pion_pt->at(pionLeadingIndex)) pionLeadingIndex = i;
        if (TREEMCs[s]->reco_pion_pt->at(i) < TREEMCs[s]->reco_pion_pt->at(pionSubSubLeadingIndex)) pionSubSubLeadingIndex = i;
      }
      for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
        if (i != pionLeadingIndex && i != pionSubSubLeadingIndex) pionSubLeadingIndex = i;
      }
    }*/
    
    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    tau_weight = 1;
    candCounter = 0;
    temp_tau_rho1 = 0.; temp_tau_rho2 = 0.;
    if (temp_nch == 3){
      for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++) {
        int binx = xaxis_eff->FindBin(TREEMCs[s]->reco_pion_eta->at(i));
        int biny = yaxis_eff->FindBin(TREEMCs[s]->reco_pion_pt->at(i));
        float eff = pion_eff->GetBinContent(binx,biny);
        //if (eff) tau_weight /= eff;
      }
    }
    
    //if(!TREEMCs[s]->BsTauTau_tau_pt) continue;
    float vprob = -1;
    float temp_pt = 0;
    for (int i=0; i<(int)TREEMCs[s]->BsTauTau_tau_pt->size(); i++) {
      if (muon_charge*TREEMCs[s]->BsTauTau_tau_q->at(i) != MuTauCharge) continue;
      //if (TREEMCs[s]->BsTauTau_tau_pt->at(i) < 2) continue;
      //if (temp_nch != 1 && TREEMCs[s]->BsTauTau_tau_pt->at(i) < 3.0) continue;
      //if(TMath::Abs(TREEMCs[s]->BsTauTau_tau_eta->at(i)) > 2) continue;
      candCounter += 1;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_pt) {
        vprob = 100*TREEMCs[s]->BsTauTau_tau_vprob->at(i);
        if (temp_nch == 1) vprob = 100*TREEMCs[s]->BsTauTau_B_vprob->at(i);
        temp_pt = TREEMCs[s]->BsTauTau_tau_pt->at(i);
      }
      if (TREEMCs[s]->BsTauTau_tau_vprob->at(i) < tau_hadron_vertexprob) continue;
      //if (temp_nch == 1 && TREEMCs[s]->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      passedtau = true;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) / TREEMCs[s]->BsTauTau_tau_matched_gentaupt->at(i) > 0.85) matchedpT_rank += 1;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        tauh_charge = TREEMCs[s]->BsTauTau_tau_q->at(i);
        //if (s == 0) tau_weight *= 0.985;
        if (s == 0) tau_weight *= TREEMCs[s]->BsTauTau_tau_rEff->at(i);
        // Mixing tau and muon SF
        tau_weight *= muonSF;
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        //tau_hadron.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_tau_pt->at(i), TREEMCs[s]->BsTauTau_tau_eta->at(i), TREEMCs[s]->BsTauTau_tau_phi->at(i), muon_mass);
        tau_hadron.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_tau_pt->at(i), TREEMCs[s]->BsTauTau_tau_eta->at(i), TREEMCs[s]->BsTauTau_tau_phi->at(i), TREEMCs[s]->BsTauTau_tau_mass->at(i));
        if (threeProng[s+1]){
          temp_tau_rho1 = TREEMCs[s]->BsTauTau_tau_rhomass1->at(i);
          if (temp_nch>=3) temp_tau_rho2 = TREEMCs[s]->BsTauTau_tau_rhomass2->at(i);
        }
      }
    } // loop over the size of the reco taus
    
    
    //if ((tau_muon+tau_hadron).Pt() > MET_cut) passedMET = false;
    
    double sumHFp = 0;
    double maxHFp = 0;
    int sizeHFp = 0;
    double sumHFm = 0;
    double maxHFm = 0;
    int sizeHFm = 0;
    double maxHF = 0;
    double leadingHFeta = 0;
    double averageHFeta = 0;

#ifdef nch_cut
    if (triggered && passedmu && passedtau && temp_nch>=min_nch && temp_nch<=max_nch)
#else
    if (triggered && passedmu && passedtau && (temp_nch<min_nch || temp_nch>max_nch))
#endif
    {
      for (int i=0; i<(int)TREEMCs[s]->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREEMCs[s]->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREEMCs[s]->BsTauTau_calo_energyHFm->at(i);
        double etaCal = TREEMCs[s]->BsTauTau_calo_eta->at(i);
        if (eHFp != -1){
          if (eHFp > maxHFp) maxHFp = eHFp;
          if (maxHFp > maxHF) {maxHF = maxHFp; leadingHFeta = etaCal;};
          averageHFeta += eHFp * etaCal;
          if (passedNch) h_calo_energyHFp[s+1]->Fill(eHFp); sizeHFp++; sumHFp += eHFp;
          if (passedNch) h_calo_HF_eta[s+1]->Fill(etaCal);
          if (passedNch) h_calo_HF_energy_eta[s+1]->Fill(etaCal,eHFp);        }
        if (eHFm != -1){
          if (eHFm > maxHFm) maxHFm = eHFm;
          if (maxHFm > maxHF) {maxHF = maxHFm; leadingHFeta = etaCal;};
          averageHFeta += eHFm * etaCal;
          if (passedNch) h_calo_energyHFm[s+1]->Fill(eHFm); sizeHFm++; sumHFm += eHFm;
          if (passedNch) h_calo_HF_eta[s+1]->Fill(etaCal);
          if (passedNch) h_calo_HF_energy_eta[s+1]->Fill(etaCal,eHFm);        }
      }
      averageHFeta /= (sumHFp+sumHFm);
      if (maxHFp > HFpLeading_high || maxHFm > HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcalo = false;
      if (passedcalo && passedNch) h_calo_energyHFp_sum[s+1]->Fill(sumHFp,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHFm_sum[s+1]->Fill(sumHFm,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHF_pm[s+1]->Fill(sumHFm,sumHFp,tau_weight);
      //h_calo_energyHFp_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),maxHFp); h_calo_energyHFm_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),maxHFm,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHFp_size[s+1]->Fill(sizeHFp,tau_weight);
      if (passedcalo && passedNch) h_calo_energyHFm_size[s+1]->Fill(sizeHFm,tau_weight);
      if (passedNch) h_calo_leadingHFp[s+1]->Fill(maxHFp,tau_weight);
      if (passedNch) h_calo_leadingHFm[s+1]->Fill(maxHFm,tau_weight);
      if (passedNch) h_calo_leadingHF_pm[s+1]->Fill(maxHFm,maxHFp,tau_weight);
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    }
    
    float gen_tau_hadron_pt = -1000;
    float gen_tau_hadron_visible_pt = 0;
    TLorentzVector gen_tau_hadron_visible;
    int nAccGenPions = 0;
    int nMatchedPions = 0;
    if (gen_tautau_to_mu3prong){
      for (int i=0; i<(int)TREEMCs[s]->gen_tau_daughter_pdgId->size(); i++) {
        if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==13){
          if((TREEMCs[s]->gen_tau_daughter_pt->at(i) < 3.5 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i)) < 1.2) || (TREEMCs[s]->gen_tau_daughter_pt->at(i) < 2.5 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i)) >= 1.2)) continue;
          if(TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))>2.4) continue;
        }
        if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<2.5){
        //if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<1){
          nAccGenPions++;
          h_N_genPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
          h_N_genPionEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i));
          h_N_genPionPtEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i),TREEMCs[s]->gen_tau_daughter_pt->at(i));
          if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<1) h_N_genCentralPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
          if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > 2) h_N_genHighPtPionEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i));
          gen_tau_hadron_visible_pt += TREEMCs[s]->gen_tau_daughter_pt->at(i);
          TLorentzVector temp_gen;
          temp_gen.SetPtEtaPhiM(TREEMCs[s]->gen_tau_daughter_pt->at(i),TREEMCs[s]->gen_tau_daughter_eta->at(i),TREEMCs[s]->gen_tau_daughter_phi->at(i),0.13957);
          gen_tau_hadron_visible += temp_gen;
          for (int j=0; j<(int)TREEMCs[s]->reco_pion_pt->size(); j++) {
            TLorentzVector temp_reco;
            temp_reco.SetPtEtaPhiM(TREEMCs[s]->reco_pion_pt->at(j),TREEMCs[s]->reco_pion_eta->at(j),TREEMCs[s]->reco_pion_phi->at(j),0.13957);
            if (TMath::Abs(temp_reco.Pt()/temp_gen.Pt()-1) < 0.07 && temp_gen.DeltaR(temp_reco) < 0.015){
              nMatchedPions++;
              h_resPt_matchedPionReco[s+1]->Fill(temp_reco.Pt()/temp_gen.Pt());
              h_resR_matchedPionReco[s+1]->Fill(temp_gen.DeltaR(temp_reco));
              h_resPhi_matchedPionReco[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_phi->at(i),TREEMCs[s]->reco_pion_phi->at(j));
              h_eff_matchedPionReco_genPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
              if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<1) h_eff_matchedCentralPionReco_genPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
              h_eff_matchedPionReco_genPionEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i));
              h_eff_matchedPionReco_genPionPtEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i),TREEMCs[s]->gen_tau_daughter_pt->at(i));
              h_eff_matchedPionReco_genPionPtEtaZoomed[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i),TREEMCs[s]->gen_tau_daughter_pt->at(i));
              if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > 2) h_eff_matchedHighPtPionReco_genPionEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i));
            }
          }
        }
      } // loop on gen daughters
      gen_tau_hadron_pt = TREEMCs[s]->gen_tau_pt->at(0);
      if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(1))==13) gen_tau_hadron_pt = TREEMCs[s]->gen_tau_pt->at(1);
      if (nAccGenPions == 3){
        h_N_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
        acceptedGen++;
      }
      if(nAccGenPions == 3 && nMatchedPions == 3 && temp_nch>=3){
        h_eff_3prongReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
        matchedTauHad = true;
      }
    }
    
    TLorentzVector recoPiZero, genPiZero;
    int Ngamma = 0;
    int NgenPiZero = 0;
    float minDeltaR = 1000;
    float minDeltaM = 100000; //MeV
    if (!threeProng[s+1]){
      for (int i=0; i<(int)TREEMCs[s]->gen_neutral_pion_daughter_pdgId->size(); i++) {
        if (TMath::Abs(TREEMCs[s]->gen_neutral_pion_daughter_pdgId->at(i))==22){
          gamma.SetPtEtaPhiM (TREEMCs[s]->gen_neutral_pion_daughter_pt->at(i),TREEMCs[s]->gen_neutral_pion_daughter_eta->at(i),TREEMCs[s]->gen_neutral_pion_daughter_phi->at(i),0);
          for (int j=i+1; j<(int)TREEMCs[s]->gen_neutral_pion_daughter_pdgId->size(); j++) {
            if (TMath::Abs(TREEMCs[s]->gen_neutral_pion_daughter_pdgId->at(j))==22){
              tempGamma.SetPtEtaPhiM (TREEMCs[s]->gen_neutral_pion_daughter_pt->at(j),TREEMCs[s]->gen_neutral_pion_daughter_eta->at(j),TREEMCs[s]->gen_neutral_pion_daughter_phi->at(j),0);
              float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
              if (deltaR < minDeltaR) minDeltaR = deltaR;
              genPiZero = gamma+tempGamma;
              double pionMass = 1000*genPiZero.M(); // MeV
              if (TMath::Abs(pionMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = pionMass-PiZeroMass;
              h_genPiZeroDeltaM[s+1]->Fill(pionMass-PiZeroMass);
              if (TMath::Abs(pionMass-PiZeroMass) < 50){
                h_genPiZero_pt[s+1]->Fill(genPiZero.Pt());
                h_genPiZero_eta[s+1]->Fill(genPiZero.Eta());
                h_genPiZero_phi[s+1]->Fill(genPiZero.Phi());
                h_genPiZero_deltaphi_muon[s+1]->Fill(TMath::Abs(genPiZero.DeltaPhi(tau_muon)));
                h_genPiZero_deltaphi_pion[s+1]->Fill(TMath::Abs(genPiZero.DeltaPhi(tau_hadron)));
                NgenPiZero += 1;
              }
            }
          }
          Ngamma += 1;
          h_genGamma_pt[s+1]->Fill(TREEMCs[s]->gen_neutral_pion_daughter_pt->at(i));
          h_genGamma_eta[s+1]->Fill(TREEMCs[s]->gen_neutral_pion_daughter_eta->at(i));
          h_genGamma_phi[s+1]->Fill(TREEMCs[s]->gen_neutral_pion_daughter_phi->at(i));
          h_genGamma_deltaphi_muon[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_muon)));
          h_genGamma_deltaphi_pion[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_hadron)));
        }
      }
      h_NgenGamma[s+1]->Fill(Ngamma);
      h_NgenPiZero[s+1]->Fill(NgenPiZero);
      h_genGammasMinDeltaR[s+1]->Fill(minDeltaR);
      h_genPiZeroMinDeltaM[s+1]->Fill(minDeltaM);
    }
    
    int NrecoPiZero = 0;
    minDeltaR = 1000;
    minDeltaM = 100000; //MeV
    if (!threeProng[s+1]){
    Ngamma = TREEMCs[s]->BsTauTau_nGammas->at(0);
    h_NrecoGamma[s+1]->Fill(Ngamma);
    full_tau_hadron = tau_hadron;
    for (int i=0; i<Ngamma; i++){
      gamma.SetPtEtaPhiM (TREEMCs[s]->reco_gamma_pt->at(i),TREEMCs[s]->reco_gamma_eta->at(i),TREEMCs[s]->reco_gamma_phi->at(i),0);
      //if (i == 0) recoNeutralPion = gamma;
      //else recoNeutralPion += gamma;
      h_recoGamma_pt[s+1]->Fill(TREEMCs[s]->reco_gamma_pt->at(i));
      h_recoGamma_eta[s+1]->Fill(TREEMCs[s]->reco_gamma_eta->at(i));
      //if (TMath::Abs(TREEMCs[s]->reco_gamma_eta->at(i)) > 1.7) passedGamma = false;
      h_recoGamma_phi[s+1]->Fill(TREEMCs[s]->reco_gamma_phi->at(i));
      h_recoGamma_deltaphi_muon[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_muon)));
      h_recoGamma_deltaphi_pion[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_hadron)));
      for (int j=i+1; j<Ngamma; j++){
        tempGamma.SetPtEtaPhiM (TREEMCs[s]->reco_gamma_pt->at(j),TREEMCs[s]->reco_gamma_eta->at(j),TREEMCs[s]->reco_gamma_phi->at(j),0);
        float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
        if (deltaR < minDeltaR) minDeltaR = deltaR;
        recoPiZero = gamma+tempGamma;
        double pionMass = 1000*recoPiZero.M(); // MeV
        if (TMath::Abs(pionMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = pionMass-PiZeroMass;
        h_recoPiZeroDeltaM[s+1]->Fill(pionMass-PiZeroMass);
        if (TMath::Abs(pionMass-PiZeroMass) < 50){
          h_recoPiZero_pt[s+1]->Fill(recoPiZero.Pt());
          h_recoPiZero_eta[s+1]->Fill(recoPiZero.Eta());
          h_recoPiZero_phi[s+1]->Fill(recoPiZero.Phi());
          h_recoPiZero_deltaphi_muon[s+1]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_muon)));
          h_recoPiZero_deltaphi_pion[s+1]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)));
          NrecoPiZero += 1;
          full_tau_hadron += recoPiZero;
        }
      }
    }
    h_NrecoPiZero[s+1]->Fill(NrecoPiZero);
    h_recoGammasMinDeltaR[s+1]->Fill(minDeltaR);
    h_recoPiZeroMinDeltaM[s+1]->Fill(minDeltaM);
    }
    
    //cout << passedmu << endl << (temp_nch>=min_nch && temp_nch<=max_nch) << endl << passedcalo << endl << passedMET << endl << " ---- " << endl;
    delta_phi = TMath::Abs(tau_muon.DeltaPhi(tau_hadron));
    full_delta_phi = TMath::Abs(tau_muon.DeltaPhi(full_tau_hadron));
    //float reco_pions_sum_pt = TREEMCs[s]->BsTauTau_tau_pi1_pt->at(0) + TREEMCs[s]->BsTauTau_tau_pi2_pt->at(0) + TREEMCs[s]->BsTauTau_tau_pi3_pt->at(0);
    
    /*
    h_cutflow->Fill(0); //input from Ntuplizer
    if () h_cutflow->Fill(1); //triggered
    if (passedmu) h_cutflow->Fill(2); //tau mu
    //if () h_cutflow->Fill(3); //tau pT
    //if () h_cutflow->Fill(4); //vertex prob
    if (temp_nch>=min_nch && temp_nch<=max_nch) h_cutflow->Fill(5); //nCh
    if (passedcalo) h_cutflow->Fill(6); //HF
    if (passedtau) h_cutflow->Fill(7); //tau hadron */
    
    //if (s==0) passedcalo = true;
    
    //if (s > 2) gen_tautau_to_mu3prong = true;
    
#ifdef nch_cut
    if (passedmu && passedNch && passedcalo && passedGamma && triggered /*&& gen_tautau_to_mu3prong*/)
#else
    if (passedmu && passedcalo && passedGamma && triggered /*&& gen_tautau_to_mu3prong*/ && (temp_nch<min_nch || temp_nch>max_nch))
#endif
    {
      //if (!TREEMCs[s]->BsTauTau_tau_isRight->at(0)) continue;
      h_tau_hadron_vprob[s+1]->Fill(vprob,tau_weight);
      if (!passedtau) continue;
      charge_counter[s+1][muon_charge*tauh_charge + 1] += 1;
      if (muon_charge*tauh_charge == 0) continue;
      if(gen_tautau_to_mu3prong && nAccGenPions == 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_eff_tauReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
      h_deltaphi_tau_mu_tau_hadron[s+1]->Fill(delta_phi,tau_weight);
      h_deltaphi_tau_mu_full_tau_hadron[s+1]->Fill(full_delta_phi,tau_weight);
      h_deltaphi_tau_mu_tau_hadron_mueta[s+1]->Fill(delta_phi,tau_muon.Eta(),tau_weight);
      h_deltaphi_tau_mu_tau_hadron_deltaeta[s+1]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()),tau_weight);
      if(delta_phi < deltaPhi_cut) continue;
      //if (!TREEMCs[s]->BsTauTau_nPions->empty()) temp_nch = TREEMCs[s]->BsTauTau_nPions->at(0);
      h_tau_hadron_nch[s+1]->Fill(temp_nch,tau_weight);
      //h_tau_hadron_nch[s+1]->Fill(temp_nch,tau_weight);
      h_tau_hadron_ncand_final[s+1]->Fill(candCounter,tau_weight);
      h_tau_hadron_matched_pt_index[s+1]->Fill(matchedpT_rank,tau_weight);
      h_tau_mu_p[s+1]->Fill(tau_muon.P(),tau_weight);
      h_tau_mu_pz[s+1]->Fill(tau_muon.Pz(),tau_weight);
      h_tau_mu_pt[s+1]->Fill(tau_muon.Pt(),tau_weight);
      h_tau_mu_eta[s+1]->Fill(tau_muon.Eta(),tau_weight);
      h_tau_mu_phi[s+1]->Fill(tau_muon.Phi(),tau_weight);
      h_tau_hadron_p[s+1]->Fill(tau_hadron.P(),tau_weight);
      h_tau_hadron_pz[s+1]->Fill(tau_hadron.Pz(),tau_weight);
      h_tau_hadron_pt[s+1]->Fill(tau_hadron.Pt(),tau_weight);
      h_gen_tau_hadron_visible_pt[s+1]->Fill(gen_tau_hadron_visible_pt,tau_weight);
      //if (gen_tautau_to_mu3prong/* && TREEMCs[s]->BsTauTau_tau_isRight->at(0)*/) h_resVisTauPt[s+1]->Fill(tau_hadron.Pt() - gen_tau_hadron_visible_pt);
      h_tau_hadron_eta[s+1]->Fill(tau_hadron.Eta(),tau_weight);
      h_tau_hadron_phi[s+1]->Fill(tau_hadron.Phi(),tau_weight);
      if (gen_tautau_to_mu3prong/* && !TREEMCs[s]->BsTauTau_tau_isRight->empty()*/){
        //if (!TREEMCs[s]->BsTauTau_tau_isRight->at(0)){
          h_resVisTauPt[s+1]->Fill(tau_hadron.Pt() - gen_tau_hadron_visible.Pt());
          h_resVisTauEta[s+1]->Fill(tau_hadron.Eta() - gen_tau_hadron_visible.Eta());
          h_resVisTauPhi[s+1]->Fill(tau_hadron.Phi() - gen_tau_hadron_visible.Phi());
        //}
      }
      h_tau_hadron_rhomass[s+1][0]->Fill(temp_tau_rho1,tau_weight);
      h_tau_hadron_rhomass[s+1][1]->Fill(temp_tau_rho2,tau_weight);
      h_tau_hadron_rhomass2D[s+1]->Fill(temp_tau_rho1,temp_tau_rho2,tau_weight);
      h_tau_hadron_mass[s+1]->Fill(tau_hadron.M(),tau_weight);
//      h_tau_hadron_track_pvz[s+1][0]->Fill(tau_track1.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][1]->Fill(tau_track2.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][2]->Fill(tau_track3.Z() - temp_pvz);
      h_calo_muEta_leadingHFeta[s+1]->Fill(tau_muon.Eta(),leadingHFeta);
      h_calo_tauEta_leadingHFeta[s+1]->Fill(tau_hadron.Eta(),leadingHFeta);
      h_calo_muEta_averageHFeta[s+1]->Fill(tau_muon.Eta(),averageHFeta);
      h_calo_tauEta_averageHFeta[s+1]->Fill(tau_hadron.Eta(),averageHFeta);
      h_mueta_taueta[s+1]->Fill(tau_muon.Eta(),tau_hadron.Eta(),tau_weight);
      h_PV_N[s+1]->Fill(temp_npv,tau_weight);
      TLorentzVector MET;
      MET.SetPtEtaPhiM((tau_muon+tau_hadron).Pt(),(tau_muon+tau_hadron).Eta(),-(tau_muon+tau_hadron).Phi(),0);
      h_MET[s+1]->Fill(MET.Pt(),tau_weight);
      //h_ditau_mass[s+1]->Fill((tau_muon+tau_hadron+MET).M(),tau_weight);
      h_ditau_mass[s+1]->Fill((tau_muon+tau_hadron).M(),tau_weight);
      h_ditau_p[s+1]->Fill((tau_muon+tau_hadron).P(),tau_weight);
      h_ditau_pz[s+1]->Fill((tau_muon+tau_hadron).Pz(),tau_weight);
      
      if (temp_nch == 3){
        h_pion_leading_pt[s+1]->Fill(TREEMCs[s]->reco_pion_pt->at(pionLeadingIndex),tau_weight);
        h_pion_subleading_pt[s+1]->Fill(TREEMCs[s]->reco_pion_pt->at(pionSubLeadingIndex),tau_weight);
        h_pion_subsubleading_pt[s+1]->Fill(TREEMCs[s]->reco_pion_pt->at(pionSubSubLeadingIndex),tau_weight);
        h_pion_leading_eta[s+1]->Fill(TREEMCs[s]->reco_pion_eta->at(pionLeadingIndex),tau_weight);
        h_pion_subleading_eta[s+1]->Fill(TREEMCs[s]->reco_pion_eta->at(pionSubLeadingIndex),tau_weight);
        h_pion_subsubleading_eta[s+1]->Fill(TREEMCs[s]->reco_pion_eta->at(pionSubSubLeadingIndex),tau_weight);
        h_pion_leading_phi[s+1]->Fill(TREEMCs[s]->reco_pion_phi->at(pionLeadingIndex),tau_weight);
        h_pion_subleading_phi[s+1]->Fill(TREEMCs[s]->reco_pion_phi->at(pionSubLeadingIndex),tau_weight);
        h_pion_subsubleading_phi[s+1]->Fill(TREEMCs[s]->reco_pion_phi->at(pionSubSubLeadingIndex),tau_weight);
      }
      
      h_AP[s+1]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2,tau_weight);
      //h_MET[s+1]->Fill(TREEMCs[s]->MET_sumEt->at(0));
      if (matchedTauHad){
        acceptedReco += tau_weight;
        acceptedRecoRaw++;
      }
    }
    
    
    if (passedmu && TREE->triggered->at(0) && passedtau && muon_charge*tauh_charge == -1){
      if (temp_nch >= 5 && !passedcalo) A_highNch_highHF[s+1]->Fill(delta_phi,tau_weight);
      if (temp_nch == 3 && !passedcalo) B_lowNch_highHF[s+1]->Fill(delta_phi,tau_weight);
      if (temp_nch >= 5 && passedcalo) C_highNch_lowHF[s+1]->Fill(delta_phi,tau_weight);
      if (temp_nch == 3 && passedcalo) D_lowNch_lowHF[s+1]->Fill(delta_phi,tau_weight);
    }
    
  } // loop on entries
  
  cout << "Muon charge times hadronic tau charge for " << tag[s+1] << ":\n -1: " << charge_counter[s+1][0] << "\n  0: " << charge_counter[s+1][1] << "\n +1: " << charge_counter[s+1][2] << endl;
  } // loop on samples

  // *******************************************
  // *******************************************
  // Plotting starts here
  // *******************************************
  // *******************************************
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0,1);
  
  
  string print = "Found ";
  for (int i = 0; i < nSamples; i++){
    if (i == nSamples-1) print += "and ";
    print += to_string(int(h_tau_mu_pt[i]->GetEntries())) + " " + tag[i];
    if (i != nSamples-1) print += ", ";
  }
  print += " events.\n";
  cout << print;
  cout << "Scale Factor / scaled events" << endl;
  for (int i = 0; i < nSamples; i++){
    if (i == nSamples-1) cout << "and ";
    cout << SF[i] << " / " << SF[i]*h_tau_mu_pt[i]->Integral();
    if (i != nSamples-1) cout << ", ";
  }
  cout << endl;

  for (int i = 0; i < nSamples; i++){
    //if (i != 0 && h_tau_hadron_nch[i]->GetMaximum()!=0) SF[i] = 0.8 * h_tau_hadron_nch[0]->GetMaximum()/h_tau_hadron_nch[i]->GetMaximum();
    //SF[i] = 1./h_tau_mu_pt[i]->Integral();
    //cout << "Temporary SF: " << SF[i] << endl;
    h_tau_mu_p[i]->Scale(SF[i]); h_tau_mu_pz[i]->Scale(SF[i]); h_tau_mu_pt[i]->Scale(SF[i]); h_tau_mu_eta[i]->Scale(SF[i]); h_tau_mu_phi[i]->Scale(SF[i]);
    h_tau_hadron_p[i]->Scale(SF[i]); h_tau_hadron_pz[i]->Scale(SF[i]); h_tau_hadron_pt[i]->Scale(SF[i]); h_gen_tau_hadron_visible_pt[i]->Scale(SF[i]);
    h_tau_hadron_eta[i]->Scale(SF[i]); h_tau_hadron_phi[i]->Scale(SF[i]);
    for (int j = 0; j < 2; j++) h_tau_hadron_rhomass[i][j]->Scale(SF[i]);
    h_tau_hadron_nch[i]->Scale(SF[i]);
    h_tau_hadron_mass[i]->Scale(SF[i]);
    //h_tau_hadron_mass[i]->Scale(1./h_tau_hadron_mass[i]->GetMaximum());
    h_tau_hadron_ncand_final[i]->Scale(SF[i]);
    h_tau_hadron_vprob[i]->Scale(1./h_tau_hadron_vprob[i]->GetMaximum());
    for (int j = 0; j < 3; j++) h_tau_hadron_track_pvz[i][j]->Scale(SF[i]);
    h_deltaphi_tau_mu_tau_hadron[i]->Scale(SF[i]);
    h_deltaphi_tau_mu_full_tau_hadron[i]->Scale(SF[i]);
    
    //if (i == 0) h_deltaphi_tau_mu_tau_hadron[i]->Scale(0.1);
    
    h_PV_N[i]->Scale(SF[i]);
    h_calo_energyHFp[i]->Scale(SF[i]);
    //h_calo_energyHFp_sum[i]->Scale(1./h_calo_energyHFp_sum[i]->GetMaximum());
    h_calo_energyHFp_sum[i]->Scale(SF[i]);
    h_sumZDCplus[i]->Scale(SF[i]);
    h_sumZDCminus[i]->Scale(SF[i]);
    h_calo_energyHFp_size[i]->Scale(1./h_calo_energyHFp_size[i]->GetMaximum());
    h_calo_energyHFm[i]->Scale(SF[i]);
    //h_calo_energyHFm_sum[i]->Scale(1./h_calo_energyHFm_sum[i]->GetMaximum());
    h_calo_energyHFm_sum[i]->Scale(SF[i]);
    h_calo_energyHFm_size[i]->Scale(1./h_calo_energyHFm_size[i]->GetMaximum());
    h_calo_HF_eta[i]->Scale(SF[i]);
    h_calo_leadingHFp[i]->Scale(SF[i]);
    h_calo_leadingHFm[i]->Scale(SF[i]);
    h_MET[i]->Scale(SF[i]);
    //h_ditau_mass[i]->Scale(SF[i]);
    h_ditau_p[i]->Scale(SF[i]);
    h_ditau_pz[i]->Scale(SF[i]);
    //h_ditau_mass[i]->Scale(1./h_ditau_mass[i]->GetMaximum());
    h_ditau_mass[i]->Scale(SF[i]);
    //for (int bin=1; bin <=4; bin++) cutflow[i]->SetBinContent(bin,SF[i]*cutflow[i]->GetBinContent(bin));
    h_resVisTauPt[i]->Scale(SF[i]);
    
    h_pion_leading_pt[i]->Scale(SF[i]);
    h_pion_subleading_pt[i]->Scale(SF[i]);
    h_pion_subsubleading_pt[i]->Scale(SF[i]);
    h_pion_leading_eta[i]->Scale(SF[i]);
    h_pion_subleading_eta[i]->Scale(SF[i]);
    h_pion_subsubleading_eta[i]->Scale(SF[i]);
    h_pion_leading_phi[i]->Scale(SF[i]);
    h_pion_subleading_phi[i]->Scale(SF[i]);
    h_pion_subsubleading_phi[i]->Scale(SF[i]);
    
    A_highNch_highHF[i]->Scale(SF[i]);
    B_lowNch_highHF[i]->Scale(SF[i]);
    C_highNch_lowHF[i]->Scale(SF[i]);
    D_lowNch_lowHF[i]->Scale(SF[i]);
  }
  
  // Cross section and errors
  
  float nA = A_highNch_highHF[0]->GetEntries();
  float nB = B_lowNch_highHF[0]->GetEntries();
  float nC = C_highNch_lowHF[0]->GetEntries();
  float countBackgroundABCD = nB * nC / nA;
  
  TH1F *background = (TH1F*)B_lowNch_highHF[0]->Clone("background");
  background->SetTitle("Estimated background from ABCD");
  background->Multiply(C_highNch_lowHF[0]);
  background->Divide(A_highNch_highHF[0]);
  
  float analysisEff = acceptedReco*1.0/acceptedGen;
  //float analysisAccEff = acceptedReco*1.0/nMu3prong;
  float analysisAccEff = acceptedReco*1.0/cutflow[1]->GetBinContent(1);
  cout << "#gen tau-tau events in all decay modes (with pT cut): " << cutflow[1]->GetBinContent(1) << endl;
  
  float BR = 2 * 0.1741 * 0.1457; // errors: 0.0004 and 0.0007
  int nData = h_tau_mu_pt[0]->GetEntries();
  double nRecoErr;
  //cout << "Integral and error: " << h_MET[1]->IntegralAndError(1,h_MET[1]->GetNbinsX(),nRecoErr,"");
  //cout << " +- " << nRecoErr << endl;
  float sigmaFiducial = (nData - countBackgroundABCD) / (dataLumi*analysisEff*BR);
  
  float countBackErr2 = countBackgroundABCD*countBackgroundABCD*(1./nA + 1./nB + 1./nC);
  
  cout << endl << "#events in" << endl << "A_highNch_highHF: " << nA << endl << "B_lowNch_highHF: " << nB << endl << "C_highNch_lowHF: " << nC << endl << "D_lowNch_lowHF: " << D_lowNch_lowHF[0]->GetEntries() << endl << "Estimated number of background with ABCD count: " << countBackgroundABCD << " +- " << countBackErr2 << endl << endl;
  cout << "Estimated number of background events from ABCD shape: " << background->Integral() << endl;
  float relEstatFid = sqrt( (nData+countBackErr2)/((nData-countBackgroundABCD)*(nData-countBackgroundABCD)) + 1./acceptedRecoRaw + 1./acceptedGen );
  float statisticalErrorOnFiducial = sigmaFiducial * relEstatFid;
  
  cout << "\nmu+3prong: "<< nMu3prong << " , accepted gen: " << acceptedGen << " , reconstructed: " << acceptedReco << " (" << acceptedRecoRaw << " raw) , efficiency: " << analysisEff << " , acceptance x efficiency: " << analysisAccEff << endl;
  
  cout << endl << " **** **** **** **** " << endl;
  cout << "Fiducial cross section: " << sigmaFiducial << " ub +- " << statisticalErrorOnFiducial << " (" << 100*statisticalErrorOnFiducial/sigmaFiducial << "%) (stat)" << endl;
  
  float sigmaInclusive = (nData - countBackgroundABCD) / (dataLumi*analysisAccEff*pT_SF);
  float relEstatInc = sqrt( (nData+countBackErr2)/((nData-countBackgroundABCD)*(nData-countBackgroundABCD)) + 1./acceptedRecoRaw + 1./cutflow[1]->GetBinContent(1));// +  1./above3GeVSuperChic + 1./nAllEventsSuperChic);
  float statisticalErrorOnInclusive = sigmaInclusive * relEstatInc;
  
  cout << "Inclusive cross section: " << sigmaInclusive << " ub +- " << statisticalErrorOnInclusive << " (" << 100*statisticalErrorOnInclusive/sigmaInclusive << "%) (stat)" << endl;
  cout << " **** **** **** **** " << endl << endl;
  //cout << "Inclusive cross section (Method 2): " << (nData - countBackgroundABCD) * crossSectionMC[0] / h_tau_mu_pt[1]->Integral() << " nb" << endl;
  

  
  TCanvas *c_tau_had_pv = new TCanvas("c_tau_had_pv", "c_tau_had_pv", 1500, 500); c_tau_had_pv->Divide(3,1);
  for (int j = 0; j < 3; j++){
    for (int i = 1; i < nSamples; i++) hs_tau_hadron_track_pvz[j]->Add(h_tau_hadron_track_pvz[i][j]);
    c_tau_had_pv->cd(j+1);
    //hs_tau_hadron_track_pvz[j]->Draw("he");
    h_tau_hadron_track_pvz[0][j]->Draw("ex0");
    //c_tau_had_pv->cd(j+1); h_tau_hadron_track_pvz[0][j]->Draw("ex0"); h_tau_hadron_track_pvz[0][j]->GetYaxis()->SetRangeUser(0., 0.6);
    for (int i = 1; i < nSamples; i++){
      h_tau_hadron_track_pvz[i][j]->Draw("hesame");
      h_tau_hadron_track_pvz[i][j]->GetYaxis()->SetRangeUser(0., 0.1);
    }
  }
#ifdef nch_cut
  c_tau_had_pv->SaveAs((basePlotDir+"/plots/nch_cut/tau_had_pv."+plotFormat).c_str());
#else
  c_tau_had_pv->SaveAs((basePlotDir+"/plots/tau_had_pv."+plotFormat).c_str());
#endif
    
  TCanvas *c_delta_phi = new TCanvas("c_delta_phi", "c_delta_phi", 800, 800);
  //for (int i = 1; i < nSamples; i++) hs_deltaphi_tau_mu_tau_hadron->Add(h_deltaphi_tau_mu_tau_hadron[i]);
  //if (hs_deltaphi_tau_mu_tau_hadron->GetMaximum() < h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum()) hs_deltaphi_tau_mu_tau_hadron->SetMaximum(h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum());
  //c_delta_phi->cd(1); hs_deltaphi_tau_mu_tau_hadron->Draw("he"); h_deltaphi_tau_mu_tau_hadron[0]->Draw("e1same");
  for (int i = 1; i < nSamples; i++){
    hs_deltaphi_tau_mu_tau_hadron->Add(h_deltaphi_tau_mu_tau_hadron[i]);
    if (h_deltaphi_tau_mu_tau_hadron[i]->GetMaximum() > h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum()) h_deltaphi_tau_mu_tau_hadron[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron[i]->GetMaximum());
    if (h_deltaphi_tau_mu_tau_hadron[i]->GetMinimum() < h_deltaphi_tau_mu_tau_hadron[0]->GetMinimum() && h_deltaphi_tau_mu_tau_hadron[i]->GetMinimum()  > 0) h_deltaphi_tau_mu_tau_hadron[0]->SetMinimum(0.8*h_deltaphi_tau_mu_tau_hadron[i]->GetMinimum());
    h_deltaphi_tau_mu_tau_hadron[i]->SetLineWidth(3);
  }
  
  h_deltaphi_tau_mu_tau_hadron[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum());
  h_deltaphi_tau_mu_tau_hadron[0]->SetMinimum(0.005);
  
  string binSize = to_string(deltaphi_bins).substr(0, 4);
  string t = "Events / (#pi/";
  h_deltaphi_tau_mu_tau_hadron[0]->GetYaxis()->SetTitle((t+binSize+")").c_str());
  h_deltaphi_tau_mu_tau_hadron[0]->Draw("ex0");
  //c_delta_phi->SetGrid();
  //c_delta_phi->SetLogy(1);
  
  for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_tau_hadron[i]->Draw("hesame"); 
  //legend = new TLegend(0.5,0.6,0.75,0.88);
  TLegend *legend = new TLegend(0.75,0.80,0.9,0.93);
  legend->SetFillStyle(0);
  legend->SetFillColor(0); legend->SetLineColor(0); legend->SetShadowColor(0); legend->SetTextSize(0.03);
  //legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    (tag[0]).c_str(), "ep");
  legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data 2015", "ep");
  for (int i = 1; i < nSamples; i++) legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tag[i]).c_str(), "l");
  //legend->Draw();
  
  CMS_lumi( c_delta_phi, iPeriod, iPos );
  
  background->SetMarkerStyle(kFullTriangleDown);
  background->SetLineWidth(3);
  background->SetLineColor(colors[nSamples]);
  background->Draw("hesame");
  
  TLegend *tempGend = new TLegend(0.75,0.80,0.9,0.93);
  tempGend->SetFillStyle(0);
  tempGend->SetFillColor(0); tempGend->SetLineColor(0); tempGend->SetShadowColor(0); tempGend->SetTextSize(0.03);
  //tempGend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    (tag[0]).c_str(), "ep");
  tempGend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data 2015", "ep");
  for (int i = 1; i < nSamples; i++) tempGend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tag[i]).c_str(), "l");
  tempGend->AddEntry(background,    "background", "lp");
  tempGend->Draw();
  
  c_delta_phi->SaveAs((basePlotDir+"/singlePlot/delta_phi."+plotFormat).c_str());
  
#ifdef nch_cut
  c_delta_phi->SaveAs((basePlotDir+"/plots/nch_cut/delta_phi."+plotFormat).c_str());
#else
  c_delta_phi->SaveAs((basePlotDir+"/plots/delta_phi."+plotFormat).c_str());
#endif
  
  TCanvas *c_SB = new TCanvas("c_SB", "c_SB", 800, 800);
  for (int i = 1; i < deltaphi_bins+1; i++){
    double S_error;
    double signal = D_lowNch_lowHF[1]->IntegralAndError(i, deltaphi_bins, S_error, "");
    //cout << "signal: " << signal << " +- " << S_error << endl;
    double B_error;
    double bckgrd = background->IntegralAndError(i, deltaphi_bins, B_error, "");
    //cout << "background: " << bckgrd << " +- " << B_error << endl;
    float S_B = signal / sqrt(signal + bckgrd);
    float S_B_error = S_B * sqrt(pow(S_error/signal,2) + 0.25*((pow(S_error,2)+pow(B_error,2))/pow(signal+bckgrd,2)));
    //cout << "signal / sqrt(signal+background): " << S_B << " +- " << S_B_error << endl;
    h_SB_deltaphi->SetBinContent(i,S_B);
    h_SB_deltaphi->SetBinError(i,S_B_error);
  }
  h_SB_deltaphi->SetMaximum(1.2*h_SB_deltaphi->GetMaximum());
  h_SB_deltaphi->Draw("he");
  //c_SB->SetLogy(1);
  CMS_lumi( c_SB, iPeriod, iPos );
  
  c_SB->SaveAs((basePlotDir+"/singlePlot/S_sqrt_S+B_delta_phi."+plotFormat).c_str());
#ifdef nch_cut
  c_SB->SaveAs((basePlotDir+"/plots/nch_cut/S_sqrt_S+B_delta_phi."+plotFormat).c_str());
#else
  c_SB->SaveAs((basePlotDir+"/plots/S_sqrt_S+B_delta_phi."+plotFormat).c_str());
#endif
  
  TCanvas *c_full_delta_phi = new TCanvas("c_full_delta_phi", "c_full_delta_phi", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_deltaphi_tau_mu_full_tau_hadron[i]->GetMaximum() > h_deltaphi_tau_mu_full_tau_hadron[0]->GetMaximum()) h_deltaphi_tau_mu_full_tau_hadron[0]->SetMaximum(1.2*h_deltaphi_tau_mu_full_tau_hadron[i]->GetMaximum());
    if (h_deltaphi_tau_mu_full_tau_hadron[i]->GetMinimum() < h_deltaphi_tau_mu_full_tau_hadron[0]->GetMinimum() && h_deltaphi_tau_mu_full_tau_hadron[i]->GetMinimum()  > 0) h_deltaphi_tau_mu_full_tau_hadron[0]->SetMinimum(0.8*h_deltaphi_tau_mu_full_tau_hadron[i]->GetMinimum());
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetLineWidth(3);
  }
  h_deltaphi_tau_mu_full_tau_hadron[0]->Draw("ex0");
  //c_full_delta_phi->SetGrid();
  c_full_delta_phi->SetLogy(1);
  
  for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_full_tau_hadron[i]->Draw("hesame"); 
  legend->Draw();
#ifdef nch_cut
  c_full_delta_phi->SaveAs((basePlotDir+"/plots/nch_cut/full_delta_phi."+plotFormat).c_str());
#else
  c_full_delta_phi->SaveAs((basePlotDir+"/plots/full_delta_phi."+plotFormat).c_str());
#endif
  
  TCanvas *c_ABCD = new TCanvas("c_ABCD", "c_ABCD", 1600, 2400); c_ABCD->Divide(2,3);
  
  c_ABCD->cd(1); for (int i = 1; i < nSamples; i++){
    if (A_highNch_highHF[i]->GetMaximum() > A_highNch_highHF[0]->GetMaximum()) A_highNch_highHF[0]->SetMaximum(1.2*A_highNch_highHF[i]->GetMaximum());
    if (A_highNch_highHF[i]->GetMinimum() < A_highNch_highHF[0]->GetMinimum() && A_highNch_highHF[i]->GetMinimum() > 0) A_highNch_highHF[0]->SetMinimum(0.8*A_highNch_highHF[i]->GetMinimum());
  }
  A_highNch_highHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) A_highNch_highHF[i]->Draw("hesame");
  
  TCanvas *c_A_highNch_highHF = new TCanvas("c_A_highNch_highHF", "c_A_highNch_highHF", 800, 800);
  A_highNch_highHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) A_highNch_highHF[i]->Draw("hesame");
  c_A_highNch_highHF->SetLogy(1);
  CMS_lumi( c_A_highNch_highHF, iPeriod, iPos );
  c_A_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/A_highNch_highHF."+plotFormat).c_str());
  
  c_ABCD->cd(2); for (int i = 1; i < nSamples; i++){
    if (B_lowNch_highHF[i]->GetMaximum() > B_lowNch_highHF[0]->GetMaximum()) B_lowNch_highHF[0]->SetMaximum(1.2*B_lowNch_highHF[i]->GetMaximum());
    if (B_lowNch_highHF[i]->GetMinimum() < B_lowNch_highHF[0]->GetMinimum() && B_lowNch_highHF[i]->GetMinimum() > 0) B_lowNch_highHF[0]->SetMinimum(0.8*B_lowNch_highHF[i]->GetMinimum());
  }
  B_lowNch_highHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) B_lowNch_highHF[i]->Draw("hesame");
  
  TCanvas *c_B_lowNch_highHF = new TCanvas("c_B_lowNch_highHF", "c_B_lowNch_highHF", 800, 800);
  B_lowNch_highHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) B_lowNch_highHF[i]->Draw("hesame");
  c_B_lowNch_highHF->SetLogy(1);
  CMS_lumi( c_B_lowNch_highHF, iPeriod, iPos );
  c_B_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/B_lowNch_highHF."+plotFormat).c_str());
  
  c_ABCD->cd(3); for (int i = 1; i < nSamples; i++){
    if (C_highNch_lowHF[i]->GetMaximum() > C_highNch_lowHF[0]->GetMaximum()) C_highNch_lowHF[0]->SetMaximum(1.2*C_highNch_lowHF[i]->GetMaximum());
    if (C_highNch_lowHF[i]->GetMinimum() < C_highNch_lowHF[0]->GetMinimum() && C_highNch_lowHF[i]->GetMinimum() > 0) C_highNch_lowHF[0]->SetMinimum(0.8*C_highNch_lowHF[i]->GetMinimum());
  }
  C_highNch_lowHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) C_highNch_lowHF[i]->Draw("hesame");
  
  TCanvas *c_C_highNch_lowHF = new TCanvas("c_C_highNch_lowHF", "c_C_highNch_lowHF", 800, 800);
  C_highNch_lowHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) C_highNch_lowHF[i]->Draw("hesame");
  c_C_highNch_lowHF->SetLogy(1);
  CMS_lumi( c_C_highNch_lowHF, iPeriod, iPos );
  c_C_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/C_highNch_lowHF."+plotFormat).c_str());
  
  c_ABCD->cd(4); for (int i = 1; i < nSamples; i++){
    if (D_lowNch_lowHF[i]->GetMaximum() > D_lowNch_lowHF[0]->GetMaximum()) D_lowNch_lowHF[0]->SetMaximum(1.2*D_lowNch_lowHF[i]->GetMaximum());
    if (D_lowNch_lowHF[i]->GetMinimum() < D_lowNch_lowHF[0]->GetMinimum() && D_lowNch_lowHF[i]->GetMinimum() > 0) D_lowNch_lowHF[0]->SetMinimum(0.8*D_lowNch_lowHF[i]->GetMinimum());
  }
  D_lowNch_lowHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) D_lowNch_lowHF[i]->Draw("hesame");
  
  TCanvas *c_D_lowNch_lowHF = new TCanvas("c_D_lowNch_lowHF", "c_D_lowNch_lowHF", 800, 800);
  D_lowNch_lowHF[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) D_lowNch_lowHF[i]->Draw("hesame");
  c_D_lowNch_lowHF->SetLogy(1);
  CMS_lumi( c_D_lowNch_lowHF, iPeriod, iPos );
  c_D_lowNch_lowHF->SaveAs((basePlotDir+"/singlePlot/D_lowNch_lowHF."+plotFormat).c_str());
  
  c_ABCD->cd(5); background->Draw("ex0");
  
#ifdef nch_cut
  c_ABCD->SaveAs((basePlotDir+"/plots/nch_cut/ABCD."+plotFormat).c_str());
#else
  c_ABCD->SaveAs((basePlotDir+"/plots/ABCD."+plotFormat).c_str());
#endif
  
  TCanvas *c_PV = new TCanvas("c_PV", "c_PV", 1000, 500); c_PV->Divide(2,1);
  c_PV->cd(1); h_deltaphi_tau_mu_tau_hadron_nch->Draw("COLZ");
  for (int i = 1; i < nSamples; i++){
    if (h_PV_N[i]->GetMaximum() > h_PV_N[0]->GetMaximum()) h_PV_N[0]->SetMaximum(1.2*h_PV_N[i]->GetMaximum());
  }
  h_PV_N[0]->GetYaxis()->SetRangeUser(0.1, 1.2*h_PV_N[0]->GetMaximum());
  c_PV->cd(2); h_PV_N[0]->Draw("ex0");
  for (int i = 1; i < nSamples; i++) h_PV_N[i]->Draw("hesame"); legend->Draw();
  c_PV->cd(2)->SetLogy(1);
#ifdef nch_cut
  c_PV->SaveAs((basePlotDir+"/plots/nch_cut/PV."+plotFormat).c_str());
#else
  c_PV->SaveAs((basePlotDir+"/plots/PV."+plotFormat).c_str());
#endif
  
  int nGenPlots = 15;
  TCanvas *c_neutral_pion = new TCanvas("c_neutral_pion", "c_neutral_pion", nGenPlots*600, 1000); c_neutral_pion->Divide(nGenPlots,2);
  c_neutral_pion->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_genGamma_pt[i]->GetMaximum() > h_genGamma_pt[1]->GetMaximum()) h_genGamma_pt[1]->SetMaximum(1.2*h_genGamma_pt[i]->GetMaximum());
    h_genGamma_pt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_genGamma_eta[i]->GetMaximum() > h_genGamma_eta[1]->GetMaximum()) h_genGamma_eta[1]->SetMaximum(1.2*h_genGamma_eta[i]->GetMaximum());
    if (h_genGamma_eta[i]->GetMinimum() < h_genGamma_eta[1]->GetMinimum()) h_genGamma_eta[1]->SetMinimum(0.8*h_genGamma_eta[i]->GetMinimum());
    h_genGamma_eta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(3); for (int i = 1; i < nSamples; i++){
    if (h_genGamma_phi[i]->GetMaximum() > h_genGamma_phi[1]->GetMaximum()) h_genGamma_phi[1]->SetMaximum(1.2*h_genGamma_phi[i]->GetMaximum());
    if (h_genGamma_phi[i]->GetMinimum() < h_genGamma_phi[1]->GetMinimum()) h_genGamma_phi[1]->SetMinimum(0.8*h_genGamma_phi[i]->GetMinimum());
    h_genGamma_phi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(4); for (int i = 1; i < nSamples; i++){
    if (h_genGamma_deltaphi_muon[i]->GetMaximum() > h_genGamma_deltaphi_muon[1]->GetMaximum()) h_genGamma_deltaphi_muon[1]->SetMaximum(1.2*h_genGamma_deltaphi_muon[i]->GetMaximum());
    if (h_genGamma_deltaphi_muon[i]->GetMinimum() < h_genGamma_deltaphi_muon[1]->GetMinimum() && h_genGamma_deltaphi_muon[i]->GetMinimum() > 0) h_genGamma_deltaphi_muon[1]->SetMinimum(0.8*h_genGamma_deltaphi_muon[i]->GetMinimum());
    h_genGamma_deltaphi_muon[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(4)->SetLogy(1);
  
  c_neutral_pion->cd(5); for (int i = 1; i < nSamples; i++){
    if (h_genGamma_deltaphi_pion[i]->GetMaximum() > h_genGamma_deltaphi_pion[1]->GetMaximum()) h_genGamma_deltaphi_pion[1]->SetMaximum(1.2*h_genGamma_deltaphi_pion[i]->GetMaximum());
    if (h_genGamma_deltaphi_pion[i]->GetMinimum() < h_genGamma_deltaphi_pion[1]->GetMinimum() && h_genGamma_deltaphi_pion[i]->GetMinimum() > 0) h_genGamma_deltaphi_pion[1]->SetMinimum(0.8*h_genGamma_deltaphi_pion[i]->GetMinimum());
    h_genGamma_deltaphi_pion[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(5)->SetLogy(1);
  
  c_neutral_pion->cd(6); for (int i = 1; i < nSamples; i++){
    if (h_NgenGamma[i]->GetMaximum() > h_NgenGamma[1]->GetMaximum()) h_NgenGamma[1]->SetMaximum(1.2*h_NgenGamma[i]->GetMaximum());
    h_NgenGamma[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(6)->SetLogy(1);
  
  c_neutral_pion->cd(7); for (int i = 1; i < nSamples; i++){
    if (h_genGammasMinDeltaR[i]->GetMaximum() > h_genGammasMinDeltaR[1]->GetMaximum()) h_genGammasMinDeltaR[1]->SetMaximum(1.2*h_genGammasMinDeltaR[i]->GetMaximum());
    h_genGammasMinDeltaR[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(8); for (int i = 1; i < nSamples; i++){
    if (h_genPiZeroMinDeltaM[i]->GetMaximum() > h_genPiZeroMinDeltaM[1]->GetMaximum()) h_genPiZeroMinDeltaM[1]->SetMaximum(1.2*h_genPiZeroMinDeltaM[i]->GetMaximum());
    h_genPiZeroMinDeltaM[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(9); for (int i = 1; i < nSamples; i++){
    if (h_genPiZeroDeltaM[i]->GetMaximum() > h_genPiZeroDeltaM[1]->GetMaximum()) h_genPiZeroDeltaM[1]->SetMaximum(1.2*h_genPiZeroDeltaM[i]->GetMaximum());
    h_genPiZeroDeltaM[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(10); for (int i = 1; i < nSamples; i++){
    if (h_genPiZero_pt[i]->GetMaximum() > h_genPiZero_pt[1]->GetMaximum()) h_genPiZero_pt[1]->SetMaximum(1.2*h_genPiZero_pt[i]->GetMaximum());
    h_genPiZero_pt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(11); for (int i = 1; i < nSamples; i++){
    if (h_genPiZero_eta[i]->GetMaximum() > h_genPiZero_eta[1]->GetMaximum()) h_genPiZero_eta[1]->SetMaximum(1.2*h_genPiZero_eta[i]->GetMaximum());
    if (h_genPiZero_eta[i]->GetMinimum() < h_genPiZero_eta[1]->GetMinimum()) h_genPiZero_eta[1]->SetMinimum(0.8*h_genPiZero_eta[i]->GetMinimum());
    h_genPiZero_eta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(12); for (int i = 1; i < nSamples; i++){
    if (h_genPiZero_phi[i]->GetMaximum() > h_genPiZero_phi[1]->GetMaximum()) h_genPiZero_phi[1]->SetMaximum(1.2*h_genPiZero_phi[i]->GetMaximum());
    if (h_genPiZero_phi[i]->GetMinimum() < h_genPiZero_phi[1]->GetMinimum()) h_genPiZero_phi[1]->SetMinimum(0.8*h_genPiZero_phi[i]->GetMinimum());
    h_genPiZero_phi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(13); for (int i = 1; i < nSamples; i++){
    if (h_genPiZero_deltaphi_muon[i]->GetMaximum() > h_genPiZero_deltaphi_muon[1]->GetMaximum()) h_genPiZero_deltaphi_muon[1]->SetMaximum(1.2*h_genPiZero_deltaphi_muon[i]->GetMaximum());
    if (h_genPiZero_deltaphi_muon[i]->GetMinimum() < h_genPiZero_deltaphi_muon[1]->GetMinimum() && h_genPiZero_deltaphi_muon[i]->GetMinimum() > 0) h_genPiZero_deltaphi_muon[1]->SetMinimum(0.8*h_genPiZero_deltaphi_muon[i]->GetMinimum());
    h_genPiZero_deltaphi_muon[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(13)->SetLogy(1);
  
  c_neutral_pion->cd(14); for (int i = 1; i < nSamples; i++){
    if (h_genPiZero_deltaphi_pion[i]->GetMaximum() > h_genPiZero_deltaphi_pion[1]->GetMaximum()) h_genPiZero_deltaphi_pion[1]->SetMaximum(1.2*h_genPiZero_deltaphi_pion[i]->GetMaximum());
    if (h_genPiZero_deltaphi_pion[i]->GetMinimum() < h_genPiZero_deltaphi_pion[1]->GetMinimum() && h_genPiZero_deltaphi_pion[i]->GetMinimum() > 0) h_genPiZero_deltaphi_pion[1]->SetMinimum(0.8*h_genPiZero_deltaphi_pion[i]->GetMinimum());
    h_genPiZero_deltaphi_pion[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(14)->SetLogy(1);
  
  c_neutral_pion->cd(15); for (int i = 1; i < nSamples; i++){
    if (h_NgenPiZero[i]->GetMaximum() > h_NgenPiZero[1]->GetMaximum()) h_NgenPiZero[1]->SetMaximum(1.2*h_NgenPiZero[i]->GetMaximum());
    h_NgenPiZero[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(15)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+1); h_recoGamma_pt[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_pt[i]->GetMaximum() > h_recoGamma_pt[0]->GetMaximum()) h_recoGamma_pt[0]->SetMaximum(1.2*h_recoGamma_pt[i]->GetMaximum());
    h_recoGamma_pt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+2); h_recoGamma_eta[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_eta[i]->GetMaximum() > h_recoGamma_eta[0]->GetMaximum()) h_recoGamma_eta[0]->SetMaximum(1.2*h_recoGamma_eta[i]->GetMaximum());
    if (h_recoGamma_eta[i]->GetMinimum() < h_recoGamma_eta[0]->GetMinimum()) h_recoGamma_eta[0]->SetMinimum(0.8*h_recoGamma_eta[i]->GetMinimum());
    h_recoGamma_eta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+3); h_recoGamma_phi[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_phi[i]->GetMaximum() > h_recoGamma_phi[0]->GetMaximum()) h_recoGamma_phi[0]->SetMaximum(1.2*h_recoGamma_phi[i]->GetMaximum());
    if (h_recoGamma_phi[i]->GetMinimum() < h_recoGamma_phi[0]->GetMinimum()) h_recoGamma_phi[0]->SetMinimum(0.8*h_recoGamma_phi[i]->GetMinimum());
    h_recoGamma_phi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+4); h_recoGamma_deltaphi_muon[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_deltaphi_muon[i]->GetMaximum() > h_recoGamma_deltaphi_muon[0]->GetMaximum()) h_recoGamma_deltaphi_muon[0]->SetMaximum(1.2*h_recoGamma_deltaphi_muon[i]->GetMaximum());
    if (h_recoGamma_deltaphi_muon[i]->GetMinimum() < h_recoGamma_deltaphi_muon[0]->GetMinimum() && h_recoGamma_deltaphi_muon[i]->GetMinimum() > 0) h_recoGamma_deltaphi_muon[0]->SetMinimum(0.8*h_recoGamma_deltaphi_muon[i]->GetMinimum());
    h_recoGamma_deltaphi_muon[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+4)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+5); h_recoGamma_deltaphi_pion[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_deltaphi_pion[i]->GetMaximum() > h_recoGamma_deltaphi_pion[0]->GetMaximum()) h_recoGamma_deltaphi_pion[0]->SetMaximum(1.2*h_recoGamma_deltaphi_pion[i]->GetMaximum());
    if (h_recoGamma_deltaphi_pion[i]->GetMinimum() < h_recoGamma_deltaphi_pion[0]->GetMinimum() && h_recoGamma_deltaphi_pion[i]->GetMinimum() > 0) h_recoGamma_deltaphi_pion[0]->SetMinimum(0.8*h_recoGamma_deltaphi_pion[i]->GetMinimum());
    h_recoGamma_deltaphi_pion[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+5)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+6); h_NrecoGamma[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_NrecoGamma[i]->GetMaximum() > h_NrecoGamma[0]->GetMaximum()) h_NrecoGamma[0]->SetMaximum(1.2*h_NrecoGamma[i]->GetMaximum());
    h_NrecoGamma[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+6)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+7); h_recoGammasMinDeltaR[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGammasMinDeltaR[i]->GetMaximum() > h_recoGammasMinDeltaR[0]->GetMaximum()) h_recoGammasMinDeltaR[0]->SetMaximum(1.2*h_recoGammasMinDeltaR[i]->GetMaximum());
    h_recoGammasMinDeltaR[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+8); h_recoPiZeroMinDeltaM[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZeroMinDeltaM[i]->GetMaximum() > h_recoPiZeroMinDeltaM[0]->GetMaximum()) h_recoPiZeroMinDeltaM[0]->SetMaximum(1.2*h_recoPiZeroMinDeltaM[i]->GetMaximum());
    h_recoPiZeroMinDeltaM[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+9); h_recoPiZeroDeltaM[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZeroDeltaM[i]->GetMaximum() > h_recoPiZeroDeltaM[0]->GetMaximum()) h_recoPiZeroDeltaM[0]->SetMaximum(1.2*h_recoPiZeroDeltaM[i]->GetMaximum());
    h_recoPiZeroDeltaM[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+10); h_recoPiZero_pt[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_pt[i]->GetMaximum() > h_recoPiZero_pt[0]->GetMaximum()) h_recoPiZero_pt[0]->SetMaximum(1.2*h_recoPiZero_pt[i]->GetMaximum());
    h_recoPiZero_pt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+11); h_recoPiZero_eta[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_eta[i]->GetMaximum() > h_recoPiZero_eta[0]->GetMaximum()) h_recoPiZero_eta[0]->SetMaximum(1.2*h_recoPiZero_eta[i]->GetMaximum());
    if (h_recoPiZero_eta[i]->GetMinimum() < h_recoPiZero_eta[0]->GetMinimum()) h_recoPiZero_eta[0]->SetMinimum(0.8*h_recoPiZero_eta[i]->GetMinimum());
    h_recoPiZero_eta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+12); h_recoPiZero_phi[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_phi[i]->GetMaximum() > h_recoPiZero_phi[0]->GetMaximum()) h_recoPiZero_phi[0]->SetMaximum(1.2*h_recoPiZero_phi[i]->GetMaximum());
    if (h_recoPiZero_phi[i]->GetMinimum() < h_recoPiZero_phi[0]->GetMinimum()) h_recoPiZero_phi[0]->SetMinimum(0.8*h_recoPiZero_phi[i]->GetMinimum());
    h_recoPiZero_phi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+13); h_recoPiZero_deltaphi_muon[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_deltaphi_muon[i]->GetMaximum() > h_recoPiZero_deltaphi_muon[0]->GetMaximum()) h_recoPiZero_deltaphi_muon[0]->SetMaximum(1.2*h_recoPiZero_deltaphi_muon[i]->GetMaximum());
    if (h_recoPiZero_deltaphi_muon[i]->GetMinimum() < h_recoPiZero_deltaphi_muon[0]->GetMinimum() && h_recoPiZero_deltaphi_muon[i]->GetMinimum() > 0) h_recoPiZero_deltaphi_muon[0]->SetMinimum(0.8*h_recoPiZero_deltaphi_muon[i]->GetMinimum());
    h_recoPiZero_deltaphi_muon[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+13)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+14); h_recoPiZero_deltaphi_pion[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_deltaphi_pion[i]->GetMaximum() > h_recoPiZero_deltaphi_pion[0]->GetMaximum()) h_recoPiZero_deltaphi_pion[0]->SetMaximum(1.2*h_recoPiZero_deltaphi_pion[i]->GetMaximum());
    if (h_recoPiZero_deltaphi_pion[i]->GetMinimum() < h_recoPiZero_deltaphi_pion[0]->GetMinimum() && h_recoPiZero_deltaphi_pion[i]->GetMinimum() > 0) h_recoPiZero_deltaphi_pion[0]->SetMinimum(0.8*h_recoPiZero_deltaphi_pion[i]->GetMinimum());
    h_recoPiZero_deltaphi_pion[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+14)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+15); h_NrecoPiZero[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++){
    if (h_NrecoPiZero[i]->GetMaximum() > h_NrecoPiZero[0]->GetMaximum()) h_NrecoPiZero[0]->SetMaximum(1.2*h_NrecoPiZero[i]->GetMaximum());
    h_NrecoPiZero[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+15)->SetLogy(1);
  
#ifdef nch_cut
  c_neutral_pion->SaveAs((basePlotDir+"/plots/nch_cut/neutral_pion."+plotFormat).c_str());
#else
  c_neutral_pion->SaveAs((basePlotDir+"/plots/neutral_pion."+plotFormat).c_str());
#endif


  TCanvas *c_efficiency = new TCanvas("c_efficiency", "c_efficiency", 1500, 1000+nSamples*500); c_efficiency->Divide(3,2+nSamples);
  c_efficiency->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_N_genTauHadpt[i]->GetMaximum() > h_N_genTauHadpt[1]->GetMaximum()) h_N_genTauHadpt[1]->SetMaximum(1.2*h_N_genTauHadpt[i]->GetMaximum());
    h_N_genTauHadpt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_efficiency->cd(2); for (int i = 1; i < nSamples; i++){
    h_eff_3prongReco_genTauHadpt[i]->Divide(h_N_genTauHadpt[i]);
    if (h_eff_3prongReco_genTauHadpt[i]->GetMaximum() > h_eff_3prongReco_genTauHadpt[1]->GetMaximum()) h_eff_3prongReco_genTauHadpt[1]->SetMaximum(1.2*h_eff_3prongReco_genTauHadpt[i]->GetMaximum());
    h_eff_3prongReco_genTauHadpt[i]->Draw("hist same");
  }
  //legend->Draw();
  
  c_efficiency->cd(3); for (int i = 1; i < nSamples; i++){
    h_eff_tauReco_genTauHadpt[i]->Divide(h_N_genTauHadpt[i]);
    if (h_eff_tauReco_genTauHadpt[i]->GetMaximum() > h_eff_tauReco_genTauHadpt[1]->GetMaximum()) h_eff_tauReco_genTauHadpt[1]->SetMaximum(1.2*h_eff_tauReco_genTauHadpt[i]->GetMaximum());
    h_eff_tauReco_genTauHadpt[i]->Draw("hist same");
  }
  legend->Draw();
  
  c_efficiency->cd(4); for (int i = 1; i < nSamples; i++){
    if (h_N_genPionPt[i]->GetMaximum() > h_N_genPionPt[1]->GetMaximum()) h_N_genPionPt[1]->SetMaximum(1.2*h_N_genPionPt[i]->GetMaximum());
    h_N_genPionPt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_efficiency->cd(5); for (int i = 1; i < nSamples; i++){
    h_eff_matchedPionReco_genPionPt[i]->Divide(h_N_genPionPt[i]);
    if (h_eff_matchedPionReco_genPionPt[i]->GetMaximum() > h_eff_matchedPionReco_genPionPt[1]->GetMaximum()) h_eff_matchedPionReco_genPionPt[1]->SetMaximum(1.2*h_eff_matchedPionReco_genPionPt[i]->GetMaximum());
    h_eff_matchedPionReco_genPionPt[i]->Draw("hist same");
  }
  //legend->Draw();
  
  c_efficiency->cd(6); for (int i = 1; i < nSamples; i++){
    h_eff_matchedCentralPionReco_genPionPt[i]->Divide(h_N_genCentralPionPt[i]);
    if (h_eff_matchedCentralPionReco_genPionPt[i]->GetMaximum() > h_eff_matchedCentralPionReco_genPionPt[1]->GetMaximum()) h_eff_matchedCentralPionReco_genPionPt[1]->SetMaximum(1.2*h_eff_matchedCentralPionReco_genPionPt[i]->GetMaximum());
    h_eff_matchedCentralPionReco_genPionPt[i]->Draw("hist same");
  }
  //legend->Draw();
  
  c_efficiency->cd(7); for (int i = 1; i < nSamples; i++){
    if (h_N_genPionEta[i]->GetMaximum() > h_N_genPionEta[1]->GetMaximum()) h_N_genPionEta[1]->SetMaximum(1.2*h_N_genPionEta[i]->GetMaximum());
    h_N_genPionEta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_efficiency->cd(8); for (int i = 1; i < nSamples; i++){
    h_eff_matchedPionReco_genPionEta[i]->Divide(h_N_genPionEta[i]);
    if (h_eff_matchedPionReco_genPionEta[i]->GetMaximum() > h_eff_matchedPionReco_genPionEta[1]->GetMaximum()) h_eff_matchedPionReco_genPionEta[1]->SetMaximum(1.2*h_eff_matchedPionReco_genPionEta[i]->GetMaximum());
    h_eff_matchedPionReco_genPionEta[i]->Draw("hist same");
  }
  legend->Draw();
  
  c_efficiency->cd(9); for (int i = 1; i < nSamples; i++){
    h_eff_matchedHighPtPionReco_genPionEta[i]->Divide(h_N_genHighPtPionEta[i]);
    if (h_eff_matchedHighPtPionReco_genPionEta[i]->GetMaximum() > h_eff_matchedHighPtPionReco_genPionEta[1]->GetMaximum()) h_eff_matchedHighPtPionReco_genPionEta[1]->SetMaximum(1.2*h_eff_matchedHighPtPionReco_genPionEta[i]->GetMaximum());
    h_eff_matchedHighPtPionReco_genPionEta[i]->Draw("hist same");
  }
  //legend->Draw();
  
  
  for (int i = 1; i < nSamples; i++){
    c_efficiency->cd(7+3*i); 
    h_N_genPionPtEta[i]->Draw("colz");
    c_efficiency->cd(8+3*i); 
    h_eff_matchedPionReco_genPionPtEta[i]->Divide(h_N_genPionPtEta[i]);
    h_eff_matchedPionReco_genPionPtEta[i]->Draw("colz");
    c_efficiency->cd(9+3*i);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->Divide(h_N_genPionPtEta[i]);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->GetZaxis()->SetRangeUser(0.9, 1.1);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->Draw("colz");
  }
  
#ifdef nch_cut
  c_efficiency->SaveAs((basePlotDir+"/plots/nch_cut/efficiency."+plotFormat).c_str());
#else
  c_efficiency->SaveAs((basePlotDir+"/plots/efficiency."+plotFormat).c_str());
#endif


  TCanvas *c_tau_muon_kinematics = new TCanvas("c_tau_muon_kinematics", "c_tau_muon_kinematics", 1500, 1000); c_tau_muon_kinematics->Divide(3,2);
  
  c_tau_muon_kinematics->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_tau_mu_p[i]->GetMaximum() > h_tau_mu_p[0]->GetMaximum()) h_tau_mu_p[0]->SetMaximum(1.2*h_tau_mu_p[i]->GetMaximum());
    if (h_tau_mu_p[i]->GetMinimum() < h_tau_mu_p[0]->GetMinimum() && h_tau_mu_p[i]->GetMinimum() > 0) h_tau_mu_p[0]->SetMinimum(0.8*h_tau_mu_p[i]->GetMinimum());
  }
  h_tau_mu_p[0]->Draw("ex0"); legend->Draw();
  c_tau_muon_kinematics->cd(1)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_mu_p[i]->Draw("hesame");
  
  c_tau_muon_kinematics->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_tau_mu_pz[i]->GetMaximum() > h_tau_mu_pz[0]->GetMaximum()) h_tau_mu_pz[0]->SetMaximum(1.2*h_tau_mu_pz[i]->GetMaximum());
    if (h_tau_mu_pz[i]->GetMinimum() < h_tau_mu_pz[0]->GetMinimum() && h_tau_mu_pz[i]->GetMinimum() > 0) h_tau_mu_pz[0]->SetMinimum(0.8*h_tau_mu_pz[i]->GetMinimum());
  }
  h_tau_mu_pz[0]->Draw("ex0"); legend->Draw();
  c_tau_muon_kinematics->cd(2)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_mu_pz[i]->Draw("hesame");
  
  c_tau_muon_kinematics->cd(3); for (int i = 1; i < nSamples; i++){
    hs_tau_mu_pt->Add(h_tau_mu_pt[i]);
    if (h_tau_mu_pt[i]->GetMaximum() > h_tau_mu_pt[0]->GetMaximum()) h_tau_mu_pt[0]->SetMaximum(1.2*h_tau_mu_pt[i]->GetMaximum());
    if (h_tau_mu_pt[i]->GetMinimum() < h_tau_mu_pt[0]->GetMinimum() && h_tau_mu_pt[i]->GetMinimum() > 0) h_tau_mu_pt[0]->SetMinimum(0.8*h_tau_mu_pt[i]->GetMinimum());
  }
  h_tau_mu_pt[0]->Draw("ex0"); legend->Draw();
  c_tau_muon_kinematics->cd(3)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_mu_pt[i]->Draw("hesame");
  
  TCanvas *c_tau_muon_pt = new TCanvas("c_tau_muon_pt", "c_tau_muon_pt", 800, 800);
  h_tau_mu_pt[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_pt[i]->Draw("hesame");
  CMS_lumi( c_tau_muon_pt, iPeriod, iPos );
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt."+plotFormat).c_str());
  
  c_tau_muon_kinematics->cd(4); for (int i = 1; i < nSamples; i++){
    hs_tau_mu_eta->Add(h_tau_mu_eta[i]);
    if (h_tau_mu_eta[i]->GetMaximum() > h_tau_mu_eta[0]->GetMaximum()) h_tau_mu_eta[0]->SetMaximum(1.2*h_tau_mu_eta[i]->GetMaximum());
  }
  h_tau_mu_eta[0]->SetMaximum(1.4*h_tau_mu_eta[0]->GetMaximum());
  h_tau_mu_eta[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_eta[i]->Draw("hesame");
  
  
  TCanvas *c_tau_muon_eta = new TCanvas("c_tau_muon_eta", "c_tau_muon_eta", 800, 800);
  h_tau_mu_eta[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_eta[i]->Draw("hesame");
  CMS_lumi( c_tau_muon_eta, iPeriod, iPos );
  c_tau_muon_eta->SaveAs((basePlotDir+"/singlePlot/tau_muon_eta."+plotFormat).c_str());
  
  c_tau_muon_kinematics->cd(5); for (int i = 1; i < nSamples; i++){
    hs_tau_mu_phi->Add(h_tau_mu_phi[i]);
    if (h_tau_mu_phi[i]->GetMaximum() > h_tau_mu_phi[0]->GetMaximum()) h_tau_mu_phi[0]->SetMaximum(1.2*h_tau_mu_phi[i]->GetMaximum());
  }
  h_tau_mu_phi[0]->SetMaximum(1.9*h_tau_mu_phi[0]->GetMaximum());
  h_tau_mu_phi[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_phi[i]->Draw("hesame");
  
  TCanvas *c_tau_muon_phi = new TCanvas("c_tau_muon_phi", "c_tau_muon_phi", 800, 800);
  h_tau_mu_phi[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_phi[i]->Draw("hesame");
  CMS_lumi( c_tau_muon_phi, iPeriod, iPos );
  c_tau_muon_phi->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi."+plotFormat).c_str());
  
#ifdef nch_cut
  c_tau_muon_kinematics->SaveAs((basePlotDir+"/plots/nch_cut/tau_muon_kinematics."+plotFormat).c_str());
#else
  c_tau_muon_kinematics->SaveAs((basePlotDir+"/plots/tau_muon_kinematics."+plotFormat).c_str());
#endif

TCanvas *c_tau_had = new TCanvas("c_tau_had", "c_tau_had", 2000, 3000); c_tau_had->Divide(4,6);

  c_tau_had->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_tau_hadron_p[i]->GetMaximum() > h_tau_hadron_p[0]->GetMaximum()) h_tau_hadron_p[0]->SetMaximum(1.2*h_tau_hadron_p[i]->GetMaximum());
    if (h_tau_hadron_p[i]->GetMinimum() < h_tau_hadron_p[0]->GetMinimum() && h_tau_hadron_p[i]->GetMinimum() > 0) h_tau_hadron_p[0]->SetMinimum(0.8*h_tau_hadron_p[i]->GetMinimum());
  }
  h_tau_hadron_p[0]->Draw("ex0");
  c_tau_had->cd(1)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_p[i]->Draw("hesame");
  
  c_tau_had->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_tau_hadron_pz[i]->GetMaximum() > h_tau_hadron_pz[0]->GetMaximum()) h_tau_hadron_pz[0]->SetMaximum(1.2*h_tau_hadron_pz[i]->GetMaximum());
    if (h_tau_hadron_pz[i]->GetMinimum() < h_tau_hadron_pz[0]->GetMinimum() && h_tau_hadron_pz[i]->GetMinimum() > 0) h_tau_hadron_pz[0]->SetMinimum(0.8*h_tau_hadron_pz[i]->GetMinimum());
  }
  h_tau_hadron_pz[0]->Draw("ex0");
  c_tau_had->cd(2)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_pz[i]->Draw("hesame");
  
  c_tau_had->cd(3); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_pt->Add(h_tau_hadron_pt[i]);
    if (h_tau_hadron_pt[i]->GetMaximum() > h_tau_hadron_pt[0]->GetMaximum()) h_tau_hadron_pt[0]->SetMaximum(1.2*h_tau_hadron_pt[i]->GetMaximum());
  }
  h_tau_hadron_pt[0]->Draw("ex0");
  c_tau_had->cd(3)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_pt[i]->Draw("hesame");
  legend->Draw();
  
  TCanvas *c_tau_hadron_pt = new TCanvas("c_tau_hadron_pt", "c_tau_hadron_pt", 800, 800);
  h_tau_hadron_pt[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_pt[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_pt, iPeriod, iPos );
  c_tau_hadron_pt->SaveAs((basePlotDir+"/singlePlot/tau_hadron_pt."+plotFormat).c_str());
  
  c_tau_had->cd(4); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_eta->Add(h_tau_hadron_eta[i]);
    if (h_tau_hadron_eta[i]->GetMaximum() > h_tau_hadron_eta[0]->GetMaximum()) h_tau_hadron_eta[0]->SetMaximum(1.2*h_tau_hadron_eta[0]->GetMaximum());
  }
  h_tau_hadron_eta[0]->Draw("ex0");
  c_tau_had->cd(4)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_eta[i]->Draw("hesame");
  
  TCanvas *c_tau_hadron_eta = new TCanvas("c_tau_hadron_eta", "c_tau_hadron_eta", 800, 800);
  h_tau_hadron_eta[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_eta[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_eta, iPeriod, iPos );
  c_tau_hadron_eta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_eta."+plotFormat).c_str());
  
  c_tau_had->cd(5); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_phi->Add(h_tau_hadron_phi[i]);
    if (h_tau_hadron_phi[i]->GetMaximum() > h_tau_hadron_phi[0]->GetMaximum()) h_tau_hadron_phi[0]->SetMaximum(1.2*h_tau_hadron_phi[0]->GetMaximum());
  }
  h_tau_hadron_phi[0]->Draw("ex0");
  c_tau_had->cd(5)->SetLogy(1);
  h_tau_hadron_phi[0]->GetYaxis()->SetRangeUser(0.1,1.2*h_tau_hadron_phi[0]->GetMaximum());
  for (int i = 1; i < nSamples; i++) h_tau_hadron_phi[i]->Draw("hesame");
  
  TCanvas *c_tau_hadron_phi = new TCanvas("c_tau_hadron_phi", "c_tau_hadron_phi", 800, 800);
  h_tau_hadron_phi[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_phi[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_phi, iPeriod, iPos );
  c_tau_hadron_phi->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi."+plotFormat).c_str());
  
  for (int j = 0; j < 2; j++){
    c_tau_had->cd(j+6);
    for (int i = 1; i < nSamples; i++){
      hs_tau_hadron_rhomass[j]->Add(h_tau_hadron_rhomass[i][j]);
      if (h_tau_hadron_rhomass[i][j]->GetMaximum() > h_tau_hadron_rhomass[0][j]->GetMaximum()) h_tau_hadron_rhomass[0][j]->SetMaximum(h_tau_hadron_rhomass[i][j]->GetMaximum());
    }
    h_tau_hadron_rhomass[0][j]->Draw("ex0");
    for (int i = 1; i < nSamples; i++) h_tau_hadron_rhomass[i][j]->Draw("hesame");
  }
  c_tau_had->cd(8);
  for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_nch->Add(h_tau_hadron_nch[i]);
    if (h_tau_hadron_nch[i]->GetMaximum() > h_tau_hadron_nch[0]->GetMaximum()) h_tau_hadron_nch[0]->SetMaximum(1.2*h_tau_hadron_nch[i]->GetMaximum());
  }
  //if (hs_tau_hadron_nch->GetMaximum() < h_tau_hadron_nch[0]->GetMaximum()) hs_tau_hadron_nch->SetMaximum(h_tau_hadron_nch[0]->GetMaximum());
  /*hs_tau_hadron_nch->Draw("he");*/ h_tau_hadron_nch[0]->Draw("e1same"); for (int i = 1; i < nSamples; i++) h_tau_hadron_nch[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_tau_hadron_nch = new TCanvas("c_tau_hadron_nch", "c_tau_hadron_nch", 800, 800);
  h_tau_hadron_nch[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_nch[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_nch, iPeriod, iPos );
  c_tau_hadron_nch->SaveAs((basePlotDir+"/singlePlot/tau_hadron_nch."+plotFormat).c_str());
  
  c_tau_had->cd(9); h_tau_hadron_ncand_final[0]->Draw("e1same"); for (int i = 1; i < nSamples; i++) h_tau_hadron_ncand_final[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(10); h_tau_hadron_vprob[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++) h_tau_hadron_vprob[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(11); for (int i = 1; i < nSamples; i++) h_tau_hadron_matched_pt_index[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(12); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_mass->Add(h_tau_hadron_mass[i]);
    if (h_tau_hadron_mass[i]->GetMaximum() > h_tau_hadron_mass[0]->GetMaximum()) h_tau_hadron_mass[0]->SetMaximum(1.2*h_tau_hadron_mass[i]->GetMaximum());
  }
  //h_tau_hadron_mass[0]->SetTitle("#tau_{hadron} visible mass");
  h_tau_hadron_mass[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_mass[i]->Draw("hesame");
  
  TCanvas *c_tau_hadron_mass = new TCanvas("c_tau_hadron_mass", "c_tau_hadron_mass", 800, 800);
  h_tau_hadron_mass[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_mass[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_mass, iPeriod, iPos );
  c_tau_hadron_mass->SaveAs((basePlotDir+"/singlePlot/tau_hadron_mass."+plotFormat).c_str());
  
  c_tau_had->cd(13); for (int i = 1; i < nSamples; i++){
    hs_ditau_mass->Add(h_ditau_mass[i]);
    if (h_ditau_mass[i]->GetMaximum() > h_ditau_mass[0]->GetMaximum()) h_ditau_mass[0]->SetMaximum(1.2*h_ditau_mass[i]->GetMaximum());
  }
  //h_ditau_mass[0]->SetTitle("#tau#tau visible mass");
  h_ditau_mass[0]->Draw("e1same");
  for (int i = 1; i < nSamples; i++) h_ditau_mass[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_ditau_mass = new TCanvas("c_ditau_mass", "c_ditau_mass", 800, 800);
  h_ditau_mass[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_ditau_mass[i]->Draw("hesame");
  CMS_lumi( c_ditau_mass, iPeriod, iPos );
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass."+plotFormat).c_str());
  
  c_tau_had->cd(14); for (int i = 1; i < nSamples; i++){
    if (h_ditau_p[i]->GetMaximum() > h_ditau_p[0]->GetMaximum()) h_ditau_p[0]->SetMaximum(1.2*h_ditau_p[i]->GetMaximum());
    if (h_ditau_p[i]->GetMinimum() < h_ditau_p[0]->GetMinimum() && h_ditau_p[i]->GetMinimum() > 0) h_ditau_p[0]->SetMinimum(0.8*h_ditau_p[i]->GetMinimum());
  }
  h_ditau_p[0]->Draw("e1same");
  c_tau_had->cd(14)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_ditau_p[i]->Draw("hesame"); legend->Draw();
  
  c_tau_had->cd(15); for (int i = 1; i < nSamples; i++){
    if (h_ditau_pz[i]->GetMaximum() > h_ditau_pz[0]->GetMaximum()) h_ditau_pz[0]->SetMaximum(1.2*h_ditau_pz[i]->GetMaximum());
    if (h_ditau_pz[i]->GetMinimum() < h_ditau_pz[0]->GetMinimum() && h_ditau_pz[i]->GetMinimum() > 0) h_ditau_pz[0]->SetMinimum(0.8*h_ditau_pz[i]->GetMinimum());
  }
  h_ditau_pz[0]->Draw("e1same");
  c_tau_had->cd(15)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_ditau_pz[i]->Draw("hesame"); legend->Draw();
  
  c_tau_had->cd(16); h_tau_hadron_nPions[0]->Scale(1./h_tau_hadron_nPions[0]->GetMaximum()); h_tau_hadron_nPions[0]->Draw("ex0");
  for (int i = 1; i < nSamples; i++){
    h_tau_hadron_nPions[i]->Scale(1./h_tau_hadron_nPions[i]->GetMaximum());
    h_tau_hadron_nPions[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_tau_had->cd(17); for (int i = 1; i < nSamples; i++){
    h_resVisTauPt[i]->Scale(1./h_resVisTauPt[i]->GetMaximum());
    h_resVisTauPt[i]->Draw("hesame");
  }
  legend->Draw();
  c_tau_had->cd(18); for (int i = 1; i < nSamples; i++){
    h_resVisTauEta[i]->Scale(1./h_resVisTauEta[i]->GetMaximum());
    h_resVisTauEta[i]->Draw("hesame");
  }
  legend->Draw();
  c_tau_had->cd(19); for (int i = 1; i < nSamples; i++){
    h_resVisTauPhi[i]->Scale(1./h_resVisTauPhi[i]->GetMaximum());
    h_resVisTauPhi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_tau_had->cd(20); for (int i = 1; i < nSamples; i++){
    if (h_resPt_matchedPionReco[i]->GetMaximum() > h_resPt_matchedPionReco[1]->GetMaximum()) h_resPt_matchedPionReco[1]->SetMaximum(1.2*h_resPt_matchedPionReco[i]->GetMaximum());
    h_resPt_matchedPionReco[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_tau_had->cd(21); for (int i = 1; i < nSamples; i++){
    if (h_resR_matchedPionReco[i]->GetMaximum() > h_resR_matchedPionReco[1]->GetMaximum()) h_resR_matchedPionReco[1]->SetMaximum(1.2*h_resR_matchedPionReco[i]->GetMaximum());
    h_resR_matchedPionReco[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_tau_had->cd(22); h_resPhi_matchedPionReco[1]->Draw("colz");
  
  c_tau_had->cd(23); h_reco_pion_energy_HCAL_ECAL[0]->Draw("colz");
  
  c_tau_had->cd(24); for (int i = 1; i < nSamples; i++){
    if (h_gen_tau_hadron_visible_pt[i]->GetMaximum() > h_gen_tau_hadron_visible_pt[1]->GetMaximum()) h_gen_tau_hadron_visible_pt[1]->SetMaximum(1.2*h_gen_tau_hadron_visible_pt[i]->GetMaximum());
  }
  c_tau_had->cd(24)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_gen_tau_hadron_visible_pt[i]->Draw("hesame");
  legend->Draw();
  
  
#ifdef nch_cut
  c_tau_had->SaveAs((basePlotDir+"/plots/nch_cut/tau_had."+plotFormat).c_str());
#else
  c_tau_had->SaveAs((basePlotDir+"/plots/tau_had."+plotFormat).c_str());
#endif

  TCanvas *c_recoPion = new TCanvas("c_recoPion", "c_recoPion", 3*800, 3*800); c_recoPion->Divide(3,3);
  
  c_recoPion->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_pion_leading_pt[i]->GetMaximum() > h_pion_leading_pt[0]->GetMaximum()) h_pion_leading_pt[0]->SetMaximum(1.2*h_pion_leading_pt[i]->GetMaximum());
  }
  h_pion_leading_pt[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_leading_pt[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_leading_pt = new TCanvas("c_pion_leading_pt", "c_pion_leading_pt", 800, 800);
  h_pion_leading_pt[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_leading_pt[i]->Draw("hesame");
  CMS_lumi( c_pion_leading_pt, iPeriod, iPos );
  c_pion_leading_pt->SaveAs((basePlotDir+"/singlePlot/pion_leading_pt."+plotFormat).c_str());
  
  c_recoPion->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_pion_subleading_pt[i]->GetMaximum() > h_pion_subleading_pt[0]->GetMaximum()) h_pion_subleading_pt[0]->SetMaximum(1.2*h_pion_subleading_pt[i]->GetMaximum());
  }
  h_pion_subleading_pt[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subleading_pt[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subleading_pt = new TCanvas("c_pion_subleading_pt", "c_pion_subleading_pt", 800, 800);
  h_pion_subleading_pt[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_subleading_pt[i]->Draw("hesame");
  CMS_lumi( c_pion_subleading_pt, iPeriod, iPos );
  c_pion_subleading_pt->SaveAs((basePlotDir+"/singlePlot/pion_subleading_pt."+plotFormat).c_str());
  
  c_recoPion->cd(3); for (int i = 1; i < nSamples; i++){
    if (h_pion_subsubleading_pt[i]->GetMaximum() > h_pion_subsubleading_pt[0]->GetMaximum()) h_pion_subsubleading_pt[0]->SetMaximum(1.2*h_pion_subsubleading_pt[i]->GetMaximum());
  }
  h_pion_subsubleading_pt[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_pt[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subsubleading_pt = new TCanvas("c_pion_subsubleading_pt", "c_pion_subsubleading_pt", 800, 800);
  h_pion_subsubleading_pt[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_pt[i]->Draw("hesame");
  CMS_lumi( c_pion_subsubleading_pt, iPeriod, iPos );
  c_pion_subsubleading_pt->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_pt."+plotFormat).c_str());
  
  c_recoPion->cd(4); for (int i = 1; i < nSamples; i++){
    if (h_pion_leading_eta[i]->GetMaximum() > h_pion_leading_eta[0]->GetMaximum()) h_pion_leading_eta[0]->SetMaximum(1.2*h_pion_leading_eta[i]->GetMaximum());
  }
  h_pion_leading_eta[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_leading_eta[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_leading_eta = new TCanvas("c_pion_leading_eta", "c_pion_leading_eta", 800, 800);
  h_pion_leading_eta[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_leading_eta[i]->Draw("hesame");
  CMS_lumi( c_pion_leading_eta, iPeriod, iPos );
  c_pion_leading_eta->SaveAs((basePlotDir+"/singlePlot/pion_leading_eta."+plotFormat).c_str());
  
  c_recoPion->cd(5); for (int i = 1; i < nSamples; i++){
    if (h_pion_subleading_eta[i]->GetMaximum() > h_pion_subleading_eta[0]->GetMaximum()) h_pion_subleading_eta[0]->SetMaximum(1.2*h_pion_subleading_eta[i]->GetMaximum());
  }
  h_pion_subleading_eta[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subleading_eta[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subleading_eta = new TCanvas("c_pion_subleading_eta", "c_pion_subleading_eta", 800, 800);
  h_pion_subleading_eta[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_subleading_eta[i]->Draw("hesame");
  CMS_lumi( c_pion_subleading_eta, iPeriod, iPos );
  c_pion_subleading_eta->SaveAs((basePlotDir+"/singlePlot/pion_subleading_eta."+plotFormat).c_str());
  
  c_recoPion->cd(6); for (int i = 1; i < nSamples; i++){
    if (h_pion_subsubleading_eta[i]->GetMaximum() > h_pion_subsubleading_eta[0]->GetMaximum()) h_pion_subsubleading_eta[0]->SetMaximum(1.2*h_pion_subsubleading_eta[i]->GetMaximum());
  }
  h_pion_subsubleading_eta[0]->SetMaximum(1.3*h_pion_subsubleading_eta[0]->GetMaximum());
  h_pion_subsubleading_eta[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_eta[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subsubleading_eta = new TCanvas("c_pion_subsubleading_eta", "c_pion_subsubleading_eta", 800, 800);
  h_pion_subsubleading_eta[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_eta[i]->Draw("hesame");
  CMS_lumi( c_pion_subsubleading_eta, iPeriod, iPos );
  c_pion_subsubleading_eta->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_eta."+plotFormat).c_str());
  
  c_recoPion->cd(7); for (int i = 1; i < nSamples; i++){
    if (h_pion_leading_phi[i]->GetMaximum() > h_pion_leading_phi[0]->GetMaximum()) h_pion_leading_phi[0]->SetMaximum(1.2*h_pion_leading_phi[i]->GetMaximum());
  }
  h_pion_leading_phi[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_leading_phi[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_leading_phi = new TCanvas("c_pion_leading_phi", "c_pion_leading_phi", 800, 800);
  h_pion_leading_phi[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_leading_phi[i]->Draw("hesame");
  CMS_lumi( c_pion_leading_phi, iPeriod, iPos );
  c_pion_leading_phi->SaveAs((basePlotDir+"/singlePlot/pion_leading_phi."+plotFormat).c_str());
  
  c_recoPion->cd(8); for (int i = 1; i < nSamples; i++){
    if (h_pion_subleading_phi[i]->GetMaximum() > h_pion_subleading_phi[0]->GetMaximum()) h_pion_subleading_phi[0]->SetMaximum(1.2*h_pion_subleading_phi[i]->GetMaximum());
  }
  h_pion_subleading_phi[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subleading_phi[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subleading_phi = new TCanvas("c_pion_subleading_phi", "c_pion_subleading_phi", 800, 800);
  h_pion_subleading_phi[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_subleading_phi[i]->Draw("hesame");
  CMS_lumi( c_pion_subleading_phi, iPeriod, iPos );
  c_pion_subleading_phi->SaveAs((basePlotDir+"/singlePlot/pion_subleading_phi."+plotFormat).c_str());
  
  c_recoPion->cd(9); for (int i = 1; i < nSamples; i++){
    if (h_pion_subsubleading_phi[i]->GetMaximum() > h_pion_subsubleading_phi[0]->GetMaximum()) h_pion_subsubleading_phi[0]->SetMaximum(1.2*h_pion_subsubleading_phi[i]->GetMaximum());
  }
  h_pion_subsubleading_phi[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_phi[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subsubleading_phi = new TCanvas("c_pion_subsubleading_phi", "c_pion_subsubleading_phi", 800, 800);
  h_pion_subsubleading_phi[0]->Draw("ex0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_phi[i]->Draw("hesame");
  CMS_lumi( c_pion_subsubleading_phi, iPeriod, iPos );
  c_pion_subsubleading_phi->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_phi."+plotFormat).c_str());
  
  
#ifdef nch_cut
  c_recoPion->SaveAs((basePlotDir+"/plots/nch_cut/reco_pion."+plotFormat).c_str());
#else
  c_recoPion->SaveAs((basePlotDir+"/plots/reco_pion."+plotFormat).c_str());
#endif

/*  TCanvas *c_track_activity = new TCanvas("c_track_activity", "c_track_activity", 1200, 600); c_track_activity->Divide(2,1);
//  c_track_activity->cd(1); h_track_activity_pt->Draw("ex0");
  c_track_activity->cd(1); h_track_activity_pt->DrawNormalized("e1");
//  c_track_activity->cd(2); h_track_activity_pt_eta->Draw("COLZ");
  c_track_activity->cd(2); h_track_activity_pt_eta->DrawNormalized("COLZ");
#ifdef nch_cut
  c_track_activity->SaveAs((basePlotDir+"/plots/nch_cut/track_activity."+plotFormat).c_str());
#else
  c_track_activity->SaveAs((basePlotDir+"/plots/track_activity."+plotFormat).c_str());
#endif
*/
  TCanvas *c_calo = new TCanvas("c_calo", "c_calo", 1600, 5*800); c_calo->Divide(2,5);
  for (int i = 1; i < nSamples; i++) hs_calo_energyHFp->Add(h_calo_energyHFp[i]);
  c_calo->cd(1); h_calo_energyHFp[0]->Draw("ex0"); h_calo_energyHFp[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_energyHFp[0]->GetMaximum());
  c_calo->cd(1)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_energyHFp[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++) hs_calo_energyHFm->Add(h_calo_energyHFm[i]);
  c_calo->cd(2); h_calo_energyHFm[0]->Draw("ex0"); h_calo_energyHFm[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_energyHFm[0]->GetMaximum());
  c_calo->cd(2)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_energyHFm[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++){
    hs_calo_leadingHFp->Add(h_calo_leadingHFp[i]);
    if (h_calo_leadingHFp[i]->GetMaximum() > h_calo_leadingHFp[0]->GetMaximum()) h_calo_leadingHFp[0]->SetMaximum(1.2*h_calo_leadingHFp[i]->GetMaximum());
  }
  c_calo->cd(3); h_calo_leadingHFp[0]->Draw("ex0"); //h_calo_leadingHFp[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_leadingHFp[0]->GetMaximum());
  c_calo->cd(3)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_leadingHFp[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++){
    hs_calo_leadingHFm->Add(h_calo_leadingHFm[i]);
    if (h_calo_leadingHFm[i]->GetMaximum() > h_calo_leadingHFm[0]->GetMaximum()) h_calo_leadingHFm[0]->SetMaximum(1.2*h_calo_leadingHFm[i]->GetMaximum());
  }
  c_calo->cd(4); h_calo_leadingHFm[0]->Draw("ex0"); //h_calo_leadingHFm[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_leadingHFm[0]->GetMaximum());
  c_calo->cd(4)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_leadingHFm[i]->Draw("hesame"); legend->Draw();
  c_calo->cd(5); h_calo_energyHFp_nch->Draw("COLZ"); c_calo->cd(5)->SetLogz(1);
  c_calo->cd(6); h_calo_energyHFm_nch->Draw("COLZ"); c_calo->cd(6)->SetLogz(1);
  c_calo->cd(7); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFp_sum[i]->Draw("e1same");
  c_calo->cd(8); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFm_sum[i]->Draw("e1same");
  c_calo->cd(9); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFp_size[i]->Draw("e1same");
  c_calo->cd(10); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFm_size[i]->Draw("e1same");
#ifdef nch_cut
  c_calo->SaveAs((basePlotDir+"/plots/nch_cut/calo."+plotFormat).c_str());
#else
  c_calo->SaveAs((basePlotDir+"/plots/calo."+plotFormat).c_str());
#endif

  TCanvas *c_zdc = new TCanvas("c_zdc", "c_zdc", 1000, 1500); c_zdc->Divide(2,3);
  c_zdc->cd(1); for (int i = 0; i < nSamples; i++) h_sumZDCplus[i]->Draw("e1same");
  c_zdc->cd(1)->SetLogy(1);
  c_zdc->cd(2); for (int i = 0; i < nSamples; i++) h_sumZDCminus[i]->Draw("e1same");
  c_zdc->cd(2)->SetLogy(1);
  c_zdc->cd(3); h_sumZDC_pm[0]->Draw("colz");
  c_zdc->cd(3)->SetLogz(1);
  c_zdc->cd(4); h_averageZDCside_averageHFeta[0]->Draw("colz");
  c_zdc->cd(4)->SetLogz(1);
  c_zdc->cd(5); h_muEta_averageZDCside[0]->Draw("colz");
  c_zdc->cd(5)->SetLogz(1);
  c_zdc->cd(6); h_tauEta_averageZDCside[0]->Draw("colz");
  c_zdc->cd(6)->SetLogz(1);
#ifdef nch_cut
  c_zdc->SaveAs((basePlotDir+"/plots/nch_cut/zdc."+plotFormat).c_str());
#else
  c_zdc->SaveAs((basePlotDir+"/plots/zdc."+plotFormat).c_str());
#endif


  TCanvas *c_HF_eta = new TCanvas("c_HF_eta", "c_HF_eta", 500*(nSamples+1), 500); c_HF_eta->Divide(nSamples+1,1);
  for (int i = 1; i < nSamples; i++) hs_calo_HF_eta->Add(h_calo_HF_eta[i]);
  c_HF_eta->cd(1); h_calo_HF_eta[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++) h_calo_HF_eta[i]->Draw("hesame"); legend->Draw();
  for (int i = 0; i < nSamples; i++){ c_HF_eta->cd(i+2); h_calo_HF_energy_eta[i]->Draw("colz"); c_HF_eta->cd(i+2)->SetLogz(1);}
#ifdef nch_cut
  c_HF_eta->SaveAs((basePlotDir+"/plots/nch_cut/HF_eta."+plotFormat).c_str());
#else
  c_HF_eta->SaveAs((basePlotDir+"/plots/HF_eta."+plotFormat).c_str());
#endif

  TCanvas *c_HFpm = new TCanvas("c_HFpm", "c_HFpm", 500*nSamples, 3000); c_HFpm->Divide(nSamples,6);
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(i+1); h_calo_energyHF_pm[i]->Draw("colz"); c_HFpm->cd(i+1)->SetLogz(1);}
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(nSamples+i+1); h_calo_leadingHF_pm[i]->Draw("colz"); c_HFpm->cd(nSamples+i+1)->SetLogz(1);}
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(2*nSamples+i+1); h_calo_muEta_leadingHFeta[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(3*nSamples+i+1); h_calo_tauEta_leadingHFeta[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(4*nSamples+i+1); h_calo_muEta_averageHFeta[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(5*nSamples+i+1); h_calo_tauEta_averageHFeta[i]->Draw("colz"); }
#ifdef nch_cut
  c_HFpm->SaveAs((basePlotDir+"/plots/nch_cut/HFpm."+plotFormat).c_str());
#else
  c_HFpm->SaveAs((basePlotDir+"/plots/HFpm."+plotFormat).c_str());
#endif

  TCanvas *c_delta_phi_eta = new TCanvas("c_delta_phi_eta", "c_delta_phi_eta", 500*nSamples, 1500); c_delta_phi_eta->Divide(nSamples,3);
  for (int i = 0; i < nSamples; i++){
    c_delta_phi_eta->cd(i+1); h_deltaphi_tau_mu_tau_hadron_mueta[i]->Draw("colz");
    c_delta_phi_eta->cd(nSamples+i+1); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->Draw("colz");
    c_delta_phi_eta->cd(2*nSamples+i+1); h_mueta_taueta[i]->Draw("colz");
  }
#ifdef nch_cut
  c_delta_phi_eta->SaveAs((basePlotDir+"/plots/nch_cut/delta_phi_eta."+plotFormat).c_str());
#else
  c_delta_phi_eta->SaveAs((basePlotDir+"/plots/delta_phi_eta."+plotFormat).c_str());
#endif

  TCanvas *c_AP = new TCanvas("c_AP", "c_AP", 500*nSamples, 500); c_AP->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_AP->cd(i+1); h_AP[i]->Draw("colz");}
#ifdef nch_cut
  c_AP->SaveAs((basePlotDir+"/plots/nch_cut/AP."+plotFormat).c_str());
#else
  c_AP->SaveAs((basePlotDir+"/plots/AP."+plotFormat).c_str());
#endif

  TCanvas *c_rho = new TCanvas("c_rho", "c_rho", 500*nSamples, 500); c_rho->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_rho->cd(i+1); h_tau_hadron_rhomass2D[i]->Draw("colz");}
#ifdef nch_cut
  c_rho->SaveAs((basePlotDir+"/plots/nch_cut/rho."+plotFormat).c_str());
#else
  c_rho->SaveAs((basePlotDir+"/plots/rho."+plotFormat).c_str());
#endif

  TCanvas *c_MET = new TCanvas("c_MET", "c_MET", 800, 800);
  c_MET->cd(1); h_MET[0]->Draw("ex0"); for (int i = 1; i < nSamples; i++) h_MET[i]->Draw("hesame"); legend->Draw();
#ifdef nch_cut
  c_MET->SaveAs((basePlotDir+"/plots/nch_cut/MET."+plotFormat).c_str());
#else
  c_MET->SaveAs((basePlotDir+"/plots/MET."+plotFormat).c_str());
#endif

  TCanvas *c_cutflow = new TCanvas("c_cutflow", "c_cutflow", 1000, 500); c_cutflow->Divide(2,1);
  c_cutflow->cd(1);
  cutflow[0]->Draw("ex0");
  cutflow[0]->GetYaxis()->SetRangeUser(1,2*cutflow[0]->GetMaximum());
  for (int i = 1; i < nSamples; i++) cutflow[i]->Draw("hesame");
  legend->Draw();
  c_cutflow->cd(1)->SetLogy(1);
  c_cutflow->cd(2);
  h_cutflow[0]->Draw("ex0");
  c_cutflow->cd(2)->SetLogy(1);
#ifdef nch_cut
  c_cutflow->SaveAs((basePlotDir+"/plots/nch_cut/cutflow."+plotFormat).c_str());
#else
  c_cutflow->SaveAs((basePlotDir+"/plots/cutflow."+plotFormat).c_str());
#endif

#ifdef nch_cut
  TFile *outputFile = new TFile("histsNch_cut.root","RECREATE");
  outputFile->cd();
  outputFile->mkdir("nch_cut");
  outputFile->cd("nch_cut");
#else
  TFile *outputFile = new TFile("histsNOch_cut.root","RECREATE");
  outputFile->cd();
  outputFile->mkdir("noch_cut");
  outputFile->cd("noch_cut");
#endif
  TCanvas *c_output = new TCanvas("c_output", "c_output", 500*nSamples, 500); c_output->Divide(nSamples,1);
  TH1F *h_deltaphi_tau_mu_tau_hadron_new = (TH1F*)h_deltaphi_tau_mu_tau_hadron[0]->Clone("data_obs");
  c_output->cd(1); h_deltaphi_tau_mu_tau_hadron_new->DrawCopy("e1");
  h_deltaphi_tau_mu_tau_hadron_new->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron_new->Write();
  for (int i = 1; i < nSamples; i++){
    c_output->cd(i+1); h_deltaphi_tau_mu_tau_hadron[i]->DrawCopy("e1");
    h_deltaphi_tau_mu_tau_hadron[i]->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron[i]->Write();
  }
  //outputFile->Write();
  outputFile->Close();
  
  }