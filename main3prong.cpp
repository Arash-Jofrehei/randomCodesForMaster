#include <vector>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <cmath>
using namespace std;
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TF1.h"
//#include <TStyle.h>
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
//#include "TStyle.h"
#include "TREE_DATA.C"
#include "TREE_MC.C"
//#include "muon_tnp_weight.h"
#include "GK_tnp_weight.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
//#include "/Applications/root/macros/AtlasStyle.C"

#define nch_cut
//#define post
//#define ratio

#ifdef post
string anaTag = "postfit";
#else
string anaTag = "test";
//string anaTag = "prefit";
#endif

string plotFormat = "png";


// **************************
// cut definition
  bool tau_muon_isSoft = true;
  bool tau_muon_isGlobal = false;
  bool tau_muon_isTracker = false;
  double tau_hadron_vertexprob = 0.005;
  float HFpLeading_low = 0;
  float HFmLeading_low = 0;
  float HFpLeading_high = 6.0;
  float HFmLeading_high = 6.0;
  float pionLeadingCut = 0.5;
  float pionSubLeadingCut = 0.3;
  float pionSubSubLeadingCut = 0.3;
  float ditauMassCut = 0;
  float gammaPtCut = 0.;
  float minZDCp = 0;
  float minZDCm = 0;
  float maxZDCp = 42000;
  float maxZDCm = 60000;
  float ratioZDCpm = 4.2/6;
  int min_nch = 3;
  int max_nch = 3;
  float MET_cut = 20000;
  //float deltaPhi_cut = 2.78816;
  float deltaPhi_cut = 0;
  int MuTauCharge = -1;
  bool hasZDCinfo = false;
  int nNchCategories = 4;
  int firstNchCategory = 5;
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

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcgPb_AOD_SuperChic_2015_020721.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_orderedPion_041021.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_orderedPion_041021.root"};


//string inputFiles[] = {"flatTuple_AOD_data_2015_lowMuPtCut_241121.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_lowMuPtCut_241121.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_noTauPtCut_171021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_noTauPtCut_061021.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_lowMuPtCut_241121.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_noTauPtCut_171021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};


//string inputFiles[] = {"flatTuple_AOD_data_2015_noTauPtCut_061021.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_noTauPtCut_061021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_120422.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_noTauPtCut_171021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_110522.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_100522.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_110522.root","mu3prong_UPCgen_atau-10E-2_2018_270123.root","mu3prong_UPCgen_atau0E-2_2018_270123.root","mu3prong_UPCgen_atau10E-2_2018_270123.root"};

string inputFiles[] = {"flatTuple_AOD_subdata_2018_110522.root","mu3prong_UPCgen_atau0E-2_2018_270123.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_190522.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_190522.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_190522.root"};

//string inputFiles[] = {"mu3prong_AOD_subdata_2018_040922.root","ggTauTau_SuperChic_mu3prong_2018_050922.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2prong_2015_300921.root","flatTuple_mcggTauTau_AOD_pp_2pi_MadGraph_2015_280421.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_090621.root"};
  
//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_data_2015_190421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_050521.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root","flatTuple_mcggCCbar_2015_300421-1.root","flatTuple_mcggBBbar_2015_060521.root"};

//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_subdata_2018_210421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_270421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_subdata_2018_210421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_pp_2pi_MadGraph_2015_280421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_2pi_MadGraph_2018_270421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_2pi_SuperChic_2018_270421.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015_220321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_embedded.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018.root", "../ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018.root", "/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggBBbar_4f_AOD_MadGraph_2015.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_080321.root"};

//string inputFiles[] = {"../ntuples/flatTuple_data_2015.root","../ntuples/flatTuple_mcggTauTau_2015.root","../ntuples/flatTuple_mcggBBbar_2015.root","../ntuples/flatTuple_mcggCCbar_2015.root","../ntuples/flatTuple_mcggTauTau_HydjetDrumMB_2018.root"};

void main3prong(const int nSamples = 2, string files[] = inputFiles){

  string basePlotDir = "/eos/user/a/ajofrehe/www/gtau/2018/3prong/"+anaTag;
  
  string ABCDsysNames[2+nNchCategories];
  ABCDsysNames[0] = "ABCD-sys-HFDown";
  ABCDsysNames[1] = "ABCD-sys-HFUp";
  string baseABCDsysName = "highNch_";
  for (int n = 0; n < nNchCategories; n++) ABCDsysNames[n+2] = baseABCDsysName + to_string(firstNchCategory+n);

  //mkdir(directory.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
//void run(vector<string> files){
  //const int nSamples = files.size();
  TFile *root_file[nSamples];
  TREE_DATA *TREE;
  TREE_MC *TREEMCs[nSamples-1];
  TH1F *cutflow[nSamples];
  //SetAtlasStyle();
  setTDRStyle();
  int colors[] = {1,2,3,4,6,7,9,11,38,41,46};
  int styles[] = {20,1,1,1,1,1,1,1,1,1,1};

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
  
  //TFile *pionEffFile = new TFile("pionEff.root","READ");
  //TH2F *pion_eff = (TH2F*) pionEffFile->Get("pion_eff");
  //TAxis *xaxis_eff = pion_eff->GetXaxis();
  //TAxis *yaxis_eff = pion_eff->GetYaxis();
  //xaxis_eff->kCanExtend = false;
  //yaxis_eff->kCanExtend = false;
  
  TFile *hists = new TFile("hists.root","READ");
  TH1F *h_eff_matchedPionPt = (TH1F*) hists->Get("3prong/h_eff_matchedPionPt");
  
// *******************************************
  
  
  writeExtraText = true;       // if extra text
  extraText  = "work in progress";  // default extra text is "Preliminary"
  //lumi_5TeV  = "PbPb - 404 #mub^{-1} ";
  lumi_5TeV  = "PbPb - 425 #mub^{-1} ";
  //lumi_5TeV  = "PbPb - 600 #mub^{-1} ";
  lumi_sqrtS = " (#sqrt{s_{NN}} = 5.02 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

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
  
// *******************************************
// histogram declaration
  
  //float dataLumi = 404.318415215; // ub - 2015 data
  float dataLumi = 425.911750;  // ub - 2018 sub data
  //float dataLumi = 600;  // ub - 2018 sub data
  double PiZeroMass = 134.98; //MeV
  //string tag[] = {"", "Signal 2015 MadGraph", "Signal 2018 MadGraph", "Signal 2018 SuperChic", "Signal 2018 SuperChic + Embedded MB"};
  //string tag[] = {"", "MadGraph 2015", "CCbar 2015", "BBbar 2015", "BBbar 2018"};
  //string tag[] = {"","UPCgen atau-10E-2","UPCgen atau0E-2","UPCgen atau10E-2"};
  //string tags[] = {"","UPCgen atau-10E-2","UPCgen atau0E-2","UPCgen atau10E-2"};
  string tag[] = {"","signal"};
  string tags[] = {"","#gamma#gamma#rightarrow#tau#tau"};
  //string tag[] = {"", "MG signal", "SC signal"};
  //string tags[] = {"", "MG signal", "SC signal"};
  //string tags[] = {"", "#gamma#gamma#rightarrow#tau_{#mu}#tau_{3prong}", "SuperChic 2018", "CCbar 2015", "BBbar 2015"};
  //string tag[] = {"", "MC Signal", "MC 1prong", "MC 1prong + 0#pi^{0}", "MC 1prong + 1#pi^{0}", "MC 1prong + n#pi^{0}", "MC Signal 2018"};
  //string tag[] = {"", "PbPb reco 2015", "pp reco 2015", "pp reco 2018"};
  //string tag[] = {"", "MC Signal 2015", "Signal 2018", "SuperChic 2018"};
  float nEvents[] = {1,1184600,2000000,623000,30000000,1,1,1,1};
  for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1);
  //nEvents[4] /= 10;
  //for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1)/200;
  //float crossSectionMC[4] = {570000,1500,300000,570000};
  //float crossSectionMC[] = {570000,300000,1500,1500};
  //float crossSectionMC[] = {570,570,570,570,570,570};
  //float crossSectionMC[] = {680.211,847.772,1139.995};
  float crossSectionMC[] = {847.772};
  float SF[] = {1,1,1,1,1,1,1};
  //float pT_SF = 1677539 / 2000000.0;
  float above3GeVSuperChic = 322461;
  float nAllEventsSuperChic = 2000000;
  //float pT_SF = above3GeVSuperChic / nAllEventsSuperChic;
  float pT_SF = 0.21;
  float tau_pT_SF[] = {1,1,1,1,1,1};
  //float tau_pT_SF[] = {1,pT_SF,pT_SF,pT_SF,pT_SF,pT_SF,1};
  bool threeProng[] = {0,1,1,1,1,1,1};
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
  int deltaR_bins = 50;
  int deltaEta_bins = 30;
  
  vector<TH1F**> histograms = std::vector<TH1F**>();
  vector<THStack*> stacks = std::vector<THStack*>();
  vector<TH1F**> ratios = std::vector<TH1F**>();
  vector<string> plotNames = std::vector<string>();
  vector<TCanvas*> ratioCanvas = std::vector<TCanvas*>();
  
  
  TH1F *h_cutflow[nSamples];
  TH1F *h_tau_mu_p[nSamples+4];
  histograms.push_back(h_tau_mu_p);
  plotNames.push_back("tau_muon_p");
  TH1F *h_tau_mu_pz[nSamples+4];
  histograms.push_back(h_tau_mu_pz);
  plotNames.push_back("tau_muon_pz");
  TH1F *h_tau_mu_dz[nSamples+4];
  histograms.push_back(h_tau_mu_dz);
  plotNames.push_back("tau_muon_dz");
  THStack *hs_tau_mu_pt = new THStack("hs_tau_mu_pt","#tau_{#mu} p_{T} [GeV]");
  TH1F *h_tau_mu_pt[nSamples+4];
  histograms.push_back(h_tau_mu_pt);
  plotNames.push_back("tau_muon_pt");
  THStack *hs_tau_mu_eta = new THStack("hs_tau_mu_eta","#tau_{#mu} #eta");
  TH1F *h_tau_mu_eta[nSamples+4];
  histograms.push_back(h_tau_mu_eta);
  plotNames.push_back("tau_muon_eta");
  THStack *hs_tau_mu_phi = new THStack("hs_tau_mu_phi","#tau_{#mu} #phi");
  TH1F *h_tau_mu_phi[nSamples+4];
  histograms.push_back(h_tau_mu_phi);
  plotNames.push_back("tau_muon_phi");
  TH2F *h_tau_mu_tau_hadron_phi[nSamples];
  TH1F *h_tau_hadron_p[nSamples+4];
  histograms.push_back(h_tau_hadron_p);
  plotNames.push_back("tau_hadron_p");
  TH1F *h_tau_hadron_pz[nSamples+4];
  histograms.push_back(h_tau_hadron_pz);
  plotNames.push_back("tau_hadron_pz");
  TH1F *h_tau_hadron_ptScalar[nSamples+4];
  histograms.push_back(h_tau_hadron_ptScalar);
  plotNames.push_back("tau_hadron_scalar_pt");
  THStack *hs_tau_hadron_pt = new THStack("hs_tau_hadron_pt","#tau_{3prong} p_{T} [GeV]");
  TH1F *h_tau_hadron_pt[nSamples+4];
  histograms.push_back(h_tau_hadron_pt);
  plotNames.push_back("tau_hadron_pt");
  TH1F *h_gen_tau_hadron_visible_pt[nSamples];
  THStack *hs_tau_hadron_eta = new THStack("hs_tau_hadron_eta","#tau_{3prong} #eta");
  TH1F *h_tau_hadron_eta[nSamples+4];
  histograms.push_back(h_tau_hadron_eta);
  plotNames.push_back("tau_hadron_eta");
  THStack *hs_tau_hadron_phi = new THStack("hs_tau_hadron_phi","#tau_{3prong} #phi");
  TH1F *h_tau_hadron_phi[nSamples+4];
  histograms.push_back(h_tau_hadron_phi);
  plotNames.push_back("tau_hadron_phi");
  THStack *hs_tau_hadron_rhomass[2];
  for (int j=0; j < 2;j++) hs_tau_hadron_rhomass[j] = new THStack(("hs_tau_hadron_rhomass" + to_string(j)).c_str(),("#tau_{3prong} #rho mass" + to_string(j) + "[GeV]").c_str());
  TH1F *h_tau_hadron_rhomass[nSamples][2];
  TH2F *h_tau_hadron_rhomass2D[nSamples];
  THStack *hs_tau_hadron_nch = new THStack("hs_tau_hadron_nch","#tau_{3prong} nch");
  TH1F *h_tau_hadron_nch[nSamples];
  TH1F *h_tau_hadron_nch_highHF[nSamples];
  TH1F *h_tau_hadron_nPions[nSamples];
  THStack *hs_tau_hadron_ncand_final = new THStack("hs_tau_hadron_ncand_final","#tau_{3prong} cands final");
  TH1F *h_tau_hadron_ncand_final[nSamples];
  THStack *hs_tau_hadron_vprob = new THStack("hs_tau_hadron_vprob","#tau_{3prong} vprob (%)");
  TH1F *h_tau_hadron_vprob[nSamples+4];
  histograms.push_back(h_tau_hadron_vprob);
  plotNames.push_back("tau_hadron_vprob");
  TH1F *h_tau_hadron_matched_pt_index[nSamples+4];
  THStack *hs_tau_hadron_mass = new THStack("hs_tau_hadron_mass","#tau_{3prong} mass [GeV]");
  TH1F *h_tau_hadron_mass[nSamples+4];
  histograms.push_back(h_tau_hadron_mass);
  plotNames.push_back("tau_hadron_mass");
  TH1F *h_ditau_p[nSamples+4];
  histograms.push_back(h_ditau_p);
  plotNames.push_back("ditau_p");
  TH1F *h_ditau_pz[nSamples+4];
  histograms.push_back(h_ditau_pz);
  plotNames.push_back("ditau_pz");
  TH1F *h_ditau_pt[nSamples+4];
  histograms.push_back(h_ditau_pt);
  plotNames.push_back("ditau_pt");
  TH1F *h_ditau_ptScalar[nSamples+4];
  histograms.push_back(h_ditau_ptScalar);
  plotNames.push_back("ditau_scalar_pt");
  THStack *hs_ditau_mass = new THStack("hs_tau_hadron_mass","#tau#tau mass [GeV]");
  TH1F *h_ditau_mass[nSamples+4];
  histograms.push_back(h_ditau_mass);
  plotNames.push_back("ditau_mass");
  THStack *hs_resVisTauPt = new THStack("hs_resTauVis","visible #tau p_{T} resolution [GeV]");
  TH1F *h_resVisTauPt[nSamples+4];
  TH1F *h_resVisTauEta[nSamples+4];
  TH1F *h_resVisTauPhi[nSamples+4];
  THStack *hs_tau_hadron_track_pvz[3];
  for (int j=0; j < 3;j++) hs_tau_hadron_track_pvz[j] = new THStack(("h_tau_hadron_track" + to_string(j) + "_pvz").c_str(),("#tau_{3prong} track" + to_string(j) + " pvz").c_str());
  TH1F *h_tau_hadron_track_pvz[nSamples][3];
  THStack *hs_deltaphi_tau_mu_tau_hadron = new THStack("hs_deltaphi_tau_mu_tau_hadron","#Delta#phi(#tau_{#mu}, #tau_{3prong})");
  TH1F *h_deltaphi_tau_mu_tau_hadron[nSamples+4];
  histograms.push_back(h_deltaphi_tau_mu_tau_hadron);
  plotNames.push_back("delta_phi");
  TH1F *h_deltaphi_tau_mu_tau_hadron_zoomed[nSamples+4];
  histograms.push_back(h_deltaphi_tau_mu_tau_hadron_zoomed);
  plotNames.push_back("delta_phi_zoomed");
  TH1F *h_deltaphi_tau_mu_full_tau_hadron[nSamples];
  
  TH1F *h_minPionPionDeltaPhi[nSamples+4];
  histograms.push_back(h_minPionPionDeltaPhi);
  plotNames.push_back("minPionPionDeltaPhi");
  TH1F *h_maxPionPionDeltaPhi[nSamples+4];
  histograms.push_back(h_maxPionPionDeltaPhi);
  plotNames.push_back("maxPionPionDeltaPhi");
  TH1F *h_minPionPionDeltaEta[nSamples+4];
  histograms.push_back(h_minPionPionDeltaEta);
  plotNames.push_back("minPionPionDeltaEta");
  TH1F *h_maxPionPionDeltaEta[nSamples+4];
  histograms.push_back(h_maxPionPionDeltaEta);
  plotNames.push_back("maxPionPionDeltaEta");
  TH1F *h_minPionPionDeltaR[nSamples+4];
  histograms.push_back(h_minPionPionDeltaR);
  plotNames.push_back("minPionPionDeltaR");
  TH1F *h_maxPionPionDeltaR[nSamples+4];
  histograms.push_back(h_maxPionPionDeltaR);
  plotNames.push_back("maxPionPionDeltaR");
  
  TH1F *h_minPionMuDeltaPhi[nSamples+4];
  histograms.push_back(h_minPionMuDeltaPhi);
  plotNames.push_back("minPionMuDeltaPhi");
  TH1F *h_maxPionMuDeltaPhi[nSamples+4];
  histograms.push_back(h_maxPionMuDeltaPhi);
  plotNames.push_back("maxPionMuDeltaPhi");
  TH1F *h_minPionMuDeltaEta[nSamples+4];
  histograms.push_back(h_minPionMuDeltaEta);
  plotNames.push_back("minPionMuDeltaEta");
  TH1F *h_maxPionMuDeltaEta[nSamples+4];
  histograms.push_back(h_maxPionMuDeltaEta);
  plotNames.push_back("maxPionMuDeltaEta");
  TH1F *h_minPionMuDeltaR[nSamples+4];
  histograms.push_back(h_minPionMuDeltaR);
  plotNames.push_back("minPionMuDeltaR");
  TH1F *h_maxPionMuDeltaR[nSamples+4];
  histograms.push_back(h_maxPionMuDeltaR);
  plotNames.push_back("maxPionMuDeltaR");
  
  THStack *hs_PV_N = new THStack("hs_PV_N","number of PV");
  TH1F *h_PV_N[nSamples+4];
  TH1F *h_sumZDCplus[nSamples+4];
  histograms.push_back(h_sumZDCplus);
  plotNames.push_back("sumZDCplus");
  TH1F *h_sumZDCminus[nSamples+4];
  histograms.push_back(h_sumZDCminus);
  plotNames.push_back("sumZDCminus");
  TH2F *h_sumZDC_pm[nSamples];
  TH2F *h_muEta_averageZDCside[nSamples];
  TH2F *h_averageZDCside_averageHFeta[nSamples];
  TH2F *h_tauEta_averageZDCside[nSamples];
  THStack *hs_calo_energyHFp = new THStack("hs_calo_energyHFp","energy HF+ [GeV]");
  TH1F *h_calo_energyHFp[nSamples];
  THStack *hs_calo_leadingHFp = new THStack("hs_calo_leadingHFp","energy leading tower HF+ [GeV]");
  TH1F *h_calo_leadingHFp[nSamples];
  TH1F *h_calo_leadingHFp_highNch[nSamples];
  TH1F *h_nHitsPixel[nSamples+4];
  histograms.push_back(h_nHitsPixel);
  plotNames.push_back("nHitsPixel");
  //TProfile *h_nTrkHitsPixelEta[nSamples+4];
  //histograms.push_back(h_nTrkHitsPixelEta);
  plotNames.push_back("nTrkHitsPixel");
  TH1F *h_calo_Et[nSamples+4];
  histograms.push_back(h_calo_Et);
  plotNames.push_back("calo_Et");
  TH2F *h_calo_Et_eta[nSamples];
  TH1F *h_calo_sumEt[nSamples+4];
  histograms.push_back(h_calo_sumEt);
  plotNames.push_back("calo_sumEt");
  TH1F *h_calo_leadingEt[nSamples+4];
  histograms.push_back(h_calo_leadingEt);
  plotNames.push_back("calo_leadingEt");
  TH2F *h_calo_leadingEt_eta[nSamples];
  TH1F *h_calo_E[nSamples+4];
  histograms.push_back(h_calo_E);
  plotNames.push_back("calo_E");
  TH2F *h_calo_E_eta[nSamples];
  TH1F *h_calo_sumE[nSamples+4];
  histograms.push_back(h_calo_sumE);
  plotNames.push_back("calo_sumE");
  TH1F *h_calo_leadingE[nSamples+4];
  histograms.push_back(h_calo_leadingE);
  plotNames.push_back("calo_leadingE");
  TH2F *h_calo_leadingE_eta[nSamples];
  TH1F *h_nCaloTowers[nSamples+4];
  histograms.push_back(h_nCaloTowers);
  plotNames.push_back("nCaloTowers");
  TH1F *h_calo_energyHFp_sum[nSamples];
  TH1F *h_calo_energyHFp_size[nSamples];
  THStack *hs_calo_energyHFm = new THStack("hs_calo_energyHFm","energy HF- [GeV]");
  TH1F *h_calo_energyHFm[nSamples];
  THStack *hs_calo_leadingHFm = new THStack("hs_calo_leadingHFm","energy leading tower HF- [GeV]");
  TH1F *h_calo_leadingHFm[nSamples];
  TH1F *h_calo_leadingHFm_highNch[nSamples];
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
  
  TH1F *h_N_recoMuonPt[nSamples];
  TH1F *h_eff_trigger[nSamples];
  TH1F *h_N_genTauHadpt[nSamples];
  TH1F *h_eff_3prongReco_genTauHadpt[nSamples];
  TH1F *h_eff_tauReco_genTauHadpt[nSamples];
  TH1F *h_N_genPionPt[nSamples];
  TH1F *h_N_matched_genPionPt[nSamples];
  TH1F *h_N_matched_leading_genPionPt[nSamples];
  TH1F *h_N_matched_subleading_genPionPt[nSamples];
  TH1F *h_N_matched_subsubleading_genPionPt[nSamples];
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
  
  TH1F *h_recoPionPt[nSamples];
  TH1F *h_leading_recoPionPt[nSamples];
  TH1F *h_subleading_recoPionPt[nSamples];
  TH1F *h_subsubleading_recoPionPt[nSamples];
  
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
  
  TH1F *h_pion_leading_pt[nSamples+4];
  histograms.push_back(h_pion_leading_pt);
  plotNames.push_back("pion_leading_pt");
  TH1F *h_pion_subleading_pt[nSamples+4];
  histograms.push_back(h_pion_subleading_pt);
  plotNames.push_back("pion_subleading_pt");
  TH1F *h_pion_subsubleading_pt[nSamples+4];
  histograms.push_back(h_pion_subsubleading_pt);
  plotNames.push_back("pion_subsubleading_pt");
  
  TH1F *h_pion_leading_eta[nSamples+4];
  histograms.push_back(h_pion_leading_eta);
  plotNames.push_back("pion_leading_eta");
  TH1F *h_pion_subleading_eta[nSamples+4];
  histograms.push_back(h_pion_subleading_eta);
  plotNames.push_back("pion_subleading_eta");
  TH1F *h_pion_subsubleading_eta[nSamples+4];
  histograms.push_back(h_pion_subsubleading_eta);
  plotNames.push_back("pion_subsubleading_eta");
  
  TH1F *h_pion_leading_phi[nSamples+4];
  histograms.push_back(h_pion_leading_phi);
  plotNames.push_back("pion_leading_phi");
  TH1F *h_pion_subleading_phi[nSamples+4];
  histograms.push_back(h_pion_subleading_phi);
  plotNames.push_back("pion_subleading_phi");
  TH1F *h_pion_subsubleading_phi[nSamples+4];
  histograms.push_back(h_pion_subsubleading_phi);
  plotNames.push_back("pion_subsubleading_phi");
  
  TH1F *h_NrecoGamma[nSamples+4];
  //histograms.push_back(h_NrecoGamma);
  TH1F *h_recoGamma_pt[nSamples+4];
  //histograms.push_back(h_recoGamma_pt);
  TH1F *h_recoGamma_eta[nSamples+4];
  //histograms.push_back(h_recoGamma_eta);
  TH1F *h_recoGamma_phi[nSamples+4];
  //histograms.push_back(h_recoGamma_phi);
  TH1F *h_recoGamma_deltaphi_muon[nSamples+4];
  //histograms.push_back(h_recoGamma_deltaphi_muon);
  TH1F *h_recoGamma_deltaphi_pion[nSamples+4];
  //histograms.push_back(h_recoGamma_deltaphi_pion);
  TH1F *h_recoGammasMinDeltaR[nSamples+4];
  //histograms.push_back(h_recoGammasMinDeltaR);
  TH1F *h_recoPiZeroMinDeltaM[nSamples+4];
  //histograms.push_back(h_recoPiZeroMinDeltaM);
  TH1F *h_recoPiZeroDeltaM[nSamples+4];
  //histograms.push_back(h_recoPiZeroDeltaM);
  TH1F *h_NrecoPiZero[nSamples+4];
  //histograms.push_back(h_NrecoPiZero);
  TH1F *h_recoPiZero_pt[nSamples+4];
  //histograms.push_back(h_recoPiZero_pt);
  TH1F *h_recoPiZero_eta[nSamples+4];
  //histograms.push_back(h_recoPiZero_eta);
  TH1F *h_recoPiZero_phi[nSamples+4];
  //histograms.push_back(h_recoPiZero_phi);
  TH1F *h_recoPiZero_deltaphi_muon[nSamples+4];
  //histograms.push_back(h_recoPiZero_deltaphi_muon);
  TH1F *h_recoPiZero_deltaphi_pion[nSamples+4];
  //histograms.push_back(h_recoPiZero_deltaphi_pion);
  TH2F *h_reco_pion_energy_HCAL_ECAL[nSamples];
  
  TH1F *A_highNch_highHF[nSamples+2+nNchCategories];
  TH1F *B_lowNch_highHF[nSamples+2+nNchCategories];
  TH1F *C_highNch_lowHF[nSamples+2+nNchCategories];
  TH1F *D_lowNch_lowHF[nSamples+2+nNchCategories];
  
  // ABCD validation
  TH1F *ABCD_validation[6];
  ABCD_validation[0] = new TH1F("ABCD_validation_highNch_highHF","ABCD validation highNch highHF;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[1] = new TH1F("ABCD_validation_lowNch_highHF","ABCD validation lowNch highHF;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[2] = new TH1F("ABCD_validation_highNch_lowHF","ABCD validation highNch lowHF;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[3] = new TH1F("ABCD_validation_lowNch_lowHF","ABCD validation lowNch lowHF;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[4] = new TH1F("ABCD_validation_highNch_HF_ratio","ABCD validation highNch HF ratio;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[5] = new TH1F("ABCD_validation_lowNch_HF_ratio","ABCD validation lowNch HF ratio;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  
  for (int cat = 0; cat < 4; cat++){
    ABCD_validation[cat]->SetLineColor(colors[cat+1]); ABCD_validation[cat]->SetMarkerColor(colors[cat+1]);
    ABCD_validation[cat]->Sumw2();
  }
  for (int cat = 4; cat < 6; cat++){
    ABCD_validation[cat]->SetLineColor(colors[cat+1-4]); ABCD_validation[cat]->SetMarkerColor(colors[cat+1-4]);
    ABCD_validation[cat]->Sumw2();
  }
  
  // SF effect
  TH1F *h_muon_pt_before_SF = new TH1F("h_muon_pt_before_SF","visible #tau_{#mu} p_{T} before SF; reco. muon p_{T}",16, -0.5, 15.5);
  h_muon_pt_before_SF->Sumw2(); h_muon_pt_before_SF->SetLineColor(kBlue); h_muon_pt_before_SF->SetMarkerStyle(1);
  TH1F *h_muon_pt_after_SF = new TH1F("h_muon_pt_after_SF","visible #tau_{#mu} p_{T} after SF; reco. muon p_{T}",16, -0.5, 15.5);
  h_muon_pt_after_SF->Sumw2(); h_muon_pt_after_SF->SetLineColor(kGreen); h_muon_pt_after_SF->SetMarkerStyle(1);
  
  
  // MET
  TH1F *h_MET[nSamples+4];
  histograms.push_back(h_MET);
  plotNames.push_back("MET");
  
  TH1F *h_sys_muon_SF_Down = new TH1F("muon-SFDown","muon SF uncertainty down;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  TH1F *h_sys_muon_SF_Up = new TH1F("muon-SFUp","muon SF uncertainty up;#Delta#phi(#tau_{#mu}, #tau_{3prong})", deltaphi_bins, 0, TMath::Pi());
  h_sys_muon_SF_Down->Sumw2(); h_sys_muon_SF_Up->Sumw2();
  
  
  for (int i = 0; i < nSamples; i++){
    
    h_cutflow[i] = new TH1F(("h_cutflow_" + tag[i]).c_str(),("Analysis cutflow - " + tag[i]).c_str(),8, 0, 8);
    std::string cutflow_bins_string[] = {"input from Ntuplizer", "#tau #mu", "nCh", "#tau_{3prong}", "HF", "HF & #tau_{3prong}", "HF & #tau_{3prong} & nch", "..."};
    for(size_t j=0; j< 8; j++){
      h_cutflow[0]->GetXaxis()->SetBinLabel(j+1, (cutflow_bins_string[j]).c_str());
    }
    h_cutflow[i]->Sumw2();
    /*if (i != 0) {*/h_cutflow[i]->SetLineColor(colors[i]); h_cutflow[i]->SetMarkerStyle(styles[i]);
  
    A_highNch_highHF[i] = new TH1F(("A_highNch_highHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - high HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    A_highNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - high HF"); A_highNch_highHF[i]->Sumw2();
    A_highNch_highHF[i]->SetLineColor(colors[i]); A_highNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    B_lowNch_highHF[i] = new TH1F(("B_lowNch_highHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - high HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    B_lowNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - high HF"); B_lowNch_highHF[i]->Sumw2();
    B_lowNch_highHF[i]->SetLineColor(colors[i]); B_lowNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    C_highNch_lowHF[i] = new TH1F(("C_highNch_lowHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - low HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    C_highNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - low HF"); C_highNch_lowHF[i]->Sumw2();
    C_highNch_lowHF[i]->SetLineColor(colors[i]); C_highNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
    D_lowNch_lowHF[i] = new TH1F(("D_lowNch_lowHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - low HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    D_lowNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - low HF"); D_lowNch_lowHF[i]->Sumw2();
    D_lowNch_lowHF[i]->SetLineColor(colors[i]); D_lowNch_lowHF[i]->SetMarkerStyle(styles[i]);
    
    h_pion_leading_pt[i] = new TH1F(("h_pion_leading_pt_" + tag[i]).c_str(),("leading pion pt " + tag[i]).c_str(),10, 0, 10*pionLeadingCut);
    h_pion_leading_pt[i]->SetXTitle("leading pion pt [GeV]"); h_pion_leading_pt[i]->Sumw2();
    /*if (i != 0) {*/h_pion_leading_pt[i]->SetLineColor(colors[i]); h_pion_leading_pt[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_pt[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subleading_pt[i] = new TH1F(("h_pion_subleading_pt_" + tag[i]).c_str(),("subleading pion pt " + tag[i]).c_str(),10, 0, 10*pionSubLeadingCut);
    h_pion_subleading_pt[i]->SetXTitle("subleading pion pt [GeV]"); h_pion_subleading_pt[i]->Sumw2();
    /*if (i != 0) {*/h_pion_subleading_pt[i]->SetLineColor(colors[i]); h_pion_subleading_pt[i]->SetMarkerStyle(styles[i]);
    h_pion_subleading_pt[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_subsubleading_pt[i] = new TH1F(("h_pion_subsubleading_pt_" + tag[i]).c_str(),("subsubleading pion pt " + tag[i]).c_str(),10, 0, 10*pionSubSubLeadingCut);
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
    
    h_NrecoGamma[i] = new TH1F(("h_NrecoGamma_" + tag[i]).c_str(),("total number of reco gammas " + tag[i]).c_str(),30, -0.5, 29.5);
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
    
    h_NrecoPiZero[i] = new TH1F(("h_NrecoPiZero_" + tag[i]).c_str(),("number of reco PiZero candidates " + tag[i]).c_str(),20, -0.5, 19.5);
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
    
    h_NgenGamma[i] = new TH1F(("h_NgenGamma_" + tag[i]).c_str(),("total number of gen gammas " + tag[i]).c_str(),30, -0.5, 29.5);
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
    
    h_NgenPiZero[i] = new TH1F(("h_NgenPiZero_" + tag[i]).c_str(),("number of gen PiZero candidates " + tag[i]).c_str(),20, -0.5, 19.5);
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
    
    h_recoPionPt[i] = new TH1F(("h_recoPionPt_" + tag[i]).c_str(),("pion p_{T} - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_recoPionPt[i]->SetXTitle("pion p_{T} [GeV]"); h_recoPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_recoPionPt[i]->SetLineColor(colors[i]); h_recoPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_leading_recoPionPt[i] = new TH1F(("h_leading_recoPionPt_" + tag[i]).c_str(),("leading pion p_{T} - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_leading_recoPionPt[i]->SetXTitle("leading pion p_{T} [GeV]"); h_leading_recoPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_leading_recoPionPt[i]->SetLineColor(colors[i]); h_leading_recoPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_subleading_recoPionPt[i] = new TH1F(("h_subleading_recoPionPt_" + tag[i]).c_str(),("subleading pion p_{T} - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_subleading_recoPionPt[i]->SetXTitle("subleading pion p_{T} [GeV]"); h_subleading_recoPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_subleading_recoPionPt[i]->SetLineColor(colors[i]); h_subleading_recoPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_subsubleading_recoPionPt[i] = new TH1F(("h_subsubleading_recoPionPt_" + tag[i]).c_str(),("subsubleading pion p_{T} - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_subsubleading_recoPionPt[i]->SetXTitle("subsubleading pion p_{T} [GeV]"); h_subsubleading_recoPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_subsubleading_recoPionPt[i]->SetLineColor(colors[i]); h_subsubleading_recoPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_matched_genPionPt[i] = new TH1F(("h_N_matched_genPionPt_" + tag[i]).c_str(),("Number of matched gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_matched_genPionPt[i]->SetXTitle("matched pion gen p_{T} [GeV]"); h_N_matched_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_matched_genPionPt[i]->SetLineColor(colors[i]); h_N_matched_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_matched_leading_genPionPt[i] = new TH1F(("h_N_matched_leading_genPionPt_" + tag[i]).c_str(),("Number of matched leading gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_matched_leading_genPionPt[i]->SetXTitle("leading matched pion gen p_{T} [GeV]"); h_N_matched_leading_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_matched_leading_genPionPt[i]->SetLineColor(colors[i]); h_N_matched_leading_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_matched_subleading_genPionPt[i] = new TH1F(("h_N_matched_subleading_genPionPt_" + tag[i]).c_str(),("Number of matched subleading gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_matched_subleading_genPionPt[i]->SetXTitle("subleading matched pion gen p_{T} [GeV]"); h_N_matched_subleading_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_matched_subleading_genPionPt[i]->SetLineColor(colors[i]); h_N_matched_subleading_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_matched_subsubleading_genPionPt[i] = new TH1F(("h_N_matched_subsubleading_genPionPt_" + tag[i]).c_str(),("Number of matched subsubleading gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_matched_subsubleading_genPionPt[i]->SetXTitle("subsubleading matched pion gen p_{T} [GeV]"); h_N_matched_subsubleading_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_matched_subsubleading_genPionPt[i]->SetLineColor(colors[i]); h_N_matched_subsubleading_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_genPionPt[i] = new TH1F(("h_N_genPionPt_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_genPionPt[i]->SetXTitle("pion gen p_{T} [GeV]"); h_N_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genPionPt[i]->SetLineColor(colors[i]); h_N_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_genCentralPionPt[i] = new TH1F(("h_N_genCentralPionPt_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_genCentralPionPt[i]->SetXTitle("pion gen p_{T} [GeV]"); h_N_genCentralPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genCentralPionPt[i]->SetLineColor(colors[i]); h_N_genCentralPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedPionReco_genPionPt[i] = new TH1F(("h_eff_matchedPionReco_genPionPt_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_eff_matchedPionReco_genPionPt[i]->SetXTitle("p_{T} [GeV]");
    h_eff_matchedPionReco_genPionPt[i]->SetYTitle("pion efficiency");
    h_eff_matchedPionReco_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedPionReco_genPionPt[i]->SetLineColor(colors[i]); h_eff_matchedPionReco_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedCentralPionReco_genPionPt[i] = new TH1F(("h_eff_matchedCentralPionReco_genPionPt_" + tag[i]).c_str(),("Efficiency of reconstructing central matched pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_eff_matchedCentralPionReco_genPionPt[i]->SetXTitle("pion gen p_{T} [GeV]");
    h_eff_matchedCentralPionReco_genPionPt[i]->SetYTitle("central matched pion efficiency");
    h_eff_matchedCentralPionReco_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedCentralPionReco_genPionPt[i]->SetLineColor(colors[i]); h_eff_matchedCentralPionReco_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_genPionEta[i] = new TH1F(("h_N_genPionEta_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),27, -2.7, 2.7);
    h_N_genPionEta[i]->SetXTitle("pion gen #eta"); h_N_genPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_N_genPionEta[i]->SetLineColor(colors[i]); h_N_genPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_N_genHighPtPionEta[i] = new TH1F(("h_N_genHighPtPionEta_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),27, -2.7, 2.7);
    h_N_genHighPtPionEta[i]->SetXTitle("pion gen #eta"); h_N_genHighPtPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_N_genHighPtPionEta[i]->SetLineColor(colors[i]); h_N_genHighPtPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedPionReco_genPionEta[i] = new TH1F(("h_eff_matchedPionReco_genPionEta_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),27, -2.7, 2.7);
    h_eff_matchedPionReco_genPionEta[i]->SetXTitle("#eta");
    h_eff_matchedPionReco_genPionEta[i]->SetYTitle("pion efficiency");
    h_eff_matchedPionReco_genPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedPionReco_genPionEta[i]->SetLineColor(colors[i]); h_eff_matchedPionReco_genPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_eff_matchedHighPtPionReco_genPionEta[i] = new TH1F(("h_eff_matchedHighPtPionReco_genPionEta_" + tag[i]).c_str(),("Efficiency of reconstructing high P_{T} matched pions - " + tag[i]).c_str(),27, -2.7, 2.7);
    h_eff_matchedHighPtPionReco_genPionEta[i]->SetXTitle("pion gen #eta");
    h_eff_matchedHighPtPionReco_genPionEta[i]->SetYTitle("high p_{T} matched pion efficiency");
    h_eff_matchedHighPtPionReco_genPionEta[i]->Sumw2();
    /*if (i != 0) {*/h_eff_matchedHighPtPionReco_genPionEta[i]->SetLineColor(colors[i]); h_eff_matchedHighPtPionReco_genPionEta[i]->SetMarkerStyle(styles[i]);
    
    h_N_genPionPtEta[i] = new TH2F(("h_N_genPionPtEta_" + tag[i]).c_str(),("Number of gen pions - " + tag[i]).c_str(),30, -3, 3, 40, -0.1, 7.9);
    h_N_genPionPtEta[i]->SetXTitle("pion gen #eta"); h_N_genPionPtEta[i]->SetYTitle("pion gen p_{T} [GeV]");
    
    h_eff_matchedPionReco_genPionPtEta[i] = new TH2F(("h_eff_matchedPionReco_genPionPtEta_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),30, -3, 3, 40, -0.1, 7.9);
    h_eff_matchedPionReco_genPionPtEta[i]->SetXTitle("pion gen #eta"); h_eff_matchedPionReco_genPionPtEta[i]->SetYTitle("pion gen p_{T} [GeV]");
    
    h_eff_matchedPionReco_genPionPtEtaZoomed[i] = new TH2F(("h_eff_matchedPionReco_genPionPtEtaZoomed_" + tag[i]).c_str(),("Efficiency of reconstructing matched pions - " + tag[i]).c_str(),30, -3, 3, 40, -0.1, 7.9);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->SetXTitle("pion gen #eta"); h_eff_matchedPionReco_genPionPtEtaZoomed[i]->SetYTitle("pion gen p_{T} [GeV]");
    
    h_resPt_matchedPionReco[i] = new TH1F(("h_resPt_matchedPionReco_" + tag[i]).c_str(),("relative p_{T} resolution of reconstructed matched pions - " + tag[i]).c_str(),31, 0.8, 1.3);
    h_resPt_matchedPionReco[i]->SetXTitle("reco pion p_{T} / gen pion p_{T} [GeV]"); h_resPt_matchedPionReco[i]->Sumw2();
    /*if (i != 0) {*/h_resPt_matchedPionReco[i]->SetLineColor(colors[i]); h_resPt_matchedPionReco[i]->SetMarkerStyle(styles[i]);
    
    h_resR_matchedPionReco[i] = new TH1F(("h_resR_matchedPionReco_" + tag[i]).c_str(),("#DeltaR resolution of reconstructed matched pions - " + tag[i]).c_str(),20, 0, 2);
    h_resR_matchedPionReco[i]->SetXTitle("100 #times #DeltaR matched gen & reco pions [GeV]"); h_resR_matchedPionReco[i]->Sumw2();
    /*if (i != 0) {*/h_resR_matchedPionReco[i]->SetLineColor(colors[i]); h_resR_matchedPionReco[i]->SetMarkerStyle(styles[i]);
    
    h_resPhi_matchedPionReco[i] = new TH2F(("h_resPhi_matchedPionReco_" + tag[i]).c_str(),("#phi resolution of reconstructed matched pions - " + tag[i]).c_str(),20, -0.02, 0.02, 20, -0.02, 0.02);
    h_resPhi_matchedPionReco[i]->SetXTitle("#phi matched gen pions [GeV]");
    h_resPhi_matchedPionReco[i]->SetYTitle("#phi matched reco pions [GeV]");
    
    h_N_genTauHadpt[i] = new TH1F(("h_N_genTauHadpt_" + tag[i]).c_str(),("Number of events with a reconstructed muon - " + tag[i]).c_str(),20, 0, 20);
    h_N_genTauHadpt[i]->SetXTitle("#tau_{3prong} gen p_{T} [GeV]"); h_N_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genTauHadpt[i]->SetLineColor(colors[i]); h_N_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_3prongReco_genTauHadpt[i] = new TH1F(("h_eff_3prongReco_genTauHadpt_" + tag[i]).c_str(),("Efficiency of reconstructing 3 pions - " + tag[i]).c_str(),20, 0, 20);
    h_eff_3prongReco_genTauHadpt[i]->SetXTitle("#tau_{3prong} gen p_{T} [GeV]");
    h_eff_3prongReco_genTauHadpt[i]->SetYTitle("Efficiency of reconstructing 3 pions");
    h_eff_3prongReco_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_3prongReco_genTauHadpt[i]->SetLineColor(colors[i]); h_eff_3prongReco_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_tauReco_genTauHadpt[i] = new TH1F(("h_eff_tauReco_genTauHadpt_" + tag[i]).c_str(),("Efficiency of reconstructing a hadronic tau - " + tag[i]).c_str(),20, 0, 20);
    h_eff_tauReco_genTauHadpt[i]->SetXTitle("#tau_{3prong} gen p_{T} [GeV]");
    h_eff_tauReco_genTauHadpt[i]->SetYTitle("#tau_{3prong} reconstruction efficiency");
    h_eff_tauReco_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_tauReco_genTauHadpt[i]->SetLineColor(colors[i]); h_eff_tauReco_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_N_recoMuonPt[i] = new TH1F(("h_N_recoMuonPt_" + tag[i]).c_str(),("Events passing all analysis cuts " + tag[i]).c_str(),15, 0.5, 15.5);
    h_N_recoMuonPt[i]->SetXTitle("#mu reco p_{T} [GeV]"); h_N_recoMuonPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_recoMuonPt[i]->SetLineColor(colors[i]); h_N_recoMuonPt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_trigger[i] = new TH1F(("h_eff_trigger_" + tag[i]).c_str(),("Trigger Efficiency " + tag[i]).c_str(),15, 0.5, 15.5);
    h_eff_trigger[i]->SetXTitle("#mu reco p_{T} [GeV]"); h_eff_trigger[i]->SetYTitle("trigger efficiency"); h_eff_trigger[i]->Sumw2();
    /*if (i != 0) {*/h_eff_trigger[i]->SetLineColor(colors[i]); h_eff_trigger[i]->SetMarkerStyle(styles[i]);
    
    h_tau_mu_p[i] = new TH1F(("h_tau_mu_p_" + tag[i]).c_str(),("visible #tau_{#mu} total p " + tag[i]).c_str(),30, 0, 30);
    h_tau_mu_p[i]->SetXTitle("visible #tau_{#mu} total p [GeV]"); h_tau_mu_p[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_p[i]->SetLineColor(colors[i]); h_tau_mu_p[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_p[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_pz[i] = new TH1F(("h_tau_mu_pz_" + tag[i]).c_str(),("visible #tau_{#mu} p_{z} " + tag[i]).c_str(),30, -30, 30);
    h_tau_mu_pz[i]->SetXTitle("visible #tau_{#mu} p_{z} [GeV]"); h_tau_mu_pz[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_pz[i]->SetLineColor(colors[i]); h_tau_mu_pz[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_pz[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_dz[i] = new TH1F(("h_tau_mu_dz_" + tag[i]).c_str(),("#Deltaz(#mu,PV) " + tag[i]).c_str(),19, -0.95, 0.95);
    h_tau_mu_dz[i]->SetXTitle("#Deltaz(#mu,PV) [mm]"); h_tau_mu_dz[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_dz[i]->SetLineColor(colors[i]); h_tau_mu_dz[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_dz[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_pt[i] = new TH1F(("h_tau_mu_pt_" + tag[i]).c_str(),("visible #tau_{#mu} p_{T} " + tag[i]).c_str(),20, -0.5, 19.5);
    h_tau_mu_pt[i]->SetXTitle("visible #tau_{#mu} p_{T} [GeV]"); h_tau_mu_pt[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_pt[i]->SetLineColor(colors[i]); h_tau_mu_pt[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_pt[i]->SetBinErrorOption(TH1::kPoisson);
        
    h_tau_mu_eta[i] = new TH1F(("h_tau_mu_eta_" + tag[i]).c_str(),("visible #tau_{#mu} #eta " + tag[i]).c_str(),tau_mu_eta_bins, -2.5, 2.5);
    h_tau_mu_eta[i]->SetXTitle("visible #tau_{#mu} #eta"); h_tau_mu_eta[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_eta[i]->SetLineColor(colors[i]); h_tau_mu_eta[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_eta[i]->SetBinErrorOption(TH1::kPoisson);
        
    h_tau_mu_phi[i] = new TH1F(("h_tau_mu_phi_" + tag[i]).c_str(),("visible #tau_{#mu} #phi " + tag[i]).c_str(),tau_mu_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_mu_phi[i]->SetXTitle("visible #tau_{#mu} #phi"); h_tau_mu_phi[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_phi[i]->SetLineColor(colors[i]); h_tau_mu_phi[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_phi[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_tau_hadron_phi[i] = new TH2F(("h_tau_mu_tau_hadron_phi_" + tag[i]).c_str(),("visible #tau_{#mu}-#tau_{3prong} #phi " + tag[i] +";visible #tau_{#mu} #phi;visible #tau_{3prong} #phi").c_str(),tau_mu_phi_bins, -TMath::Pi(), TMath::Pi(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    
    // hadron related
    
    h_tau_hadron_p[i] = new TH1F(("h_tau_hadron_p_" + tag[i]).c_str(),("#tau_{3prong} vector sum p " + tag[i]).c_str(),10, 0, 20);
    h_tau_hadron_p[i]->SetXTitle("#tau_{3prong} vector sum p [GeV]"); h_tau_hadron_p[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_p[i]->SetLineColor(colors[i]); h_tau_hadron_p[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_pz[i] = new TH1F(("h_tau_hadron_pz_" + tag[i]).c_str(),("#tau_{3prong} p_{z} " + tag[i]).c_str(),20, -10, 20);
    h_tau_hadron_pz[i]->SetXTitle("#tau_{3prong} p_{z} [GeV]"); h_tau_hadron_pz[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_pz[i]->SetLineColor(colors[i]); h_tau_hadron_pz[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_ptScalar[i] = new TH1F(("h_tau_hadron_ptScalar_" + tag[i]).c_str(),("3prong scalar sum p_{T} " + tag[i]).c_str(),10, 0, 20);
    h_tau_hadron_ptScalar[i]->SetXTitle("3prong scalar sum p_{T} [GeV]"); h_tau_hadron_ptScalar[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_ptScalar[i]->SetLineColor(colors[i]); h_tau_hadron_ptScalar[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_pt[i] = new TH1F(("h_tau_hadron_pt_" + tag[i]).c_str(),("#tau_{3prong} p_{T} " + tag[i]).c_str(),20, 0, 20);
    h_tau_hadron_pt[i]->SetXTitle("#tau_{3prong} p_{T} [GeV]"); h_tau_hadron_pt[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_pt[i]->SetLineColor(colors[i]); h_tau_hadron_pt[i]->SetMarkerStyle(styles[i]);
    
    h_gen_tau_hadron_visible_pt[i] = new TH1F(("h_gen_tau_hadron_visible_pt_" + tag[i]).c_str(),("gen #tau_{3prong} visible p_{T} " + tag[i]).c_str(),20, 0, 20);
    h_gen_tau_hadron_visible_pt[i]->SetXTitle("gen #tau_{3prong} p_{T} [GeV]"); h_gen_tau_hadron_visible_pt[i]->Sumw2();
    /*if (i != 0) {*/h_gen_tau_hadron_visible_pt[i]->SetLineColor(colors[i]); h_gen_tau_hadron_visible_pt[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_eta[i] = new TH1F(("h_tau_hadron_eta_" + tag[i]).c_str(),("#tau_{3prong} #eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_tau_hadron_eta[i]->SetXTitle("#tau_{3prong} #eta"); h_tau_hadron_eta[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_eta[i]->SetLineColor(colors[i]); h_tau_hadron_eta[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_phi[i] = new TH1F(("h_tau_hadron_phi_" + tag[i]).c_str(),("#tau_{3prong} #phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_hadron_phi[i]->SetXTitle("#tau_{3prong} #phi"); h_tau_hadron_phi[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_phi[i]->SetLineColor(colors[i]); h_tau_hadron_phi[i]->SetMarkerStyle(styles[i]);
        
    for (int j=0; j < 2;j++){
      h_tau_hadron_rhomass[i][j] = new TH1F(("h_tau_hadron_rhomass" + to_string(j) + "_" + tag[i]).c_str(),("#tau_{3prong} #rho_{"+to_string(j+1)+"} mass [GeV] " + tag[i]).c_str(),tau_hadron_rhomass_bins, 0.2, 1.4);
      h_tau_hadron_rhomass[i][j]->SetXTitle(("#rho_{"+to_string(j+1)+"} mass [GeV]").c_str()); h_tau_hadron_rhomass[i][j]->Sumw2();
      /*if (i != 0) {*/h_tau_hadron_rhomass[i][j]->SetLineColor(colors[i]); h_tau_hadron_rhomass[i][j]->SetMarkerStyle(styles[i]);
    }
    
    h_tau_hadron_rhomass2D[i] = new TH2F(("h_tau_hadron_rhomass2D_" + tag[i]).c_str(),("#tau_{3prong} #rho mass [GeV] " + tag[i]).c_str(),2*tau_hadron_rhomass_bins, 0.2, 1.4, 2*tau_hadron_rhomass_bins, 0.2, 1.4);
    h_tau_hadron_rhomass2D[i]->SetXTitle("#rho_{1} mass [GeV]"); h_tau_hadron_rhomass2D[i]->SetYTitle("#rho_{2} mass [GeV]");
    
    h_tau_hadron_nch[i] = new TH1F(("h_tau_hadron_nch_" + tag[i]).c_str(),("h_tau_hadron_nch_" + tag[i]).c_str(),16, -0.5, 15.5);
    h_tau_hadron_nch[i]->SetXTitle("number of charged pions"); h_tau_hadron_nch[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_nch[i]->SetLineColor(colors[i]); h_tau_hadron_nch[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_nch_highHF[i] = new TH1F(("h_tau_hadron_nch_highHF_" + tag[i]).c_str(),("h_tau_hadron_nch_highHF_" + tag[i]).c_str(),16, -0.5, 15.5);
    h_tau_hadron_nch_highHF[i]->SetXTitle("number of charged pions"); h_tau_hadron_nch_highHF[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_nch_highHF[i]->SetLineColor(colors[i]); h_tau_hadron_nch_highHF[i]->SetMarkerStyle(styles[i]);
    
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
    
    h_tau_hadron_mass[i] = new TH1F(("h_tau_hadron_mass_" + tag[i]).c_str(),("#tau_{3prong} visible mass " + tag[i]).c_str(), 16, 0.1, 1.7);
    h_tau_hadron_mass[i]->SetXTitle("visible #tau_{3prong} mass [GeV]"); h_tau_hadron_mass[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_mass[i]->SetLineColor(colors[i]); h_tau_hadron_mass[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_mass[i] = new TH1F(("h_ditau_mass_" + tag[i]).c_str(),("#tau#tau visible invariant mass " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_mass[i]->SetXTitle("visible #tau#tau invariant mass [GeV]"); h_ditau_mass[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_mass[i]->SetLineColor(colors[i]); h_ditau_mass[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_p[i] = new TH1F(("h_ditau_p_" + tag[i]).c_str(),("#tau#tau visible total p " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_p[i]->SetXTitle("#tau#tau total p [GeV]"); h_ditau_p[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_p[i]->SetLineColor(colors[i]); h_ditau_p[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_pz[i] = new TH1F(("h_ditau_pz_" + tag[i]).c_str(),("#tau#tau visible p_{z} " + tag[i]).c_str(), 40, -40, 40);
    h_ditau_pz[i]->SetXTitle("visible #tau#tau p_{z} [GeV]"); h_ditau_pz[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_pz[i]->SetLineColor(colors[i]); h_ditau_pz[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_pt[i] = new TH1F(("h_ditau_pt_" + tag[i]).c_str(),("#tau#tau visible p_{T} " + tag[i]).c_str(), 10, 0, 10);
    h_ditau_pt[i]->SetXTitle("visible #tau#tau p_{T} [GeV]"); h_ditau_pt[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_pt[i]->SetLineColor(colors[i]); h_ditau_pt[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_ptScalar[i] = new TH1F(("h_ditau_ptScalar_" + tag[i]).c_str(),("#mu + 3prong scalar sum p_{T} " + tag[i]).c_str(), 15, 0, 30);
    h_ditau_ptScalar[i]->SetXTitle("#mu + 3prong scalar sum p_{T} [GeV]"); h_ditau_ptScalar[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_ptScalar[i]->SetLineColor(colors[i]); h_ditau_ptScalar[i]->SetMarkerStyle(styles[i]);
    
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
        
    //h_deltaphi_tau_mu_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) " + tag[i]).c_str(),deltaphi_bins/4, 0.75*TMath::Pi(), TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) " + tag[i]).c_str(),deltaphi_bins, 0.0*TMath::Pi(), TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong})"); h_deltaphi_tau_mu_tau_hadron[i]->Sumw2();
    h_deltaphi_tau_mu_tau_hadron[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_tau_hadron[i]->SetMarkerStyle(styles[i]);
    h_deltaphi_tau_mu_tau_hadron[i]->SetBinErrorOption(TH1::kPoisson);
    /*h_deltaphi_tau_mu_tau_hadron[i]->SetFillColor(i);*/
        
    h_deltaphi_tau_mu_tau_hadron_zoomed[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_zoomed_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) " + tag[i]).c_str(),9, 71*TMath::Pi()/80, TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong})"); h_deltaphi_tau_mu_tau_hadron_zoomed[i]->Sumw2();
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetMarkerStyle(styles[i]);
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetBinErrorOption(TH1::kPoisson);
    /*h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetFillColor(i);*/
        
    h_deltaphi_tau_mu_full_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_full_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}+#pi^{0}(s)) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}+#pi^{0}(s))"); h_deltaphi_tau_mu_full_tau_hadron[i]->Sumw2();
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_full_tau_hadron[i]->SetMarkerStyle(styles[i]);
    /*h_deltaphi_tau_mu_full_tau_hadron[i]->SetFillColor(i);*/
        
    h_minPionPionDeltaPhi[i] = new TH1F(("h_minPionPionDeltaPhi_" + tag[i]).c_str(),("min #Delta#phi(#pi^{#pm},#pi^{#pm}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_minPionPionDeltaPhi[i]->SetXTitle("min #Delta#phi(#pi^{#pm},#pi^{#pm})"); h_minPionPionDeltaPhi[i]->Sumw2();
    h_minPionPionDeltaPhi[i]->SetLineColor(colors[i]); h_minPionPionDeltaPhi[i]->SetMarkerStyle(styles[i]);
    /*h_minPionPionDeltaPhi[i]->SetFillColor(i);*/
        
    h_maxPionPionDeltaPhi[i] = new TH1F(("h_maxPionPionDeltaPhi_" + tag[i]).c_str(),("max #Delta#phi(#pi^{#pm},#pi^{#pm}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_maxPionPionDeltaPhi[i]->SetXTitle("max #Delta#phi(#pi^{#pm},#pi^{#pm})"); h_maxPionPionDeltaPhi[i]->Sumw2();
    h_maxPionPionDeltaPhi[i]->SetLineColor(colors[i]); h_maxPionPionDeltaPhi[i]->SetMarkerStyle(styles[i]);
    /*h_maxPionPionDeltaPhi[i]->SetFillColor(i);*/
        
    h_minPionPionDeltaEta[i] = new TH1F(("h_minPionPionDeltaEta_" + tag[i]).c_str(),("min #Delta#eta(#pi^{#pm},#pi^{#pm}) " + tag[i]).c_str(),deltaEta_bins, 0, 5);
    h_minPionPionDeltaEta[i]->SetXTitle("min #Delta#eta(#pi^{#pm},#pi^{#pm})"); h_minPionPionDeltaEta[i]->Sumw2();
    h_minPionPionDeltaEta[i]->SetLineColor(colors[i]); h_minPionPionDeltaEta[i]->SetMarkerStyle(styles[i]);
    /*h_minPionPionDeltaEta[i]->SetFillColor(i);*/
        
    h_maxPionPionDeltaEta[i] = new TH1F(("h_maxPionPionDeltaEta_" + tag[i]).c_str(),("max #Delta#eta(#pi^{#pm},#pi^{#pm}) " + tag[i]).c_str(),deltaEta_bins, 0, 5);
    h_maxPionPionDeltaEta[i]->SetXTitle("max #Delta#eta(#pi^{#pm},#pi^{#pm})"); h_maxPionPionDeltaEta[i]->Sumw2();
    h_maxPionPionDeltaEta[i]->SetLineColor(colors[i]); h_maxPionPionDeltaEta[i]->SetMarkerStyle(styles[i]);
    /*h_maxPionPionDeltaEta[i]->SetFillColor(i);*/
        
    h_minPionPionDeltaR[i] = new TH1F(("h_minPionPionDeltaR_" + tag[i]).c_str(),("min #DeltaR(#pi^{#pm},#pi^{#pm}) " + tag[i]).c_str(),deltaR_bins, 0, 5);
    h_minPionPionDeltaR[i]->SetXTitle("min #DeltaR(#pi^{#pm},#pi^{#pm})"); h_minPionPionDeltaR[i]->Sumw2();
    h_minPionPionDeltaR[i]->SetLineColor(colors[i]); h_minPionPionDeltaR[i]->SetMarkerStyle(styles[i]);
    /*h_minPionPionDeltaR[i]->SetFillColor(i);*/
        
    h_maxPionPionDeltaR[i] = new TH1F(("h_maxPionPionDeltaR_" + tag[i]).c_str(),("max #DeltaR(#pi^{#pm},#pi^{#pm}) " + tag[i]).c_str(),deltaR_bins, 0, 5);
    h_maxPionPionDeltaR[i]->SetXTitle("max #DeltaR(#pi^{#pm},#pi^{#pm})"); h_maxPionPionDeltaR[i]->Sumw2();
    h_maxPionPionDeltaR[i]->SetLineColor(colors[i]); h_maxPionPionDeltaR[i]->SetMarkerStyle(styles[i]);
    /*h_maxPionPionDeltaR[i]->SetFillColor(i);*/
        
    h_minPionMuDeltaPhi[i] = new TH1F(("h_minPionMuDeltaPhi_" + tag[i]).c_str(),("min #Delta#phi(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_minPionMuDeltaPhi[i]->SetXTitle("min #Delta#phi(#mu,#pi^{#pm})"); h_minPionMuDeltaPhi[i]->Sumw2();
    h_minPionMuDeltaPhi[i]->SetLineColor(colors[i]); h_minPionMuDeltaPhi[i]->SetMarkerStyle(styles[i]);
    /*h_minPionMuDeltaPhi[i]->SetFillColor(i);*/
        
    h_maxPionMuDeltaPhi[i] = new TH1F(("h_maxPionMuDeltaPhi_" + tag[i]).c_str(),("max #Delta#phi(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_maxPionMuDeltaPhi[i]->SetXTitle("max #Delta#phi(#mu,#pi^{#pm})"); h_maxPionMuDeltaPhi[i]->Sumw2();
    h_maxPionMuDeltaPhi[i]->SetLineColor(colors[i]); h_maxPionMuDeltaPhi[i]->SetMarkerStyle(styles[i]);
    /*h_maxPionMuDeltaPhi[i]->SetFillColor(i);*/
        
    h_minPionMuDeltaEta[i] = new TH1F(("h_minPionMuDeltaEta_" + tag[i]).c_str(),("min #Delta#eta(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaEta_bins, 0, 5);
    h_minPionMuDeltaEta[i]->SetXTitle("min #Delta#eta(#mu,#pi^{#pm})"); h_minPionMuDeltaEta[i]->Sumw2();
    h_minPionMuDeltaEta[i]->SetLineColor(colors[i]); h_minPionMuDeltaEta[i]->SetMarkerStyle(styles[i]);
    /*h_minPionMuDeltaEta[i]->SetFillColor(i);*/
        
    h_maxPionMuDeltaEta[i] = new TH1F(("h_maxPionMuDeltaEta_" + tag[i]).c_str(),("max #Delta#eta(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaEta_bins, 0, 5);
    h_maxPionMuDeltaEta[i]->SetXTitle("max #Delta#eta(#mu,#pi^{#pm})"); h_maxPionMuDeltaEta[i]->Sumw2();
    h_maxPionMuDeltaEta[i]->SetLineColor(colors[i]); h_maxPionMuDeltaEta[i]->SetMarkerStyle(styles[i]);
    /*h_maxPionMuDeltaEta[i]->SetFillColor(i);*/
        
    h_minPionMuDeltaR[i] = new TH1F(("h_minPionMuDeltaR_" + tag[i]).c_str(),("min #DeltaR(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaR_bins, 0, 5);
    h_minPionMuDeltaR[i]->SetXTitle("min #DeltaR(#mu,#pi^{#pm})"); h_minPionMuDeltaR[i]->Sumw2();
    h_minPionMuDeltaR[i]->SetLineColor(colors[i]); h_minPionMuDeltaR[i]->SetMarkerStyle(styles[i]);
    /*h_minPionMuDeltaR[i]->SetFillColor(i);*/
        
    h_maxPionMuDeltaR[i] = new TH1F(("h_maxPionMuDeltaR_" + tag[i]).c_str(),("max #DeltaR(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaR_bins, 0, 5);
    h_maxPionMuDeltaR[i]->SetXTitle("max #DeltaR(#mu,#pi^{#pm})"); h_maxPionMuDeltaR[i]->Sumw2();
    h_maxPionMuDeltaR[i]->SetLineColor(colors[i]); h_maxPionMuDeltaR[i]->SetMarkerStyle(styles[i]);
    /*h_maxPionMuDeltaR[i]->SetFillColor(i);*/
    
    h_deltaphi_tau_mu_tau_hadron_mueta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_mueta_" + tag[i]).c_str(),("#eta_{#mu} vs #Delta#phi(#tau_{#mu}, #tau_{3prong}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong})"); h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetYTitle("#eta_{#mu}");
    
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_deltaeta_" + tag[i]).c_str(),("#Delta(abs(#eta))_{#mu_#tau} vs #Delta#phi(#tau_{#mu}, #tau_{3prong}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong})"); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetYTitle("abs(#tau_{3prong} #eta) - abs(#mu #eta)");
    
    h_mueta_taueta[i] = new TH2F(("h_mueta_taueta_" + tag[i]).c_str(),("#tau_{3prong} #eta vs #tau_{#mu} #eta " + tag[i]).c_str(),18,-3,3, 18,-3,3);
    h_mueta_taueta[i]->SetXTitle("#tau_{3prong} #eta"); h_mueta_taueta[i]->SetYTitle("#tau_{#mu} #eta");
    
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
    
    h_tauEta_averageZDCside[i] = new TH2F(("h_tauEta_averageZDCside_" + tag[i]).c_str(),("average ZDC side vs #tau_{3prong} #eta " + tag[i]).c_str(),19, -2.5, 2.5, 21, -1, 1);
    h_tauEta_averageZDCside[i]->SetXTitle("#tau_{3prong} #eta"); h_tauEta_averageZDCside[i]->SetYTitle("average ZDC side");
    
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
    
    h_calo_leadingHFp_highNch[i] = new TH1F(("h_calo_leadingHFp_highNch_" + tag[i]).c_str(),("energy leading tower HF+ " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHFp_highNch[i]->SetXTitle("leading tower energy HF+ [GeV]"); h_calo_leadingHFp_highNch[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFp_highNch[i]->SetLineColor(colors[i]); h_calo_leadingHFp_highNch[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFm[i] = new TH1F(("h_calo_leadingHFm_" + tag[i]).c_str(),("energy leading tower HF- " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHFm[i]->SetXTitle("leading tower energy HF- [GeV]"); h_calo_leadingHFm[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFm[i]->SetLineColor(colors[i]); h_calo_leadingHFm[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFm_highNch[i] = new TH1F(("h_calo_leadingHFm_highNch_" + tag[i]).c_str(),("energy leading tower HF- " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHFm_highNch[i]->SetXTitle("leading tower energy HF- [GeV]"); h_calo_leadingHFm_highNch[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFm_highNch[i]->SetLineColor(colors[i]); h_calo_leadingHFm_highNch[i]->SetMarkerStyle(styles[i]);
    
    h_nHitsPixel[i] = new TH1F(("h_nHitsPixel_" + tag[i]).c_str(),("pixel hits in the event " + tag[i]).c_str(),30, 0, 60);
    h_nHitsPixel[i]->SetXTitle("pixel hits in the event"); h_nHitsPixel[i]->Sumw2();
    h_nHitsPixel[i]->SetLineColor(colors[i]); h_nHitsPixel[i]->SetMarkerStyle(styles[i]);
    
    /*h_nTrkHitsPixelEta[i] = new TProfile(("h_nTrkHitsPixelEta_" + tag[i]).c_str(),("average pixel hits per track vs eta " + tag[i]).c_str(),30, -3, 3, 0,25);
    h_nTrkHitsPixelEta[i]->SetXTitle("#eta"); h_nTrkHitsPixelEta[i]->SetYTitle("average pixel hits per track");
    h_nTrkHitsPixelEta[i]->Sumw2();
    h_nTrkHitsPixelEta[i]->SetLineColor(colors[i]); h_nTrkHitsPixelEta[i]->SetMarkerStyle(styles[i]);*/
    
    h_calo_Et_eta[i] = new TH2F(("h_calo_Et_eta_" + tag[i]).c_str(),("calo E_{T} vs #eta " + tag[i] +";#eta;calo E_{T} [GeV]").c_str(),50, -5, 5, 36, -1, 71);
    
    h_calo_Et[i] = new TH1F(("h_calo_Et_" + tag[i]).c_str(),("calo E_{T} " + tag[i]).c_str(),50, -1, 99);
    h_calo_Et[i]->SetXTitle("calo E_{T} [GeV]"); h_calo_Et[i]->Sumw2();
    /*if (i != 0) {*/h_calo_Et[i]->SetLineColor(colors[i]); h_calo_Et[i]->SetMarkerStyle(styles[i]);
    
    h_calo_Et_eta[i] = new TH2F(("h_calo_Et_eta_" + tag[i]).c_str(),("calo E_{T} vs #eta " + tag[i] +";#eta;calo E_{T} [GeV]").c_str(),50, -5, 5, 36, -1, 71);
    
    h_calo_leadingEt_eta[i] = new TH2F(("h_calo_leadingEt_eta_" + tag[i]).c_str(),("calo leading E_{T} vs #eta " + tag[i] +";#eta of leading calo tower;E_{T} of leading calo tower [GeV]").c_str(),50, -5, 5, 26, -1, 51);
    
    h_calo_sumEt[i] = new TH1F(("h_calo_sumEt_" + tag[i]).c_str(),("calo E_{T} sum " + tag[i]).c_str(),70, -1, 139);
    h_calo_sumEt[i]->SetXTitle("calo E_{T} sum [GeV]"); h_calo_sumEt[i]->Sumw2();
    /*if (i != 0) {*/h_calo_sumEt[i]->SetLineColor(colors[i]); h_calo_sumEt[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingEt[i] = new TH1F(("h_calo_leadingEt_" + tag[i]).c_str(),("leading calo E_{T} " + tag[i]).c_str(),15, -2, 58);
    h_calo_leadingEt[i]->SetXTitle("leading calo E_{T} [GeV]"); h_calo_leadingEt[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingEt[i]->SetLineColor(colors[i]); h_calo_leadingEt[i]->SetMarkerStyle(styles[i]);
    
    h_calo_E[i] = new TH1F(("h_calo_E_" + tag[i]).c_str(),("calo energy " + tag[i]).c_str(),50, -1, 99);
    h_calo_E[i]->SetXTitle("calo energy [GeV]"); h_calo_E[i]->Sumw2();
    /*if (i != 0) {*/h_calo_E[i]->SetLineColor(colors[i]); h_calo_E[i]->SetMarkerStyle(styles[i]);
    
    h_calo_E_eta[i] = new TH2F(("h_calo_E_eta_" + tag[i]).c_str(),("calo energy vs #eta " + tag[i] +";#eta;calo energy [GeV]").c_str(),50, -5, 5, 36, -1, 71);
    
    h_calo_leadingE_eta[i] = new TH2F(("h_calo_leadingE_eta_" + tag[i]).c_str(),("calo leading energy vs #eta " + tag[i] +";#eta of leading calo tower;energy of leading calo tower [GeV]").c_str(),50, -5, 5, 26, -1, 51);
    
    h_calo_sumE[i] = new TH1F(("h_calo_sumE_" + tag[i]).c_str(),("calo energy sum " + tag[i]).c_str(),100, -7, 1393);
    h_calo_sumE[i]->SetXTitle("calo energy sum [GeV]"); h_calo_sumE[i]->Sumw2();
    /*if (i != 0) {*/h_calo_sumE[i]->SetLineColor(colors[i]); h_calo_sumE[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingE[i] = new TH1F(("h_calo_leadingE_" + tag[i]).c_str(),("leading calo energy " + tag[i]).c_str(),15, -2, 58);
    h_calo_leadingE[i]->SetXTitle("leading calo energy [GeV]"); h_calo_leadingE[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingE[i]->SetLineColor(colors[i]); h_calo_leadingE[i]->SetMarkerStyle(styles[i]);
    
    h_nCaloTowers[i] = new TH1F(("h_nCaloTowers_" + tag[i]).c_str(),("number of calo towers " + tag[i]).c_str(),60, -5, 595);
    h_nCaloTowers[i]->SetXTitle("number of calo towers"); h_nCaloTowers[i]->Sumw2();
    /*if (i != 0) {*/h_nCaloTowers[i]->SetLineColor(colors[i]); h_nCaloTowers[i]->SetMarkerStyle(styles[i]);
    
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
    
    h_calo_tauEta_leadingHFeta[i] = new TH2F(("h_calo_tauEta_leadingHFeta_" + tag[i]).c_str(),("leading HF #eta vs #tau_{3prong} #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_tauEta_leadingHFeta[i]->SetXTitle("#tau_{3prong} #eta"); h_calo_tauEta_leadingHFeta[i]->SetYTitle("leading HF #eta");
    
    h_calo_muEta_averageHFeta[i] = new TH2F(("h_calo_muEta_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs muon #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_muEta_averageHFeta[i]->SetXTitle("muon #eta"); h_calo_muEta_averageHFeta[i]->SetYTitle("average HF #eta");
    
    h_calo_tauEta_averageHFeta[i] = new TH2F(("h_calo_tauEta_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs #tau_{3prong} #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_tauEta_averageHFeta[i]->SetXTitle("#tau_{3prong} #eta"); h_calo_tauEta_averageHFeta[i]->SetYTitle("average HF #eta");
    
    // MET
    
    h_MET[i] = new TH1F(("h_MET_" + tag[i]).c_str(),("MET " + tag[i]).c_str(),35, -0.5, 34.5);
    h_MET[i]->SetXTitle("MET"); h_MET[i]->Sumw2();
    /*if (i != 0) {*/h_MET[i]->SetLineColor(colors[i]); h_MET[i]->SetMarkerStyle(styles[i]);
    
  } //loop on nSamples
  
  for (int i = nSamples; i < nSamples+2+nNchCategories; i++){
    
    A_highNch_highHF[i] = new TH1F(("A_highNch_highHF_" + to_string(i)).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - high HF " + to_string(i)).c_str(), deltaphi_bins, 0, TMath::Pi());
    A_highNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - high HF"); A_highNch_highHF[i]->Sumw2();
    A_highNch_highHF[i]->SetLineColor(colors[i]); A_highNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    B_lowNch_highHF[i] = new TH1F((ABCDsysNames[i-nSamples]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) - " + ABCDsysNames[i-nSamples]).c_str(), deltaphi_bins, 0, TMath::Pi());
    B_lowNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - high HF"); B_lowNch_highHF[i]->Sumw2();
    B_lowNch_highHF[i]->SetLineColor(colors[i]); B_lowNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    C_highNch_lowHF[i] = new TH1F(("C_highNch_lowHF_" + to_string(i)).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - low HF " + to_string(i)).c_str(), deltaphi_bins, 0, TMath::Pi());
    C_highNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) high Nch - low HF"); C_highNch_lowHF[i]->Sumw2();
    C_highNch_lowHF[i]->SetLineColor(colors[i]); C_highNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
    D_lowNch_lowHF[i] = new TH1F(("D_lowNch_lowHF_" + to_string(i)).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - low HF " + to_string(i)).c_str(), deltaphi_bins, 0, TMath::Pi());
    D_lowNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong}) low Nch - low HF"); D_lowNch_lowHF[i]->Sumw2();
    D_lowNch_lowHF[i]->SetLineColor(colors[i]); D_lowNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
  }
  
  for (int h = 0; h < histograms.size(); h++){
    string histogramsName = histograms.at(h)[0]->GetName();
    histograms.at(h)[nSamples] = (TH1F*)histograms.at(h)[0]->Clone(("A - "+histogramsName).c_str());
    histograms.at(h)[nSamples+1] = (TH1F*)histograms.at(h)[0]->Clone(("B - "+histogramsName).c_str());
    histograms.at(h)[nSamples+2] = (TH1F*)histograms.at(h)[0]->Clone(("C - "+histogramsName).c_str());
  }

  TH2F *h_calo_energyHFp_nch = new TH2F("h_calo_energyHFp_nch","leading tower HF+ vs nch",9 , 2.5, 11.5, 12, 1, 5);
  h_calo_energyHFp_nch->SetXTitle("nch"); h_calo_energyHFp_nch->SetYTitle("leading tower HF+ [GeV]");
  
  TH2F *h_calo_energyHFm_nch = new TH2F("h_calo_energyHFm_nch","leading tower HF- vs nch",9 , 2.5, 11.5, 12, 1, 5);
  h_calo_energyHFm_nch->SetXTitle("nch"); h_calo_energyHFm_nch->SetYTitle("leading tower HF- [GeV]");

    
  TH1F *h_SB_deltaphi = new TH1F("h_SB_deltaphi","S/sqrt(S+B) for #Delta#phi(#tau_{#mu}, #tau_{3prong})",deltaphi_bins, 0, TMath::Pi());
  h_SB_deltaphi->SetXTitle("S/sqrt(S+B) for #Delta#phi(#tau_{#mu}, #tau_{3prong})"); h_SB_deltaphi->Sumw2();
  h_SB_deltaphi->SetLineColor(colors[0]); h_SB_deltaphi->SetMarkerStyle(styles[1]);
  h_SB_deltaphi->SetBinErrorOption(TH1::kPoisson);
    
  // other
  TH1F *h_track_activity_pt     = new TH1F("h_track_activity_pt",     "h_track_activity_pt",     25, 0, 10); h_track_activity_pt->Sumw2(); h_track_activity_pt->SetXTitle("track p_{T} [GeV]");
  TH2F *h_track_activity_pt_eta = new TH2F("h_track_activity_pt_eta", "h_track_activity_pt_eta", 25, 0, 10, 10, -2.5, 2.5);
  h_track_activity_pt_eta->SetXTitle("track p_{T} [GeV]"); h_track_activity_pt_eta->SetYTitle("track #eta");
  
  TH2F *h_deltaphi_tau_mu_tau_hadron_nch = new TH2F("h_deltaphi_tau_mu_tau_hadron_nch", "nch vs #Delta#phi(#tau_{#mu}, #tau_{3prong})", 8, 0, TMath::Pi(), 5 , 0.5, 5.5);
  h_deltaphi_tau_mu_tau_hadron_nch->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{3prong})"); h_deltaphi_tau_mu_tau_hadron_nch->SetYTitle("nch");


// end of histogram declaration
// *******************************************

//  double mass1, mass2, mass3;
//  TFile *ntuple = new TFile("ntuple.root", "RECREATE");
//  TTree *aux;
//  aux = new TTree("tree", "tree");
//  aux->Branch("mass1", &mass1);
//  aux->Branch("mass2", &mass2);
//  aux->Branch("mass3", &mass3);

  double muon_mass = 105.6583 / 1000.; //GeV
  double pionMass = 139.570 / 1000.; //GeV
  int entries = (TREE->fChain)->GetEntries();
  
  TLorentzVector tau_muon, tau_hadron, full_tau_hadron, gamma, tempGamma, pion[3];
  int charge_counter[nSamples][3];
  for (int s = 0; s < nSamples; s++){
    for (int i = 0; i < 3; i++) charge_counter[s][i] = 0;
  }
  
  TFile *eventFile = new TFile("events.root","recreate");
  TTree *eventID = new TTree ("eventID","eventID");
  int event;
  int run;
  int lumiBlock;
  float deltaPhi;
  eventID->Branch("event",&event,"event/I");
  eventID->Branch("run",&run,"run/I");
  eventID->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  eventID->Branch("deltaPhi",&deltaPhi,"deltaPhi/F");
  
  int pionLeadingIndex = -1;
  int pionSubLeadingIndex = -1;
  int pionSubSubLeadingIndex = -1;
  float pionLeadingPt = -1;
  float pionSubLeadingPt = -1;
  float pionSubSubLeadingPt = -1;
  //float muonSF = tnp_weight_trk_pbpb(0);
  //float muonSFerror = TMath::Max(TMath::Abs(tnp_weight_trk_pbpb(-1)-tnp_weight_trk_pbpb(0)),TMath::Abs(tnp_weight_trk_pbpb(0)-tnp_weight_trk_pbpb(-2)));
  //cout << "muon inner tracker SF: " << muonSF << endl;
  //cout << "uncertainty on muon inner tracker SF: " << muonSFerror << endl;
  int muon_charge = 0;
  int tauh_charge = 0;
  float muon_weight = 1;
  float muon_weight_error = 0;
  bool passedNch = false;
  int category = -1;
  int candCounter = 0;
  float tau_hadron_ptScalar = 0;
  float ditau_ptScalar = 0;
  float minPionPionDeltaPhi = 1000;
  float maxPionPionDeltaPhi = 0;
  float minPionPionDeltaEta = 1000;
  float maxPionPionDeltaEta = 0;
  float minPionPionDeltaR = 1000;
  float maxPionPionDeltaR = 0;
  float minPionMuDeltaPhi = 1000;
  float maxPionMuDeltaPhi = 0;
  float minPionMuDeltaEta = 1000;
  float maxPionMuDeltaEta = 0;
  float minPionMuDeltaR = 1000;
  float maxPionMuDeltaR = 0;
  float delta_phi = 0;
  float full_delta_phi = 0;
  float maxSB = 0;
  int maxSBbin = 1;
  int nBatches = 20;
  int printing_every = entries/nBatches;
  cout << "Running on Data ..." << endl;
  for(int iEntry=0; iEntry<entries; iEntry++){
    (TREE->fChain)->GetEntry(iEntry);
    
    if (!(iEntry%printing_every)) cout << "\r" << int((100.0/nBatches)*(iEntry/printing_every)) << "%" << flush;
    //if (TREE->BsTauTau_mu1_pt->at(0) < 2.5) continue; // fix me
    //if (TREE->BsTauTau_nch->at(0) > 8) continue; // fix me
    
    TLorentzVector tau_hadron_pi1, tau_hadron_pi2, tau_hadron_pi3;
    
    bool triggered = TREE->triggered->at(0);
    if (!triggered) continue; // fix me
    
    pionLeadingIndex = -1;
    pionSubLeadingIndex = -1;
    pionSubSubLeadingIndex = -1;
    pionLeadingPt = -1;
    pionSubLeadingPt = -1;
    pionSubSubLeadingPt = -1;
    bool passedmu  = false;
    bool passedtau = false;
    bool passedcalo = true;
    bool passedcaloDown = true;
    bool passedcaloUp = true;
    bool passedZDC = true;
    bool passedMET = true;
    bool passedGamma = true;
    passedNch = false;
    category = -1;
    tau_hadron_ptScalar = 0;
    ditau_ptScalar = 0;
    minPionPionDeltaPhi = 1000;
    maxPionPionDeltaPhi = 0;
    minPionPionDeltaEta = 1000;
    maxPionPionDeltaEta = 0;
    minPionPionDeltaR = 1000;
    maxPionPionDeltaR = 0;
    minPionMuDeltaPhi = 1000;
    maxPionMuDeltaPhi = 0;
    minPionMuDeltaEta = 1000;
    maxPionMuDeltaEta = 0;
    minPionMuDeltaR = 1000;
    maxPionMuDeltaR = 0;
    double temp_tau_pt_comparison = 0.;
    muon_charge = 0;
    int temp_nch = 0.;
    double temp_npv = 0.;
    double temp_pvz = 0.;
    double temp_pv_trk_1 = 0.;
    double temp_pv_trk_2 = 0.;
    double temp_pv_trk_3 = 0.;
    for (int i=0; i<(int)TREE->BsTauTau_mu1_pt->size(); i++) {
      if (tau_muon_isSoft){
        if (!TREE->BsTauTau_mu1_isSoft->at(i)) continue;
      }
      if (tau_muon_isGlobal){
        if (!TREE->BsTauTau_mu1_isGlobal->at(i)) continue;
      }
      if (tau_muon_isTracker){
        if (!TREE->BsTauTau_mu1_isTracker->at(i)) continue;
      }
if( (TREE->BsTauTau_mu1_pt->at(i) < 3.5 && TMath::Abs(TREE->BsTauTau_mu1_eta->at(i)) < 1.2) || (TREE->BsTauTau_mu1_pt->at(i) < 2.5 && TMath::Abs(TREE->BsTauTau_mu1_eta->at(i)) > 1.2))  continue;
      //if (TREE->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }
      if (TREE->BsTauTau_mu1_pt->at(i) > 2.5) passedmu = true;
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
    
    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    candCounter = 0;
    int chosenTauIndex = -1;
    float vprob = -1;
    float temp_pt = 0;
    double temp_tau_rho1 = 0.; double temp_tau_rho2 = 0.;
    for (int i=0; i<(int)TREE->BsTauTau_tau_pt->size(); i++) {
      if (TREE->cand_leadingPionPt->at(i) > pionLeadingCut && TREE->cand_subleadingPionPt->at(i) > pionSubLeadingCut && TREE->cand_subsubleadingPionPt->at(i) > pionSubSubLeadingCut){
        candCounter++;
        chosenTauIndex = i;
        tau_hadron_ptScalar = TREE->cand_leadingPionPt->at(i) + TREE->cand_subleadingPionPt->at(i) + TREE->cand_subsubleadingPionPt->at(i);
      }
      else continue;
      if (muon_charge*TREE->BsTauTau_tau_q->at(i) != MuTauCharge) continue;
      if (TREE->BsTauTau_tau_pt->at(i) < 2) continue;
      //if (temp_nch != 1 && TREE->BsTauTau_tau_pt->at(i) < 3.0) continue; //temporary commented
      //if(TMath::Abs(TREE->BsTauTau_tau_eta->at(i)) > 2) continue;
      if (TREE->BsTauTau_tau_pt->at(i) > temp_pt) {
        vprob = 100*TREE->BsTauTau_tau_vprob->at(i);
        if (temp_nch == 1) vprob = 100*TREE->BsTauTau_B_vprob->at(i);
        temp_pt = TREE->BsTauTau_tau_pt->at(i);
      }
      if (TREE->BsTauTau_tau_vprob->at(i) < tau_hadron_vertexprob) continue;
      passedtau = true;
      pion[0].SetPtEtaPhiM (TREE->cand_leadingPionPt->at(i),TREE->cand_leadingPionEta->at(i),TREE->cand_leadingPionPhi->at(i),pionMass);
      pion[1].SetPtEtaPhiM (TREE->cand_subleadingPionPt->at(i),TREE->cand_subleadingPionEta->at(i),TREE->cand_subleadingPionPhi->at(i),pionMass);
      pion[2].SetPtEtaPhiM (TREE->cand_subsubleadingPionPt->at(i),TREE->cand_subsubleadingPionEta->at(i),TREE->cand_subsubleadingPionPhi->at(i),pionMass);
      //if (temp_nch == 1 && TREE->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      
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
        //temp_pv_trk_1 = TREE->BsTauTau_tau_pi1_z->at(i); temp_pv_trk_2 = TREE->BsTauTau_tau_pi2_z->at(i); temp_pv_trk_3 = TREE->BsTauTau_tau_pi3_z->at(i);
        }
      }
    } // loop over the size of the tau candidates
    if (candCounter != 1) passedNch = false;
    
    ditau_ptScalar = tau_hadron_ptScalar + tau_muon.Pt();
    bool passedDitau = (tau_muon+tau_hadron).M() > ditauMassCut;
    
    for (int p1 = 0; p1 < 3; p1++){
      float tempDeltaPhiPionMu = TMath::Abs(tau_muon.DeltaPhi(pion[p1]));
      float tempDeltaEtaPionMu = TMath::Abs(tau_muon.Eta()-pion[p1].Eta());
      float tempDeltaRPionMu = TMath::Abs(tau_muon.DeltaR(pion[p1]));
      if (minPionMuDeltaPhi > tempDeltaPhiPionMu) minPionMuDeltaPhi = tempDeltaPhiPionMu;
      if (maxPionMuDeltaPhi < tempDeltaPhiPionMu) maxPionMuDeltaPhi = tempDeltaPhiPionMu;
      if (minPionMuDeltaEta > tempDeltaEtaPionMu) minPionMuDeltaEta = tempDeltaEtaPionMu;
      if (maxPionMuDeltaEta < tempDeltaEtaPionMu) maxPionMuDeltaEta = tempDeltaEtaPionMu;
      if (minPionMuDeltaR > tempDeltaRPionMu) minPionMuDeltaR = tempDeltaRPionMu;
      if (maxPionMuDeltaR < tempDeltaRPionMu) maxPionMuDeltaR = tempDeltaRPionMu;
      for (int p2 = p1+1; p2 < 3; p2++){
        float tempDeltaPhiPionPion = TMath::Abs(pion[p2].DeltaPhi(pion[p1]));
        float tempDeltaEtaPionPion = TMath::Abs(pion[p2].Eta()-pion[p1].Eta());
        float tempDeltaRPionPion = TMath::Abs(pion[p2].DeltaR(pion[p1]));
        if (minPionPionDeltaPhi > tempDeltaPhiPionPion) minPionPionDeltaPhi = tempDeltaPhiPionPion;
        if (maxPionPionDeltaPhi < tempDeltaPhiPionPion) maxPionPionDeltaPhi = tempDeltaPhiPionPion;
        if (minPionPionDeltaEta > tempDeltaEtaPionPion) minPionPionDeltaEta = tempDeltaEtaPionPion;
        if (maxPionPionDeltaEta < tempDeltaEtaPionPion) maxPionPionDeltaEta = tempDeltaEtaPionPion;
        if (minPionPionDeltaR > tempDeltaRPionPion) minPionPionDeltaR = tempDeltaRPionPion;
        if (maxPionPionDeltaR < tempDeltaRPionPion) maxPionPionDeltaR = tempDeltaRPionPion;
      }
    }
    
    //if ((tau_muon+tau_hadron).Pt() > MET_cut) passedMET = false;

    bool found_mu_track = false;
    bool found_pi1_track = false;
    bool found_pi2_track = false;
    bool found_pi3_track = false;
    bool high_activity = false;

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
    double sumEt = 0;
    double sumE = 0;
    double leadingE = 0;
    double leadingE_Eta = -100;
    double leadingEt = 0;
    double leadingEt_Eta = -100;

    if (triggered && passedmu && passedtau)
    {
      for (int i=0; i<(int)TREE->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREE->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREE->BsTauTau_calo_energyHFm->at(i);
        /*if (passedNch) h_calo_E[0]->Fill(TREE->BsTauTau_calo_energy->at(i));
        if (passedNch) h_calo_E_eta[0]->Fill(TREE->BsTauTau_calo_eta->at(i),TREE->BsTauTau_calo_energy->at(i));
        if (passedNch) h_calo_Et[0]->Fill(TREE->BsTauTau_calo_eT->at(i));
        if (passedNch) h_calo_Et_eta[0]->Fill(TREE->BsTauTau_calo_eta->at(i),TREE->BsTauTau_calo_eT->at(i));
        sumE += TREE->BsTauTau_calo_energy->at(i);
        sumEt += TREE->BsTauTau_calo_eT->at(i);*/
        double etaCal = TREE->BsTauTau_calo_eta->at(i);
        /*if (TREE->BsTauTau_calo_energy->at(i) > leadingE){
          leadingE = TREE->BsTauTau_calo_energy->at(i);
          leadingE_Eta = etaCal;
        }
        if (TREE->BsTauTau_calo_eT->at(i) > leadingEt){
          leadingEt = TREE->BsTauTau_calo_eT->at(i);
          leadingEt_Eta = etaCal;
        }*/
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
      if (maxHFp > 0.9*HFpLeading_high || maxHFm > 0.9*HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcaloDown = false;
      if (maxHFp > 1.1*HFpLeading_high || maxHFm > 1.1*HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcaloUp = false;
      if (passedcalo && passedNch) h_calo_sumE[0]->Fill(sumE);
      if (passedcalo && passedNch) h_calo_leadingE[0]->Fill(leadingE);
      if (passedNch) h_calo_leadingE_eta[0]->Fill(leadingE_Eta,leadingE);
      if (passedcalo && passedNch) h_calo_sumEt[0]->Fill(sumEt);
      if (passedcalo && passedNch) h_calo_leadingEt[0]->Fill(leadingEt);
      if (passedNch) h_calo_leadingEt_eta[0]->Fill(leadingEt_Eta,leadingEt);
      //if (passedcalo && passedNch) h_nCaloTowers[0]->Fill(TREE->BsTauTau_calo_eta->size());
      if (passedcalo && passedNch) h_calo_energyHFp_sum[0]->Fill(sumHFp);
      if (passedcalo && passedNch) h_calo_energyHFm_sum[0]->Fill(sumHFm);
      if (passedcalo && passedNch) h_calo_energyHF_pm[0]->Fill(sumHFm,sumHFp);
      if (passedcalo && passedNch) h_calo_energyHFp_size[0]->Fill(sizeHFp);
      if (passedcalo && passedNch) h_calo_energyHFm_size[0]->Fill(sizeHFm);
      if (passedNch && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFp[0]->Fill(maxHFp);
      if (passedNch && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low) h_calo_leadingHFm[0]->Fill(maxHFm);
      if (temp_nch>=5 && temp_nch<=8 && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFp_highNch[0]->Fill(maxHFp);
      if (temp_nch>=5 && temp_nch<=8 && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low) h_calo_leadingHFm_highNch[0]->Fill(maxHFm);
      if (passedNch) h_calo_leadingHF_pm[0]->Fill(maxHFm,maxHFp);
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    }
    

    if (passedmu && passedNch && passedcalo && triggered)
    {
    TLorentzVector recoPiZero;
    int Ngamma = 0;
    float minDeltaR = 1000;
    float minDeltaM = 100000; //MeV
    int NrecoPiZero = 0;
    if (!threeProng[0]){
    Ngamma = TREE->BsTauTau_nGammas->at(0);
    full_tau_hadron = tau_hadron;
    for (int i=0; i<Ngamma; i++){
      if (TREE->reco_gamma_pt->at(i)<gammaPtCut || TMath::Abs(TREE->reco_gamma_eta->at(i)) > 2.2){
        Ngamma -= 1;
        continue;
      }
      gamma.SetPtEtaPhiM (TREE->reco_gamma_pt->at(i),TREE->reco_gamma_eta->at(i),TREE->reco_gamma_phi->at(i),0);
      //if (i == 0) recoNeutralPion = gamma;
      //else recoNeutralPion += gamma;
      h_recoGamma_pt[0]->Fill(TREE->reco_gamma_pt->at(i));
      h_recoGamma_eta[0]->Fill(TREE->reco_gamma_eta->at(i));
      //if (TMath::Abs(TREE->reco_gamma_eta->at(i)) > 2.2) passedGamma = false;
      h_recoGamma_phi[0]->Fill(TREE->reco_gamma_phi->at(i));
      h_recoGamma_deltaphi_muon[0]->Fill(TMath::Abs(gamma.DeltaPhi(tau_muon)));
      h_recoGamma_deltaphi_pion[0]->Fill(TMath::Abs(gamma.DeltaPhi(tau_hadron)));
      for (int j=i+1; j<Ngamma; j++){
        if (TREE->reco_gamma_pt->at(j)<gammaPtCut || TMath::Abs(TREE->reco_gamma_eta->at(j)) > 2.2) continue;
        tempGamma.SetPtEtaPhiM (TREE->reco_gamma_pt->at(j),TREE->reco_gamma_eta->at(j),TREE->reco_gamma_phi->at(j),0);
        float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
        if (deltaR < minDeltaR) minDeltaR = deltaR;
        recoPiZero = gamma+tempGamma;
        double recoPiZeroMass = 1000*recoPiZero.M(); // MeV
        if (TMath::Abs(recoPiZeroMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = recoPiZeroMass-PiZeroMass;
        h_recoPiZeroDeltaM[0]->Fill(recoPiZeroMass-PiZeroMass);
        if (TMath::Abs(recoPiZeroMass-PiZeroMass) < 50){
          h_recoPiZero_pt[0]->Fill(recoPiZero.Pt());
          h_recoPiZero_eta[0]->Fill(recoPiZero.Eta());
          h_recoPiZero_phi[0]->Fill(recoPiZero.Phi());
          h_recoPiZero_deltaphi_muon[0]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_muon)));
          h_recoPiZero_deltaphi_pion[0]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)));
          NrecoPiZero += 1;
          if (gamma.Pt()>gammaPtCut && tempGamma.Pt()>gammaPtCut) full_tau_hadron += recoPiZero;
        }
      }
    }
    h_NrecoPiZero[0]->Fill(NrecoPiZero);
    h_recoGammasMinDeltaR[0]->Fill(minDeltaR);
    h_recoPiZeroMinDeltaM[0]->Fill(minDeltaM);
    h_NrecoGamma[0]->Fill(Ngamma);
    } // if not three prong
    } // if in the signal region
    

    delta_phi = TMath::Abs(tau_muon.DeltaPhi(tau_hadron));
    full_delta_phi = TMath::Abs(tau_muon.DeltaPhi(full_tau_hadron));
    
    
    h_cutflow[0]->Fill(0); //input from Ntuplizer
    if (passedmu) h_cutflow[0]->Fill(1); //tau mu
    if (passedNch) h_cutflow[0]->Fill(2); //nCh
    if (passedtau) h_cutflow[0]->Fill(3); //tau hadron
    if (passedcalo) h_cutflow[0]->Fill(4); //HF
    if (passedcalo && passedtau) h_cutflow[0]->Fill(5); //HF and tau hadron
    if (passedcalo && passedtau && passedNch) h_cutflow[0]->Fill(6); //HF and tau hadron and nch
    //if () h_cutflow[0]->Fill(4); //vertex prob
    
    if (temp_nch >= 5 && temp_nch <= 8 && !passedcalo) category = nSamples;
    if (passedNch && !passedcalo) category = nSamples+1;
    if (temp_nch >= 5 && temp_nch <= 8 && passedcalo) category = nSamples+2;
    if (passedNch && passedcalo) category = 0;
    
    bool keepEvent = false;
    
    if (passedmu && passedMET && passedGamma && triggered)
    {
      keepEvent = true;
    }
    if (keepEvent && category != -1){
      h_tau_hadron_vprob[category]->Fill(vprob);
    }
    if (!passedtau) keepEvent = false;
    if (keepEvent && category == 0){
      charge_counter[0][muon_charge*tauh_charge + 1] += 1;
    }
    if (muon_charge*tauh_charge == 0) keepEvent = false;
    if (keepEvent && category != -1){
      h_sumZDCplus[category]->Fill(sumZDCp);
      h_sumZDCminus[category]->Fill(sumZDCm);
      h_sumZDC_pm[0]->Fill(sumZDCm,sumZDCp);
      float averageZDC = (sumZDCp-ratioZDCpm*sumZDCm)/(sumZDCp+ratioZDCpm*sumZDCm);
      h_muEta_averageZDCside[0]->Fill(tau_muon.Eta(),averageZDC);
      h_tauEta_averageZDCside[0]->Fill(tau_hadron.Eta(),averageZDC);
      h_averageZDCside_averageHFeta[0]->Fill(averageZDC,averageHFeta);
    }
    if(!passedZDC) keepEvent = false;
    if (keepEvent && category != -1){
      h_deltaphi_tau_mu_tau_hadron[category]->Fill(delta_phi);
      //h_deltaphi_tau_mu_tau_hadron[category]->Fill(delta_phi);
      h_deltaphi_tau_mu_tau_hadron_zoomed[category]->Fill(delta_phi);
      if (!category) h_deltaphi_tau_mu_full_tau_hadron[0]->Fill(full_delta_phi);
      h_deltaphi_tau_mu_tau_hadron_nch->Fill(delta_phi,temp_nch);
      if (!category) h_deltaphi_tau_mu_tau_hadron_mueta[0]->Fill(delta_phi,tau_muon.Eta());
      if (!category) h_deltaphi_tau_mu_tau_hadron_deltaeta[0]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()));
    }
    if(delta_phi < deltaPhi_cut) keepEvent = false;
    if (keepEvent){
      if (passedcalo) h_tau_hadron_nch[0]->Fill(temp_nch);
      if (!passedcalo) h_tau_hadron_nch_highHF[0]->Fill(temp_nch);
      if (passedcalo) h_tau_hadron_ncand_final[0]->Fill(candCounter);
      h_calo_energyHFp_nch->Fill(temp_nch,maxHFp); h_calo_energyHFm_nch->Fill(temp_nch,maxHFm);
    }
    if (category == -1) keepEvent = false;
    if (keepEvent){
      h_minPionPionDeltaPhi[category]->Fill(minPionPionDeltaPhi);
      h_maxPionPionDeltaPhi[category]->Fill(maxPionPionDeltaPhi);
      h_minPionPionDeltaEta[category]->Fill(minPionPionDeltaEta);
      h_maxPionPionDeltaEta[category]->Fill(maxPionPionDeltaEta);
      h_minPionPionDeltaR[category]->Fill(minPionPionDeltaR);
      h_maxPionPionDeltaR[category]->Fill(maxPionPionDeltaR);
      h_minPionMuDeltaPhi[category]->Fill(minPionMuDeltaPhi);
      h_maxPionMuDeltaPhi[category]->Fill(maxPionMuDeltaPhi);
      h_minPionMuDeltaEta[category]->Fill(minPionMuDeltaEta);
      h_maxPionMuDeltaEta[category]->Fill(maxPionMuDeltaEta);
      h_minPionMuDeltaR[category]->Fill(minPionMuDeltaR);
      h_maxPionMuDeltaR[category]->Fill(maxPionMuDeltaR);
      for (int i = 0; i < (int)TREE->reco_pion_ecalEnergy->size(); i++) h_reco_pion_energy_HCAL_ECAL[0]->Fill(TREE->reco_pion_ecalEnergy->at(i),TREE->reco_pion_hcalEnergy->at(i));
      h_tau_mu_p[category]->Fill(tau_muon.P());
      h_tau_mu_pz[category]->Fill(tau_muon.Pz());
      h_tau_mu_dz[category]->Fill(10*(TREE->BsTauTau_mu1_vz->at(0)-TREE->BsTauTau_PV_vz->at(0)));
      h_tau_mu_pt[category]->Fill(tau_muon.Pt());
      h_tau_mu_eta[category]->Fill(tau_muon.Eta());
      h_tau_mu_phi[category]->Fill(tau_muon.Phi());
      h_tau_hadron_p[category]->Fill(tau_hadron.P());
      h_tau_hadron_pz[category]->Fill(tau_hadron.Pz());
      h_tau_hadron_ptScalar[category]->Fill(tau_hadron_ptScalar);
      h_tau_hadron_pt[category]->Fill(tau_hadron.Pt());
      h_tau_hadron_eta[category]->Fill(tau_hadron.Eta());
      h_tau_hadron_phi[category]->Fill(tau_hadron.Phi());
      h_tau_hadron_mass[category]->Fill(tau_hadron.M());
      if (!category){
        h_tau_mu_tau_hadron_phi[0]->Fill(tau_muon.Phi(),tau_hadron.Phi());
        h_tau_hadron_rhomass[0][0]->Fill(temp_tau_rho1);
        h_tau_hadron_rhomass[0][1]->Fill(temp_tau_rho2);
        h_tau_hadron_rhomass2D[0]->Fill(temp_tau_rho1,temp_tau_rho2);
        h_calo_muEta_leadingHFeta[0]->Fill(tau_muon.Eta(),leadingHFeta);
        h_calo_tauEta_leadingHFeta[0]->Fill(tau_hadron.Eta(),leadingHFeta);
        h_calo_muEta_averageHFeta[0]->Fill(tau_muon.Eta(),averageHFeta);
        h_calo_tauEta_averageHFeta[0]->Fill(tau_hadron.Eta(),averageHFeta);
        h_mueta_taueta[0]->Fill(tau_muon.Eta(),tau_hadron.Eta());
        h_tau_hadron_track_pvz[0][0]->Fill(temp_pvz - temp_pv_trk_1);
        h_tau_hadron_track_pvz[0][1]->Fill(temp_pvz - temp_pv_trk_2);
        h_tau_hadron_track_pvz[0][2]->Fill(temp_pvz - temp_pv_trk_3);
        h_PV_N[0]->Fill(temp_npv);
        h_AP[0]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2);
      }
      TLorentzVector MET;
      MET.SetPtEtaPhiM((tau_muon+tau_hadron).Pt(),(tau_muon+tau_hadron).Eta(),-(tau_muon+tau_hadron).Phi(),0);
      h_MET[category]->Fill(MET.Pt());
      //h_ditau_mass[0]->Fill((tau_muon+tau_hadron+MET).M());
      h_ditau_mass[category]->Fill((tau_muon+tau_hadron).M());
      h_ditau_p[category]->Fill((tau_muon+tau_hadron).P());
      h_ditau_pz[category]->Fill((tau_muon+tau_hadron).Pz());
      h_ditau_pt[category]->Fill((tau_muon+tau_hadron).Pt());
      h_ditau_ptScalar[category]->Fill(ditau_ptScalar);
      h_pion_leading_pt[category]->Fill(TREE->cand_leadingPionPt->at(chosenTauIndex));
      h_pion_subleading_pt[category]->Fill(TREE->cand_subleadingPionPt->at(chosenTauIndex));
      h_pion_subsubleading_pt[category]->Fill(TREE->cand_subsubleadingPionPt->at(chosenTauIndex));
      h_pion_leading_eta[category]->Fill(TREE->cand_leadingPionEta->at(chosenTauIndex));
      h_pion_subleading_eta[category]->Fill(TREE->cand_subleadingPionEta->at(chosenTauIndex));
      h_pion_subsubleading_eta[category]->Fill(TREE->cand_subsubleadingPionEta->at(chosenTauIndex));
      h_pion_leading_phi[category]->Fill(TREE->cand_leadingPionPhi->at(chosenTauIndex));
      h_pion_subleading_phi[category]->Fill(TREE->cand_subleadingPionPhi->at(chosenTauIndex));
      h_pion_subsubleading_phi[category]->Fill(TREE->cand_subsubleadingPionPhi->at(chosenTauIndex));
      
      //h_nHitsPixel[category]->Fill(TREE->nevtPixelHits->at(0));
      //for () h_nTrkHitsPixelEta[category]->Fill(##ETA##,TREE->ntrkPixelHits->at(i));
      //cout << TREE->BsTauTau_calo_zdcSumPlus->at(0) << "," << TREE->BsTauTau_calo_zdcSumMinus->at(0) << endl;
      
      //cout << (tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()) << " ---- " << tau_muon.Pt() << endl;
      //h_MET[0]->Fill(TREE->MET_sumEt->at(0));
      
      
      if (!category && delta_phi > 2.78816){
        event = TREE->EVENT_event;
        run = TREE->EVENT_run;
        lumiBlock = TREE->EVENT_lumiBlock;
        deltaPhi = delta_phi;
        //if (deltaPhi>3.07 && deltaPhi<3.08) cout << "event: " << event << " run: " << run << " lumiBlock: " << lumiBlock << endl;
        /*if (TREE->cand_subsubleadingPionPt->at(chosenTauIndex) > 0.7 && (tau_muon+tau_hadron).M() > 20 && category == 0){
          cout << "\n\nevent: " << event << " run: " << run << " lumiBlock: " << lumiBlock << endl;
          cout << "delta phi: " << deltaPhi << " subsubleading pion pT: " << TREE->cand_subsubleadingPionPt->at(chosenTauIndex) << " muon pT: " << tau_muon.Pt() << endl;
          cout << "ditau mass: " << (tau_muon+tau_hadron).M() << " ditau pT: " << (tau_muon+tau_hadron).Pt() << endl;
        }*/
        eventID->Fill();
      }
    } //if (keepEvent)
    
    if (passedmu && triggered && passedtau && muon_charge*tauh_charge == -1){// && delta_phi >= deltaPhi_cut){
      if (temp_nch >= firstNchCategory && temp_nch < firstNchCategory+nNchCategories && !passedcalo) A_highNch_highHF[0]->Fill(delta_phi);
      if (temp_nch == 3 && passedNch && !passedcalo) B_lowNch_highHF[0]->Fill(delta_phi);
      if (temp_nch >= firstNchCategory && temp_nch < firstNchCategory+nNchCategories && passedcalo) C_highNch_lowHF[0]->Fill(delta_phi);
      if (temp_nch == 3 && passedNch && passedcalo) D_lowNch_lowHF[0]->Fill(delta_phi);
      
      if (temp_nch >= 5 && temp_nch < firstNchCategory+nNchCategories && !passedcaloDown) A_highNch_highHF[nSamples]->Fill(delta_phi);
      if (temp_nch == 3 && passedNch && !passedcaloDown) B_lowNch_highHF[nSamples]->Fill(delta_phi);
      if (temp_nch >= firstNchCategory && temp_nch < firstNchCategory+nNchCategories && passedcaloDown) C_highNch_lowHF[nSamples]->Fill(delta_phi);
      if (temp_nch == 3 && passedNch && passedcaloDown) D_lowNch_lowHF[nSamples]->Fill(delta_phi);
      
      if (temp_nch >= firstNchCategory && temp_nch < firstNchCategory+nNchCategories && !passedcaloUp) A_highNch_highHF[nSamples+1]->Fill(delta_phi);
      if (temp_nch == 3 && passedNch && !passedcaloUp) B_lowNch_highHF[nSamples+1]->Fill(delta_phi);
      if (temp_nch >= firstNchCategory && temp_nch < firstNchCategory+nNchCategories && passedcaloUp) C_highNch_lowHF[nSamples+1]->Fill(delta_phi);
      if (temp_nch == 3 && passedNch && passedcaloUp) D_lowNch_lowHF[nSamples+1]->Fill(delta_phi);
      
      for (int cat = 0; cat < nNchCategories; cat++){      
        if (temp_nch == firstNchCategory+cat && !passedcalo) A_highNch_highHF[nSamples+2+cat]->Fill(delta_phi);
        if (temp_nch == 3 && passedNch && !passedcalo) B_lowNch_highHF[nSamples+2+cat]->Fill(delta_phi);
        if (temp_nch == firstNchCategory+cat && passedcalo) C_highNch_lowHF[nSamples+2+cat]->Fill(delta_phi);
        if (temp_nch == 3 && passedNch && passedcalo) D_lowNch_lowHF[nSamples+2+cat]->Fill(delta_phi);
      }
      
      
      // ABCD validation
      if (category == nSamples){
        bool ultraHighHF = (maxHFp > 20 || maxHFm > 20);
        if (temp_nch > 6 && ultraHighHF) ABCD_validation[0]->Fill(delta_phi);
        if (temp_nch <= 6 && ultraHighHF) ABCD_validation[1]->Fill(delta_phi);
        if (temp_nch > 6 && !ultraHighHF) ABCD_validation[2]->Fill(delta_phi);
        if (temp_nch > 6 && !ultraHighHF) ABCD_validation[4]->Fill(delta_phi);
        if (temp_nch <= 6 && !ultraHighHF) ABCD_validation[3]->Fill(delta_phi);
        if (temp_nch <= 6 && !ultraHighHF) ABCD_validation[5]->Fill(delta_phi);
      }
    }

//      aux->Fill();

  } // loop over the entries
  
  
  eventFile->Write();
  eventFile->Close();
  
  
  TCanvas *c_ABCD_subcategories = new TCanvas("c_ABCD_subcategories", "c_ABCD_subcategories", 800, 800);
  ABCD_validation[0]->SetMinimum(0.9);
  for (int cat = 0; cat < 4; cat++) ABCD_validation[cat]->Draw("he0e1x0same");
  CMS_lumi( c_ABCD_subcategories, iPeriod, iPos );
  c_ABCD_subcategories->SetLogy(1);
  
  TLegend *legend_ABCD_subcategories = new TLegend(0.4,0.75,0.55,0.9);
  legend_ABCD_subcategories->SetFillStyle(0); legend_ABCD_subcategories->SetFillColor(0);
  legend_ABCD_subcategories->SetLineColor(0); legend_ABCD_subcategories->SetShadowColor(0); legend_ABCD_subcategories->SetTextSize(0.035);
  legend_ABCD_subcategories->AddEntry(ABCD_validation[0],    "N_{ch} = 7,8 & HF > 20", "l");
  legend_ABCD_subcategories->AddEntry(ABCD_validation[1],    "N_{ch} = 5,6 & HF > 20", "l");
  legend_ABCD_subcategories->AddEntry(ABCD_validation[2],    "N_{ch} = 7,8 & 4 #leq HF < 20", "l");
  legend_ABCD_subcategories->AddEntry(ABCD_validation[3],    "N_{ch} = 5,6 & 4 #leq HF < 20", "l");
  legend_ABCD_subcategories->Draw();
  
  c_ABCD_subcategories->SaveAs((basePlotDir+"/singlePlot/ABCD_subcategories.png").c_str());
  c_ABCD_subcategories->SaveAs((basePlotDir+"/singlePlot/ABCD_subcategories.pdf").c_str());
  
  
  TCanvas *c_ABCD_validation = new TCanvas("c_ABCD_validation", "c_ABCD_validation", 800, 800);
  ABCD_validation[4]->Divide(ABCD_validation[0]);
  ABCD_validation[5]->Divide(ABCD_validation[1]);
  ABCD_validation[4]->SetMaximum(1.2 * ABCD_validation[5]->GetMaximum());
  for (int cat = 4; cat < 6; cat++) ABCD_validation[cat]->Draw("he0e1x0same");
  CMS_lumi( c_ABCD_validation, iPeriod, iPos );
  
  TLegend *legend_ABCD_validation = new TLegend(0.65,0.75,0.8,0.9);
  legend_ABCD_validation->SetFillStyle(0); legend_ABCD_validation->SetFillColor(0);
  legend_ABCD_validation->SetLineColor(0); legend_ABCD_validation->SetShadowColor(0); legend_ABCD_validation->SetTextSize(0.035);
  legend_ABCD_validation->AddEntry((TObject*)0, "#frac{4 #leq HF < 20}{HF > 20}", "");
  legend_ABCD_validation->AddEntry((TObject*)0, "", "");
  legend_ABCD_validation->AddEntry(ABCD_validation[4],    "N_{ch} = 7,8", "l");
  legend_ABCD_validation->AddEntry(ABCD_validation[5],    "N_{ch} = 5,6", "l");
  legend_ABCD_validation->Draw();
  
  c_ABCD_validation->SaveAs((basePlotDir+"/singlePlot/ABCD_validation.png").c_str());
  c_ABCD_validation->SaveAs((basePlotDir+"/singlePlot/ABCD_validation.pdf").c_str());
  
  cout << "\nMuon charge times hadronic tau charge for Data:\n -1: " << charge_counter[0][0] << "\n  0: " << charge_counter[0][1] << "\n +1: " << charge_counter[0][2] << endl;

  // *********************************************************
  // *********************************************************
  // MC STARTS HERE
  // *********************************************************
  // *********************************************************
  
  //Pion SF table:
  float pt_bins[14] = {0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1};
  float pion_SF_hist_bins[15] = {0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 5};
  float eta_bins[3] = {0, 0.8, 1.5};
  float pion_rEff[3][14] = {
  {1.00, 1.01, 1.02, 1.03, 1.04, 1.07, 1.08, 1.05, 1.03, 1.02, 1.01, 1.00, 0.97, 0.96},
  {1.03, 1.03, 1.01, 1.04, 1.07, 1.06, 1.06, 1.12, 1.09, 1.10, 1.07, 1.02, 1.10, 1.13},
  {1.00, 0.99, 1.00, 0.99, 1.05, 1.05, 1.10, 1.10, 1.09, 1.11, 1.12, 1.10, 1.05, 1.04}};
  
  float pion_rEff_error[3][14] = {
  {0.03, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.04, 0.03, 0.04, 0.05, 0.05, 0.06, 0.06},
  {0.05, 0.05, 0.04, 0.08, 0.05, 0.04, 0.09, 0.05, 0.07, 0.06, 0.07, 0.08, 0.10, 0.18},
  {0.00, 0.13, 0.12, 0.07, 0.08, 0.08, 0.09, 0.12, 0.11, 0.13, 0.16, 0.13, 0.18, 0.24}};
  
  TH1F *h_pion_SF[3];
  h_pion_SF[0] = new TH1F("h_pion_SF_1","pion SF;pion p_{T};pion SF (#epsilon_{data}/#epsilon_{MC})",14,pion_SF_hist_bins);
  h_pion_SF[1] = new TH1F("h_pion_SF_2","pion SF;pion p_{T};pion SF (#epsilon_{data}/#epsilon_{MC})",14,pion_SF_hist_bins);
  h_pion_SF[2] = new TH1F("h_pion_SF_3","pion SF;pion p_{T};pion SF (#epsilon_{data}/#epsilon_{MC})",14,pion_SF_hist_bins);
  for (int ebin = 0; ebin < 3; ebin++){
    h_pion_SF[ebin]->SetLineColor(colors[ebin+1]); h_pion_SF[ebin]->SetMarkerColor(colors[ebin+1]);
    for (int pbin = 0; pbin < 14; pbin++){
      h_pion_SF[ebin]->SetBinContent(pbin+1,pion_rEff[ebin][pbin]);
      h_pion_SF[ebin]->SetBinError(pbin+1,pion_rEff_error[ebin][pbin]);
    }
  }
  
  TCanvas *c_pionSF = new TCanvas("c_pionSF", "c_pionSF", 800, 800);
  for (int ebin = 0; ebin < 3; ebin++) h_pion_SF[ebin]->Draw("he0e1x0same");
  h_pion_SF[0]->GetYaxis()->SetRangeUser(0.77,1.35);
  //CMS_lumi( c_pionSF, iPeriod, iPos );
  
  TLegend *legend_pionSF = new TLegend(0.39,0.18,0.54,0.33);
  legend_pionSF->SetFillStyle(0);
  legend_pionSF->SetFillColor(0); legend_pionSF->SetLineColor(0); legend_pionSF->SetShadowColor(0); legend_pionSF->SetTextSize(0.035);
  legend_pionSF->AddEntry(h_pion_SF[0],    "|#eta| < 0.8", "l");
  legend_pionSF->AddEntry(h_pion_SF[1],    "0.8 #leq |#eta| < 1.5", "l");
  legend_pionSF->AddEntry(h_pion_SF[2],    "1.5 #leq |#eta|", "l");
  legend_pionSF->Draw();
  
  c_pionSF->SaveAs((basePlotDir+"/singlePlot/pion_SF.png").c_str());
  c_pionSF->SaveAs((basePlotDir+"/singlePlot/pion_SF.pdf").c_str());
  
  // 3prong tau SF ups and down per pion pt and eta bins - treat each bin uncorrelated - 3*14 nuisance parameters
  TH1F *h_tauSFup[3][14];
  TH1F *h_tauSFdown[3][14];
  TH1F *h_tauSFmean = new TH1F("pionSF_mean","pionSF_mean",deltaphi_bins, 0, TMath::Pi());
  h_tauSFmean->Sumw2();
  for (int ebin = 0; ebin < 3; ebin++){
    for (int pbin = 0; pbin < 14; pbin++){
      h_tauSFup[ebin][pbin] = new TH1F(("pionSF_pt"+to_string(int(10*pt_bins[pbin]))+"_eta"+to_string(int(10*eta_bins[ebin]))+"Up").c_str(),("pionSF_pt"+to_string(int(10*pt_bins[pbin]))+"_eta"+to_string(int(10*eta_bins[ebin]))+"Up").c_str(),deltaphi_bins, 0, TMath::Pi());
      h_tauSFdown[ebin][pbin] = new TH1F(("pionSF_pt"+to_string(int(10*pt_bins[pbin]))+"_eta"+to_string(int(10*eta_bins[ebin]))+"Down").c_str(),("pionSF_pt"+to_string(int(10*pt_bins[pbin]))+"_eta"+to_string(int(10*eta_bins[ebin]))+"Down").c_str(),deltaphi_bins, 0, TMath::Pi());
      h_tauSFup[ebin][pbin]->Sumw2(); h_tauSFdown[ebin][pbin]->Sumw2();
    }
  }
  
  int entriesMC;
  //TLorentzVector tau_muon_mc, tau_hadron_mc;
  double temp_tau_rho1 = 0.; double temp_tau_rho2 = 0.;
  bool passedmu  = false;
  bool passedtau = false;
  bool passedcalo = true;
  bool passedGamma = true;
  double temp_tau_pt_comparison = 0.;
  int temp_nch = 0.;
  double temp_npv = 0.;
  double temp_pvz = 0.;
  int matchedpT_rank = 0;
  int acceptedGen = 0;
  float acceptedReco = 0;
  float acceptedRecoRaw = 0;
  float acceptedRecoSquare = 0;
  bool matchedTauHad = false;
  int nMu3prong = 0;
  float gen_muon_pt_cut = 2.5;
  
  for (int s = 0; s < nSamples-1; s++){
  cout << "Starting " << tag[s+1] << endl;
  entriesMC = (TREEMCs[s]->fChain)->GetEntries();
  //if (s == 3) entriesMC /= 10;
  nMu3prong = cutflow[1]->GetBinContent(8);
  acceptedReco = 0;
  acceptedRecoRaw = 0;
  acceptedRecoSquare = 0;
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
    tau_hadron_ptScalar = 0;
    ditau_ptScalar = 0;
    minPionPionDeltaPhi = 1000;
    maxPionPionDeltaPhi = 0;
    minPionPionDeltaEta = 1000;
    maxPionPionDeltaEta = 0;
    minPionPionDeltaR = 1000;
    maxPionPionDeltaR = 0;
    minPionMuDeltaPhi = 1000;
    maxPionMuDeltaPhi = 0;
    minPionMuDeltaEta = 1000;
    maxPionMuDeltaEta = 0;
    minPionMuDeltaR = 1000;
    maxPionMuDeltaR = 0;
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
      if (tau_muon_isSoft){
        if (!TREEMCs[s]->BsTauTau_mu1_isSoft->at(i)) continue;
      }
      if (tau_muon_isGlobal){
        if (!TREEMCs[s]->BsTauTau_mu1_isGlobal->at(i)) continue;
      }
      if (tau_muon_isTracker){
        if (!TREEMCs[s]->BsTauTau_mu1_isTracker->at(i)) continue;
      }
      if( (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < 3.5 && TMath::Abs(TREEMCs[s]->BsTauTau_mu1_eta->at(i)) < 1.2) || (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < 2.5 && TMath::Abs(TREEMCs[s]->BsTauTau_mu1_eta->at(i)) > 1.2))  continue;
      //if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) > 2.5) passedmu = true;
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_mu1_pt->at(i);
        tau_muon.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_mu1_pt->at(i), TREEMCs[s]->BsTauTau_mu1_eta->at(i), TREEMCs[s]->BsTauTau_mu1_phi->at(i), muon_mass);
        muon_charge = TREEMCs[s]->BsTauTau_mu1_q->at(i);
      }
    }    
    
    
    if (!TREEMCs[s]->BsTauTau_nch->empty()) temp_nch = TREEMCs[s]->BsTauTau_nch->at(0);
    if (temp_nch>=min_nch && temp_nch<=max_nch) passedNch = true;
    
    if (!TREEMCs[s]->BsTauTau_bbPV_vz->empty()) temp_pvz = TREEMCs[s]->BsTauTau_bbPV_vz->at(0);
    temp_npv = TREEMCs[s]->PV_N;
    h_tau_hadron_nPions[s+1]->Fill(TREEMCs[s]->BsTauTau_nPions->at(0));
    
    /*
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
    */

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
    float muon_weight_trg = tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),0);
    float muon_weight_trk = tnp_weight_trk_pbpb(0);
    muon_weight = muon_weight_trg*muon_weight_trk;
    muon_weight_error = 0;
    candCounter = 0;
    int chosenTauIndex = -1;
    temp_tau_rho1 = 0.; temp_tau_rho2 = 0.;
    
    //if(!TREEMCs[s]->BsTauTau_tau_pt) continue;
    float vprob = -1;
    float temp_pt = 0;
    for (int i=0; i<(int)TREEMCs[s]->BsTauTau_tau_pt->size(); i++) {
      if (TREEMCs[s]->cand_leadingPionPt->at(i) > pionLeadingCut && TREEMCs[s]->cand_subleadingPionPt->at(i) > pionSubLeadingCut && TREEMCs[s]->cand_subsubleadingPionPt->at(i) > pionSubSubLeadingCut){
        candCounter++;
        chosenTauIndex = i;
        tau_hadron_ptScalar = TREEMCs[s]->cand_leadingPionPt->at(i) + TREEMCs[s]->cand_subleadingPionPt->at(i) + TREEMCs[s]->cand_subsubleadingPionPt->at(i);
      }
      else continue;
      if (muon_charge*TREEMCs[s]->BsTauTau_tau_q->at(i) != MuTauCharge) continue;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) < 2) continue;
      //if (temp_nch != 1 && TREEMCs[s]->BsTauTau_tau_pt->at(i) < 2.0) continue;
      //if(TMath::Abs(TREEMCs[s]->BsTauTau_tau_eta->at(i)) > 2) continue;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_pt) {
        vprob = 100*TREEMCs[s]->BsTauTau_tau_vprob->at(i);
        if (temp_nch == 1) vprob = 100*TREEMCs[s]->BsTauTau_B_vprob->at(i);
        temp_pt = TREEMCs[s]->BsTauTau_tau_pt->at(i);
      }
      if (TREEMCs[s]->BsTauTau_tau_vprob->at(i) < tau_hadron_vertexprob) continue;
      passedtau = true;
      pion[0].SetPtEtaPhiM (TREEMCs[s]->cand_leadingPionPt->at(i),TREEMCs[s]->cand_leadingPionEta->at(i),TREEMCs[s]->cand_leadingPionPhi->at(i),pionMass);
      pion[1].SetPtEtaPhiM (TREEMCs[s]->cand_subleadingPionPt->at(i),TREEMCs[s]->cand_subleadingPionEta->at(i),TREEMCs[s]->cand_subleadingPionPhi->at(i),pionMass);
      pion[2].SetPtEtaPhiM (TREEMCs[s]->cand_subsubleadingPionPt->at(i),TREEMCs[s]->cand_subsubleadingPionEta->at(i),TREEMCs[s]->cand_subsubleadingPionPhi->at(i),pionMass);
      //if (temp_nch == 1 && TREEMCs[s]->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) / TREEMCs[s]->BsTauTau_tau_matched_gentaupt->at(i) > 0.85) matchedpT_rank += 1;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        tauh_charge = TREEMCs[s]->BsTauTau_tau_q->at(i);
        //if (s == 0) tau_weight *= 0.985;
        //if (s == 0) tau_weight *= TREEMCs[s]->BsTauTau_tau_rEff->at(i);
        // Mixing tau and muon SF
        //tau_weight *= muonSF;
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        tau_hadron.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_tau_pt->at(i), TREEMCs[s]->BsTauTau_tau_eta->at(i), TREEMCs[s]->BsTauTau_tau_phi->at(i), TREEMCs[s]->BsTauTau_tau_mass->at(i));
        if (threeProng[s+1]){
          temp_tau_rho1 = TREEMCs[s]->BsTauTau_tau_rhomass1->at(i);
          if (temp_nch>=3) temp_tau_rho2 = TREEMCs[s]->BsTauTau_tau_rhomass2->at(i);
        }
      }
    } // loop over the size of the tau candidates
    if (candCounter != 1) passedNch = false;
    bool passedDitau = (tau_muon+tau_hadron).M() > ditauMassCut;
    
    // Applying pion SF and up/down variations
    float tauSFup[3][14];
    float tauSFdown[3][14];
    float tauSF = 1;
    float tauSFmean = 1; //should be equal to tauSF. This is a crosscheck
    if (passedtau){
      int etaIndex[3] = {0,0,0};
      int ptIndex[3] = {0,0,0};
      for (int etabin = 0; etabin < 3; etabin++){
        if (TMath::Abs(TREEMCs[s]->cand_leadingPionEta->at(chosenTauIndex)) > eta_bins[etabin]) etaIndex[0] = etabin;
        if (TMath::Abs(TREEMCs[s]->cand_subleadingPionEta->at(chosenTauIndex)) > eta_bins[etabin]) etaIndex[1] = etabin;
        if (TMath::Abs(TREEMCs[s]->cand_subsubleadingPionEta->at(chosenTauIndex)) > eta_bins[etabin]) etaIndex[2] = etabin;
      }
      for (int pTbin = 0; pTbin < 14; pTbin++){
        if (TREEMCs[s]->cand_leadingPionPt->at(chosenTauIndex) > pt_bins[pTbin]) ptIndex[0] = pTbin;
        if (TREEMCs[s]->cand_subleadingPionPt->at(chosenTauIndex) > pt_bins[pTbin]) ptIndex[1] = pTbin;
        if (TREEMCs[s]->cand_subsubleadingPionPt->at(chosenTauIndex) > pt_bins[pTbin]) ptIndex[2] = pTbin;
      }
      
      for (int p = 0; p < 3; p++) tauSF *= pion_rEff[etaIndex[p]][ptIndex[p]];
      
      for (int etabin = 0; etabin < 3; etabin++){
        for (int pTbin = 0; pTbin < 14; pTbin++){
          tauSFup[etabin][pTbin] = 1;
          tauSFdown[etabin][pTbin] = 1;
          tauSFmean = 1;
          //if (etabin==2 && pTbin==13) cout << "**** **** **** ****" << endl;
          for (int p = 0; p < 3; p++){
            if (etabin == etaIndex[p] && pTbin == ptIndex[p]){
              tauSFup[etabin][pTbin] *= pion_rEff[etabin][pTbin] + pion_rEff_error[etabin][pTbin];
              tauSFdown[etabin][pTbin] *= pion_rEff[etabin][pTbin] - pion_rEff_error[etabin][pTbin];
              tauSFmean *= pion_rEff[etabin][pTbin];
            }
            else{
              tauSFup[etabin][pTbin] *= pion_rEff[etaIndex[p]][ptIndex[p]];
              tauSFdown[etabin][pTbin] *= pion_rEff[etaIndex[p]][ptIndex[p]];
              tauSFmean *= pion_rEff[etaIndex[p]][ptIndex[p]];
            }
            //if (etabin==2 && pTbin==13) cout << "tauSF: " << tauSF << " tauSFmean: " << tauSFmean << " etaIndex: " << etaIndex[p] << " ptIndex: " << ptIndex[p] << endl;
          } // loop on pions
        } // loop on pt bins
      } // loop on eta bins
    } // if passedtau
    
    
    
    ditau_ptScalar = tau_hadron_ptScalar + tau_muon.Pt();
    
    for (int p1 = 0; p1 < 3; p1++){
      float tempDeltaPhiPionMu = TMath::Abs(tau_muon.DeltaPhi(pion[p1]));
      float tempDeltaEtaPionMu = TMath::Abs(tau_muon.Eta()-pion[p1].Eta());
      float tempDeltaRPionMu = TMath::Abs(tau_muon.DeltaR(pion[p1]));
      if (minPionMuDeltaPhi > tempDeltaPhiPionMu) minPionMuDeltaPhi = tempDeltaPhiPionMu;
      if (maxPionMuDeltaPhi < tempDeltaPhiPionMu) maxPionMuDeltaPhi = tempDeltaPhiPionMu;
      if (minPionMuDeltaEta > tempDeltaEtaPionMu) minPionMuDeltaEta = tempDeltaEtaPionMu;
      if (maxPionMuDeltaEta < tempDeltaEtaPionMu) maxPionMuDeltaEta = tempDeltaEtaPionMu;
      if (minPionMuDeltaR > tempDeltaRPionMu) minPionMuDeltaR = tempDeltaRPionMu;
      if (maxPionMuDeltaR < tempDeltaRPionMu) maxPionMuDeltaR = tempDeltaRPionMu;
      for (int p2 = p1+1; p2 < 3; p2++){
        float tempDeltaPhiPionPion = TMath::Abs(pion[p2].DeltaPhi(pion[p1]));
        float tempDeltaEtaPionPion = TMath::Abs(pion[p2].Eta()-pion[p1].Eta());
        float tempDeltaRPionPion = TMath::Abs(pion[p2].DeltaR(pion[p1]));
        if (minPionPionDeltaPhi > tempDeltaPhiPionPion) minPionPionDeltaPhi = tempDeltaPhiPionPion;
        if (maxPionPionDeltaPhi < tempDeltaPhiPionPion) maxPionPionDeltaPhi = tempDeltaPhiPionPion;
        if (minPionPionDeltaEta > tempDeltaEtaPionPion) minPionPionDeltaEta = tempDeltaEtaPionPion;
        if (maxPionPionDeltaEta < tempDeltaEtaPionPion) maxPionPionDeltaEta = tempDeltaEtaPionPion;
        if (minPionPionDeltaR > tempDeltaRPionPion) minPionPionDeltaR = tempDeltaRPionPion;
        if (maxPionPionDeltaR < tempDeltaRPionPion) maxPionPionDeltaR = tempDeltaRPionPion;
      }
    }
    
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
    double sumE = 0;
    double leadingE = 0;
    double leadingE_Eta = -100;
    double sumEt = 0;
    double leadingEt = 0;
    double leadingEt_Eta = -100;
    
    if (triggered && passedmu && passedtau && passedNch)
    {
      for (int i=0; i<(int)TREEMCs[s]->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREEMCs[s]->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREEMCs[s]->BsTauTau_calo_energyHFm->at(i);
        /*if (passedNch) h_calo_E[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_energy->at(i));
        if (passedNch) h_calo_E_eta[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eta->at(i),TREEMCs[s]->BsTauTau_calo_energy->at(i));
        //if (passedNch) h_calo_Et[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eT->at(i));
        //if (passedNch) h_calo_Et_eta[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eta->at(i),TREEMCs[s]->BsTauTau_calo_eT->at(i));
        sumE += TREEMCs[s]->BsTauTau_calo_energy->at(i);
        sumEt += TREEMCs[s]->BsTauTau_calo_eT->at(i);*/
        double etaCal = TREEMCs[s]->BsTauTau_calo_eta->at(i);
        /*if (TREEMCs[s]->BsTauTau_calo_energy->at(i) > leadingE){
          leadingE = TREEMCs[s]->BsTauTau_calo_energy->at(i);
          leadingE_Eta = etaCal;
        }
        if (TREEMCs[s]->BsTauTau_calo_eT->at(i) > leadingEt){
          leadingEt = TREEMCs[s]->BsTauTau_calo_eT->at(i);
          leadingEt_Eta = etaCal;
        }*/
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
      if (passedcalo && passedNch) h_calo_sumE[s+1]->Fill(sumE);
      if (passedcalo && passedNch) h_calo_leadingE[s+1]->Fill(leadingE);
      if (passedNch) h_calo_leadingE_eta[s+1]->Fill(leadingE_Eta,leadingE);
      if (passedcalo && passedNch) h_calo_sumEt[s+1]->Fill(sumEt);
      if (passedcalo && passedNch) h_calo_leadingEt[s+1]->Fill(leadingEt);
      if (passedNch) h_calo_leadingEt_eta[s+1]->Fill(leadingEt_Eta,leadingEt);
      //if (passedcalo && passedNch) h_nCaloTowers[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eta->size());
      if (passedcalo && passedNch) h_calo_energyHFp_sum[s+1]->Fill(sumHFp,tauSF*muon_weight);
      if (passedcalo && passedNch) h_calo_energyHFm_sum[s+1]->Fill(sumHFm,tauSF*muon_weight);
      if (passedcalo && passedNch) h_calo_energyHF_pm[s+1]->Fill(sumHFm,sumHFp,tauSF*muon_weight);
      //h_calo_energyHFp_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),maxHFp); h_calo_energyHFm_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),maxHFm,tauSF*muon_weight);
      if (passedcalo && passedNch) h_calo_energyHFp_size[s+1]->Fill(sizeHFp,tauSF*muon_weight);
      if (passedcalo && passedNch) h_calo_energyHFm_size[s+1]->Fill(sizeHFm,tauSF*muon_weight);
      if (passedNch && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFp[s+1]->Fill(maxHFp,tauSF*muon_weight);
      if (passedNch && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFm[s+1]->Fill(maxHFm,tauSF*muon_weight);
      if (temp_nch>=5 && temp_nch<=8 && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFp_highNch[s+1]->Fill(maxHFp,tauSF*muon_weight);
      if (temp_nch>=5 && temp_nch<=8 && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low) h_calo_leadingHFm_highNch[s+1]->Fill(maxHFm,tauSF*muon_weight);
      if (passedNch) h_calo_leadingHF_pm[s+1]->Fill(maxHFm,maxHFp,tauSF*muon_weight);
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    }
    
    float gen_tau_hadron_pt = -1000;
    float gen_tau_hadron_visible_pt = 0;
    TLorentzVector gen_tau_hadron_visible;
    int nAccGenPions = 0;
    int nLeadingPion = 0;
    int nSubLeadingPion = 0;
    int nSubSubLeadingPion = 0;
    int nMatchedPions = 0;
    float gen_muon_pt = -1;
    if (gen_tautau_to_mu3prong){
      float leadingGenPt = 0;
      float subleadingGenPt = 0;
      float subsubleadingGenPt = 1000;
      for (int i=0; i<(int)TREEMCs[s]->gen_tau_daughter_pdgId->size(); i++) { // loop on gen daughters
        if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==13){
          if((TREEMCs[s]->gen_tau_daughter_pt->at(i) < 3.5 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i)) < 1.2) || (TREEMCs[s]->gen_tau_daughter_pt->at(i) < 2.5 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i)) >= 1.2)) continue;
          if(TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))>2.4) continue;
          gen_muon_pt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
        }
        if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<2.5){
        //if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<1){
          if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > pionLeadingCut) nLeadingPion++;
          if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > pionSubLeadingCut) nSubLeadingPion++;
          if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > pionSubSubLeadingCut) nSubSubLeadingPion++;
          nAccGenPions++;
          if (leadingGenPt < TREEMCs[s]->gen_tau_daughter_pt->at(i)) leadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
          if (subsubleadingGenPt > TREEMCs[s]->gen_tau_daughter_pt->at(i)) subsubleadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
          h_N_genPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
          h_N_genPionEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i));
          h_N_genPionPtEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i),TREEMCs[s]->gen_tau_daughter_pt->at(i));
          if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<1) h_N_genCentralPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
          if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > 2) h_N_genHighPtPionEta[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_eta->at(i));
          gen_tau_hadron_visible_pt += TREEMCs[s]->gen_tau_daughter_pt->at(i); //wrong. They should be added as vectors. 
          TLorentzVector temp_gen;
          temp_gen.SetPtEtaPhiM(TREEMCs[s]->gen_tau_daughter_pt->at(i),TREEMCs[s]->gen_tau_daughter_eta->at(i),TREEMCs[s]->gen_tau_daughter_phi->at(i),0.13957);
          gen_tau_hadron_visible += temp_gen;
          bool found_pion_match = false;
          for (int j=0; j<(int)TREEMCs[s]->reco_pion_pt->size(); j++){
            if(found_pion_match) continue;
            TLorentzVector temp_reco;
            temp_reco.SetPtEtaPhiM(TREEMCs[s]->reco_pion_pt->at(j),TREEMCs[s]->reco_pion_eta->at(j),TREEMCs[s]->reco_pion_phi->at(j),0.13957);
            if (TMath::Abs(temp_reco.Pt()/temp_gen.Pt()-1) < 0.07 && temp_gen.DeltaR(temp_reco) < 0.015){
              found_pion_match = true;
              nMatchedPions++;
              h_resPt_matchedPionReco[s+1]->Fill(temp_reco.Pt()/temp_gen.Pt());
              h_resR_matchedPionReco[s+1]->Fill(100*temp_gen.DeltaR(temp_reco));
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
      if (gen_tautau_to_mu3prong && nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && gen_tau_hadron_visible.Pt() >= 2 && gen_muon_pt > gen_muon_pt_cut/*2.5*/){
      //if (nAccGenPions == 3){
        h_N_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
        acceptedGen++;
      }
      if(gen_tautau_to_mu3prong && nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && nMatchedPions == 3 && temp_nch>=3){
      //if(nAccGenPions == 3 && nMatchedPions == 3 && temp_nch>=3){
        h_eff_3prongReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
        matchedTauHad = true;
      }
      if (passedtau && matchedTauHad && passedmu){
        for (int i=0; i<(int)TREEMCs[s]->reco_pion_pt->size(); i++){
          if (TREEMCs[s]->reco_pion_pt->size() != 3) continue;
          if (TMath::Abs(TREEMCs[s]->reco_pion_eta->at(i)) > 2.5) continue;
          float pionEff = h_eff_matchedPionPt->GetBinContent(h_eff_matchedPionPt->FindBin(TREEMCs[s]->reco_pion_pt->at(i)));
          h_recoPionPt[s+1]->Fill(TREEMCs[s]->reco_pion_pt->at(i),1.0/pionEff);
        }
        h_leading_recoPionPt[s+1]->Fill(TREEMCs[s]->cand_leadingPionPt->at(chosenTauIndex),1.0/h_eff_matchedPionPt->GetBinContent(h_eff_matchedPionPt->FindBin(TREEMCs[s]->cand_leadingPionPt->at(chosenTauIndex))));
        h_subleading_recoPionPt[s+1]->Fill(TREEMCs[s]->cand_subleadingPionPt->at(chosenTauIndex),1.0/h_eff_matchedPionPt->GetBinContent(h_eff_matchedPionPt->FindBin(TREEMCs[s]->cand_subleadingPionPt->at(chosenTauIndex))));
        h_subsubleading_recoPionPt[s+1]->Fill(TREEMCs[s]->cand_subsubleadingPionPt->at(chosenTauIndex),1.0/h_eff_matchedPionPt->GetBinContent(h_eff_matchedPionPt->FindBin(TREEMCs[s]->cand_subsubleadingPionPt->at(chosenTauIndex))));
        for (int i=0; i<(int)TREEMCs[s]->gen_tau_daughter_pdgId->size(); i++) {
          if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<2.5){
            h_N_matched_genPionPt[s+1]->Fill(TREEMCs[s]->gen_tau_daughter_pt->at(i));
            if (TREEMCs[s]->gen_tau_daughter_pt->at(i) != leadingGenPt && TREEMCs[s]->gen_tau_daughter_pt->at(i) != subsubleadingGenPt) subleadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
            /*if (TREEMCs[s]->gen_tau_daughter_pt->at(i) >= TREEMCs[s]->gen_tau_daughter_pt->at((i+1)%3) && TREEMCs[s]->gen_tau_daughter_pt->at(i) >= TREEMCs[s]->gen_tau_daughter_pt->at((i+2)%3)) leadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
            else if (TREEMCs[s]->gen_tau_daughter_pt->at(i) <= TREEMCs[s]->gen_tau_daughter_pt->at((i+1)%3) && TREEMCs[s]->gen_tau_daughter_pt->at(i) <= TREEMCs[s]->gen_tau_daughter_pt->at((i+2)%3)) subsubleadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
            else subleadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);*/
          }
        }
        h_N_matched_leading_genPionPt[s+1]->Fill(leadingGenPt);
        h_N_matched_subleading_genPionPt[s+1]->Fill(subleadingGenPt);
        h_N_matched_subsubleading_genPionPt[s+1]->Fill(subsubleadingGenPt);
      } // if passedtau and matchedTauHad and passedmu
    }
    
    TLorentzVector recoPiZero, genPiZero;
    int Ngamma = 0;
    int NgenPiZero = 0;
    float minDeltaR = 1000;
    float minDeltaM = 100000; //MeV
    if (!threeProng[s+1]){
      for (int i=0; i<(int)TREEMCs[s]->gen_neutral_pion_daughter_pdgId->size(); i++) {
        if (TMath::Abs(TREEMCs[s]->gen_neutral_pion_daughter_pdgId->at(i))==22){
          if (TREEMCs[s]->gen_neutral_pion_daughter_pt->at(i) < gammaPtCut || TMath::Abs(TREEMCs[s]->gen_neutral_pion_daughter_eta->at(i)) > 2.2) continue;
          gamma.SetPtEtaPhiM (TREEMCs[s]->gen_neutral_pion_daughter_pt->at(i),TREEMCs[s]->gen_neutral_pion_daughter_eta->at(i),TREEMCs[s]->gen_neutral_pion_daughter_phi->at(i),0);
          for (int j=i+1; j<(int)TREEMCs[s]->gen_neutral_pion_daughter_pdgId->size(); j++) {
            if (TMath::Abs(TREEMCs[s]->gen_neutral_pion_daughter_pdgId->at(j))==22){
              if (TREEMCs[s]->gen_neutral_pion_daughter_pt->at(j) < gammaPtCut || TMath::Abs(TREEMCs[s]->gen_neutral_pion_daughter_eta->at(j)) > 2.2) continue;
              tempGamma.SetPtEtaPhiM (TREEMCs[s]->gen_neutral_pion_daughter_pt->at(j),TREEMCs[s]->gen_neutral_pion_daughter_eta->at(j),TREEMCs[s]->gen_neutral_pion_daughter_phi->at(j),0);
              float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
              if (deltaR < minDeltaR) minDeltaR = deltaR;
              genPiZero = gamma+tempGamma;
              double recoPiZeroMass = 1000*genPiZero.M(); // MeV
              if (TMath::Abs(recoPiZeroMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = recoPiZeroMass-PiZeroMass;
              h_genPiZeroDeltaM[s+1]->Fill(recoPiZeroMass-PiZeroMass);
              if (TMath::Abs(recoPiZeroMass-PiZeroMass) < 50){
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
    full_tau_hadron = tau_hadron;
    for (int i=0; i<Ngamma; i++){
      if (TREEMCs[s]->reco_gamma_pt->at(i)<gammaPtCut || TMath::Abs(TREEMCs[s]->reco_gamma_eta->at(i)) > 2.2){
        Ngamma -= 1;
        continue;
      }
      gamma.SetPtEtaPhiM (TREEMCs[s]->reco_gamma_pt->at(i),TREEMCs[s]->reco_gamma_eta->at(i),TREEMCs[s]->reco_gamma_phi->at(i),0);
      //if (i == 0) recoNeutralPion = gamma;
      //else recoNeutralPion += gamma;
      h_recoGamma_pt[s+1]->Fill(TREEMCs[s]->reco_gamma_pt->at(i));
      h_recoGamma_eta[s+1]->Fill(TREEMCs[s]->reco_gamma_eta->at(i));
      //if (TMath::Abs(TREEMCs[s]->reco_gamma_eta->at(i)) > 2.2) passedGamma = false;
      h_recoGamma_phi[s+1]->Fill(TREEMCs[s]->reco_gamma_phi->at(i));
      h_recoGamma_deltaphi_muon[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_muon)));
      h_recoGamma_deltaphi_pion[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_hadron)));
      for (int j=i+1; j<Ngamma; j++){
        if (TREEMCs[s]->reco_gamma_pt->at(j)<gammaPtCut || TMath::Abs(TREEMCs[s]->reco_gamma_eta->at(j)) > 2.2) continue;
        tempGamma.SetPtEtaPhiM (TREEMCs[s]->reco_gamma_pt->at(j),TREEMCs[s]->reco_gamma_eta->at(j),TREEMCs[s]->reco_gamma_phi->at(j),0);
        float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
        if (deltaR < minDeltaR) minDeltaR = deltaR;
        recoPiZero = gamma+tempGamma;
        double recoPiZeroMass = 1000*recoPiZero.M(); // MeV
        if (TMath::Abs(recoPiZeroMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = recoPiZeroMass-PiZeroMass;
        h_recoPiZeroDeltaM[s+1]->Fill(recoPiZeroMass-PiZeroMass);
        if (TMath::Abs(recoPiZeroMass-PiZeroMass) < 50){
          h_recoPiZero_pt[s+1]->Fill(recoPiZero.Pt());
          h_recoPiZero_eta[s+1]->Fill(recoPiZero.Eta());
          h_recoPiZero_phi[s+1]->Fill(recoPiZero.Phi());
          h_recoPiZero_deltaphi_muon[s+1]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_muon)));
          h_recoPiZero_deltaphi_pion[s+1]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)));
          NrecoPiZero += 1;
          if (gamma.Pt()>gammaPtCut && tempGamma.Pt()>gammaPtCut) full_tau_hadron += recoPiZero;
        }
      }
    }
    h_NrecoPiZero[s+1]->Fill(NrecoPiZero);
    h_recoGammasMinDeltaR[s+1]->Fill(minDeltaR);
    h_recoPiZeroMinDeltaM[s+1]->Fill(minDeltaM);
    h_NrecoGamma[s+1]->Fill(Ngamma);
    }
    
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
    
    if (passedmu && passedcalo && passedGamma && gen_tautau_to_mu3prong && passedtau && muon_charge*tauh_charge != 0 && delta_phi >= deltaPhi_cut && passedNch && nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_N_recoMuonPt[s+1]->Fill(tau_muon.Pt());
    
    //if (passedmu && passedcalo && passedGamma && gen_tautau_to_mu3prong && passedtau && muon_charge*tauh_charge != 0 && delta_phi >= deltaPhi_cut && passedNch && nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_N_recoMuonPt[s+1]->Fill(gen_muon_pt);
    
    
    if (passedmu && passedGamma && triggered/* && gen_tautau_to_mu3prong && matchedTauHad && TREEMCs[s]->gen_tau_pt->at(0) > 0 && TREEMCs[s]->gen_tau_pt->at(1) > 0 && gen_muon_pt > gen_muon_pt_cut*/)
    {
      //if (!TREEMCs[s]->BsTauTau_tau_isRight->at(0)) continue;
      if (passedNch && passedcalo) h_tau_hadron_vprob[s+1]->Fill(vprob,tauSF*muon_weight);
      if (!passedtau && passedcalo) continue;
      if (passedNch && passedcalo) charge_counter[s+1][muon_charge*tauh_charge + 1] += 1;
      if (muon_charge*tauh_charge == 0 && passedcalo) continue;
      if(passedNch && gen_tautau_to_mu3prong && nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0) && passedcalo) h_eff_tauReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
      //if(gen_tautau_to_mu3prong && nAccGenPions == 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_eff_tauReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
      if (passedNch && passedcalo){
        h_deltaphi_tau_mu_tau_hadron_zoomed[s+1]->Fill(delta_phi,tauSF*muon_weight);
        h_deltaphi_tau_mu_full_tau_hadron[s+1]->Fill(full_delta_phi,tauSF*muon_weight);
        h_deltaphi_tau_mu_tau_hadron_mueta[s+1]->Fill(delta_phi,tau_muon.Eta(),tauSF*muon_weight);
        h_deltaphi_tau_mu_tau_hadron_deltaeta[s+1]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()),tauSF*muon_weight);
      }
      
      
      if (s==0 && passedNch && passedcalo){
        
        // #### Muon SF systematic uncertainties ####
        float muon_weight_sys_trigger = muon_weight_trk*max(abs(tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),-1)-muon_weight_trg),abs(tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),-2)-muon_weight_trg));
        
        float muon_weight_sys_tracking = muon_weight_trg*max(abs(tnp_weight_trk_pbpb(-1)-muon_weight_trk),abs(tnp_weight_trk_pbpb(-2)-muon_weight_trk));
        
        float muon_weight_sys_muid = muon_weight*max(abs(tnp_weight_muid_pbpb(tau_muon.Pt(),tau_muon.Eta(),-1)-1),abs(tnp_weight_muid_pbpb(tau_muon.Pt(),tau_muon.Eta(),-2)-1));
        
        float muon_weight_sys_sta = muon_weight*max(abs(tnp_weight_sta_pbpb(tau_muon.Pt(),tau_muon.Eta(),-1)-1),abs(tnp_weight_sta_pbpb(tau_muon.Pt(),tau_muon.Eta(),-2)-1));
        
        float muon_weight_sys_binned = tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),-10)*muon_weight_trk-muon_weight;
        
        
        // #### Muon SF statistical uncertainties (squared)####
        int nToysTrigger = 100;
        int nToys = 100;
        float muon_weight_stat_trigger = 0; // squared
        float muon_weight_stat_muid = 0; // squared
        float muon_weight_stat_sta = 0; // squared
        for (int toy = 1; toy <= nToys; toy++){
          muon_weight_stat_muid += pow(tnp_weight_muid_pbpb(tau_muon.Pt(),tau_muon.Eta(),toy)-1,2);
          muon_weight_stat_sta += pow(tnp_weight_sta_pbpb(tau_muon.Pt(),tau_muon.Eta(),toy)-1,2);
        }
        for (int toy = 1; toy <= nToysTrigger; toy++){
          muon_weight_stat_trigger += pow(tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),toy)-muon_weight_trg,2);
          //if (abs(tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),toy)-muon_weight_trg) > 0.4) cout << "SF: " << muon_weight_trg << " +- " << tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),toy)-muon_weight_trg << " pt: " << tau_muon.Pt() << " eta: " << tau_muon.Eta() << endl;
        }
        muon_weight_stat_trigger *= muon_weight_trk*muon_weight_trk/nToysTrigger;
        muon_weight_stat_muid *= muon_weight*muon_weight/nToys;
        muon_weight_stat_sta *= muon_weight*muon_weight/nToys;
        
        
        // #### Adding up all Muon SF uncertainties ####
        muon_weight_error = sqrt(pow(muon_weight_sys_trigger,2)+pow(muon_weight_sys_tracking,2)+pow(muon_weight_sys_muid,2)/*+pow(muon_weight_sys_sta,2)*/+pow(muon_weight_sys_binned,2)+muon_weight_stat_trigger+muon_weight_stat_muid/*+muon_weight_stat_sta*/);
        
        
        //float muon_weight_stat = sqrt(muon_weight_stat_trigger+muon_weight_stat_muid+muon_weight_stat_sta);
        //float muon_weight_sys = sqrt(pow(muon_weight_sys_trigger,2)+pow(muon_weight_sys_tracking,2)+pow(muon_weight_sys_muid,2)+pow(muon_weight_sys_sta,2)+pow(muon_weight_sys_binned,2));
        //cout << "muon SF: " << muon_weight << " +- " << muon_weight_stat << "(stat) +- " << muon_weight_sys << "(sys)" << endl;
        
        h_deltaphi_tau_mu_tau_hadron[s+1]->Fill(delta_phi,tauSF*muon_weight);
        h_sys_muon_SF_Down->Fill(delta_phi,tauSF*(muon_weight-muon_weight_error));
        h_sys_muon_SF_Up->Fill(delta_phi,tauSF*(muon_weight+muon_weight_error));
        
        // #### 3prong tau SF systematic uncertainties ####
        for (int ebin = 0; ebin < 3; ebin++){
          for (int pbin = 0; pbin < 14; pbin++){
            h_tauSFup[ebin][pbin]->Fill(delta_phi,tauSFup[ebin][pbin]*muon_weight);
            h_tauSFdown[ebin][pbin]->Fill(delta_phi,tauSFdown[ebin][pbin]*muon_weight);
          }
        }
        
        tauSF = 1;
        muon_weight = 1;
        
        h_tauSFmean->Fill(delta_phi,tauSFmean*muon_weight);
        
      } // if(s==0) for muon and pion SF uncertainties
      
        tauSF = 1;
        muon_weight = 1;
      
      
      if(delta_phi < deltaPhi_cut && passedcalo) continue;
      if (!passedcalo) h_tau_hadron_nch_highHF[s+1]->Fill(temp_nch,tauSF*muon_weight);
      if (!passedcalo) continue;
      //if (!TREEMCs[s]->BsTauTau_nPions->empty()) temp_nch = TREEMCs[s]->BsTauTau_nPions->at(0);
      h_tau_hadron_nch[s+1]->Fill(temp_nch,tauSF*muon_weight);
      //h_tau_hadron_nch[s+1]->Fill(temp_nch,tauSF*muon_weight);
      h_tau_hadron_ncand_final[s+1]->Fill(candCounter,tauSF*muon_weight);
      if (!passedNch) continue;
      h_minPionPionDeltaPhi[s+1]->Fill(minPionPionDeltaPhi,tauSF*muon_weight);
      h_maxPionPionDeltaPhi[s+1]->Fill(maxPionPionDeltaPhi,tauSF*muon_weight);
      h_minPionPionDeltaEta[s+1]->Fill(minPionPionDeltaEta,tauSF*muon_weight);
      h_maxPionPionDeltaEta[s+1]->Fill(maxPionPionDeltaEta,tauSF*muon_weight);
      h_minPionPionDeltaR[s+1]->Fill(minPionPionDeltaR,tauSF*muon_weight);
      h_maxPionPionDeltaR[s+1]->Fill(maxPionPionDeltaR,tauSF*muon_weight);
      h_minPionMuDeltaPhi[s+1]->Fill(minPionMuDeltaPhi,tauSF*muon_weight);
      h_maxPionMuDeltaPhi[s+1]->Fill(maxPionMuDeltaPhi,tauSF*muon_weight);
      h_minPionMuDeltaEta[s+1]->Fill(minPionMuDeltaEta,tauSF*muon_weight);
      h_maxPionMuDeltaEta[s+1]->Fill(maxPionMuDeltaEta,tauSF*muon_weight);
      h_minPionMuDeltaR[s+1]->Fill(minPionMuDeltaR,tauSF*muon_weight);
      h_maxPionMuDeltaR[s+1]->Fill(maxPionMuDeltaR,tauSF*muon_weight);
      //if(nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_eff_trigger[s+1]->Fill(gen_muon_pt);
      if(nLeadingPion >= 1 && nSubLeadingPion >= 2 && nSubSubLeadingPion >= 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_eff_trigger[s+1]->Fill(tau_muon.Pt());
      
      h_tau_hadron_matched_pt_index[s+1]->Fill(matchedpT_rank,tauSF*muon_weight);
      h_tau_mu_p[s+1]->Fill(tau_muon.P(),tauSF*muon_weight);
      h_tau_mu_pz[s+1]->Fill(tau_muon.Pz(),tauSF*muon_weight);
      h_tau_mu_dz[s+1]->Fill(10*(TREEMCs[s]->BsTauTau_mu1_vz->at(0)-TREEMCs[s]->BsTauTau_PV_vz->at(0)),tauSF*muon_weight);
      h_tau_mu_pt[s+1]->Fill(tau_muon.Pt(),tauSF*muon_weight);
      if (s==0) h_muon_pt_before_SF->Fill(tau_muon.Pt(),tauSF);
      if (s==0) h_muon_pt_after_SF->Fill(tau_muon.Pt(),tauSF*muon_weight);
      h_tau_mu_eta[s+1]->Fill(tau_muon.Eta(),tauSF*muon_weight);
      h_tau_mu_phi[s+1]->Fill(tau_muon.Phi(),tauSF*muon_weight);
      h_tau_mu_tau_hadron_phi[s+1]->Fill(tau_muon.Phi(),tau_hadron.Phi());
      h_tau_hadron_p[s+1]->Fill(tau_hadron.P(),tauSF*muon_weight);
      h_tau_hadron_pz[s+1]->Fill(tau_hadron.Pz(),tauSF*muon_weight);
      h_tau_hadron_ptScalar[s+1]->Fill(tau_hadron_ptScalar,tauSF*muon_weight);
      h_tau_hadron_pt[s+1]->Fill(tau_hadron.Pt(),tauSF*muon_weight);
      h_gen_tau_hadron_visible_pt[s+1]->Fill(gen_tau_hadron_visible.Pt(),tauSF*muon_weight);
      //if (gen_tautau_to_mu3prong/* && TREEMCs[s]->BsTauTau_tau_isRight->at(0)*/) h_resVisTauPt[s+1]->Fill(tau_hadron.Pt() - gen_tau_hadron_visible_pt);
      h_tau_hadron_eta[s+1]->Fill(tau_hadron.Eta(),tauSF*muon_weight);
      h_tau_hadron_phi[s+1]->Fill(tau_hadron.Phi(),tauSF*muon_weight);
      if (gen_tautau_to_mu3prong/* && !TREEMCs[s]->BsTauTau_tau_isRight->empty()*/){
        //if (!TREEMCs[s]->BsTauTau_tau_isRight->at(0)){
          h_resVisTauPt[s+1]->Fill(tau_hadron.Pt() - gen_tau_hadron_visible.Pt());
          h_resVisTauEta[s+1]->Fill(tau_hadron.Eta() - gen_tau_hadron_visible.Eta());
          h_resVisTauPhi[s+1]->Fill(tau_hadron.Phi() - gen_tau_hadron_visible.Phi());
        //}
      }
      h_tau_hadron_rhomass[s+1][0]->Fill(temp_tau_rho1,tauSF*muon_weight);
      h_tau_hadron_rhomass[s+1][1]->Fill(temp_tau_rho2,tauSF*muon_weight);
      h_tau_hadron_rhomass2D[s+1]->Fill(temp_tau_rho1,temp_tau_rho2,tauSF*muon_weight);
      h_tau_hadron_mass[s+1]->Fill(tau_hadron.M(),tauSF*muon_weight);
//      h_tau_hadron_track_pvz[s+1][0]->Fill(tau_track1.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][1]->Fill(tau_track2.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][2]->Fill(tau_track3.Z() - temp_pvz);
      h_calo_muEta_leadingHFeta[s+1]->Fill(tau_muon.Eta(),leadingHFeta);
      h_calo_tauEta_leadingHFeta[s+1]->Fill(tau_hadron.Eta(),leadingHFeta);
      h_calo_muEta_averageHFeta[s+1]->Fill(tau_muon.Eta(),averageHFeta);
      h_calo_tauEta_averageHFeta[s+1]->Fill(tau_hadron.Eta(),averageHFeta);
      h_mueta_taueta[s+1]->Fill(tau_muon.Eta(),tau_hadron.Eta(),tauSF*muon_weight);
      h_PV_N[s+1]->Fill(temp_npv,tauSF*muon_weight);
      TLorentzVector MET;
      MET.SetPtEtaPhiM((tau_muon+tau_hadron).Pt(),(tau_muon+tau_hadron).Eta(),-(tau_muon+tau_hadron).Phi(),0);
      h_MET[s+1]->Fill(MET.Pt(),tauSF*muon_weight);
      //h_ditau_mass[s+1]->Fill((tau_muon+tau_hadron+MET).M(),tauSF*muon_weight);
      h_ditau_mass[s+1]->Fill((tau_muon+tau_hadron).M(),tauSF*muon_weight);
      h_ditau_p[s+1]->Fill((tau_muon+tau_hadron).P(),tauSF*muon_weight);
      h_ditau_pz[s+1]->Fill((tau_muon+tau_hadron).Pz(),tauSF*muon_weight);
      h_ditau_pt[s+1]->Fill((tau_muon+tau_hadron).Pt(),tauSF*muon_weight);
      h_ditau_ptScalar[s+1]->Fill(ditau_ptScalar,tauSF*muon_weight);
      
      if (passedNch){
        h_pion_leading_pt[s+1]->Fill(TREEMCs[s]->cand_leadingPionPt->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_subleading_pt[s+1]->Fill(TREEMCs[s]->cand_subleadingPionPt->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_subsubleading_pt[s+1]->Fill(TREEMCs[s]->cand_subsubleadingPionPt->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_leading_eta[s+1]->Fill(TREEMCs[s]->cand_leadingPionEta->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_subleading_eta[s+1]->Fill(TREEMCs[s]->cand_subleadingPionEta->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_subsubleading_eta[s+1]->Fill(TREEMCs[s]->cand_subsubleadingPionEta->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_leading_phi[s+1]->Fill(TREEMCs[s]->cand_leadingPionPhi->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_subleading_phi[s+1]->Fill(TREEMCs[s]->cand_subleadingPionPhi->at(chosenTauIndex),tauSF*muon_weight);
        h_pion_subsubleading_phi[s+1]->Fill(TREEMCs[s]->cand_subsubleadingPionPhi->at(chosenTauIndex),tauSF*muon_weight);
      }
      
      h_AP[s+1]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2,tauSF*muon_weight);
      //h_MET[s+1]->Fill(TREEMCs[s]->MET_sumEt->at(0));
      
      //h_nHitsPixel[s+1]->Fill(TREEMCs[s]->nevtPixelHits->at(0));
      
      //if (matchedTauHad){
      if (true){ // Wrong: For efficiency, we make sure gen is solid but reco can be a fake just as in data efficiency.
        acceptedReco += tauSF*muon_weight;
        acceptedRecoSquare += pow(tauSF*muon_weight,2);
        acceptedRecoRaw++;
      }
    }
    
    
    if (passedmu && triggered && passedtau && muon_charge*tauh_charge == -1){
      if (temp_nch >= 5 && temp_nch <= 8 && !passedcalo) A_highNch_highHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
      if (temp_nch == 3 && passedNch  && !passedcalo) B_lowNch_highHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
      if (temp_nch >= 5 && temp_nch <= 8 && passedcalo) C_highNch_lowHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
      if (temp_nch == 3 && passedNch  && passedcalo) D_lowNch_lowHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
    }
    
  } // loop on entries
  
  cout << "\nMuon charge times hadronic tau charge for " << tag[s+1] << ":\n -1: " << charge_counter[s+1][0] << "\n  0: " << charge_counter[s+1][1] << "\n +1: " << charge_counter[s+1][2] << endl;
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
  
  if (nSamples > 0) h_muon_pt_before_SF->Scale(SF[1]);
  if (nSamples > 0) h_muon_pt_after_SF->Scale(SF[1]);
  
  for (int i = 0; i < nSamples; i++){
    //if (i != 0 && h_tau_hadron_nch[i]->GetMaximum()!=0) SF[i] = 0.8 * h_tau_hadron_nch[0]->GetMaximum()/h_tau_hadron_nch[i]->GetMaximum();
    //SF[i] = 1./h_tau_mu_pt[i]->Integral();
    //cout << "Temporary SF: " << SF[i] << endl;
    h_tau_mu_p[i]->Scale(SF[i]); h_tau_mu_pz[i]->Scale(SF[i]); h_tau_mu_pt[i]->Scale(SF[i]); h_tau_mu_eta[i]->Scale(SF[i]); h_tau_mu_phi[i]->Scale(SF[i]);
    h_tau_mu_dz[i]->Scale(SF[i]);
    h_tau_hadron_p[i]->Scale(SF[i]); h_tau_hadron_pz[i]->Scale(SF[i]); h_tau_hadron_pt[i]->Scale(SF[i]);
    //h_gen_tau_hadron_visible_pt[i]->Scale(SF[i]);
    h_gen_tau_hadron_visible_pt[i]->Scale(1./h_gen_tau_hadron_visible_pt[i]->GetMaximum());
    h_tau_hadron_eta[i]->Scale(SF[i]); h_tau_hadron_phi[i]->Scale(SF[i]);
    for (int j = 0; j < 2; j++) h_tau_hadron_rhomass[i][j]->Scale(SF[i]);
    h_tau_hadron_nch[i]->Scale(SF[i]);
    h_tau_hadron_nch_highHF[i]->Scale(SF[i]);
    h_tau_hadron_mass[i]->Scale(SF[i]);
    //h_tau_hadron_mass[i]->Scale(1./h_tau_hadron_mass[i]->GetMaximum());
    h_tau_hadron_ncand_final[i]->Scale(SF[i]);
    //h_tau_hadron_vprob[i]->Scale(1./h_tau_hadron_vprob[i]->GetMaximum());
    h_tau_hadron_vprob[i]->Scale(SF[i]);
    for (int j = 0; j < 3; j++) h_tau_hadron_track_pvz[i][j]->Scale(SF[i]);
    h_deltaphi_tau_mu_tau_hadron[i]->Scale(SF[i]);
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->Scale(SF[i]);
    h_deltaphi_tau_mu_full_tau_hadron[i]->Scale(SF[i]);
    h_minPionPionDeltaPhi[i]->Scale(SF[i]);
    h_maxPionPionDeltaPhi[i]->Scale(SF[i]);
    h_minPionPionDeltaEta[i]->Scale(SF[i]);
    h_maxPionPionDeltaEta[i]->Scale(SF[i]);
    h_minPionPionDeltaR[i]->Scale(SF[i]);
    h_maxPionPionDeltaR[i]->Scale(SF[i]);
    h_minPionMuDeltaPhi[i]->Scale(SF[i]);
    h_maxPionMuDeltaPhi[i]->Scale(SF[i]);
    h_minPionMuDeltaEta[i]->Scale(SF[i]);
    h_maxPionMuDeltaEta[i]->Scale(SF[i]);
    h_minPionMuDeltaR[i]->Scale(SF[i]);
    h_maxPionMuDeltaR[i]->Scale(SF[i]);
    
    //if (i == 0) h_deltaphi_tau_mu_tau_hadron[i]->Scale(0.1);
    
    h_PV_N[i]->Scale(SF[i]);
    h_calo_energyHFp[i]->Scale(SF[i]);
    //h_calo_energyHFp_sum[i]->Scale(1./h_calo_energyHFp_sum[i]->GetMaximum());
    h_calo_Et[i]->Scale(SF[i]);
    h_calo_sumEt[i]->Scale(SF[i]);
    h_calo_leadingEt[i]->Scale(SF[i]);
    h_nCaloTowers[i]->Scale(SF[i]);
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
    h_calo_leadingHFp_highNch[i]->Scale(SF[i]);
    h_calo_leadingHFm_highNch[i]->Scale(SF[i]);
    h_MET[i]->Scale(SF[i]);
    //h_ditau_mass[i]->Scale(SF[i]);
    h_ditau_p[i]->Scale(SF[i]);
    h_ditau_pz[i]->Scale(SF[i]);
    h_ditau_pt[i]->Scale(SF[i]);
    //h_ditau_mass[i]->Scale(1./h_ditau_mass[i]->GetMaximum());
    h_ditau_mass[i]->Scale(SF[i]);
    //for (int bin=1; bin <=4; bin++) cutflow[i]->SetBinContent(bin,SF[i]*cutflow[i]->GetBinContent(bin));
    h_resVisTauPt[i]->Scale(SF[i]);
    h_tau_hadron_ptScalar[i]->Scale(SF[i]);
    h_ditau_ptScalar[i]->Scale(SF[i]);
    
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
    //h_NrecoGamma[i]->Scale(SF[i]);
    
    h_nHitsPixel[i]->Scale(SF[i]);
  }
  
  h_sys_muon_SF_Down->Scale(SF[1]);
  h_sys_muon_SF_Up->Scale(SF[1]);
  
  cout << "Signal: " << h_deltaphi_tau_mu_tau_hadron[1]->Integral() << " muon-SF down: " << h_sys_muon_SF_Down->Integral() << " muon-SF up: " << h_sys_muon_SF_Up->Integral() << " -> " << 100.0*(h_sys_muon_SF_Up->Integral()-h_sys_muon_SF_Down->Integral())/(2*h_deltaphi_tau_mu_tau_hadron[1]->Integral()) << " % variation" << endl;
  
  cout << "\npion SF variations:\n";
  for (int ebin = 0; ebin < 3; ebin++){
    for (int pbin = 0; pbin < 14; pbin++){
      h_tauSFup[ebin][pbin]->Scale(SF[1]);
      h_tauSFdown[ebin][pbin]->Scale(SF[1]);
      //cout << "pt: " << pt_bins[pbin] << " eta: " << eta_bins[ebin] << "\nSignal: " << h_deltaphi_tau_mu_tau_hadron[1]->Integral() << " tau-SF down: " << h_tauSFdown[ebin][pbin]->Integral() << " tau-SF up: " << h_tauSFup[ebin][pbin]->Integral() << " -> " << 100.0*(h_tauSFup[ebin][pbin]->Integral()-h_tauSFdown[ebin][pbin]->Integral())/(2*h_deltaphi_tau_mu_tau_hadron[1]->Integral()) << " % variation" << endl;
    }
  }
  h_tauSFmean->Scale(SF[1]);
  
  
  int shift = 0;
  float nPostSignal = 0;
  float nPostBackground = 0;
  
#ifdef post
  TFile *fitDiagnostics = new TFile("fitDiagnostics/fitDiagnostics.root","READ");
  
  TH1F *temp_post_background = (TH1F*) fitDiagnostics->Get("shapes_fit_s/ggtautau_3prong_1_2015/background");
  TH1F *temp_post_signal = (TH1F*) fitDiagnostics->Get("shapes_fit_s/ggtautau_3prong_1_2015/MC-signal");
  TH1F *post_background = new TH1F("post_background","#Delta#phi(#tau_{#mu}, #tau_{3prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  TH1F *post_signal = new TH1F("post_signal","#Delta#phi(#tau_{#mu}, #tau_{3prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  
  if ((1+deltaphi_bins-maxSBbin) != temp_post_signal->GetNbinsX()) shift = maxSBbin-1;
  for (int i = 1; i <= 1+deltaphi_bins-maxSBbin; i++){
    post_background->SetBinContent(i,temp_post_background->GetBinContent(i+shift));
    post_signal->SetBinContent(i,temp_post_signal->GetBinContent(i+shift));
  }
  nPostSignal = post_signal->Integral();
  nPostBackground = post_background->Integral();
#endif

  TH1F *background = (TH1F*)B_lowNch_highHF[0]->Clone("background");
  background->SetTitle("Estimated background from ABCD;#Delta#phi(#tau_{#mu}, #tau_{3prong});background in the signal region");
  background->Multiply(C_highNch_lowHF[0]);
  background->Divide(A_highNch_highHF[0]);
  background->SetMarkerStyle(1);
  //background->SetMarkerStyle(kFullTriangleDown);
  background->SetLineWidth(3);
  background->SetLineColor(colors[nSamples]);
  
  
  // Cross section and errors
  
  for (int i = 1; i < deltaphi_bins+1; i++){
    double S_error;
    double signal = D_lowNch_lowHF[1]->IntegralAndError(i, deltaphi_bins, S_error, "");
    //cout << "signal: " << signal << " +- " << S_error << endl;
    double B_error;
    double bckgrd = background->IntegralAndError(i, deltaphi_bins, B_error, "");
    //cout << "background: " << bckgrd << " +- " << B_error << endl;
    float S_B = signal / sqrt(signal + bckgrd);
    if (S_B > maxSB){
      maxSBbin = i;
      maxSB = S_B;
    }
    float S_B_error = S_B * sqrt(pow(S_error/signal,2) + 0.25*((pow(S_error,2)+pow(B_error,2))/pow(signal+bckgrd,2)));
    //cout << "signal / sqrt(signal+background): " << S_B << " +- " << S_B_error << endl;
    h_SB_deltaphi->SetBinContent(i,S_B);
    h_SB_deltaphi->SetBinError(i,S_B_error);
  }
  cout << "\nmaximum S/sqrt(S+B) at bin number " << maxSBbin << " corresponding to delta phi of " << maxSBbin * TMath::Pi()/deltaphi_bins << endl;
  
  
  maxSBbin = 1;
  
  
  TH1F *pre_background = new TH1F("pre_background","#Delta#phi(#tau_{#mu}, #tau_{3prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  TH1F *pre_signal = new TH1F("pre_signal","#Delta#phi(#tau_{#mu}, #tau_{3prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  TH1F *pre_data = new TH1F("pre_data","#Delta#phi(#tau_{#mu}, #tau_{3prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  for (int i = 1; i <= 1+deltaphi_bins-maxSBbin; i++){
    pre_background->SetBinContent(i,background->GetBinContent(i+maxSBbin-1));
    pre_background->SetBinError(i,background->GetBinError(i+maxSBbin-1));
    pre_signal->SetBinContent(i,h_deltaphi_tau_mu_tau_hadron[1]->GetBinContent(i+maxSBbin-1));
    pre_signal->SetBinError(i,h_deltaphi_tau_mu_tau_hadron[1]->GetBinError(i+maxSBbin-1));
    pre_data->SetBinContent(i,h_deltaphi_tau_mu_tau_hadron[0]->GetBinContent(i+maxSBbin-1));
    pre_data->SetBinError(i,h_deltaphi_tau_mu_tau_hadron[0]->GetBinError(i+maxSBbin-1));
  }
  pre_data->SetBinErrorOption(TH1::kPoisson);
  
  float nA = A_highNch_highHF[0]->Integral(maxSBbin, deltaphi_bins);
  float nB = B_lowNch_highHF[0]->Integral(maxSBbin, deltaphi_bins);
  float nC = C_highNch_lowHF[0]->Integral(maxSBbin, deltaphi_bins);
  float countBackgroundABCD = nB * nC / nA;
  
  float analysisEff = acceptedReco*1.0/acceptedGen;
  float analysisEffError = sqrt(acceptedRecoSquare * analysisEff * (1-analysisEff)) / acceptedGen;
  //float analysisAccEff = acceptedReco*1.0/nMu3prong;
  float analysisAccEff = acceptedReco*1.0/cutflow[1]->GetBinContent(1);
  cout << "#gen tau-tau events in all decay modes (with pT cut): " << cutflow[1]->GetBinContent(1) << endl;
  
  float BR = 2 * 0.1739 * 0.1455; // errors: 0.0004 and 0.0006
  int nData = pre_data->Integral();
  double nRecoErr;
  //cout << "Integral and error: " << h_MET[1]->IntegralAndError(1,h_MET[1]->GetNbinsX(),nRecoErr,"");
  //cout << " +- " << nRecoErr << endl;
  float sigmaFiducial = (nData - countBackgroundABCD) / (dataLumi*analysisEff*BR);
  
  float countBackErr2 = countBackgroundABCD*countBackgroundABCD*(1./nA + 1./nB + 1./nC);
  
  cout << endl << "#events in" << endl << "A_highNch_highHF: " << nA << endl << "B_lowNch_highHF: " << nB << endl << "C_highNch_lowHF: " << nC << endl << "D_lowNch_lowHF: " << nData << endl << "Estimated number of background with ABCD count: " << countBackgroundABCD << " +- " << sqrt(countBackErr2) << endl << endl;
  cout << "Estimated number of background events from ABCD shape: " << pre_background->Integral() << endl;
  float relEstatFid = sqrt( (nData+countBackErr2)/((nData-countBackgroundABCD)*(nData-countBackgroundABCD)));// + 1./acceptedRecoRaw + 1./acceptedGen );
  float statisticalErrorOnFiducial = sigmaFiducial * relEstatFid;
  
  cout << "\nmu+3prong: "<< nMu3prong << " , accepted gen: " << acceptedGen << " , reconstructed: " << acceptedReco << " (" << acceptedRecoRaw << " raw) , efficiency: " << analysisEff << " +- " << analysisEffError << " (" << analysisEffError/analysisEff << "% relative) , acceptance x efficiency: " << analysisAccEff << endl;
  
  cout << endl << " **** **** **** **** " << endl;
  cout << "Fiducial cross section: " << sigmaFiducial << " ub +- " << statisticalErrorOnFiducial << " (" << 100*statisticalErrorOnFiducial/sigmaFiducial << "%) (stat)" << endl;
  
  float sigmaInclusive = (nData - countBackgroundABCD) / (dataLumi*analysisAccEff*pT_SF);
  float relEstatInc = sqrt( (nData+countBackErr2)/((nData-countBackgroundABCD)*(nData-countBackgroundABCD)));// + 1./acceptedRecoRaw + 1./cutflow[1]->GetBinContent(1));// +  1./above3GeVSuperChic + 1./nAllEventsSuperChic);
  float statisticalErrorOnInclusive = sigmaInclusive * relEstatInc;
  
  cout << "Inclusive cross section: " << sigmaInclusive << " ub +- " << statisticalErrorOnInclusive << " (" << 100*statisticalErrorOnInclusive/sigmaInclusive << "%) (stat)" << endl;
  cout << " **** **** **** **** " << endl << endl;
  //cout << "Inclusive cross section (Method 2): " << (nData - countBackgroundABCD) * crossSectionMC[0] / h_tau_mu_pt[1]->Integral() << " nb" << endl;
  
  cout << "ABCD sys prefit HF variations: ";
  float variation_up = -1000;
  float variation_down = 1000;
  int varHFup = 0;
  int varHFdown = 0;
  for (int variation = 0; variation < 2; variation++){
    B_lowNch_highHF[nSamples+variation]->Multiply(C_highNch_lowHF[nSamples+variation]);
    B_lowNch_highHF[nSamples+variation]->Divide(A_highNch_highHF[nSamples+variation]);
    float var = B_lowNch_highHF[nSamples+variation]->Integral(maxSBbin, deltaphi_bins)-background->Integral(maxSBbin, deltaphi_bins);
    cout << var << " ";
    if (var > variation_up){
      variation_up = var;
      varHFup = variation;
    }
    if (var < variation_down){
      variation_down = var;
      varHFdown = variation;
    }
  }
  cout << "\nHF variation UP is " << variation_up << " at " << ABCDsysNames[varHFup] << endl;
  cout << "HF variation DOWN is " << variation_down << " at " << ABCDsysNames[varHFdown] << endl;
  
  
  cout << "\nABCD sys prefit nch variations: ";
  variation_up = -1000;
  variation_down = 1000;
  int varNchUp = 0;
  int varNchDown = 0;
  for (int variation = 2; variation < 2+nNchCategories; variation++){
    B_lowNch_highHF[nSamples+variation]->Multiply(C_highNch_lowHF[nSamples+variation]);
    B_lowNch_highHF[nSamples+variation]->Divide(A_highNch_highHF[nSamples+variation]);
    float var = B_lowNch_highHF[nSamples+variation]->Integral(maxSBbin, deltaphi_bins)-background->Integral(maxSBbin, deltaphi_bins);
    cout << var << " ";
    if (var > variation_up){
      variation_up = var;
      varNchUp = variation-2+firstNchCategory;
    }
    if (var < variation_down){
      variation_down = var;
      varNchDown = variation-2+firstNchCategory;
    }
  }
  cout << "\nnch variation UP is " << variation_up << " with nch = " << varNchUp << endl;
  cout << "nch variation DOWN is " << variation_down << " with nch = " << varNchDown << endl << endl;
  
  TH1F *background_sysNchDown = (TH1F*)B_lowNch_highHF[nSamples+2]->Clone("ABCD-sys-nchDown");
  background_sysNchDown->SetTitle("Highest variation from ABCD background from nch - Down");
  //background_sysDown->Multiply(C_highNch_lowHF[nSamples]);
  //background_sysDown->Divide(A_highNch_highHF[nSamples]);
  
  TH1F *background_sysNchUp = (TH1F*)B_lowNch_highHF[nSamples+3]->Clone("ABCD-sys-nchUp");
  background_sysNchUp->SetTitle("Highest variation from ABCD background from nch - Up");
  //background_sysUp->Multiply(C_highNch_lowHF[nSamples+1]);
  //background_sysUp->Divide(A_highNch_highHF[nSamples+1]);

  
  TCanvas *c_tau_had_pv = new TCanvas("c_tau_had_pv", "c_tau_had_pv", 1500, 500); c_tau_had_pv->Divide(3,1);
  for (int j = 0; j < 3; j++){
    for (int i = 1; i < nSamples; i++) hs_tau_hadron_track_pvz[j]->Add(h_tau_hadron_track_pvz[i][j]);
    c_tau_had_pv->cd(j+1);
    //hs_tau_hadron_track_pvz[j]->Draw("he");
    h_tau_hadron_track_pvz[0][j]->Draw("e0e1x0");
    //c_tau_had_pv->cd(j+1); h_tau_hadron_track_pvz[0][j]->Draw("e0e1x0"); h_tau_hadron_track_pvz[0][j]->GetYaxis()->SetRangeUser(0., 0.6);
    for (int i = 1; i < nSamples; i++){
      h_tau_hadron_track_pvz[i][j]->Draw("hesame");
      h_tau_hadron_track_pvz[i][j]->GetYaxis()->SetRangeUser(0., 0.1);
    }
  }

  c_tau_had_pv->SaveAs((basePlotDir+"/collective/tau_had_pv."+plotFormat).c_str());
    
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
  
  h_deltaphi_tau_mu_tau_hadron[0]->SetBinErrorOption(TH1::kPoisson);
  
  h_deltaphi_tau_mu_tau_hadron[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum());
  h_deltaphi_tau_mu_tau_hadron[0]->SetMinimum(0.005);
  
  string binSize = to_string(deltaphi_bins).substr(0, 4);
  string t = "Events / (#pi/";
  h_deltaphi_tau_mu_tau_hadron[0]->GetYaxis()->SetTitle((t+binSize+")").c_str());
  //h_deltaphi_tau_mu_tau_hadron[0]->SetMarkerSize(0);
  h_deltaphi_tau_mu_tau_hadron[0]->Draw("e0e1x0");
  //c_delta_phi->SetGrid();
  //c_delta_phi->SetLogy(1);
  
  for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_tau_hadron[i]->Draw("hesame");   
  
  CMS_lumi( c_delta_phi, iPeriod, iPos );
  
  background->Draw("hesame");
  
  TLegend *tempGend = new TLegend(0.6,0.75,0.9,0.93);
  tempGend->SetFillStyle(0);
  tempGend->SetFillColor(0); tempGend->SetLineColor(0); tempGend->SetShadowColor(0); tempGend->SetTextSize(0.03);
  tempGend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data", "ep");
  for (int i = 1; i < nSamples; i++) tempGend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tags[i]).c_str(), "l");
  tempGend->AddEntry(background,    "ABCD background", "l");
  tempGend->Draw();
  
  c_delta_phi->SaveAs((basePlotDir+"/singlePlot/delta_phi.png").c_str());
  c_delta_phi->SaveAs((basePlotDir+"/singlePlot/delta_phi.pdf").c_str());
  
  c_delta_phi->SaveAs((basePlotDir+"/collective/delta_phi."+plotFormat).c_str());



  //gStyle->SetErrorX(0.5);
  TLegend *legend = new TLegend(0.6,0.76,0.9,0.90);
  legend->SetFillStyle(0);
  legend->SetFillColor(0); legend->SetLineColor(0); legend->SetShadowColor(0); legend->SetTextSize(0.035);
  //legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    (tag[0]).c_str(), "ep");
  legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data", "ep");
  for (int i = 1; i < nSamples; i++) legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tags[i]).c_str(), "l");
  legend->AddEntry(background,    "ABCD background", "l");
  //legend->Draw();
  //gStyle->SetErrorX(0.);
  
  TLegend *nobackgend = new TLegend(0.7,0.80,0.9,0.93);
  nobackgend->SetFillStyle(0);
  nobackgend->SetFillColor(0); nobackgend->SetLineColor(0); nobackgend->SetShadowColor(0); nobackgend->SetTextSize(0.03);
  nobackgend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data", "ep");
  for (int i = 1; i < nSamples; i++) nobackgend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tags[i]).c_str(), "l");
  
  TCanvas *tempCanvas;
  
  // loop on histograms
  for (int h = 0; h < histograms.size(); h++){
    string histogramsName = histograms.at(h)[0]->GetName();
    histograms.at(h)[nSamples+3] = (TH1F*)histograms.at(h)[nSamples+1]->Clone(("background - "+histogramsName).c_str());
    histograms.at(h)[nSamples+3]->Multiply(histograms.at(h)[nSamples+2]);
    histograms.at(h)[nSamples+3]->Divide(histograms.at(h)[nSamples]);
    histograms.at(h)[nSamples+3]->SetLineColor(colors[nSamples]);
    histograms.at(h)[nSamples+3]->SetMarkerStyle(styles[nSamples]);
    histograms.at(h)[nSamples+1]->SetBinContent(1,h); //hack to save the index of histograms. Needed later to plot the stacked histogram.
    histograms.at(h)[0]->SetBinErrorOption(TH1::kPoisson);
#ifdef post
    histograms.at(h)[nSamples+3]->Scale(nPostBackground/histograms.at(h)[nSamples+3]->Integral());
    if (nSamples > 1) histograms.at(h)[1]->Scale(nPostSignal/histograms.at(h)[1]->Integral());
#endif
    THStack *tempStack = new THStack(("stacked - "+histogramsName).c_str(),histogramsName.c_str());
    tempStack->Add(histograms.at(h)[nSamples+3]);
    if (nSamples > 1) tempStack->Add(histograms.at(h)[1]);
    if (tempStack->GetMaximum() > histograms.at(h)[0]->GetMaximum()) histograms.at(h)[0]->SetMaximum(1.2*tempStack->GetMaximum());
    stacks.push_back(tempStack);
    
#ifdef ratio
    int nbins = histograms.at(h)[0]->GetNbinsX();
    TH1F *tempRatio[2];
    tempRatio[0] = (TH1F*)histograms.at(h)[0]->Clone(("numerator - "+histogramsName).c_str());
    tempRatio[1] = new TH1F(*(TH1F*)((tempStack->GetStack()->Last())));
    
    float val_n, err_n, val_d, err_d;
    for (int b = 1; b <= nbins; b++){
      if (tempRatio[1]->GetBinContent(b) == 0){
        val_n = 0;
        err_n = 0;
        val_d = 1;
        err_d = 0;
      } else{
        val_n = histograms.at(h)[0]->GetBinContent(b) / tempRatio[1]->GetBinContent(b);
        err_n = histograms.at(h)[0]->GetBinError(b) / tempRatio[1]->GetBinContent(b);
        val_d = 1;
        err_d = tempRatio[1]->GetBinError(b) / tempRatio[1]->GetBinContent(b);
      }
      tempRatio[0]->SetBinContent(b,val_n);
      tempRatio[0]->SetBinError(b,err_n);
      tempRatio[1]->SetBinContent(b,val_d);
      tempRatio[1]->SetBinError(b,err_d);
      if (histogramsName == "h_deltaphi_tau_mu_tau_hadron_zoomed_"){
        cout << val_n << " " << err_n << " " << val_d << " " << err_d << endl;
      }
    }
    ratios.push_back(tempRatio);
    
    //plot the histogram and ratio
    //gStyle->SetErrorX(0.5);
    tempCanvas = new TCanvas(("c_ratio_"+plotNames[h]).c_str(), ("c_ratio_"+plotNames[h]).c_str(), 800, 800);
    
    // setup the pads
    float yplot = 0.66;
    float yratio = 0.28;
    float yspace = 0.01;
    
    float padTop_scale = 1./yplot;
    float padRatio_scale = 1./yratio;
    
    // set up the coordinates of the two pads:    //
    float y1, y2, y3, y4;                         //  y4 +-------------+
    y4 = 0.99;                                    //     |             |
    y3 = y4 - yplot;                              //     |     pad1    |
    y2 = y3 - yspace;                             //  y3 |-------------|
    y1 = y2 - yratio;                             //  y2 |-------------|
    float x1, x2;                                 //     |     rp1     |
    x1 = 0.01;                                    //  y1 +-------------+
    x2 = 0.99;                                    //     x1            x2
  
    TPad* pad_ratio = new TPad("rp1", "Ratio", x1, y1, x2, y2);
  
    pad_ratio->SetTopMargin(0.035);
    pad_ratio->SetBottomMargin(0.35);
    pad_ratio->SetLeftMargin(0.13);
    pad_ratio->SetRightMargin(0.025);
    
    TPad* pad_top = new TPad("pad1", "Control Plots", x1, y3, x2, y4);
  
    pad_top->SetTopMargin(0.07);
    pad_top->SetBottomMargin(0.03);
    pad_top->SetLeftMargin(0.13);
    pad_top->SetRightMargin(0.025);
    
    pad_ratio->Draw();
    pad_top->Draw();
    
    pad_top->cd();
    histograms.at(h)[0]->GetXaxis()->SetLabelSize(0.05*padTop_scale);
    histograms.at(h)[0]->GetYaxis()->SetLabelSize(0.05*padTop_scale);
    histograms.at(h)[0]->GetXaxis()->SetLabelOffset(0.017*padTop_scale);
    histograms.at(h)[0]->GetYaxis()->SetLabelOffset(0.007*padTop_scale);
    histograms.at(h)[0]->GetXaxis()->SetTitleSize(0.06*padTop_scale);
    histograms.at(h)[0]->GetYaxis()->SetTitleSize(0.06*padTop_scale);
    histograms.at(h)[0]->Draw("e0e1x0");
    tempStack->Draw("noclearhesame");
    
    TLegend *loopgend = new TLegend(0.6,0.76,0.9,0.90);
    loopgend->SetFillStyle(0);
    loopgend->SetFillColor(0); loopgend->SetLineColor(0); loopgend->SetShadowColor(0); loopgend->SetTextSize(0.035);
    loopgend->AddEntry(histograms.at(h)[0],    "Data", "ep");
    for (int i = 1; i < nSamples; i++) loopgend->AddEntry(histograms.at(h)[i], (tags[i]).c_str(), "l");
    loopgend->AddEntry(background,    "ABCD background", "l");
    loopgend->Draw();
    
    
    pad_ratio->Clear();
    pad_ratio->cd();
    tempRatio[0]->GetXaxis()->SetLabelSize(0.04*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetLabelSize(0.04*padRatio_scale);
    tempRatio[0]->GetXaxis()->SetLabelOffset(0.007*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetLabelOffset(0.007*padRatio_scale);
    tempRatio[0]->GetXaxis()->SetTitleSize(0.05*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetTitleSize(0.04*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetTitleOffset(0.11*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetNdivisions(504);
    
    //tempRatio[0]->GetYaxis()->SetRangeUser(0., 1.2*tempRatio[0]->GetMaximum());
    tempRatio[0]->GetYaxis()->SetRangeUser(-0.01, 2.01);
    //tempRatio[0]->Draw("AXIS");
    tempRatio[1]->SetLineColor(kBlack); tempRatio[1]->SetMarkerSize(0); //tempRatio[1]->SetMarkerStyle(1);
    tempRatio[1]->SetFillColor(kGray+1); tempRatio[1]->SetFillStyle(3013);
    tempRatio[0]->GetYaxis()->SetTitle("Data / Pred.");
    tempRatio[0]->Draw("e0e1x0");
    tempRatio[1]->Draw("e2same");
    TLine *line = new TLine(tempRatio[0]->GetXaxis()->GetXmin(), 1, tempRatio[0]->GetXaxis()->GetXmax(), 1);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    gPad->Modified(); gPad->Update();
    
    CMS_lumi( tempCanvas, iPeriod, iPos );
    gStyle->SetErrorX(0.5);
    tempCanvas->SaveAs((basePlotDir+"/ratioPlots/"+plotNames.at(h)+".png").c_str());
    tempCanvas->SaveAs((basePlotDir+"/ratioPlots/"+plotNames.at(h)+".pdf").c_str());
    tempCanvas->SaveAs((basePlotDir+"/ratioPlots/"+plotNames.at(h)+".C").c_str());
    ratioCanvas.push_back(tempCanvas);
    gStyle->SetErrorX(0.);
#endif
  } // loop on histograms
  
  


  
  TCanvas *c_delta_phi_zoomed = new TCanvas("c_delta_phi_zoomed", "c_delta_phi_zoomed", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMaximum() > h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetMaximum()) h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMaximum());
    if (h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMinimum() < h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetMinimum() && h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMinimum()  > 0) h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMinimum(0.8*h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMinimum());
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetLineWidth(3);
  }
  
  h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetMaximum());
  h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMinimum(0.005);
  
  h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetYaxis()->SetTitle((t+binSize+")").c_str());
  h_deltaphi_tau_mu_tau_hadron_zoomed[0]->Draw("e0e1x0");
  //c_delta_phi_zoomed->SetGrid();
  //c_delta_phi_zoomed->SetLogy(1);
  //for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_tau_hadron_zoomed[i]->Draw("hesame"); 
  stacks.at(h_deltaphi_tau_mu_tau_hadron_zoomed[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  
  CMS_lumi( c_delta_phi_zoomed, iPeriod, iPos );
  
  background->Draw("hesame");
  
  tempGend->Draw();
  
  c_delta_phi_zoomed->SaveAs((basePlotDir+"/singlePlot/delta_phi_zoomed.png").c_str());
  c_delta_phi_zoomed->SaveAs((basePlotDir+"/singlePlot/delta_phi_zoomed.pdf").c_str());
  
  c_delta_phi_zoomed->SaveAs((basePlotDir+"/collective/delta_phi_zoomed."+plotFormat).c_str());
  
  TCanvas *c_SB = new TCanvas("c_SB", "c_SB", 800, 800);
  h_SB_deltaphi->SetMaximum(1.2*h_SB_deltaphi->GetMaximum());
  h_SB_deltaphi->Draw("he");
  //c_SB->SetLogy(1);
  CMS_lumi( c_SB, iPeriod, iPos );
  
  c_SB->SaveAs((basePlotDir+"/singlePlot/S_sqrt_S+B_delta_phi.png").c_str());
  c_SB->SaveAs((basePlotDir+"/singlePlot/S_sqrt_S+B_delta_phi.pdf").c_str());
  
  TCanvas *c_full_delta_phi = new TCanvas("c_full_delta_phi", "c_full_delta_phi", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_deltaphi_tau_mu_full_tau_hadron[i]->GetMaximum() > h_deltaphi_tau_mu_full_tau_hadron[0]->GetMaximum()) h_deltaphi_tau_mu_full_tau_hadron[0]->SetMaximum(1.2*h_deltaphi_tau_mu_full_tau_hadron[i]->GetMaximum());
    if (h_deltaphi_tau_mu_full_tau_hadron[i]->GetMinimum() < h_deltaphi_tau_mu_full_tau_hadron[0]->GetMinimum() && h_deltaphi_tau_mu_full_tau_hadron[i]->GetMinimum()  > 0) h_deltaphi_tau_mu_full_tau_hadron[0]->SetMinimum(0.8*h_deltaphi_tau_mu_full_tau_hadron[i]->GetMinimum());
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetLineWidth(3);
  }
  h_deltaphi_tau_mu_full_tau_hadron[0]->Draw("e0e1x0");
  //c_full_delta_phi->SetGrid();
  c_full_delta_phi->SetLogy(1);
  
  for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_full_tau_hadron[i]->Draw("hesame"); 
  legend->Draw();
  c_full_delta_phi->SaveAs((basePlotDir+"/collective/full_delta_phi."+plotFormat).c_str());
  
  TCanvas *c_nHitsPixel = new TCanvas("c_nHitsPixel", "c_nHitsPixel", 800, 800);
  h_nHitsPixel[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_nHitsPixel[i]->Draw("hesame");
  //stacks.at(h_nHitsPixel[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_nHitsPixel, iPeriod, iPos );
  c_nHitsPixel->SaveAs((basePlotDir+"/singlePlot/nHitsPixel.png").c_str());
  c_nHitsPixel->SaveAs((basePlotDir+"/singlePlot/nHitsPixel.pdf").c_str());
  
  TCanvas *c_minPionPionDeltaPhi = new TCanvas("c_minPionPionDeltaPhi", "c_minPionPionDeltaPhi", 800, 800);
  h_minPionPionDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_minPionPionDeltaPhi[i]->Draw("hesame");
  stacks.at(h_minPionPionDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_minPionPionDeltaPhi, iPeriod, iPos );
  c_minPionPionDeltaPhi->SaveAs((basePlotDir+"/singlePlot/minPionPionDeltaPhi.png").c_str());
  c_minPionPionDeltaPhi->SaveAs((basePlotDir+"/singlePlot/minPionPionDeltaPhi.pdf").c_str());
  
  TCanvas *c_maxPionPionDeltaPhi = new TCanvas("c_maxPionPionDeltaPhi", "c_maxPionPionDeltaPhi", 800, 800);
  h_maxPionPionDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_maxPionPionDeltaPhi[i]->Draw("hesame");
  stacks.at(h_maxPionPionDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_maxPionPionDeltaPhi, iPeriod, iPos );
  c_maxPionPionDeltaPhi->SaveAs((basePlotDir+"/singlePlot/maxPionPionDeltaPhi.png").c_str());
  c_maxPionPionDeltaPhi->SaveAs((basePlotDir+"/singlePlot/maxPionPionDeltaPhi.pdf").c_str());
  
  TCanvas *c_minPionPionDeltaEta = new TCanvas("c_minPionPionDeltaEta", "c_minPionPionDeltaEta", 800, 800);
  h_minPionPionDeltaEta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_minPionPionDeltaEta[i]->Draw("hesame");
  stacks.at(h_minPionPionDeltaEta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_minPionPionDeltaEta, iPeriod, iPos );
  c_minPionPionDeltaEta->SaveAs((basePlotDir+"/singlePlot/minPionPionDeltaEta.png").c_str());
  c_minPionPionDeltaEta->SaveAs((basePlotDir+"/singlePlot/minPionPionDeltaEta.pdf").c_str());
  
  TCanvas *c_maxPionPionDeltaEta = new TCanvas("c_maxPionPionDeltaEta", "c_maxPionPionDeltaEta", 800, 800);
  h_maxPionPionDeltaEta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_maxPionPionDeltaEta[i]->Draw("hesame");
  stacks.at(h_maxPionPionDeltaEta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_maxPionPionDeltaEta, iPeriod, iPos );
  c_maxPionPionDeltaEta->SaveAs((basePlotDir+"/singlePlot/maxPionPionDeltaEta.png").c_str());
  c_maxPionPionDeltaEta->SaveAs((basePlotDir+"/singlePlot/maxPionPionDeltaEta.pdf").c_str());
  
  TCanvas *c_minPionPionDeltaR = new TCanvas("c_minPionPionDeltaR", "c_minPionPionDeltaR", 800, 800);
  h_minPionPionDeltaR[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_minPionPionDeltaR[i]->Draw("hesame");
  stacks.at(h_minPionPionDeltaR[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_minPionPionDeltaR, iPeriod, iPos );
  c_minPionPionDeltaR->SaveAs((basePlotDir+"/singlePlot/minPionPionDeltaR.png").c_str());
  c_minPionPionDeltaR->SaveAs((basePlotDir+"/singlePlot/minPionPionDeltaR.pdf").c_str());
  
  TCanvas *c_maxPionPionDeltaR = new TCanvas("c_maxPionPionDeltaR", "c_maxPionPionDeltaR", 800, 800);
  h_maxPionPionDeltaR[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_maxPionPionDeltaR[i]->Draw("hesame");
  stacks.at(h_maxPionPionDeltaR[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_maxPionPionDeltaR, iPeriod, iPos );
  c_maxPionPionDeltaR->SaveAs((basePlotDir+"/singlePlot/maxPionPionDeltaR.png").c_str());
  c_maxPionPionDeltaR->SaveAs((basePlotDir+"/singlePlot/maxPionPionDeltaR.pdf").c_str());
  
  TCanvas *c_minPionMuDeltaPhi = new TCanvas("c_minPionMuDeltaPhi", "c_minPionMuDeltaPhi", 800, 800);
  h_minPionMuDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_minPionMuDeltaPhi[i]->Draw("hesame");
  stacks.at(h_minPionMuDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_minPionMuDeltaPhi, iPeriod, iPos );
  c_minPionMuDeltaPhi->SaveAs((basePlotDir+"/singlePlot/minPionMuDeltaPhi.png").c_str());
  c_minPionMuDeltaPhi->SaveAs((basePlotDir+"/singlePlot/minPionMuDeltaPhi.pdf").c_str());
  
  TCanvas *c_maxPionMuDeltaPhi = new TCanvas("c_maxPionMuDeltaPhi", "c_maxPionMuDeltaPhi", 800, 800);
  h_maxPionMuDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_maxPionMuDeltaPhi[i]->Draw("hesame");
  stacks.at(h_maxPionMuDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_maxPionMuDeltaPhi, iPeriod, iPos );
  c_maxPionMuDeltaPhi->SaveAs((basePlotDir+"/singlePlot/maxPionMuDeltaPhi.png").c_str());
  c_maxPionMuDeltaPhi->SaveAs((basePlotDir+"/singlePlot/maxPionMuDeltaPhi.pdf").c_str());
  
  TCanvas *c_minPionMuDeltaEta = new TCanvas("c_minPionMuDeltaEta", "c_minPionMuDeltaEta", 800, 800);
  h_minPionMuDeltaEta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_minPionMuDeltaEta[i]->Draw("hesame");
  stacks.at(h_minPionMuDeltaEta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_minPionMuDeltaEta, iPeriod, iPos );
  c_minPionMuDeltaEta->SaveAs((basePlotDir+"/singlePlot/minPionMuDeltaEta.png").c_str());
  c_minPionMuDeltaEta->SaveAs((basePlotDir+"/singlePlot/minPionMuDeltaEta.pdf").c_str());
  
  TCanvas *c_maxPionMuDeltaEta = new TCanvas("c_maxPionMuDeltaEta", "c_maxPionMuDeltaEta", 800, 800);
  h_maxPionMuDeltaEta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_maxPionMuDeltaEta[i]->Draw("hesame");
  stacks.at(h_maxPionMuDeltaEta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_maxPionMuDeltaEta, iPeriod, iPos );
  c_maxPionMuDeltaEta->SaveAs((basePlotDir+"/singlePlot/maxPionMuDeltaEta.png").c_str());
  c_maxPionMuDeltaEta->SaveAs((basePlotDir+"/singlePlot/maxPionMuDeltaEta.pdf").c_str());
  
  TCanvas *c_minPionMuDeltaR = new TCanvas("c_minPionMuDeltaR", "c_minPionMuDeltaR", 800, 800);
  h_minPionMuDeltaR[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_minPionMuDeltaR[i]->Draw("hesame");
  stacks.at(h_minPionMuDeltaR[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_minPionMuDeltaR, iPeriod, iPos );
  c_minPionMuDeltaR->SaveAs((basePlotDir+"/singlePlot/minPionMuDeltaR.png").c_str());
  c_minPionMuDeltaR->SaveAs((basePlotDir+"/singlePlot/minPionMuDeltaR.pdf").c_str());
  
  TCanvas *c_maxPionMuDeltaR = new TCanvas("c_maxPionMuDeltaR", "c_maxPionMuDeltaR", 800, 800);
  h_maxPionMuDeltaR[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_maxPionMuDeltaR[i]->Draw("hesame");
  stacks.at(h_maxPionMuDeltaR[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_maxPionMuDeltaR, iPeriod, iPos );
  c_maxPionMuDeltaR->SaveAs((basePlotDir+"/singlePlot/maxPionMuDeltaR.png").c_str());
  c_maxPionMuDeltaR->SaveAs((basePlotDir+"/singlePlot/maxPionMuDeltaR.pdf").c_str());
  
  
  TCanvas *c_ABCD = new TCanvas("c_ABCD", "c_ABCD", 1600, 2400); c_ABCD->Divide(2,3);
  
  c_ABCD->cd(1); for (int i = 1; i < nSamples; i++){
    if (A_highNch_highHF[i]->GetMaximum() > A_highNch_highHF[0]->GetMaximum()) A_highNch_highHF[0]->SetMaximum(1.2*A_highNch_highHF[i]->GetMaximum());
    if (A_highNch_highHF[i]->GetMinimum() < A_highNch_highHF[0]->GetMinimum() && A_highNch_highHF[i]->GetMinimum() > 0) A_highNch_highHF[0]->SetMinimum(0.8*A_highNch_highHF[i]->GetMinimum());
  }
  A_highNch_highHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) A_highNch_highHF[i]->Draw("hesame");
  
  TCanvas *c_A_highNch_highHF = new TCanvas("c_A_highNch_highHF", "c_A_highNch_highHF", 800, 800);
  A_highNch_highHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) A_highNch_highHF[i]->Draw("hesame");
  c_A_highNch_highHF->SetLogy(1);
  CMS_lumi( c_A_highNch_highHF, iPeriod, iPos );
  c_A_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/A_highNch_highHF.png").c_str());
  c_A_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/A_highNch_highHF.pdf").c_str());
  
  c_ABCD->cd(2); for (int i = 1; i < nSamples; i++){
    if (B_lowNch_highHF[i]->GetMaximum() > B_lowNch_highHF[0]->GetMaximum()) B_lowNch_highHF[0]->SetMaximum(1.2*B_lowNch_highHF[i]->GetMaximum());
    if (B_lowNch_highHF[i]->GetMinimum() < B_lowNch_highHF[0]->GetMinimum() && B_lowNch_highHF[i]->GetMinimum() > 0) B_lowNch_highHF[0]->SetMinimum(0.8*B_lowNch_highHF[i]->GetMinimum());
  }
  B_lowNch_highHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) B_lowNch_highHF[i]->Draw("hesame");
  
  TCanvas *c_B_lowNch_highHF = new TCanvas("c_B_lowNch_highHF", "c_B_lowNch_highHF", 800, 800);
  B_lowNch_highHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) B_lowNch_highHF[i]->Draw("hesame");
  c_B_lowNch_highHF->SetLogy(1);
  CMS_lumi( c_B_lowNch_highHF, iPeriod, iPos );
  c_B_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/B_lowNch_highHF.png").c_str());
  c_B_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/B_lowNch_highHF.pdf").c_str());
  
  c_ABCD->cd(3); for (int i = 1; i < nSamples; i++){
    if (C_highNch_lowHF[i]->GetMaximum() > C_highNch_lowHF[0]->GetMaximum()) C_highNch_lowHF[0]->SetMaximum(1.2*C_highNch_lowHF[i]->GetMaximum());
    if (C_highNch_lowHF[i]->GetMinimum() < C_highNch_lowHF[0]->GetMinimum() && C_highNch_lowHF[i]->GetMinimum() > 0) C_highNch_lowHF[0]->SetMinimum(0.8*C_highNch_lowHF[i]->GetMinimum());
  }
  C_highNch_lowHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) C_highNch_lowHF[i]->Draw("hesame");
  
  TCanvas *c_C_highNch_lowHF = new TCanvas("c_C_highNch_lowHF", "c_C_highNch_lowHF", 800, 800);
  C_highNch_lowHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) C_highNch_lowHF[i]->Draw("hesame");
  c_C_highNch_lowHF->SetLogy(1);
  CMS_lumi( c_C_highNch_lowHF, iPeriod, iPos );
  c_C_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/C_highNch_lowHF.png").c_str());
  c_C_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/C_highNch_lowHF.pdf").c_str());
  
  c_ABCD->cd(4); for (int i = 1; i < nSamples; i++){
    if (D_lowNch_lowHF[i]->GetMaximum() > D_lowNch_lowHF[0]->GetMaximum()) D_lowNch_lowHF[0]->SetMaximum(1.2*D_lowNch_lowHF[i]->GetMaximum());
    if (D_lowNch_lowHF[i]->GetMinimum() < D_lowNch_lowHF[0]->GetMinimum() && D_lowNch_lowHF[i]->GetMinimum() > 0) D_lowNch_lowHF[0]->SetMinimum(0.8*D_lowNch_lowHF[i]->GetMinimum());
  }
  D_lowNch_lowHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) D_lowNch_lowHF[i]->Draw("hesame");
  
  TCanvas *c_D_lowNch_lowHF = new TCanvas("c_D_lowNch_lowHF", "c_D_lowNch_lowHF", 800, 800);
  D_lowNch_lowHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) D_lowNch_lowHF[i]->Draw("hesame");
  c_D_lowNch_lowHF->SetLogy(1);
  CMS_lumi( c_D_lowNch_lowHF, iPeriod, iPos );
  c_D_lowNch_lowHF->SaveAs((basePlotDir+"/singlePlot/D_lowNch_lowHF.png").c_str());
  c_D_lowNch_lowHF->SaveAs((basePlotDir+"/singlePlot/D_lowNch_lowHF.pdf").c_str());
  
  background->SetMarkerStyle(kFullTriangleDown);
  c_ABCD->cd(5); background->Draw("e0e1x0");
  
  TCanvas *c_background = new TCanvas("c_background", "c_background", 800, 800);
  background->Draw("e0e1x0");
  //c_background->SetLogy(1);
  CMS_lumi( c_background, iPeriod, iPos );
  c_background->SaveAs((basePlotDir+"/singlePlot/background.png").c_str());
  c_background->SaveAs((basePlotDir+"/singlePlot/background.pdf").c_str());
  
  
  c_ABCD->SaveAs((basePlotDir+"/collective/ABCD."+plotFormat).c_str());
  
  TCanvas *c_PV = new TCanvas("c_PV", "c_PV", 1000, 500); c_PV->Divide(2,1);
  c_PV->cd(1); h_deltaphi_tau_mu_tau_hadron_nch->Draw("COLZ");
  for (int i = 1; i < nSamples; i++){
    if (h_PV_N[i]->GetMaximum() > h_PV_N[0]->GetMaximum()) h_PV_N[0]->SetMaximum(1.2*h_PV_N[i]->GetMaximum());
  }
  h_PV_N[0]->GetYaxis()->SetRangeUser(0.1, 1.2*h_PV_N[0]->GetMaximum());
  c_PV->cd(2); h_PV_N[0]->Draw("e0e1x0");
  for (int i = 1; i < nSamples; i++) h_PV_N[i]->Draw("hesame"); legend->Draw();
  c_PV->cd(2)->SetLogy(1);
  c_PV->SaveAs((basePlotDir+"/collective/PV."+plotFormat).c_str());
  
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
  
  c_neutral_pion->cd(nGenPlots+1); h_recoGamma_pt[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_pt[i]->GetMaximum() > h_recoGamma_pt[0]->GetMaximum()) h_recoGamma_pt[0]->SetMaximum(1.2*h_recoGamma_pt[i]->GetMaximum());
    h_recoGamma_pt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+2); h_recoGamma_eta[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_eta[i]->GetMaximum() > h_recoGamma_eta[0]->GetMaximum()) h_recoGamma_eta[0]->SetMaximum(1.2*h_recoGamma_eta[i]->GetMaximum());
    if (h_recoGamma_eta[i]->GetMinimum() < h_recoGamma_eta[0]->GetMinimum()) h_recoGamma_eta[0]->SetMinimum(0.8*h_recoGamma_eta[i]->GetMinimum());
    h_recoGamma_eta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+3); h_recoGamma_phi[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_phi[i]->GetMaximum() > h_recoGamma_phi[0]->GetMaximum()) h_recoGamma_phi[0]->SetMaximum(1.2*h_recoGamma_phi[i]->GetMaximum());
    if (h_recoGamma_phi[i]->GetMinimum() < h_recoGamma_phi[0]->GetMinimum()) h_recoGamma_phi[0]->SetMinimum(0.8*h_recoGamma_phi[i]->GetMinimum());
    h_recoGamma_phi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+4); h_recoGamma_deltaphi_muon[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_deltaphi_muon[i]->GetMaximum() > h_recoGamma_deltaphi_muon[0]->GetMaximum()) h_recoGamma_deltaphi_muon[0]->SetMaximum(1.2*h_recoGamma_deltaphi_muon[i]->GetMaximum());
    if (h_recoGamma_deltaphi_muon[i]->GetMinimum() < h_recoGamma_deltaphi_muon[0]->GetMinimum() && h_recoGamma_deltaphi_muon[i]->GetMinimum() > 0) h_recoGamma_deltaphi_muon[0]->SetMinimum(0.8*h_recoGamma_deltaphi_muon[i]->GetMinimum());
    h_recoGamma_deltaphi_muon[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+4)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+5); h_recoGamma_deltaphi_pion[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGamma_deltaphi_pion[i]->GetMaximum() > h_recoGamma_deltaphi_pion[0]->GetMaximum()) h_recoGamma_deltaphi_pion[0]->SetMaximum(1.2*h_recoGamma_deltaphi_pion[i]->GetMaximum());
    if (h_recoGamma_deltaphi_pion[i]->GetMinimum() < h_recoGamma_deltaphi_pion[0]->GetMinimum() && h_recoGamma_deltaphi_pion[i]->GetMinimum() > 0) h_recoGamma_deltaphi_pion[0]->SetMinimum(0.8*h_recoGamma_deltaphi_pion[i]->GetMinimum());
    h_recoGamma_deltaphi_pion[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+5)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+6); h_NrecoGamma[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_NrecoGamma[i]->GetMaximum() > h_NrecoGamma[0]->GetMaximum()) h_NrecoGamma[0]->SetMaximum(1.2*h_NrecoGamma[i]->GetMaximum());
    h_NrecoGamma[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+6)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+7); h_recoGammasMinDeltaR[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoGammasMinDeltaR[i]->GetMaximum() > h_recoGammasMinDeltaR[0]->GetMaximum()) h_recoGammasMinDeltaR[0]->SetMaximum(1.2*h_recoGammasMinDeltaR[i]->GetMaximum());
    h_recoGammasMinDeltaR[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+8); h_recoPiZeroMinDeltaM[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZeroMinDeltaM[i]->GetMaximum() > h_recoPiZeroMinDeltaM[0]->GetMaximum()) h_recoPiZeroMinDeltaM[0]->SetMaximum(1.2*h_recoPiZeroMinDeltaM[i]->GetMaximum());
    h_recoPiZeroMinDeltaM[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+9); h_recoPiZeroDeltaM[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZeroDeltaM[i]->GetMaximum() > h_recoPiZeroDeltaM[0]->GetMaximum()) h_recoPiZeroDeltaM[0]->SetMaximum(1.2*h_recoPiZeroDeltaM[i]->GetMaximum());
    h_recoPiZeroDeltaM[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+10); h_recoPiZero_pt[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_pt[i]->GetMaximum() > h_recoPiZero_pt[0]->GetMaximum()) h_recoPiZero_pt[0]->SetMaximum(1.2*h_recoPiZero_pt[i]->GetMaximum());
    h_recoPiZero_pt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+11); h_recoPiZero_eta[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_eta[i]->GetMaximum() > h_recoPiZero_eta[0]->GetMaximum()) h_recoPiZero_eta[0]->SetMaximum(1.2*h_recoPiZero_eta[i]->GetMaximum());
    if (h_recoPiZero_eta[i]->GetMinimum() < h_recoPiZero_eta[0]->GetMinimum()) h_recoPiZero_eta[0]->SetMinimum(0.8*h_recoPiZero_eta[i]->GetMinimum());
    h_recoPiZero_eta[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+12); h_recoPiZero_phi[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_phi[i]->GetMaximum() > h_recoPiZero_phi[0]->GetMaximum()) h_recoPiZero_phi[0]->SetMaximum(1.2*h_recoPiZero_phi[i]->GetMaximum());
    if (h_recoPiZero_phi[i]->GetMinimum() < h_recoPiZero_phi[0]->GetMinimum()) h_recoPiZero_phi[0]->SetMinimum(0.8*h_recoPiZero_phi[i]->GetMinimum());
    h_recoPiZero_phi[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_neutral_pion->cd(nGenPlots+13); h_recoPiZero_deltaphi_muon[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_deltaphi_muon[i]->GetMaximum() > h_recoPiZero_deltaphi_muon[0]->GetMaximum()) h_recoPiZero_deltaphi_muon[0]->SetMaximum(1.2*h_recoPiZero_deltaphi_muon[i]->GetMaximum());
    if (h_recoPiZero_deltaphi_muon[i]->GetMinimum() < h_recoPiZero_deltaphi_muon[0]->GetMinimum() && h_recoPiZero_deltaphi_muon[i]->GetMinimum() > 0) h_recoPiZero_deltaphi_muon[0]->SetMinimum(0.8*h_recoPiZero_deltaphi_muon[i]->GetMinimum());
    h_recoPiZero_deltaphi_muon[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+13)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+14); h_recoPiZero_deltaphi_pion[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_recoPiZero_deltaphi_pion[i]->GetMaximum() > h_recoPiZero_deltaphi_pion[0]->GetMaximum()) h_recoPiZero_deltaphi_pion[0]->SetMaximum(1.2*h_recoPiZero_deltaphi_pion[i]->GetMaximum());
    if (h_recoPiZero_deltaphi_pion[i]->GetMinimum() < h_recoPiZero_deltaphi_pion[0]->GetMinimum() && h_recoPiZero_deltaphi_pion[i]->GetMinimum() > 0) h_recoPiZero_deltaphi_pion[0]->SetMinimum(0.8*h_recoPiZero_deltaphi_pion[i]->GetMinimum());
    h_recoPiZero_deltaphi_pion[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+14)->SetLogy(1);
  
  c_neutral_pion->cd(nGenPlots+15); h_NrecoPiZero[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++){
    if (h_NrecoPiZero[i]->GetMaximum() > h_NrecoPiZero[0]->GetMaximum()) h_NrecoPiZero[0]->SetMaximum(1.2*h_NrecoPiZero[i]->GetMaximum());
    h_NrecoPiZero[i]->Draw("hesame");
  }
  legend->Draw();
  c_neutral_pion->cd(nGenPlots+15)->SetLogy(1);
  
  c_neutral_pion->SaveAs((basePlotDir+"/collective/neutral_pion."+plotFormat).c_str());


  TCanvas *c_efficiency = new TCanvas("c_efficiency", "c_efficiency", 1500, 1500+nSamples*500); c_efficiency->Divide(3,3+nSamples);
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
  //nobackgend->Draw();
  
  
  TCanvas *c_tau_eff = new TCanvas("c_tau_eff", "c_tau_eff", 800, 800);
  for (int i = 1; i < nSamples; i++) h_eff_tauReco_genTauHadpt[i]->Draw("hist same");
  //nobackgend->Draw();
  //CMS_lumi( c_background, iPeriod, iPos );
  c_tau_eff->SaveAs((basePlotDir+"/singlePlot/eff_tau.png").c_str());
  c_tau_eff->SaveAs((basePlotDir+"/singlePlot/eff_tau.pdf").c_str());
  
  c_efficiency->cd(4); for (int i = 1; i < nSamples; i++){
    if (h_N_genPionPt[i]->GetMaximum() > h_N_genPionPt[1]->GetMaximum()) h_N_genPionPt[1]->SetMaximum(1.2*h_N_genPionPt[i]->GetMaximum());
    h_N_genPionPt[i]->Draw("hesame");
  }
  legend->Draw();
  
  c_efficiency->cd(5); for (int i = 1; i < nSamples; i++){
    h_eff_matchedPionReco_genPionPt[i]->Divide(h_N_genPionPt[i]);
    //if (h_eff_matchedPionReco_genPionPt[i]->GetMaximum() > h_eff_matchedPionReco_genPionPt[1]->GetMaximum()) h_eff_matchedPionReco_genPionPt[1]->SetMaximum(1.2*h_eff_matchedPionReco_genPionPt[i]->GetMaximum());
    h_eff_matchedPionReco_genPionPt[1]->SetMaximum(1.1);
    h_eff_matchedPionReco_genPionPt[i]->Draw("hist same");
  }
  //legend->Draw();
  
  
  // setup the pads
  float yplot = 0.72;
  float yratio = 0.26;
  float yspace = 0.00;
  
  float padTop_scale = 1./yplot;
  float padRatio_scale = 1./yratio;
  
  // set up the coordinates of the two pads:    //
  float y1, y2, y3, y4;                         //  y4 +-------------+
  y4 = 0.99;                                    //     |             |
  y3 = y4 - yplot;                              //     |     pad1    |
  y2 = y3 - yspace;                             //  y3 |-------------|
  y1 = y2 - yratio;                             //  y2 |-------------|
  float x1, x2;                                 //     |     rp1     |
  x1 = 0.01;                                    //  y1 +-------------+
  x2 = 0.99;                                    //     x1            x2
  
  if (nSamples>1){
      
    h_recoPionPt[1]->SetLineColor(colors[2]);
    int nbins = h_recoPionPt[1]->GetNbinsX();
    TH1F *closureRatio[2];
    closureRatio[0] = (TH1F*)h_recoPionPt[1]->Clone("numerator reco");
    closureRatio[1] = (TH1F*)h_N_matched_genPionPt[1]->Clone("denominator gen");
    
    float val_n, err_n, val_d, err_d;
    for (int b = 1; b <= nbins; b++){
      if (closureRatio[1]->GetBinContent(b) == 0){
        val_n = 0;
        err_n = 0;
        val_d = 1;
        err_d = 0;
      } else{
        val_n = h_recoPionPt[1]->GetBinContent(b) / closureRatio[1]->GetBinContent(b);
        err_n = h_recoPionPt[1]->GetBinError(b) / closureRatio[1]->GetBinContent(b);
        val_d = 1;
        err_d = closureRatio[1]->GetBinError(b) / closureRatio[1]->GetBinContent(b);
      }
      closureRatio[0]->SetBinContent(b,val_n);
      closureRatio[0]->SetBinError(b,err_n);
      closureRatio[1]->SetBinContent(b,val_d);
      closureRatio[1]->SetBinError(b,err_d);
    }
    
    TCanvas *c_MCclosure_pion_pt = new TCanvas("c_MCclosure_pion_pt", "c_MCclosure_pion_pt", 800, 800);  
    
    TLegend *closureLegend = new TLegend(0.7,0.77,0.9,0.90);
    closureLegend->SetFillStyle(0);
    closureLegend->SetFillColor(0); closureLegend->SetLineColor(0); closureLegend->SetShadowColor(0); closureLegend->SetTextSize(0.03);
    closureLegend->AddEntry(h_recoPionPt[1],"weighted reco", "l");
    closureLegend->AddEntry(h_N_matched_genPionPt[1],"gen", "l");
    
    TPad* pad_ratio_closure = new TPad("rp1 closure", "Ratio", x1, y1, x2, y2);
    pad_ratio_closure->SetTopMargin(0.035);
    pad_ratio_closure->SetBottomMargin(0.5);
    pad_ratio_closure->SetLeftMargin(0.14);
    pad_ratio_closure->SetRightMargin(0.025);
    TPad* pad_top_closure = new TPad("pad1 closure", "Control Plots", x1, y3, x2, y4);
    pad_top_closure->SetTopMargin(0.07);
    pad_top_closure->SetBottomMargin(0.03);
    pad_top_closure->SetLeftMargin(0.14);
    pad_top_closure->SetRightMargin(0.025);
    pad_ratio_closure->Draw();
    pad_top_closure->Draw();
    pad_top_closure->cd();
    h_recoPionPt[1]->GetXaxis()->SetLabelSize(0.05*padTop_scale);
    h_recoPionPt[1]->GetYaxis()->SetLabelSize(0.05*padTop_scale);
    h_recoPionPt[1]->GetXaxis()->SetLabelOffset(0.017*padTop_scale);
    h_recoPionPt[1]->GetYaxis()->SetLabelOffset(0.007*padTop_scale);
    h_recoPionPt[1]->GetXaxis()->SetTitleSize(0.06*padTop_scale);
    h_recoPionPt[1]->GetYaxis()->SetTitleSize(0.06*padTop_scale);
    h_recoPionPt[1]->SetMaximum(1.2*h_recoPionPt[1]->GetMaximum());
    h_recoPionPt[1]->Draw("hesame");
    h_N_matched_genPionPt[1]->Draw("hesame");
    closureLegend->Draw();
    
    pad_ratio_closure->Clear();
    pad_ratio_closure->cd();
    closureRatio[0]->GetXaxis()->SetLabelSize(0.04*padRatio_scale);
    closureRatio[0]->GetYaxis()->SetLabelSize(0.04*padRatio_scale);
    closureRatio[0]->GetXaxis()->SetLabelOffset(0.007*padRatio_scale);
    closureRatio[0]->GetYaxis()->SetLabelOffset(0.007*padRatio_scale);
    closureRatio[0]->GetXaxis()->SetTitleSize(0.05*padRatio_scale);
    closureRatio[0]->GetYaxis()->SetTitleSize(0.04*padRatio_scale);
    closureRatio[0]->GetYaxis()->SetTitleOffset(0.1*padRatio_scale);
    closureRatio[0]->GetYaxis()->SetNdivisions(504);
    
    closureRatio[0]->GetYaxis()->SetRangeUser(0., 1.2*closureRatio[0]->GetMaximum());
    //closureRatio[1]->SetLineColor(kBlack);
    closureRatio[1]->SetMarkerSize(0); //closureRatio[1]->SetMarkerStyle(1);
    closureRatio[1]->SetFillColor(colors[1]);
    closureRatio[1]->SetFillStyle(3013);
    closureRatio[0]->GetYaxis()->SetTitle("#frac{weighted reco}{gen}");
    closureRatio[0]->Draw("hesame");
    closureRatio[1]->Draw("e2same");
    TLine *closureLine = new TLine(closureRatio[0]->GetXaxis()->GetXmin(), 1, closureRatio[0]->GetXaxis()->GetXmax(), 1);
    closureLine->SetLineWidth(2);
    closureLine->SetLineColor(kBlack);
    closureLine->Draw("same");
    gPad->Modified(); gPad->Update();
    
    CMS_lumi( c_MCclosure_pion_pt, iPeriod, iPos );
    gStyle->SetErrorX(0.5);
    c_MCclosure_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_pion_pt.png").c_str());
    c_MCclosure_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_pion_pt.pdf").c_str());
    gStyle->SetErrorX(0.);
    
    
    // **** **** **** ****
    
    TCanvas *c_MCclosure_leading_pion_pt = new TCanvas("c_MCclosure_leading_pion_pt", "c_MCclosure_leading_pion_pt", 800, 800);
    h_leading_recoPionPt[1]->SetLineColor(colors[2]);
    h_leading_recoPionPt[1]->SetMaximum(1.2*h_leading_recoPionPt[1]->GetMaximum());
    h_leading_recoPionPt[1]->Draw("hesame");
    h_N_matched_leading_genPionPt[1]->Draw("hesame");
    closureLegend->Draw();
    CMS_lumi( c_MCclosure_leading_pion_pt, iPeriod, iPos );
    c_MCclosure_leading_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_leading_pion_pt.png").c_str());
    c_MCclosure_leading_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_leading_pion_pt.pdf").c_str());
  
    TCanvas *c_MCclosure_subleading_pion_pt = new TCanvas("c_MCclosure_subleading_pion_pt", "c_MCclosure_subleading_pion_pt", 800, 800);
    h_subleading_recoPionPt[1]->SetLineColor(colors[2]);
    h_subleading_recoPionPt[1]->SetMaximum(1.2*h_subleading_recoPionPt[1]->GetMaximum());
    h_subleading_recoPionPt[1]->Draw("hesame");
    h_N_matched_subleading_genPionPt[1]->Draw("hesame");
    closureLegend->Draw();
    CMS_lumi( c_MCclosure_subleading_pion_pt, iPeriod, iPos );
    c_MCclosure_subleading_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_subleading_pion_pt.png").c_str());
    c_MCclosure_subleading_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_subleading_pion_pt.pdf").c_str());
  
    TCanvas *c_MCclosure_subsubleading_pion_pt = new TCanvas("c_MCclosure_subsubleading_pion_pt", "c_MCclosure_subsubleading_pion_pt", 800, 800);
    h_subsubleading_recoPionPt[1]->SetLineColor(colors[2]);
    h_subsubleading_recoPionPt[1]->SetMaximum(1.2*h_subsubleading_recoPionPt[1]->GetMaximum());
    h_subsubleading_recoPionPt[1]->Draw("hesame");
    h_N_matched_subsubleading_genPionPt[1]->Draw("hesame");
    closureLegend->Draw();
    CMS_lumi( c_MCclosure_subsubleading_pion_pt, iPeriod, iPos );
    c_MCclosure_subsubleading_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_subsubleading_pion_pt.png").c_str());
    c_MCclosure_subsubleading_pion_pt->SaveAs((basePlotDir+"/singlePlot/closure_subsubleading_pion_pt.pdf").c_str());
  } // if nSamples == 1
  
  
  TCanvas *c_pion_efficiency_pt = new TCanvas("c_pion_efficiency_pt", "c_pion_efficiency_pt", 800, 800);
  for (int i = 1; i < nSamples; i++) h_eff_matchedPionReco_genPionPt[i]->Draw("he same");
  //CMS_lumi( c_pion_efficiency_pt, iPeriod, iPos );
  c_pion_efficiency_pt->SaveAs((basePlotDir+"/singlePlot/pion_efficiency_pt.png").c_str());
  c_pion_efficiency_pt->SaveAs((basePlotDir+"/singlePlot/pion_efficiency_pt.pdf").c_str());
  
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
    //if (h_eff_matchedPionReco_genPionEta[i]->GetMaximum() > h_eff_matchedPionReco_genPionEta[1]->GetMaximum()) h_eff_matchedPionReco_genPionEta[1]->SetMaximum(1.2*h_eff_matchedPionReco_genPionEta[i]->GetMaximum());
    h_eff_matchedPionReco_genPionEta[1]->SetMaximum(1.1);
    h_eff_matchedPionReco_genPionEta[i]->Draw("hist same");
  }
  //legend->Draw();
  
  TCanvas *c_pion_efficiency_eta = new TCanvas("c_pion_efficiency_eta", "c_pion_efficiency_eta", 800, 800);
  for (int i = 1; i < nSamples; i++) h_eff_matchedPionReco_genPionEta[i]->Draw("he same");
  CMS_lumi( c_pion_efficiency_eta, iPeriod, iPos );
  c_pion_efficiency_eta->SaveAs((basePlotDir+"/singlePlot/pion_efficiency_eta.png").c_str());
  c_pion_efficiency_eta->SaveAs((basePlotDir+"/singlePlot/pion_efficiency_eta.pdf").c_str());
  
  c_efficiency->cd(9); for (int i = 1; i < nSamples; i++){
    h_eff_matchedHighPtPionReco_genPionEta[i]->Divide(h_N_genHighPtPionEta[i]);
    if (h_eff_matchedHighPtPionReco_genPionEta[i]->GetMaximum() > h_eff_matchedHighPtPionReco_genPionEta[1]->GetMaximum()) h_eff_matchedHighPtPionReco_genPionEta[1]->SetMaximum(1.2*h_eff_matchedHighPtPionReco_genPionEta[i]->GetMaximum());
    h_eff_matchedHighPtPionReco_genPionEta[i]->Draw("hist same");
  }
  //legend->Draw();
  
  c_efficiency->cd(10); for (int i = 1; i < nSamples; i++){
    if (h_N_recoMuonPt[i]->GetMaximum() > h_N_recoMuonPt[1]->GetMaximum()) h_N_recoMuonPt[1]->SetMaximum(1.2*h_N_recoMuonPt[i]->GetMaximum());
    h_N_recoMuonPt[i]->Draw("hesame");
  }
  //legend->Draw();
  
  c_efficiency->cd(11); for (int i = 1; i < nSamples; i++){
    h_eff_trigger[i]->Divide(h_N_recoMuonPt[i]);
    if (h_eff_trigger[i]->GetMaximum() > h_eff_trigger[1]->GetMaximum()) h_eff_trigger[1]->SetMaximum(1.2*h_eff_trigger[i]->GetMaximum());
    h_eff_trigger[i]->Draw("hesame");
  }
  //legend->Draw();
  
  TCanvas *c_trigger_efficiency = new TCanvas("c_trigger_efficiency", "c_trigger_efficiency", 800, 800);
  for (int i = 1; i < nSamples; i++) h_eff_trigger[i]->Draw("hist same");
  CMS_lumi( c_trigger_efficiency, iPeriod, iPos );
  c_trigger_efficiency->SaveAs((basePlotDir+"/singlePlot/trigger_efficiency.png").c_str());
  c_trigger_efficiency->SaveAs((basePlotDir+"/singlePlot/trigger_efficiency.pdf").c_str());
  
  for (int i = 1; i < nSamples; i++){
    c_efficiency->cd(10+3*i); 
    h_N_genPionPtEta[i]->Draw("colz");
    c_efficiency->cd(11+3*i); 
    h_eff_matchedPionReco_genPionPtEta[i]->Divide(h_N_genPionPtEta[i]);
    h_eff_matchedPionReco_genPionPtEta[i]->Draw("colz");
    c_efficiency->cd(12+3*i);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->Divide(h_N_genPionPtEta[i]);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->GetZaxis()->SetRangeUser(0.9, 1.1);
    h_eff_matchedPionReco_genPionPtEtaZoomed[i]->Draw("colz");
  }
  
  c_efficiency->SaveAs((basePlotDir+"/collective/efficiency."+plotFormat).c_str());


  TCanvas *c_tau_muon_kinematics = new TCanvas("c_tau_muon_kinematics", "c_tau_muon_kinematics", 1500, 1000); c_tau_muon_kinematics->Divide(3,2);
  
  c_tau_muon_kinematics->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_tau_mu_p[i]->GetMaximum() > h_tau_mu_p[0]->GetMaximum()) h_tau_mu_p[0]->SetMaximum(1.2*h_tau_mu_p[i]->GetMaximum());
    if (h_tau_mu_p[i]->GetMinimum() < h_tau_mu_p[0]->GetMinimum() && h_tau_mu_p[i]->GetMinimum() > 0) h_tau_mu_p[0]->SetMinimum(0.8*h_tau_mu_p[i]->GetMinimum());
  }
  h_tau_mu_p[0]->Draw("e0e1x0"); legend->Draw();
  c_tau_muon_kinematics->cd(1)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_mu_p[i]->Draw("hesame");
  
  c_tau_muon_kinematics->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_tau_mu_pz[i]->GetMaximum() > h_tau_mu_pz[0]->GetMaximum()) h_tau_mu_pz[0]->SetMaximum(1.2*h_tau_mu_pz[i]->GetMaximum());
    if (h_tau_mu_pz[i]->GetMinimum() < h_tau_mu_pz[0]->GetMinimum() && h_tau_mu_pz[i]->GetMinimum() > 0) h_tau_mu_pz[0]->SetMinimum(0.8*h_tau_mu_pz[i]->GetMinimum());
  }
  h_tau_mu_pz[0]->Draw("e0e1x0"); legend->Draw();
  c_tau_muon_kinematics->cd(2)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_mu_pz[i]->Draw("hesame");
  
  TCanvas *c_tau_muon_pz = new TCanvas("c_tau_muon_pz", "c_tau_muon_pz", 800, 800);
  h_tau_mu_pz[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_mu_pz[i]->Draw("hesame");
  stacks.at(h_tau_mu_pz[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_pz, iPeriod, iPos );
  c_tau_muon_pz->SaveAs((basePlotDir+"/singlePlot/tau_muon_pz.png").c_str());
  c_tau_muon_pz->SaveAs((basePlotDir+"/singlePlot/tau_muon_pz.pdf").c_str());
  
  
  TCanvas *c_tau_muon_dz = new TCanvas("c_tau_muon_dz", "c_tau_muon_dz", 800, 800);
  h_tau_mu_dz[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_mu_dz[i]->Draw("hesame");
  stacks.at(h_tau_mu_dz[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_dz, iPeriod, iPos );
  c_tau_muon_dz->SaveAs((basePlotDir+"/singlePlot/tau_muon_dz.png").c_str());
  c_tau_muon_dz->SaveAs((basePlotDir+"/singlePlot/tau_muon_dz.pdf").c_str());
  
  
  c_tau_muon_kinematics->cd(3); for (int i = 1; i < nSamples; i++){
    hs_tau_mu_pt->Add(h_tau_mu_pt[i]);
    if (h_tau_mu_pt[i]->GetMaximum() > h_tau_mu_pt[0]->GetMaximum()) h_tau_mu_pt[0]->SetMaximum(1.2*h_tau_mu_pt[i]->GetMaximum());
    if (h_tau_mu_pt[i]->GetMinimum() < h_tau_mu_pt[0]->GetMinimum() && h_tau_mu_pt[i]->GetMinimum() > 0) h_tau_mu_pt[0]->SetMinimum(0.8*h_tau_mu_pt[i]->GetMinimum());
  }
  h_tau_mu_pt[0]->Draw("e0e1x0"); legend->Draw();
  c_tau_muon_kinematics->cd(3)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_mu_pt[i]->Draw("hesame");
  
  TCanvas *c_tau_muon_pt = new TCanvas("c_tau_muon_pt", "c_tau_muon_pt", 800, 800);
  h_tau_mu_pt[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_mu_pt[i]->Draw("hesame");
  //h_tau_mu_pt[1]->SetMaximum(30);
  //tempGend->Draw();
  //h_tau_mu_pt[nSamples+3]->Draw("hesame");
  stacks.at(h_tau_mu_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_pt, iPeriod, iPos );
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt.png").c_str());
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt.pdf").c_str());
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt.C").c_str());
  
  bool draw_atau_weights = true;
  if (draw_atau_weights){
    TH1F *r_tau_mu_pt[3];
    int atauList[] = {-10,0,10};
    string baseNameWeights = "r_tau_mu_pt_atau";
    for (int i = 0; i < 3; i++){
      r_tau_mu_pt[i] = (TH1F*)h_tau_mu_pt[i+1]->Clone((baseNameWeights+to_string(atauList[i])+"E-2").c_str());
      r_tau_mu_pt[i]->SetTitle((baseNameWeights+to_string(atauList[i])+"E-2").c_str());
      r_tau_mu_pt[i]->Divide(h_tau_mu_pt[2]);
    }
    TCanvas *c_tau_muon_atau_weights = new TCanvas("c_tau_muon_atau_weights", "c_tau_muon_atau_weights", 800, 800);
    for (int i = 0; i < 3; i++) r_tau_mu_pt[i]->Draw("hesame");
    r_tau_mu_pt[0]->SetMaximum(12);
    tempGend->Draw();
    c_tau_muon_atau_weights->SaveAs((basePlotDir+"/singlePlot/tau_muon_atau_weights.png").c_str());
    c_tau_muon_atau_weights->SaveAs((basePlotDir+"/singlePlot/tau_muon_atau_weights.pdf").c_str());
    c_tau_muon_atau_weights->SaveAs((basePlotDir+"/singlePlot/tau_muon_atau_weights.C").c_str());
  }
  
  if (true){
    int nbins_muon_SF_effect = h_muon_pt_after_SF->GetNbinsX();
    TH1F *tempRatio[2];
    tempRatio[0] = (TH1F*)h_muon_pt_after_SF->Clone("numerator");
    tempRatio[1] = (TH1F*)h_muon_pt_before_SF->Clone("denominator");
    
    float val_n, err_n, val_d, err_d;
    for (int b = 1; b <= nbins_muon_SF_effect; b++){
      if (tempRatio[1]->GetBinContent(b) == 0){
        val_n = 0;
        err_n = 0;
        val_d = 1;
        err_d = 0;
      } else{
        val_n = tempRatio[0]->GetBinContent(b) / tempRatio[1]->GetBinContent(b);
        err_n = tempRatio[0]->GetBinError(b) / tempRatio[1]->GetBinContent(b);
        val_d = 1;
        err_d = tempRatio[1]->GetBinError(b) / tempRatio[1]->GetBinContent(b);
      }
      tempRatio[0]->SetBinContent(b,val_n);
      tempRatio[0]->SetBinError(b,err_n);
      tempRatio[1]->SetBinContent(b,val_d);
      tempRatio[1]->SetBinError(b,err_d);
    }
    
    TCanvas *c_muon_SF_effect = new TCanvas("c_muon_SF_effect", "c_muon_SF_effect", 800, 900);
    
    // setup the pads
    float yplot = 0.72;
    float yratio = 0.26;
    float yspace = 0.00;
    
    float padTop_scale = 1./yplot;
    float padRatio_scale = 1./yratio;
    
    // set up the coordinates of the two pads:    //
    float y1, y2, y3, y4;                         //  y4 +-------------+
    y4 = 0.99;                                    //     |             |
    y3 = y4 - yplot;                              //     |     pad1    |
    y2 = y3 - yspace;                             //  y3 |-------------|
    y1 = y2 - yratio;                             //  y2 |-------------|
    float x1, x2;                                 //     |     rp1     |
    x1 = 0.01;                                    //  y1 +-------------+
    x2 = 0.99;                                    //     x1            x2
  
    TPad* pad_ratio = new TPad("rp1", "Ratio", x1, y1, x2, y2);
  
    pad_ratio->SetTopMargin(0.035);
    pad_ratio->SetBottomMargin(0.5);
    pad_ratio->SetLeftMargin(0.14);
    pad_ratio->SetRightMargin(0.025);
    
    TPad* pad_top = new TPad("pad1", "Control Plots", x1, y3, x2, y4);
  
    pad_top->SetTopMargin(0.07);
    pad_top->SetBottomMargin(0.03);
    pad_top->SetLeftMargin(0.14);
    pad_top->SetRightMargin(0.025);
    
    pad_ratio->Draw();
    pad_top->Draw();
    
    pad_top->cd();
    h_muon_pt_before_SF->GetXaxis()->SetLabelSize(0.05*padTop_scale);
    h_muon_pt_before_SF->GetYaxis()->SetLabelSize(0.05*padTop_scale);
    h_muon_pt_before_SF->GetXaxis()->SetLabelOffset(0.017*padTop_scale);
    h_muon_pt_before_SF->GetYaxis()->SetLabelOffset(0.007*padTop_scale);
    h_muon_pt_before_SF->GetXaxis()->SetTitleSize(0.06*padTop_scale);
    h_muon_pt_before_SF->GetYaxis()->SetTitleSize(0.06*padTop_scale);
    h_muon_pt_before_SF->Draw("hesame");
    h_muon_pt_after_SF->Draw("hesame");
    TLegend *legend_muonSF_effect = new TLegend(0.65,0.7,0.8,0.85);
    legend_muonSF_effect->SetFillStyle(0);
    legend_muonSF_effect->SetFillColor(0); legend_muonSF_effect->SetLineColor(0); legend_muonSF_effect->SetShadowColor(0); legend_muonSF_effect->SetTextSize(0.035);
    legend_muonSF_effect->AddEntry(h_muon_pt_before_SF,    " without muon SF", "l");
    legend_muonSF_effect->AddEntry(h_muon_pt_after_SF,    " with muon SF", "l");
    legend_muonSF_effect->Draw();
    
    pad_ratio->Clear();
    pad_ratio->cd();
    tempRatio[0]->GetXaxis()->SetLabelSize(0.04*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetLabelSize(0.04*padRatio_scale);
    tempRatio[0]->GetXaxis()->SetLabelOffset(0.007*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetLabelOffset(0.007*padRatio_scale);
    tempRatio[0]->GetXaxis()->SetTitleSize(0.05*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetTitleSize(0.04*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetTitleOffset(0.08*padRatio_scale);
    tempRatio[0]->GetYaxis()->SetNdivisions(504);
    
    tempRatio[0]->GetYaxis()->SetRangeUser(-0.05, 2.05);
    //tempRatio[0]->Draw("AXIS");
    tempRatio[1]->SetLineColor(kBlack); tempRatio[1]->SetMarkerSize(0); //tempRatio[1]->SetMarkerStyle(1);
    tempRatio[1]->SetFillColor(kBlue); tempRatio[1]->SetFillStyle(3013);
    tempRatio[0]->GetYaxis()->SetTitle("after/before SF");
    tempRatio[0]->Draw("e0e1x0");
    tempRatio[1]->Draw("e2same");
    TLine *line = new TLine(tempRatio[0]->GetXaxis()->GetXmin(), 1, tempRatio[0]->GetXaxis()->GetXmax(), 1);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    gPad->Modified(); gPad->Update();
    
    CMS_lumi( c_muon_SF_effect, iPeriod, iPos );
    c_muon_SF_effect->SaveAs((basePlotDir+"/singlePlot/muon_SF_effect.png").c_str());
    c_muon_SF_effect->SaveAs((basePlotDir+"/singlePlot/muon_SF_effect.pdf").c_str());
    gStyle->SetErrorX(0.5);
    
  }//if true
  
  
  
  
  
  
  
  
  
  
  c_tau_muon_kinematics->cd(4); for (int i = 1; i < nSamples; i++){
    hs_tau_mu_eta->Add(h_tau_mu_eta[i]);
    if (h_tau_mu_eta[i]->GetMaximum() > h_tau_mu_eta[0]->GetMaximum()) h_tau_mu_eta[0]->SetMaximum(1.2*h_tau_mu_eta[i]->GetMaximum());
  }
  h_tau_mu_eta[0]->SetMaximum(1.4*h_tau_mu_eta[0]->GetMaximum());
  h_tau_mu_eta[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_eta[i]->Draw("hesame");
  
  
  TCanvas *c_tau_muon_eta = new TCanvas("c_tau_muon_eta", "c_tau_muon_eta", 800, 800);
  h_tau_mu_eta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_mu_eta[i]->Draw("hesame");
  stacks.at(h_tau_mu_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_eta, iPeriod, iPos );
  c_tau_muon_eta->SaveAs((basePlotDir+"/singlePlot/tau_muon_eta.png").c_str());
  c_tau_muon_eta->SaveAs((basePlotDir+"/singlePlot/tau_muon_eta.pdf").c_str());
  
  c_tau_muon_kinematics->cd(5); for (int i = 1; i < nSamples; i++){
    hs_tau_mu_phi->Add(h_tau_mu_phi[i]);
    if (h_tau_mu_phi[i]->GetMaximum() > h_tau_mu_phi[0]->GetMaximum()) h_tau_mu_phi[0]->SetMaximum(1.2*h_tau_mu_phi[i]->GetMaximum());
  }
  h_tau_mu_phi[0]->SetMaximum(1.9*h_tau_mu_phi[0]->GetMaximum());
  h_tau_mu_phi[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_phi[i]->Draw("hesame");
  
  TCanvas *c_tau_muon_phi = new TCanvas("c_tau_muon_phi", "c_tau_muon_phi", 800, 800);
  h_tau_mu_phi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_mu_phi[i]->Draw("hesame");
  stacks.at(h_tau_mu_phi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  
  TH1F *stackedtau_mu_phi = new TH1F(*(TH1F*)((stacks.at(h_tau_mu_phi[nSamples+1]->GetBinContent(1))->GetStack()->Last())));
  cout << " **** **** **** ****" << endl << "KS test for tau muon phi: " << h_tau_mu_phi[0]->KolmogorovTest(stackedtau_mu_phi) << endl << endl;
  
  CMS_lumi( c_tau_muon_phi, iPeriod, iPos );
  c_tau_muon_phi->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi.png").c_str());
  c_tau_muon_phi->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi.pdf").c_str());
  
  TCanvas *c_tau_muon_phi_lowNch_highHF = new TCanvas("c_tau_muon_phi_lowNch_highHF", "c_tau_muon_phi_lowNch_highHF", 800, 800);
  h_tau_mu_phi[nSamples+1]->Draw("e0e1x0");
  
  CMS_lumi( c_tau_muon_phi_lowNch_highHF, iPeriod, iPos );
  c_tau_muon_phi_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_lowNch_highHF.png").c_str());
  c_tau_muon_phi_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_lowNch_highHF.pdf").c_str());
  c_tau_muon_phi_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_lowNch_highHF.root").c_str());
  
  TCanvas *c_tau_muon_phi_highNch_lowHF = new TCanvas("c_tau_muon_phi_highNch_lowHF", "c_tau_muon_phi_highNch_lowHF", 800, 800);
  h_tau_mu_phi[nSamples+2]->Draw("e0e1x0");
  
  CMS_lumi( c_tau_muon_phi_highNch_lowHF, iPeriod, iPos );
  c_tau_muon_phi_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_highNch_lowHF.png").c_str());
  c_tau_muon_phi_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_highNch_lowHF.pdf").c_str());
  c_tau_muon_phi_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_highNch_lowHF.root").c_str());
  
  TCanvas *c_tau_muon_phi_highNch_highHF = new TCanvas("c_tau_muon_phi_highNch_highHF", "c_tau_muon_phi_highNch_highHF", 800, 800);
  h_tau_mu_phi[nSamples]->Draw("e0e1x0");// legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_mu_phi[i]->Draw("hesame");
  
  CMS_lumi( c_tau_muon_phi_highNch_highHF, iPeriod, iPos );
  c_tau_muon_phi_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_highNch_highHF.png").c_str());
  c_tau_muon_phi_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_highNch_highHF.pdf").c_str());
  c_tau_muon_phi_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_muon_phi_highNch_highHF.root").c_str());
  
  c_tau_muon_kinematics->SaveAs((basePlotDir+"/collective/tau_muon_kinematics."+plotFormat).c_str());

  TCanvas *c_tau_muon_tau_hadron_phi_data = new TCanvas("c_tau_muon_tau_hadron_phi_data", "c_tau_muon_tau_hadron_phi_data", 800, 800);
  h_tau_mu_tau_hadron_phi[0]->Draw("colz");
  float upD = 0;
  float downD = 0;
  for (int xbin = 1; xbin <= tau_mu_phi_bins; xbin++){
    for (int ybin = 1; ybin <= tau_hadron_phi_bins; ybin++){
      if(ybin>xbin) upD += h_tau_mu_tau_hadron_phi[0]->GetBinContent(xbin,ybin);
      if(ybin<xbin) downD += h_tau_mu_tau_hadron_phi[0]->GetBinContent(xbin,ybin);
    }
  }
  cout << "tau_muon_tau_hadron_phi_data up: " << upD << " down: " << downD << endl << endl;
  //CMS_lumi( c_tau_muon_tau_hadron_phi_data, iPeriod, iPos );
  c_tau_muon_tau_hadron_phi_data->SaveAs((basePlotDir+"/singlePlot/tau_muon_tau_hadron_phi_data.png").c_str());
  c_tau_muon_tau_hadron_phi_data->SaveAs((basePlotDir+"/singlePlot/tau_muon_tau_hadron_phi_data.pdf").c_str());
  
  TCanvas *c_tau_muon_tau_hadron_phi_signal = new TCanvas("c_tau_muon_tau_hadron_phi_signal", "c_tau_muon_tau_hadron_phi_signal", 800, 800);
  if (nSamples>1){
    h_tau_mu_tau_hadron_phi[1]->Draw("colz");
    float upS = 0;
    float downS = 0;
    for (int xbin = 1; xbin <= tau_mu_phi_bins; xbin++){
      for (int ybin = 1; ybin <= tau_hadron_phi_bins; ybin++){
        if(ybin>xbin) upS += h_tau_mu_tau_hadron_phi[1]->GetBinContent(xbin,ybin);
        if(ybin<xbin) downS += h_tau_mu_tau_hadron_phi[1]->GetBinContent(xbin,ybin);
      }
    }
    cout << "tau_muon_tau_hadron_phi_signal up: " << upS << " down: " << downS << endl << endl;
  }
  c_tau_muon_tau_hadron_phi_signal->SaveAs((basePlotDir+"/singlePlot/tau_muon_tau_hadron_phi_signal.png").c_str());
  c_tau_muon_tau_hadron_phi_signal->SaveAs((basePlotDir+"/singlePlot/tau_muon_tau_hadron_phi_signal.pdf").c_str());

TCanvas *c_tau_had = new TCanvas("c_tau_had", "c_tau_had", 2000, 3000); c_tau_had->Divide(4,6);

  c_tau_had->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_tau_hadron_p[i]->GetMaximum() > h_tau_hadron_p[0]->GetMaximum()) h_tau_hadron_p[0]->SetMaximum(1.2*h_tau_hadron_p[i]->GetMaximum());
    if (h_tau_hadron_p[i]->GetMinimum() < h_tau_hadron_p[0]->GetMinimum() && h_tau_hadron_p[i]->GetMinimum() > 0) h_tau_hadron_p[0]->SetMinimum(0.8*h_tau_hadron_p[i]->GetMinimum());
  }
  h_tau_hadron_p[0]->Draw("e0e1x0");
  c_tau_had->cd(1)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_p[i]->Draw("hesame");
  
  c_tau_had->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_tau_hadron_pz[i]->GetMaximum() > h_tau_hadron_pz[0]->GetMaximum()) h_tau_hadron_pz[0]->SetMaximum(1.2*h_tau_hadron_pz[i]->GetMaximum());
    if (h_tau_hadron_pz[i]->GetMinimum() < h_tau_hadron_pz[0]->GetMinimum() && h_tau_hadron_pz[i]->GetMinimum() > 0) h_tau_hadron_pz[0]->SetMinimum(0.8*h_tau_hadron_pz[i]->GetMinimum());
  }
  h_tau_hadron_pz[0]->Draw("e0e1x0");
  c_tau_had->cd(2)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_pz[i]->Draw("hesame");
  
  c_tau_had->cd(3); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_pt->Add(h_tau_hadron_pt[i]);
    if (h_tau_hadron_pt[i]->GetMaximum() > h_tau_hadron_pt[0]->GetMaximum()) h_tau_hadron_pt[0]->SetMaximum(1.2*h_tau_hadron_pt[i]->GetMaximum());
  }
  h_tau_hadron_pt[0]->Draw("e0e1x0");
  c_tau_had->cd(3)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_pt[i]->Draw("hesame");
  legend->Draw();
  
  TCanvas *c_tau_hadron_pt = new TCanvas("c_tau_hadron_pt", "c_tau_hadron_pt", 800, 800);
  h_tau_hadron_pt[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_hadron_pt[i]->Draw("hesame");
  stacks.at(h_tau_hadron_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_pt, iPeriod, iPos );
  c_tau_hadron_pt->SaveAs((basePlotDir+"/singlePlot/tau_hadron_pt.png").c_str());
  c_tau_hadron_pt->SaveAs((basePlotDir+"/singlePlot/tau_hadron_pt.pdf").c_str());
  
  TCanvas *c_tau_hadron_ptScalar = new TCanvas("c_tau_hadron_ptScalar", "c_tau_hadron_ptScalar", 800, 800);
  h_tau_hadron_ptScalar[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_hadron_ptScalar[i]->Draw("hesame");
  stacks.at(h_tau_hadron_ptScalar[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_ptScalar, iPeriod, iPos );
  c_tau_hadron_ptScalar->SaveAs((basePlotDir+"/singlePlot/tau_hadron_scalar_pt.png").c_str());
  c_tau_hadron_ptScalar->SaveAs((basePlotDir+"/singlePlot/tau_hadron_scalar_pt.pdf").c_str());
  
  c_tau_had->cd(4); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_eta->Add(h_tau_hadron_eta[i]);
    if (h_tau_hadron_eta[i]->GetMaximum() > h_tau_hadron_eta[0]->GetMaximum()) h_tau_hadron_eta[0]->SetMaximum(1.2*h_tau_hadron_eta[0]->GetMaximum());
  }
  h_tau_hadron_eta[0]->Draw("e0e1x0");
  c_tau_had->cd(4)->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_tau_hadron_eta[i]->Draw("hesame");
  
  TCanvas *c_tau_hadron_eta = new TCanvas("c_tau_hadron_eta", "c_tau_hadron_eta", 800, 800);
  h_tau_hadron_eta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_hadron_eta[i]->Draw("hesame");
  stacks.at(h_tau_hadron_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_eta, iPeriod, iPos );
  c_tau_hadron_eta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_eta.png").c_str());
  c_tau_hadron_eta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_eta.pdf").c_str());
  
  c_tau_had->cd(5); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_phi->Add(h_tau_hadron_phi[i]);
    if (h_tau_hadron_phi[i]->GetMaximum() > h_tau_hadron_phi[0]->GetMaximum()) h_tau_hadron_phi[0]->SetMaximum(1.2*h_tau_hadron_phi[0]->GetMaximum());
  }
  h_tau_hadron_phi[0]->Draw("e0e1x0");
  c_tau_had->cd(5)->SetLogy(1);
  h_tau_hadron_phi[0]->GetYaxis()->SetRangeUser(0.1,1.2*h_tau_hadron_phi[0]->GetMaximum());
  for (int i = 1; i < nSamples; i++) h_tau_hadron_phi[i]->Draw("hesame");
  
  TCanvas *c_tau_hadron_phi = new TCanvas("c_tau_hadron_phi", "c_tau_hadron_phi", 800, 800);
  h_tau_hadron_phi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_hadron_phi[i]->Draw("hesame");
  stacks.at(h_tau_hadron_phi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_phi, iPeriod, iPos );
  c_tau_hadron_phi->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi.png").c_str());
  c_tau_hadron_phi->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi.pdf").c_str());
  
  
  TCanvas *c_tau_hadron_phi_lowNch_highHF = new TCanvas("c_tau_hadron_phi_lowNch_highHF", "c_tau_hadron_phi_lowNch_highHF", 800, 800);
  h_tau_hadron_phi[nSamples+1]->Draw("e0e1x0");
  
  CMS_lumi( c_tau_hadron_phi_lowNch_highHF, iPeriod, iPos );
  c_tau_hadron_phi_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_lowNch_highHF.png").c_str());
  c_tau_hadron_phi_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_lowNch_highHF.pdf").c_str());
  c_tau_hadron_phi_lowNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_lowNch_highHF.root").c_str());
  
  TCanvas *c_tau_hadron_phi_highNch_lowHF = new TCanvas("c_tau_hadron_phi_highNch_lowHF", "c_tau_hadron_phi_highNch_lowHF", 800, 800);
  h_tau_hadron_phi[nSamples+2]->Draw("e0e1x0");
  
  CMS_lumi( c_tau_hadron_phi_highNch_lowHF, iPeriod, iPos );
  c_tau_hadron_phi_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_highNch_lowHF.png").c_str());
  c_tau_hadron_phi_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_highNch_lowHF.pdf").c_str());
  c_tau_hadron_phi_highNch_lowHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_highNch_lowHF.root").c_str());
  
  TCanvas *c_tau_hadron_phi_highNch_highHF = new TCanvas("c_tau_hadron_phi_highNch_highHF", "c_tau_hadron_phi_highNch_highHF", 800, 800);
  h_tau_hadron_phi[nSamples]->Draw("e0e1x0");
  
  CMS_lumi( c_tau_hadron_phi_highNch_highHF, iPeriod, iPos );
  c_tau_hadron_phi_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_highNch_highHF.png").c_str());
  c_tau_hadron_phi_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_highNch_highHF.pdf").c_str());
  c_tau_hadron_phi_highNch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_phi_highNch_highHF.root").c_str());
  
  
  for (int j = 0; j < 2; j++){
    c_tau_had->cd(j+6);
    for (int i = 1; i < nSamples; i++){
      hs_tau_hadron_rhomass[j]->Add(h_tau_hadron_rhomass[i][j]);
      if (h_tau_hadron_rhomass[i][j]->GetMaximum() > h_tau_hadron_rhomass[0][j]->GetMaximum()) h_tau_hadron_rhomass[0][j]->SetMaximum(h_tau_hadron_rhomass[i][j]->GetMaximum());
    }
    h_tau_hadron_rhomass[0][j]->Draw("e0e1x0");
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
  cout << "data:\nnch = 3: " << h_tau_hadron_nch[0]->GetBinContent(4) << "\nnch = 4: " << h_tau_hadron_nch[0]->GetBinContent(5) << endl;
  cout << "signal:\nnch = 3: " << h_tau_hadron_nch[1]->GetBinContent(4) << "\nnch = 4: " << h_tau_hadron_nch[1]->GetBinContent(5) << endl;
  h_tau_hadron_nch[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_nch[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_nch, iPeriod, iPos );
  c_tau_hadron_nch->SaveAs((basePlotDir+"/singlePlot/tau_hadron_nch.png").c_str());
  c_tau_hadron_nch->SaveAs((basePlotDir+"/singlePlot/tau_hadron_nch.pdf").c_str());
  
  TCanvas *c_tau_hadron_nch_highHF = new TCanvas("c_tau_hadron_nch_highHF", "c_tau_hadron_nch_highHF", 800, 800);
  h_tau_hadron_nch_highHF[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_nch_highHF[i]->Draw("hesame");
  CMS_lumi( c_tau_hadron_nch_highHF, iPeriod, iPos );
  c_tau_hadron_nch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_nch_highHF.png").c_str());
  c_tau_hadron_nch_highHF->SaveAs((basePlotDir+"/singlePlot/tau_hadron_nch_highHF.pdf").c_str());
  
  
  c_tau_had->cd(9); h_tau_hadron_ncand_final[0]->Draw("e1same"); for (int i = 1; i < nSamples; i++) h_tau_hadron_ncand_final[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(10); h_tau_hadron_vprob[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++) h_tau_hadron_vprob[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(11); for (int i = 1; i < nSamples; i++) h_tau_hadron_matched_pt_index[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(12); for (int i = 1; i < nSamples; i++){
    hs_tau_hadron_mass->Add(h_tau_hadron_mass[i]);
    if (h_tau_hadron_mass[i]->GetMaximum() > h_tau_hadron_mass[0]->GetMaximum()) h_tau_hadron_mass[0]->SetMaximum(1.2*h_tau_hadron_mass[i]->GetMaximum());
  }
  //h_tau_hadron_mass[0]->SetTitle("#tau_{3prong} visible mass");
  h_tau_hadron_mass[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_mass[i]->Draw("hesame");
  
  TCanvas *c_tau_hadron_mass = new TCanvas("c_tau_hadron_mass", "c_tau_hadron_mass", 800, 800);
  h_tau_hadron_mass[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_tau_hadron_mass[i]->Draw("hesame");
  stacks.at(h_tau_hadron_mass[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_mass, iPeriod, iPos );
  c_tau_hadron_mass->SaveAs((basePlotDir+"/singlePlot/tau_hadron_mass.png").c_str());
  c_tau_hadron_mass->SaveAs((basePlotDir+"/singlePlot/tau_hadron_mass.pdf").c_str());
  c_tau_hadron_mass->SaveAs((basePlotDir+"/singlePlot/tau_hadron_mass.C").c_str());
  
  c_tau_had->cd(13); for (int i = 1; i < nSamples; i++){
    hs_ditau_mass->Add(h_ditau_mass[i]);
    if (h_ditau_mass[i]->GetMaximum() > h_ditau_mass[0]->GetMaximum()) h_ditau_mass[0]->SetMaximum(1.2*h_ditau_mass[i]->GetMaximum());
  }
  //h_ditau_mass[0]->SetTitle("#tau#tau visible mass");
  h_ditau_mass[0]->Draw("e1same");
  for (int i = 1; i < nSamples; i++) h_ditau_mass[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_ditau_mass = new TCanvas("c_ditau_mass", "c_ditau_mass", 800, 800);
  h_ditau_mass[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_ditau_mass[i]->Draw("hesame");
  stacks.at(h_ditau_mass[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_ditau_mass, iPeriod, iPos );
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass.png").c_str());
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass.pdf").c_str());
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass.C").c_str());
  
  TCanvas *c_ditau_ptScalar = new TCanvas("c_ditau_ptScalar", "c_ditau_ptScalar", 800, 800);
  h_ditau_ptScalar[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_ditau_ptScalar[i]->Draw("hesame");
  stacks.at(h_ditau_ptScalar[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_ditau_ptScalar, iPeriod, iPos );
  c_ditau_ptScalar->SaveAs((basePlotDir+"/singlePlot/ditau_scalar_pt.png").c_str());
  c_ditau_ptScalar->SaveAs((basePlotDir+"/singlePlot/ditau_scalar_pt.pdf").c_str());
  
  TH1F *stackedDitauScalarPt = new TH1F(*(TH1F*)((stacks.at(h_ditau_ptScalar[nSamples+1]->GetBinContent(1))->GetStack()->Last())));
  cout << " **** **** **** ****" << endl << "KS test for ditau scalar pt: " << h_ditau_ptScalar[0]->KolmogorovTest(stackedDitauScalarPt) << endl << endl;
  
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
  
  TCanvas *c_ditau_pt = new TCanvas("c_ditau_pt", "c_ditau_pt", 800, 800);
  h_ditau_pt[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_ditau_pt[i]->Draw("hesame");
  stacks.at(h_ditau_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_ditau_pt, iPeriod, iPos );
  c_ditau_pt->SaveAs((basePlotDir+"/singlePlot/ditau_pt.png").c_str());
  c_ditau_pt->SaveAs((basePlotDir+"/singlePlot/ditau_pt.pdf").c_str());
  
  c_tau_had->cd(16); h_tau_hadron_nPions[0]->Scale(1./h_tau_hadron_nPions[0]->GetMaximum()); h_tau_hadron_nPions[0]->Draw("e0e1x0");
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
  
  c_tau_had->SaveAs((basePlotDir+"/collective/tau_had."+plotFormat).c_str());

  TCanvas *c_recoPion = new TCanvas("c_recoPion", "c_recoPion", 3*800, 3*800); c_recoPion->Divide(3,3);
  
  c_recoPion->cd(1); for (int i = 1; i < nSamples; i++){
    if (h_pion_leading_pt[i]->GetMaximum() > h_pion_leading_pt[0]->GetMaximum()) h_pion_leading_pt[0]->SetMaximum(1.2*h_pion_leading_pt[i]->GetMaximum());
  }
  h_pion_leading_pt[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_leading_pt[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_leading_pt = new TCanvas("c_pion_leading_pt", "c_pion_leading_pt", 800, 800);
  h_pion_leading_pt[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_leading_pt[i]->Draw("hesame");
  stacks.at(h_pion_leading_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_leading_pt, iPeriod, iPos );
  c_pion_leading_pt->SaveAs((basePlotDir+"/singlePlot/pion_leading_pt.png").c_str());
  c_pion_leading_pt->SaveAs((basePlotDir+"/singlePlot/pion_leading_pt.pdf").c_str());
  
  c_recoPion->cd(2); for (int i = 1; i < nSamples; i++){
    if (h_pion_subleading_pt[i]->GetMaximum() > h_pion_subleading_pt[0]->GetMaximum()) h_pion_subleading_pt[0]->SetMaximum(1.2*h_pion_subleading_pt[i]->GetMaximum());
  }
  h_pion_subleading_pt[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subleading_pt[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subleading_pt = new TCanvas("c_pion_subleading_pt", "c_pion_subleading_pt", 800, 800);
  h_pion_subleading_pt[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_subleading_pt[i]->Draw("hesame");
  stacks.at(h_pion_subleading_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_subleading_pt, iPeriod, iPos );
  c_pion_subleading_pt->SaveAs((basePlotDir+"/singlePlot/pion_subleading_pt.png").c_str());
  c_pion_subleading_pt->SaveAs((basePlotDir+"/singlePlot/pion_subleading_pt.pdf").c_str());
  
  c_recoPion->cd(3); for (int i = 1; i < nSamples; i++){
    if (h_pion_subsubleading_pt[i]->GetMaximum() > h_pion_subsubleading_pt[0]->GetMaximum()) h_pion_subsubleading_pt[0]->SetMaximum(1.2*h_pion_subsubleading_pt[i]->GetMaximum());
  }
  h_pion_subsubleading_pt[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_pt[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subsubleading_pt = new TCanvas("c_pion_subsubleading_pt", "c_pion_subsubleading_pt", 800, 800);
  h_pion_subsubleading_pt[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_subsubleading_pt[i]->Draw("hesame");
  stacks.at(h_pion_subsubleading_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_subsubleading_pt, iPeriod, iPos );
  c_pion_subsubleading_pt->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_pt.png").c_str());
  c_pion_subsubleading_pt->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_pt.pdf").c_str());
  
  c_recoPion->cd(4); for (int i = 1; i < nSamples; i++){
    if (h_pion_leading_eta[i]->GetMaximum() > h_pion_leading_eta[0]->GetMaximum()) h_pion_leading_eta[0]->SetMaximum(1.2*h_pion_leading_eta[i]->GetMaximum());
  }
  h_pion_leading_eta[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_leading_eta[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_leading_eta = new TCanvas("c_pion_leading_eta", "c_pion_leading_eta", 800, 800);
  h_pion_leading_eta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_leading_eta[i]->Draw("hesame");
  stacks.at(h_pion_leading_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_leading_eta, iPeriod, iPos );
  c_pion_leading_eta->SaveAs((basePlotDir+"/singlePlot/pion_leading_eta.png").c_str());
  c_pion_leading_eta->SaveAs((basePlotDir+"/singlePlot/pion_leading_eta.pdf").c_str());
  
  c_recoPion->cd(5); for (int i = 1; i < nSamples; i++){
    if (h_pion_subleading_eta[i]->GetMaximum() > h_pion_subleading_eta[0]->GetMaximum()) h_pion_subleading_eta[0]->SetMaximum(1.2*h_pion_subleading_eta[i]->GetMaximum());
  }
  h_pion_subleading_eta[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subleading_eta[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subleading_eta = new TCanvas("c_pion_subleading_eta", "c_pion_subleading_eta", 800, 800);
  h_pion_subleading_eta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_subleading_eta[i]->Draw("hesame");
  stacks.at(h_pion_subleading_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_subleading_eta, iPeriod, iPos );
  c_pion_subleading_eta->SaveAs((basePlotDir+"/singlePlot/pion_subleading_eta.png").c_str());
  c_pion_subleading_eta->SaveAs((basePlotDir+"/singlePlot/pion_subleading_eta.pdf").c_str());
  
  c_recoPion->cd(6); for (int i = 1; i < nSamples; i++){
    if (h_pion_subsubleading_eta[i]->GetMaximum() > h_pion_subsubleading_eta[0]->GetMaximum()) h_pion_subsubleading_eta[0]->SetMaximum(1.2*h_pion_subsubleading_eta[i]->GetMaximum());
  }
  h_pion_subsubleading_eta[0]->SetMaximum(1.3*h_pion_subsubleading_eta[0]->GetMaximum());
  h_pion_subsubleading_eta[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_eta[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subsubleading_eta = new TCanvas("c_pion_subsubleading_eta", "c_pion_subsubleading_eta", 800, 800);
  h_pion_subsubleading_eta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_subsubleading_eta[i]->Draw("hesame");
  stacks.at(h_pion_subsubleading_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_subsubleading_eta, iPeriod, iPos );
  c_pion_subsubleading_eta->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_eta.png").c_str());
  c_pion_subsubleading_eta->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_eta.pdf").c_str());
  
  c_recoPion->cd(7); for (int i = 1; i < nSamples; i++){
    if (h_pion_leading_phi[i]->GetMaximum() > h_pion_leading_phi[0]->GetMaximum()) h_pion_leading_phi[0]->SetMaximum(1.2*h_pion_leading_phi[i]->GetMaximum());
  }
  h_pion_leading_phi[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_leading_phi[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_leading_phi = new TCanvas("c_pion_leading_phi", "c_pion_leading_phi", 800, 800);
  h_pion_leading_phi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_leading_phi[i]->Draw("hesame");
  stacks.at(h_pion_leading_phi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_leading_phi, iPeriod, iPos );
  c_pion_leading_phi->SaveAs((basePlotDir+"/singlePlot/pion_leading_phi.png").c_str());
  c_pion_leading_phi->SaveAs((basePlotDir+"/singlePlot/pion_leading_phi.pdf").c_str());
  
  c_recoPion->cd(8); for (int i = 1; i < nSamples; i++){
    if (h_pion_subleading_phi[i]->GetMaximum() > h_pion_subleading_phi[0]->GetMaximum()) h_pion_subleading_phi[0]->SetMaximum(1.2*h_pion_subleading_phi[i]->GetMaximum());
  }
  h_pion_subleading_phi[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subleading_phi[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subleading_phi = new TCanvas("c_pion_subleading_phi", "c_pion_subleading_phi", 800, 800);
  h_pion_subleading_phi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_subleading_phi[i]->Draw("hesame");
  stacks.at(h_pion_subleading_phi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_subleading_phi, iPeriod, iPos );
  c_pion_subleading_phi->SaveAs((basePlotDir+"/singlePlot/pion_subleading_phi.png").c_str());
  c_pion_subleading_phi->SaveAs((basePlotDir+"/singlePlot/pion_subleading_phi.pdf").c_str());
  
  c_recoPion->cd(9); for (int i = 1; i < nSamples; i++){
    if (h_pion_subsubleading_phi[i]->GetMaximum() > h_pion_subsubleading_phi[0]->GetMaximum()) h_pion_subsubleading_phi[0]->SetMaximum(1.2*h_pion_subsubleading_phi[i]->GetMaximum());
  }
  h_pion_subsubleading_phi[0]->Draw("e1");
  for (int i = 1; i < nSamples; i++) h_pion_subsubleading_phi[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_pion_subsubleading_phi = new TCanvas("c_pion_subsubleading_phi", "c_pion_subsubleading_phi", 800, 800);
  h_pion_subsubleading_phi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_pion_subsubleading_phi[i]->Draw("hesame");
  stacks.at(h_pion_subsubleading_phi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_subsubleading_phi, iPeriod, iPos );
  c_pion_subsubleading_phi->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_phi.png").c_str());
  c_pion_subsubleading_phi->SaveAs((basePlotDir+"/singlePlot/pion_subsubleading_phi.pdf").c_str());
  
  c_recoPion->SaveAs((basePlotDir+"/collective/reco_pion."+plotFormat).c_str());

  TCanvas *c_calo = new TCanvas("c_calo", "c_calo", 1600, 5*800); c_calo->Divide(2,5);
  for (int i = 1; i < nSamples; i++) hs_calo_energyHFp->Add(h_calo_energyHFp[i]);
  c_calo->cd(1); h_calo_energyHFp[0]->Draw("e0e1x0"); h_calo_energyHFp[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_energyHFp[0]->GetMaximum());
  c_calo->cd(1)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_energyHFp[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++) hs_calo_energyHFm->Add(h_calo_energyHFm[i]);
  c_calo->cd(2); h_calo_energyHFm[0]->Draw("e0e1x0"); h_calo_energyHFm[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_energyHFm[0]->GetMaximum());
  c_calo->cd(2)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_energyHFm[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++){
    hs_calo_leadingHFp->Add(h_calo_leadingHFp[i]);
    if (h_calo_leadingHFp[i]->GetMaximum() > h_calo_leadingHFp[0]->GetMaximum()) h_calo_leadingHFp[0]->SetMaximum(1.2*h_calo_leadingHFp[i]->GetMaximum());
  }
  h_calo_leadingHFp[0]->SetMaximum(2.5*h_calo_leadingHFp[0]->GetMaximum());
  c_calo->cd(3); h_calo_leadingHFp[0]->Draw("e0e1x0"); //h_calo_leadingHFp[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_leadingHFp[0]->GetMaximum());
  c_calo->cd(3)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_leadingHFp[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_calo_leadingHFp = new TCanvas("c_calo_leadingHFp", "c_calo_leadingHFp", 800, 800);
  h_calo_leadingHFp[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_calo_leadingHFp[i]->Draw("hesame");
  CMS_lumi( c_calo_leadingHFp, iPeriod, iPos );
  c_calo_leadingHFp->cd()->SetLogy(1);
  c_calo_leadingHFp->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFp.png").c_str());
  c_calo_leadingHFp->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFp.pdf").c_str());
  
  TCanvas *c_calo_leadingHFp_highNch = new TCanvas("c_calo_leadingHFp_highNch", "c_calo_leadingHFp_highNch", 800, 800);
  h_calo_leadingHFp_highNch[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_calo_leadingHFp_highNch[i]->Draw("hesame");
  CMS_lumi( c_calo_leadingHFp_highNch, iPeriod, iPos );
  c_calo_leadingHFp_highNch->cd()->SetLogy(1);
  c_calo_leadingHFp_highNch->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFp_highNch.png").c_str());
  c_calo_leadingHFp_highNch->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFp_highNch.pdf").c_str());
  
  for (int i = 1; i < nSamples; i++){
    hs_calo_leadingHFm->Add(h_calo_leadingHFm[i]);
    if (h_calo_leadingHFm[i]->GetMaximum() > h_calo_leadingHFm[0]->GetMaximum()) h_calo_leadingHFm[0]->SetMaximum(1.2*h_calo_leadingHFm[i]->GetMaximum());
  }
  h_calo_leadingHFm[0]->SetMaximum(2.5*h_calo_leadingHFm[0]->GetMaximum());
  c_calo->cd(4); h_calo_leadingHFm[0]->Draw("e0e1x0"); //h_calo_leadingHFm[0]->GetYaxis()->SetRangeUser(0.1,2*h_calo_leadingHFm[0]->GetMaximum());
  c_calo->cd(4)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_leadingHFm[i]->Draw("hesame"); legend->Draw();
  
  TCanvas *c_calo_leadingHFm = new TCanvas("c_calo_leadingHFm", "c_calo_leadingHFm", 800, 800);
  h_calo_leadingHFm[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_calo_leadingHFm[i]->Draw("hesame");
  CMS_lumi( c_calo_leadingHFm, iPeriod, iPos );
  c_calo_leadingHFm->cd()->SetLogy(1);
  c_calo_leadingHFm->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFm.png").c_str());
  c_calo_leadingHFm->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFm.pdf").c_str());
  
  TCanvas *c_calo_leadingHFm_highNch = new TCanvas("c_calo_leadingHFm_highNch", "c_calo_leadingHFm_highNch", 800, 800);
  h_calo_leadingHFm_highNch[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_calo_leadingHFm_highNch[i]->Draw("hesame");
  CMS_lumi( c_calo_leadingHFm_highNch, iPeriod, iPos );
  c_calo_leadingHFm_highNch->cd()->SetLogy(1);
  c_calo_leadingHFm_highNch->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFm_highNch.png").c_str());
  c_calo_leadingHFm_highNch->SaveAs((basePlotDir+"/singlePlot/calo_leadingHFm_highNch.pdf").c_str());
  
  c_calo->cd(5); h_calo_energyHFp_nch->Draw("COLZ"); c_calo->cd(5)->SetLogz(1);
  c_calo->cd(6); h_calo_energyHFm_nch->Draw("COLZ"); c_calo->cd(6)->SetLogz(1);
  c_calo->cd(7); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFp_sum[i]->Draw("e1same");
  c_calo->cd(8); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFm_sum[i]->Draw("e1same");
  c_calo->cd(9); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFp_size[i]->Draw("e1same");
  c_calo->cd(10); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFm_size[i]->Draw("e1same");
  c_calo->SaveAs((basePlotDir+"/collective/calo."+plotFormat).c_str());
  
  TCanvas *c_calo_E_eta = new TCanvas("c_calo_E_eta", "c_calo_E_eta", nSamples*800, 800); c_calo_E_eta->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_E_eta->cd(i+1);
    h_calo_E_eta[i]->Draw("COLZ");
    c_calo_E_eta->cd(i+1)->SetLogz(1);
  }
  c_calo_E_eta->SaveAs((basePlotDir+"/collective/calo_E_eta.png").c_str());
  c_calo_E_eta->SaveAs((basePlotDir+"/collective/calo_E_eta.pdf").c_str());
  
  TCanvas *c_calo_leadingE_eta = new TCanvas("c_calo_leadingE_eta", "c_calo_leadingE_eta", nSamples*800, 800); c_calo_leadingE_eta->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_leadingE_eta->cd(i+1);
    h_calo_leadingE_eta[i]->Draw("COLZ");
    c_calo_leadingE_eta->cd(i+1)->SetLogz(1);
  }
  c_calo_leadingE_eta->SaveAs((basePlotDir+"/collective/calo_leadingE_eta.png").c_str());
  c_calo_leadingE_eta->SaveAs((basePlotDir+"/collective/calo_leadingE_eta.pdf").c_str());
  
  TCanvas *c_calo_Et_eta = new TCanvas("c_calo_Et_eta", "c_calo_Et_eta", nSamples*800, 800); c_calo_Et_eta->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_Et_eta->cd(i+1);
    h_calo_Et_eta[i]->Draw("COLZ");
    c_calo_Et_eta->cd(i+1)->SetLogz(1);
  }
  c_calo_Et_eta->SaveAs((basePlotDir+"/collective/calo_Et_eta.png").c_str());
  c_calo_Et_eta->SaveAs((basePlotDir+"/collective/calo_Et_eta.pdf").c_str());
  
  TCanvas *c_calo_leadingEt_eta = new TCanvas("c_calo_leadingEt_eta", "c_calo_leadingEt_eta", nSamples*800, 800); c_calo_leadingEt_eta->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_leadingEt_eta->cd(i+1);
    h_calo_leadingEt_eta[i]->Draw("COLZ");
    c_calo_leadingEt_eta->cd(i+1)->SetLogz(1);
  }
  c_calo_leadingEt_eta->SaveAs((basePlotDir+"/collective/calo_leadingEt_eta.png").c_str());
  c_calo_leadingEt_eta->SaveAs((basePlotDir+"/collective/calo_leadingEt_eta.pdf").c_str());
  
  TCanvas *c_calo_E = new TCanvas("c_calo_E", "c_calo_E", 800, 800);
  h_calo_E[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_calo_E[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_E, iPeriod, iPos );
  c_calo_E->cd()->SetLogy(1);
  c_calo_E->SaveAs((basePlotDir+"/singlePlot/calo_E.png").c_str());
  c_calo_E->SaveAs((basePlotDir+"/singlePlot/calo_E.pdf").c_str());
  
  TCanvas *c_calo_sumE = new TCanvas("c_calo_sumE", "c_calo_sumE", 800, 800);
  h_calo_sumE[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_calo_sumE[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_sumE, iPeriod, iPos );
  //c_calo_sumEt->cd()->SetLogy(1);
  c_calo_sumE->SaveAs((basePlotDir+"/singlePlot/calo_sumE.png").c_str());
  c_calo_sumE->SaveAs((basePlotDir+"/singlePlot/calo_sumE.pdf").c_str());
  
  TCanvas *c_calo_leadingE = new TCanvas("c_calo_leadingE", "c_calo_leadingE", 800, 800);
  h_calo_leadingE[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_calo_leadingE[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_leadingE, iPeriod, iPos );
  //c_calo_leadingE->cd()->SetLogy(1);
  c_calo_leadingE->SaveAs((basePlotDir+"/singlePlot/calo_leadingE.png").c_str());
  c_calo_leadingE->SaveAs((basePlotDir+"/singlePlot/calo_leadingE.pdf").c_str());
  
  TCanvas *c_calo_Et = new TCanvas("c_calo_Et", "c_calo_Et", 800, 800);
  h_calo_Et[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_calo_Et[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_Et, iPeriod, iPos );
  c_calo_Et->cd()->SetLogy(1);
  c_calo_Et->SaveAs((basePlotDir+"/singlePlot/calo_Et.png").c_str());
  c_calo_Et->SaveAs((basePlotDir+"/singlePlot/calo_Et.pdf").c_str());
  
  TCanvas *c_calo_sumEt = new TCanvas("c_calo_sumEt", "c_calo_sumEt", 800, 800);
  h_calo_sumEt[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_calo_sumEt[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_sumEt, iPeriod, iPos );
  //c_calo_sumEt->cd()->SetLogy(1);
  c_calo_sumEt->SaveAs((basePlotDir+"/singlePlot/calo_sumEt.png").c_str());
  c_calo_sumEt->SaveAs((basePlotDir+"/singlePlot/calo_sumEt.pdf").c_str());
  
  TCanvas *c_calo_leadingEt = new TCanvas("c_calo_leadingEt", "c_calo_leadingEt", 800, 800);
  h_calo_leadingEt[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_calo_leadingEt[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_leadingEt, iPeriod, iPos );
  //c_calo_leadingEt->cd()->SetLogy(1);
  c_calo_leadingEt->SaveAs((basePlotDir+"/singlePlot/calo_leadingEt.png").c_str());
  c_calo_leadingEt->SaveAs((basePlotDir+"/singlePlot/calo_leadingEt.pdf").c_str());
  
  TCanvas *c_nCaloTowers = new TCanvas("c_nCaloTowers", "c_nCaloTowers", 800, 800);
  h_nCaloTowers[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    h_nCaloTowers[i]->Draw("hesame");
  }
  CMS_lumi( c_nCaloTowers, iPeriod, iPos );
  //c_nCaloTowers->cd()->SetLogy(1);
  c_nCaloTowers->SaveAs((basePlotDir+"/singlePlot/nCaloTowers.png").c_str());
  c_nCaloTowers->SaveAs((basePlotDir+"/singlePlot/nCaloTowers.pdf").c_str());

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
  c_zdc->SaveAs((basePlotDir+"/collective/zdc."+plotFormat).c_str());


  TCanvas *c_HF_eta = new TCanvas("c_HF_eta", "c_HF_eta", 500*(nSamples+1), 500); c_HF_eta->Divide(nSamples+1,1);
  for (int i = 1; i < nSamples; i++) hs_calo_HF_eta->Add(h_calo_HF_eta[i]);
  c_HF_eta->cd(1); h_calo_HF_eta[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++) h_calo_HF_eta[i]->Draw("hesame"); legend->Draw();
  for (int i = 0; i < nSamples; i++){ c_HF_eta->cd(i+2); h_calo_HF_energy_eta[i]->Draw("colz"); c_HF_eta->cd(i+2)->SetLogz(1);}
  c_HF_eta->SaveAs((basePlotDir+"/collective/HF_eta."+plotFormat).c_str());

  TCanvas *c_HFpm = new TCanvas("c_HFpm", "c_HFpm", 500*nSamples, 3000); c_HFpm->Divide(nSamples,6);
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(i+1); h_calo_energyHF_pm[i]->Draw("colz"); c_HFpm->cd(i+1)->SetLogz(1);}
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(nSamples+i+1); h_calo_leadingHF_pm[i]->Draw("colz"); c_HFpm->cd(nSamples+i+1)->SetLogz(1);}
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(2*nSamples+i+1); h_calo_muEta_leadingHFeta[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(3*nSamples+i+1); h_calo_tauEta_leadingHFeta[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(4*nSamples+i+1); h_calo_muEta_averageHFeta[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(5*nSamples+i+1); h_calo_tauEta_averageHFeta[i]->Draw("colz"); }
  c_HFpm->SaveAs((basePlotDir+"/collective/HFpm."+plotFormat).c_str());

  TCanvas *c_delta_phi_eta = new TCanvas("c_delta_phi_eta", "c_delta_phi_eta", 500*nSamples, 1500); c_delta_phi_eta->Divide(nSamples,3);
  for (int i = 0; i < nSamples; i++){
    c_delta_phi_eta->cd(i+1); h_deltaphi_tau_mu_tau_hadron_mueta[i]->Draw("colz");
    c_delta_phi_eta->cd(nSamples+i+1); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->Draw("colz");
    c_delta_phi_eta->cd(2*nSamples+i+1); h_mueta_taueta[i]->Draw("colz");
  }
  c_delta_phi_eta->SaveAs((basePlotDir+"/collective/delta_phi_eta."+plotFormat).c_str());

  TCanvas *c_AP = new TCanvas("c_AP", "c_AP", 500*nSamples, 500); c_AP->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_AP->cd(i+1); h_AP[i]->Draw("colz");}
  c_AP->SaveAs((basePlotDir+"/collective/AP."+plotFormat).c_str());

  TCanvas *c_rho = new TCanvas("c_rho", "c_rho", 500*nSamples, 500); c_rho->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_rho->cd(i+1); h_tau_hadron_rhomass2D[i]->Draw("colz");}
  c_rho->SaveAs((basePlotDir+"/collective/rho."+plotFormat).c_str());

  TCanvas *c_MET = new TCanvas("c_MET", "c_MET", 800, 800);
  c_MET->cd(1); h_MET[0]->Draw("e0e1x0"); for (int i = 1; i < nSamples; i++) h_MET[i]->Draw("hesame"); legend->Draw();
  c_MET->SaveAs((basePlotDir+"/collective/MET."+plotFormat).c_str());

  TCanvas *c_cutflow = new TCanvas("c_cutflow", "c_cutflow", 1000, 500); c_cutflow->Divide(2,1);
  c_cutflow->cd(1);
  cutflow[0]->Draw("e0e1x0");
  cutflow[0]->GetYaxis()->SetRangeUser(1,2*cutflow[0]->GetMaximum());
  for (int i = 1; i < nSamples; i++) cutflow[i]->Draw("hesame");
  legend->Draw();
  c_cutflow->cd(1)->SetLogy(1);
  c_cutflow->cd(2);
  h_cutflow[0]->Draw("e0e1x0");
  c_cutflow->cd(2)->SetLogy(1);
  c_cutflow->SaveAs((basePlotDir+"/collective/cutflow."+plotFormat).c_str());

  TFile *outputFile = new TFile("hists.root","RECREATE");
  outputFile->cd();
  outputFile->mkdir("3prong");
  outputFile->cd("3prong");
  
  TCanvas *c_output = new TCanvas("c_output", "c_output", 800*3, 800); c_output->Divide(nSamples,1);
  //TH1F *h_deltaphi_tau_mu_tau_hadron_new_data = (TH1F*)pre_data->Clone("data_obs");
  TH1F *h_deltaphi_tau_mu_tau_hadron_new_data = (TH1F*)h_deltaphi_tau_mu_tau_hadron[0]->Clone("data_obs");
  c_output->cd(1); h_deltaphi_tau_mu_tau_hadron_new_data->DrawCopy("e0e1x0");
  h_deltaphi_tau_mu_tau_hadron_new_data->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron_new_data->Write();
  
  //TH1F *h_deltaphi_tau_mu_tau_hadron_new_MC = (TH1F*)pre_signal->Clone("MC-signal");
  TH1F *h_deltaphi_tau_mu_tau_hadron_new_MC = (TH1F*)h_deltaphi_tau_mu_tau_hadron[1]->Clone("MC-signal");
  c_output->cd(2); h_deltaphi_tau_mu_tau_hadron_new_MC->DrawCopy("e0e1x0");
  h_deltaphi_tau_mu_tau_hadron_new_MC->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron_new_MC->Write();
  
  //c_output->cd(3); pre_background->DrawCopy("e0e1x0");
  //pre_background->SetDirectory(gDirectory);pre_background->Write();
  c_output->cd(3); background->DrawCopy("e0e1x0");
  background->SetDirectory(gDirectory); background->Write();
  
  background_sysNchDown->SetDirectory(gDirectory); background_sysNchDown->Write();
  background_sysNchUp->SetDirectory(gDirectory); background_sysNchUp->Write();
  for (int var = 0; var < 2+nNchCategories; var++){
    B_lowNch_highHF[nSamples+var]->SetDirectory(gDirectory); B_lowNch_highHF[nSamples+var]->Write();
  }
  
  if (nSamples>1){
    TH1F *temp_h_eff_matchedPionPt = (TH1F*)h_eff_matchedPionReco_genPionPt[1]->Clone("h_eff_matchedPionPt");
    temp_h_eff_matchedPionPt->SetDirectory(gDirectory);
    temp_h_eff_matchedPionPt->Write();
  }
  
  
  h_sys_muon_SF_Down->SetDirectory(gDirectory); h_sys_muon_SF_Down->Write();
  h_sys_muon_SF_Up->SetDirectory(gDirectory); h_sys_muon_SF_Up->Write();
  
  
  for (int ebin = 0; ebin < 3; ebin++){
    for (int pbin = 0; pbin < 14; pbin++){
      h_tauSFdown[ebin][pbin]->SetDirectory(gDirectory); h_tauSFdown[ebin][pbin]->Write();
      h_tauSFup[ebin][pbin]->SetDirectory(gDirectory); h_tauSFup[ebin][pbin]->Write();
    }
  }
  h_tauSFmean->SetDirectory(gDirectory); h_tauSFmean->Write();
  
  /*for (int i = 1; i < nSamples; i++){
    c_output->cd(i+1); h_deltaphi_tau_mu_tau_hadron[i]->DrawCopy("e1");
    h_deltaphi_tau_mu_tau_hadron[i]->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron[i]->Write();
  }*/
  //outputFile->Write();
  outputFile->Close();
    
  }