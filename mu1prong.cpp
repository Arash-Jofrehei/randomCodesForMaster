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
#define plotting
//#define post
//#define ratio
//#define ZDC0n0n

#define writeSkimData
//#define readSkimData
string tagSkimData = "skim_mu1prong_ditauPt05mass29_maxNch4_maxHF20";
//string tagSkimData = "skim_mu1prong_ditauPt05massJpsi_maxNch4_maxHF20";

#ifdef post
string anaTag = "postfit";
#else
string anaTag = "test";
//string anaTag = "prefit";
#endif

string plotFormat = "png";


// **************************
// cut definition
  float HFpLeading_low = 0;
  float HFmLeading_low = 0;
  float HFpLeading_high = 6.0;
  float HFmLeading_high = 6.0;
  float pionLeadingCut = 0.2;//0.5
  float pionExtraCut = 0.0;
  float ditau_pt_cut = 0.75;
  float aco_cut_for_low_ditau_pt = 1.0;
  float ditauMassCut = 4;
  float barrel_mu_pt_cut = 3.5;
  float endcap_mu_pt_cut = 2.5;
  //float acoGammaDitau_cut[2] = {1./100,1./50};
  //float resDitauGammaPt_cut_down[2] = {-1.0,-1.0};
  //float resDitauGammaPt_cut_up[2] = {0.3,0.3};
  float resDitauGammaPt_cut_down = -0.5;
  float resDitauGammaPt_cut_up = 0.5;
  float gammaPtCut = 0.4;
  float gammaECut = 0.;
  float gammaEtaCut = 2.2; //2.4
  float gammaDeltaRCut = 0.35;
  float FSR_ditau_acoCut = 0.05;
  float piZeroPtCut = 0.8;
  float piZeroDeltaPhiCut = 0.5;
  float piZeroDeltaRCut = 0; //0.7;
  float muon_DEta = 0.1;
  float pion_DEta = 0.1;
  float muonEB_DEta = 0.3;
  float muonEE_DEta = 0.3;
  float muonHB_DEta = 0.05;
  float muonHE_DEta = 0.05;
  float pionEB_DEta = 0.3;
  float pionEE_DEta = 0.3;
  float pionHB_DEta = 0.2;
  float pionHE_DEta = 0.2;
  float MET_EB_DR = 0.2;
  float MET_EE_DR = 0.2;
  float MET_HB_DR = 0.2;
  float MET_HE_DR = 0.2;
  float minZDCp = 0;
  float minZDCm = 0;
  float maxZDCp = 10000000;
  float maxZDCm = 10000000;
#ifdef ZDC0n0n
  maxZDCp = 4200;
  maxZDCm = 6000;
#endif
  float ratioZDCpm = 4.2/6;
  int min_nch = 1;
  int max_nch = 1;
  float MET_cut = 20000;
  //float deltaPhi_cut = 2.78816;
  float deltaPhi_cut = 0.8*TMath::Pi();
  int MuTauCharge = -1;
  bool hasZDCinfo = true;
  int nNchCategories = 2;
  int firstNchCategory = 2;
  bool tau_muon_isSoft = true;
  bool tau_muon_isGlobal = false;
  bool tau_muon_isTracker = false;
  double tau_hadron_vertexprob = 0;//0.005;
// end of cut definition
// **************************

double my_pi=TMath::Pi();
double deltaphi(double phi1, double phi2)
{
 double dphi=fabs(phi1-phi2);
 return (dphi<=my_pi)? dphi : 2.*my_pi-dphi;
}

double acop(TLorentzVector a, TLorentzVector b){
  return 1-TMath::Abs(a.DeltaPhi(b)/TMath::Pi());
}

string baseDir = "/eos/user/a/ajofrehe/gtau/ntuples/";

//void run(string fileData, string fileMCS){//, string fileMCB1, string fileMCB2){
  //string files[2] = {fileData,fileMCS};
  //const int nSamples = 2;

//string inputFiles[] = {"flatTuple_AOD_data_1prong_2015_280721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_pp_MadGraph_2015_160721.root","flatTuple_mcggTauTau_1prong_AOD_SuperChic_2018_280721.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_1prong_2018_140622.root","flatTuple_mu1prong_AOD_MadGraph_2018_140622.root","flatTuple_mumu30k_AOD_2018.root","flatTuple_mumuFSR30k_AOD_2018.root","flatTuple_mu1prong_ccbar_2018_220622.root","flatTuple_mu1prong_BBbar_2018_160622.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_070722.root","flatTuple_mu1prong_AOD_MadGraph_2018_040722.root","flatTuple_mumu200k_cut6_AOD_2018.root","flatTuple_mumuFSR100k_cut6_AOD_2018.root","flatTuple_mu1prong_AOD_MadGraph_2018_040722.root","ggJpsi_coherent_mu1prong_AOD_STARlight_2018_070722.root","ggJpsi_incoherent_mu1prong_AOD_STARlight_2018_070722.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_080722.root","ggTauTau_SuperChic_mu1prong_2018_080722.root","mumuFSR250k_cut7_2018.root","ggTauTau_SuperChic_mu1prong_2018_080722.root","ggJpsi_coherent_mu1prong_AOD_STARlight_2018_070722.root","ggJpsi_incoherent_mu1prong_AOD_STARlight_2018_070722.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_150722.root","ggTauTau_SuperChic_mu1prong_2018_150722.root","mumuFSR250k_cut7_2018.root","ggTauTau_SuperChic_mu1prong_2018_150722.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_150722.root","ggTauTau_SuperChic_mu1prong_2018_150722.root","mumuFSR250k_cut7_2018.root","ggTauTau_SuperChic_mu1prong_2018_150722.root","ggTauTau_gammaUPC_ChFF_mu1prong_2018_140922.root","ggTauTau_gammaUPC_ChFFkTSmearing_mu1prong_2018_140922.root","ggTauTau_gammaUPC_EDFF_mu1prong_2018_140922.root","ggTauTau_gammaUPC_EDFFkTSmearing_mu1prong_2018_140922.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_150722.root","ggTauTau_gammaUPC_EDFFkTSmearing_mu1prong_2018_140922.root","mumuFSR250k_cut7_2018.root","ggTauTau_SuperChic_mu1prong_2018_150722.root","mu1prong_ccbar_2018_251022.root","mu1prong_BBbar_2018_251022.root"};

//string inputFiles[] = {"mu1prong_data_full2018_100223.root","mu1prong_UPCgen_atau0E-2_2018_260123.root","mu1prong_gammaUPC_mumuFSR_2018_120223.root","mu1prong_UPCgen_atau0E-2_2018_260123.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_150722.root","mu1prong_UPCgen_atau0E-2_2018_260123.root","mu1prong_UPCgen_atau0E-2_2018_260123.root","mu1prong_gammaUPC_mumuFSR_2018_230223.root"};

string inputFiles[] = {"mu1prong_data_subsetZDC2018_280223.root","mu1prong_UPCgen_atau0E-2_2018_260123.root","mu1prong_UPCgen_atau0E-2_2018_260123.root","mu1prong_gammaUPC_mumuFSR_2018_230223.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_150722.root","ggTauTau_SuperChic_mu1prong_2018_150722.root","ggTauTau_gammaUPC_ChFF_mu1prong_2018_140922.root","ggTauTau_gammaUPC_ChFFkTSmearing_mu1prong_2018_140922.root","ggTauTau_gammaUPC_EDFF_mu1prong_2018_140922.root","ggTauTau_gammaUPC_EDFFkTSmearing_mu1prong_2018_140922.root"};

//string inputFiles[] = {"mu1prong_AOD_subdata_2018_150722.root","mu1prong_UPCgen_atau-10E-2_2018_260123.root","mu1prong_UPCgen_atau0E-2_2018_260123.root","mu1prong_UPCgen_atau10E-2_2018_260123.root","ggTauTau_SuperChic_mu1prong_2018_150722.root","ggTauTau_gammaUPC_ChFF_mu1prong_2018_140922.root","ggTauTau_gammaUPC_ChFFkTSmearing_mu1prong_2018_140922.root","ggTauTau_gammaUPC_EDFF_mu1prong_2018_140922.root","ggTauTau_gammaUPC_EDFFkTSmearing_mu1prong_2018_140922.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggCCbar_2015_210521.root","flatTuple_mcggBBbar_2015_060521.root","flatTuple_mcggBBbar_5f_2018_190521.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcgPb_AOD_SuperChic_2015_020721.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_orderedPion_041021.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_orderedPion_041021.root"};


//string inputFiles[] = {"flatTuple_AOD_data_2015_lowMuPtCut_241121.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_lowMuPtCut_241121.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_noTauPtCut_171021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_noTauPtCut_061021.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_lowMuPtCut_241121.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_noTauPtCut_171021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};


//string inputFiles[] = {"flatTuple_AOD_data_2015_noTauPtCut_061021.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_noTauPtCut_061021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_120422.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_noTauPtCut_171021.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_110522.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_100522.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_noTauPtCut_061021.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_190522.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_190522.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_190522.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2prong_2015_300921.root","flatTuple_mcggTauTau_AOD_pp_2pi_MadGraph_2015_280421.root"};

//string inputFiles[] = {"flatTuple_AOD_subdata_2018_090621.root"};
  
//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_data_2015_190421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"flatTuple_AOD_data_2015_190421.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_050521.root","flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_190521.root","flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root","flatTuple_mcggCCbar_2015_300421-1.root","flatTuple_mcggBBbar_2015_060521.root"};

//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_subdata_2018_210421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_pp_MadGraph_2015_270421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_130421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_130421.root"};

//string inputFiles[] = {"/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_AOD_subdata_2018_210421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_pp_2pi_MadGraph_2015_280421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_2pi_MadGraph_2018_270421.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_2pi_SuperChic_2018_270421.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015_220321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_embedded.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018.root", "../ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018.root", "/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggBBbar_4f_AOD_MadGraph_2015.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018_080321.root","/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018_080321.root"};

//string inputFiles[] = {"../ntuples/flatTuple_data_2015.root","../ntuples/flatTuple_mcggTauTau_2015.root","../ntuples/flatTuple_mcggBBbar_2015.root","../ntuples/flatTuple_mcggCCbar_2015.root","../ntuples/flatTuple_mcggTauTau_HydjetDrumMB_2018.root"};

void mu1prong(const int nSamples = 4, string files[] = inputFiles){

  string basePlotDir = "/eos/user/a/ajofrehe/www/gtau/2018/1prong/"+anaTag;
  
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
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_5TeV  = "PbPb - 425 #mub^{-1} ";
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
  //float dataLumi = 1200;  // ub - 2018 sub data
  double PiZeroMass = 134.98; //MeV
  //string tag[] = {"", "Signal 2015 MadGraph", "Signal 2018 MadGraph", "Signal 2018 SuperChic", "Signal 2018 SuperChic + Embedded MB"};
  //string tag[] = {"", "MadGraph 2015", "CCbar 2015", "BBbar 2015", "BBbar 2018"};
  //string tag[] = {"", "signal SC 2018", "mumuFSR", "tautau except signal","signal gUPC ChFF","signal gUPC ChFF kT","signal gUPC EDFF","signal gUPC EDFF kT"};
  //string tags[] = {"", "signal SC 2018", "#gamma#gamma#rightarrow#mu#mu + FSR", "#gamma#gamma#rightarrow#tau#tau background","signal gUPC ChFF","signal gUPC ChFF kT","signal gUPC EDFF","signal gUPC EDFF kT"};
  string tag[] = {"", "signal gUPC 2018", "tautau except signal", "mumuFSR","ccbar","bbbar"};
  string tags[] = {"", "signal gUPC 2018", "#gamma#gamma#rightarrow#tau#tau background", "#gamma#gamma#rightarrow#mu#mu + FSR","ccbar","bbbar"};
  //string tag[] = {"","UPCgen atau-10E-2","UPCgen atau0E-2","UPCgen atau10E-2","SuperChic","gUPC ChFF","gUPC ChFF kT","gUPC EDFF","gUPC EDFF kT"};
  //string tags[] = {"","UPCgen atau-10E-2","UPCgen atau0E-2","UPCgen atau10E-2","SuperChic","gUPC ChFF","gUPC ChFF kT","gUPC EDFF","gUPC EDFF kT"};
  //string tag[] = {"", "signal MG 2018", "mumu", "mumuFSR", "tautau to mu3prong", "tautau to mumu", "ccbar 2018", "bbbar 2018"};
  //string tags[] = {"", "signal MG 2018", "#gamma#gamma#rightarrow#mu#mu", "#gamma#gamma#rightarrow#mu#mu + FSR", "#gamma#gamma#rightarrow#tau#tau#rightarrow#mu+3prong", "#gamma#gamma#rightarrow#tau#tau#rightarrow#mu+#mu", "ccbar 2018", "bbbar 2018"};
  //string tags[] = {"", "#gamma#gamma#rightarrow#tau_{#mu}#tau_{1prong}", "SuperChic 2018", "CCbar 2015", "BBbar 2015"};
  //string tag[] = {"", "MC Signal", "MC 1prong", "MC 1prong + 0#pi^{0}", "MC 1prong + 1#pi^{0}", "MC 1prong + n#pi^{0}", "MC Signal 2018"};
  //string tag[] = {"", "PbPb reco 2015", "pp reco 2015", "pp reco 2018"};
  //string tag[] = {"", "MC Signal 2015", "Signal 2018", "SuperChic 2018"};
  float nEvents[] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
  for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1);
  //nEvents[4] /= 10;
  //for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1)/200;
  
  
  // cross sections (ub): 
  // gammagammatautau   MadGraph              : 570
  // gammagammatautau   SuperChic             : unknown
  // gammagammatautau   gammaUPC    EDFF kT   : 858.93
  // gammagammatautau   gammaUPC    EDFF      : 858.76
  // gammagammatautau   gammaUPC    ChFF kT   : 1060.8959
  // gammagammatautau   gammaUPC    ChFF      : 1060.968
  // gammagammatautau   UPCgen atau   0E-2    : 847.772
  // gammagammatautau   UPCgen atau -10E-2    : 680.211
  // gammagammatautau   UPCgen atau +10E-2    : 1139.995
  // gammagammamumu     gammaUPC    EDFF      : 7060
  // gammagammamumuFSR  gammaUPC    EDFF      : 139.7
  // gammagammaBBbar    MadGraph              : 1.5
  // gammagammaCCbar    MadGraph              : 300
  
  
  
  
  //float crossSectionMC[4] = {570000,1500,300000,570000};
  //float crossSectionMC[] = {570000,300000,1500,1500};
  //float crossSectionMC[] = {570,7060,9.73,300,1.5,570,570,570,570};
  //float crossSectionMC[] = {570,7060,550,570,10000,1000000,300,1.5,570,570,570,570};
  float crossSectionMC[] = {858.93,858.93,139.7,300,1.5,570,570};
  //float crossSectionMC[] = {680.211,847.772,1139.995,570,1060.968,1060.8959,858.76,858.93};
  float SF[] = {1,1,1,1,1,1,1,1,1,1,1,1};
  //float pT_SF = 1677539 / 2000000.0;
  float above3GeVSuperChic = 322461;
  float nAllEventsSuperChic = 2000000;
  //float pT_SF = above3GeVSuperChic / nAllEventsSuperChic;
  float pT_SF = 0.21;
  float tau_pT_SF[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  //float tau_pT_SF[] = {1,pT_SF,pT_SF,pT_SF,pT_SF,pT_SF,1};
  bool threeProng[] = {0,0,0,0,0,0,0,0,0,0,0};
  for (int s = 1; s < nSamples; s++) if(nEvents[s] != 0) SF[s] = dataLumi * crossSectionMC[s-1] / nEvents[s];
  //float SF2[] = {1,0,0,0,0,0,0};
  for (int s = 1; s < nSamples; s++) if(cutflow[s]->GetBinContent(1) != 0) SF[s] = tau_pT_SF[s] * dataLumi * crossSectionMC[s-1] / cutflow[s]->GetBinContent(1);
  //SF[1] = nEvents[0]/nEvents[1];
  //SF[2] = nEvents[0]/nEvents[2];

  int tau_mu_pt_bins = 40;
  //int tau_mu_eta_bins = 13;
  int tau_mu_phi_bins = 13;
  int tau_hadron_pt_bins = 40;
  int tau_hadron_eta_bins = 13;
  int tau_hadron_phi_bins = 13;
  int tau_hadron_rhomass_bins = 12;
  int tau_hadron_nch_bins = 20;
  int tau_hadron_pvz_bins = 50;
  int calo_energy_bins = 125;
  int deltaphi_bins = 60;//80;
  int deltaR_bins = 50;
  int deltaEta_bins = 30;
  
  vector<TH1F**> histograms = std::vector<TH1F**>();
  vector<THStack*> stacks = std::vector<THStack*>();
  vector<TH1F**> ratios = std::vector<TH1F**>();
  vector<string> plotNames = std::vector<string>();
  vector<TCanvas*> ratioCanvas = std::vector<TCanvas*>();
  
  
  TH1F *h_cutflow[nSamples];
  TH1F *h_tau_mu_p[nSamples+5];
  histograms.push_back(h_tau_mu_p);
  plotNames.push_back("tau_muon_p");
  TH1F *h_tau_mu_pz[nSamples+5];
  histograms.push_back(h_tau_mu_pz);
  plotNames.push_back("tau_muon_pz");
  TH1F *h_tau_mu_dz[nSamples+5];
  histograms.push_back(h_tau_mu_dz);
  plotNames.push_back("tau_muon_dz");
  THStack *hs_tau_mu_pt = new THStack("hs_tau_mu_pt","#tau_{#mu} p_{T} [GeV]");
  TH1F *h_tau_mu_pt[nSamples+5];
  histograms.push_back(h_tau_mu_pt);
  plotNames.push_back("tau_muon_pt");
  THStack *hs_tau_mu_eta = new THStack("hs_tau_mu_eta","#tau_{#mu} #eta");
  TH1F *h_tau_mu_eta[nSamples+5];
  histograms.push_back(h_tau_mu_eta);
  plotNames.push_back("tau_muon_eta");
  THStack *hs_tau_mu_theta = new THStack("hs_tau_mu_theta","#tau_{#mu} #theta");
  TH1F *h_tau_mu_theta[nSamples+5];
  histograms.push_back(h_tau_mu_theta);
  plotNames.push_back("tau_muon_theta");
  THStack *hs_tau_mu_phi = new THStack("hs_tau_mu_phi","#tau_{#mu} #phi");
  TH1F *h_tau_mu_phi[nSamples+5];
  histograms.push_back(h_tau_mu_phi);
  plotNames.push_back("tau_muon_phi");
  TH2F *h_tau_mu_tau_hadron_phi[nSamples];
  TH1F *h_tau_hadron_p[nSamples+5];
  histograms.push_back(h_tau_hadron_p);
  plotNames.push_back("tau_hadron_p");
  TH1F *h_tau_hadron_pz[nSamples+5];
  histograms.push_back(h_tau_hadron_pz);
  plotNames.push_back("tau_hadron_pz");
  TH1F *h_tau_hadron_ptScalar[nSamples+5];
  histograms.push_back(h_tau_hadron_ptScalar);
  plotNames.push_back("tau_hadron_scalar_pt");
  THStack *hs_tau_hadron_pt = new THStack("hs_tau_hadron_pt","#tau_{1prong} p_{T} [GeV]");
  TH1F *h_tau_hadron_pt[nSamples+5];
  histograms.push_back(h_tau_hadron_pt);
  plotNames.push_back("tau_hadron_pt");
  TH1F *h_gen_tau_hadron_visible_pt[nSamples];
  THStack *hs_tau_hadron_eta = new THStack("hs_tau_hadron_eta","#tau_{1prong} #eta");
  TH1F *h_tau_hadron_eta[nSamples+5];
  histograms.push_back(h_tau_hadron_eta);
  plotNames.push_back("tau_hadron_eta");
  THStack *hs_tau_hadron_theta = new THStack("hs_tau_hadron_theta","#tau_{1prong} #theta");
  TH1F *h_tau_hadron_theta[nSamples+5];
  histograms.push_back(h_tau_hadron_theta);
  plotNames.push_back("tau_hadron_theta");
  THStack *hs_tau_hadron_phi = new THStack("hs_tau_hadron_phi","#tau_{1prong} #phi");
  TH1F *h_tau_hadron_phi[nSamples+5];
  histograms.push_back(h_tau_hadron_phi);
  plotNames.push_back("tau_hadron_phi");
  THStack *hs_tau_hadron_rhomass[2];
  for (int j=0; j < 2;j++) hs_tau_hadron_rhomass[j] = new THStack(("hs_tau_hadron_rhomass" + to_string(j)).c_str(),("#tau_{1prong} #rho mass" + to_string(j) + "[GeV]").c_str());
  TH1F *h_tau_hadron_rhomass[nSamples][2];
  TH2F *h_tau_hadron_rhomass2D[nSamples];
  THStack *hs_tau_hadron_nch = new THStack("hs_tau_hadron_nch","#tau_{1prong} nch");
  TH1F *h_tau_hadron_nch[nSamples];
  TH1F *h_tau_hadron_nch_highHF[nSamples];
  TH1F *h_tau_hadron_nPions[nSamples];
  THStack *hs_tau_hadron_ncand_final = new THStack("hs_tau_hadron_ncand_final","#tau_{1prong} cands final");
  TH1F *h_tau_hadron_ncand_final[nSamples];
  THStack *hs_tau_hadron_vprob = new THStack("hs_tau_hadron_vprob","#tau_{1prong} vprob (%)");
  TH1F *h_tau_hadron_vprob[nSamples+5];
  histograms.push_back(h_tau_hadron_vprob);
  plotNames.push_back("tau_hadron_vprob");
  TH1F *h_tau_hadron_matched_pt_index[nSamples+5];
  THStack *hs_tau_hadron_mass = new THStack("hs_tau_hadron_mass","#tau_{1prong} mass [GeV]");
  TH1F *h_tau_hadron_mass[nSamples+5];
  histograms.push_back(h_tau_hadron_mass);
  plotNames.push_back("tau_hadron_mass");
  TH1F *h_muon_cosThetaStar[nSamples+5];
  histograms.push_back(h_muon_cosThetaStar);
  plotNames.push_back("muon_cosThetaStar");
  TH1F *h_pion_cosThetaStar[nSamples+5];
  histograms.push_back(h_pion_cosThetaStar);
  plotNames.push_back("pion_cosThetaStar");
  TH1F *h_muon_RF_cosDeltaPhi[nSamples+5];
  histograms.push_back(h_muon_RF_cosDeltaPhi);
  plotNames.push_back("muon_RF_cosDeltaPhi");
  TH1F *h_pion_RF_cosDeltaPhi[nSamples+5];
  histograms.push_back(h_pion_RF_cosDeltaPhi);
  plotNames.push_back("pion_RF_cosDeltaPhi");
  TH2F *h_muon_RF_pz_RF_pz[nSamples+5];
  TH2F *h_muon_pz_pion_pz[nSamples+5];
  TH1F *h_ditau_p[nSamples+5];
  histograms.push_back(h_ditau_p);
  plotNames.push_back("ditau_p");
  TH1F *h_ditau_pz[nSamples+5];
  histograms.push_back(h_ditau_pz);
  plotNames.push_back("ditau_pz");
  TH1F *h_ditau_pt[nSamples+5];
  histograms.push_back(h_ditau_pt);
  plotNames.push_back("ditau_pt");
  TH1F *h_ditau_ptScalar[nSamples+5];
  histograms.push_back(h_ditau_ptScalar);
  plotNames.push_back("ditau_scalar_pt");
  TH1F *h_ditau_HF_deltaphi[nSamples+5];
  histograms.push_back(h_ditau_HF_deltaphi);
  plotNames.push_back("h_ditau_HF_deltaphi");
  THStack *hs_ditau_mass = new THStack("hs_tau_hadron_mass","#tau#tau mass [GeV]");
  TH1F *h_ditau_mass[nSamples+5];
  histograms.push_back(h_ditau_mass);
  plotNames.push_back("ditau_mass");
  THStack *hs_resVisTauPt = new THStack("hs_resTauVis","visible #tau p_{T} resolution [GeV]");
  TH1F *h_resVisTauPt[nSamples+5];
  TH1F *h_resVisTauEta[nSamples+5];
  TH1F *h_resVisTauPhi[nSamples+5];
  THStack *hs_tau_hadron_track_pvz[3];
  for (int j=0; j < 3;j++) hs_tau_hadron_track_pvz[j] = new THStack(("h_tau_hadron_track" + to_string(j) + "_pvz").c_str(),("#tau_{1prong} track" + to_string(j) + " pvz").c_str());
  TH1F *h_tau_hadron_track_pvz[nSamples][3];
  THStack *hs_acoplanarity_tau_mu_tau_hadron = new THStack("hs_acoplanarity_tau_mu_tau_hadron","#alpha(#tau_{#mu}, #tau_{1prong})");
  TH1F *h_acoplanarity_tau_mu_tau_hadron[nSamples+5];
  histograms.push_back(h_acoplanarity_tau_mu_tau_hadron);
  plotNames.push_back("acoplanarity");
  THStack *hs_deltaphi_tau_mu_tau_hadron = new THStack("hs_deltaphi_tau_mu_tau_hadron","#Delta#phi(#tau_{#mu}, #tau_{1prong})");
  TH1F *h_deltaphi_tau_mu_tau_hadron[nSamples+5];
  histograms.push_back(h_deltaphi_tau_mu_tau_hadron);
  plotNames.push_back("delta_phi");
  TH1F *h_deltaphi_tau_mu_tau_hadron_zoomed[nSamples+5];
  histograms.push_back(h_deltaphi_tau_mu_tau_hadron_zoomed);
  plotNames.push_back("delta_phi_zoomed");
  TH1F *h_deltaphi_tau_mu_full_tau_hadron[nSamples];
  
  TH1F *h_PionMuDeltaPhi[nSamples+5];
  histograms.push_back(h_PionMuDeltaPhi);
  plotNames.push_back("PionMuDeltaPhi");
  TH1F *h_PionMuDeltaEta[nSamples+5];
  histograms.push_back(h_PionMuDeltaEta);
  plotNames.push_back("PionMuDeltaEta");
  TH1F *h_PionMuDeltaR[nSamples+5];
  histograms.push_back(h_PionMuDeltaR);
  plotNames.push_back("PionMuDeltaR");
  
  THStack *hs_PV_N = new THStack("hs_PV_N","number of PV");
  TH1F *h_PV_N[nSamples+5];
  TH1F *h_sumZDCplus[nSamples+5];
  histograms.push_back(h_sumZDCplus);
  plotNames.push_back("sumZDCplus");
  TH1F *h_sumZDCminus[nSamples+5];
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
  TH1F *h_calo_leadingHE[nSamples];
  TH1F *h_calo_leadingEE[nSamples];
  TH1F *h_calo_leadingECAL[nSamples];
  TH1F *h_calo_leadingHCAL[nSamples];
  TH1F *h_calo_Et[nSamples+5];
  histograms.push_back(h_calo_Et);
  plotNames.push_back("calo_Et");
  TH2F *h_calo_Et_eta[nSamples];
  TH1F *h_calo_sumEt[nSamples+5];
  histograms.push_back(h_calo_sumEt);
  plotNames.push_back("calo_sumEt");
  TH1F *h_calo_leadingEt[nSamples+5];
  histograms.push_back(h_calo_leadingEt);
  plotNames.push_back("calo_leadingEt");
  TH2F *h_calo_leadingEt_eta[nSamples];
  TH1F *h_calo_E[nSamples+5];
  histograms.push_back(h_calo_E);
  plotNames.push_back("calo_E");
  TH2F *h_calo_E_eta[nSamples];
  TH1F *h_calo_sumE[nSamples+5];
  histograms.push_back(h_calo_sumE);
  plotNames.push_back("calo_sumE");
  TH1F *h_calo_leadingE[nSamples+5];
  histograms.push_back(h_calo_leadingE);
  plotNames.push_back("calo_leadingE");
  TH2F *h_calo_leadingE_eta[nSamples];
  TH2F *h_calo_energy_muon_deltaR[nSamples];
  TH2F *h_calo_energy_pion_deltaR[nSamples];
  TH2F *h_calo_energy_muon_deltaEta[nSamples];
  TH2F *h_calo_energy_pion_deltaEta[nSamples];
  TH1F *h_nCaloTowers[nSamples+5];
  histograms.push_back(h_nCaloTowers);
  plotNames.push_back("nCaloTowers");
  TH1F *h_nTowersMuon[nSamples+5];
  histograms.push_back(h_nTowersMuon);
  plotNames.push_back("nTowersMuon");
  TH1F *h_nTowersPion[nSamples+5];
  histograms.push_back(h_nTowersPion);
  plotNames.push_back("nTowersPion");
  TH1F *h_nTowersECALMuon[nSamples+5];
  histograms.push_back(h_nTowersECALMuon);
  plotNames.push_back("nTowersECALMuon");
  TH1F *h_nTowersECALPion[nSamples+5];
  histograms.push_back(h_nTowersECALPion);
  plotNames.push_back("nTowersECALPion");
  TH1F *h_nTowersECALFSR[nSamples+5];
  histograms.push_back(h_nTowersECALFSR);
  plotNames.push_back("nTowersECALFSR");
  TH1F *h_nTowersHCALMuon[nSamples+5];
  histograms.push_back(h_nTowersHCALMuon);
  plotNames.push_back("nTowersHCALMuon");
  TH1F *h_nTowersHCALPion[nSamples+5];
  histograms.push_back(h_nTowersHCALPion);
  plotNames.push_back("nTowersHCALPion");
  TH1F *h_nTowersHCALFSR[nSamples+5];
  histograms.push_back(h_nTowersHCALFSR);
  plotNames.push_back("nTowersHCALFSR");
  TH1F *h_ratio_ECAL_HCAL_FSR[nSamples+5];
  histograms.push_back(h_ratio_ECAL_HCAL_FSR);
  plotNames.push_back("ratio_ECAL_HCAL_FSR");
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
  TH1F *h_eff_FSR_pt[nSamples];
  
  /*TH1F *h_pion_leading_pt[nSamples+5];
  histograms.push_back(h_pion_leading_pt);
  plotNames.push_back("pion_leading_pt");
  
  TH1F *h_pion_leading_eta[nSamples+5];
  histograms.push_back(h_pion_leading_eta);
  plotNames.push_back("pion_leading_eta");
  
  TH1F *h_pion_leading_phi[nSamples+5];
  histograms.push_back(h_pion_leading_phi);
  plotNames.push_back("pion_leading_phi");*/
  
  TH1F *h_residual_leadingGenGamma_ditau_pt[nSamples];
  //histograms.push_back(h_residual_leadingGenGamma_ditau_pt);
  TH1F *h_leadingGenGamma_ditau_deltaPhi[nSamples];
  //histograms.push_back(h_leadingGenGamma_ditau_deltaPhi);
  TH1F *h_leadingGenGamma_ditau_deltaR[nSamples];
  //histograms.push_back(h_leadingGenGamma_ditau_deltaR);
  
  TH1F *h_residual_leadingRecoGamma_ditau_pt[nSamples+5];
  histograms.push_back(h_residual_leadingRecoGamma_ditau_pt);
  plotNames.push_back("residual_leadingRecoGamma_ditau_pt");
  TH1F *h_leadingRecoGamma_ditau_deltaPhi[nSamples+5];
  //histograms.push_back(h_leadingRecoGamma_ditau_deltaPhi);
  TH1F *h_leadingRecoGamma_ditau_deltaR[nSamples+5];
  //histograms.push_back(h_leadingRecoGamma_ditau_deltaR);
  
  
  TH2F *h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma[nSamples];
  TH2F *h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma[nSamples];
  TH2F *h_residual_leadingRecoGamma_ditau_pt_aco[nSamples];
  TH2F *h_leadingGammaPt_reco_gen[nSamples];
  TH2F *h_ditau_pt_aco[nSamples];
  
  TH1F *h_NrecoGamma[nSamples+5];
  //histograms.push_back(h_NrecoGamma);
  TH1F *h_recoGamma_pt[nSamples+5];
  //histograms.push_back(h_recoGamma_pt);
  TH1F *h_recoGamma_eta[nSamples+5];
  //histograms.push_back(h_recoGamma_eta);
  TH1F *h_recoGamma_phi[nSamples+5];
  //histograms.push_back(h_recoGamma_phi);
  TH1F *h_recoGamma_deltaphi_muon[nSamples+5];
  //histograms.push_back(h_recoGamma_deltaphi_muon);
  TH1F *h_recoGamma_deltaR_muon[nSamples+5];
  //histograms.push_back(h_recoGamma_deltaR_muon);
  TH1F *h_recoGamma_deltaphi_pion[nSamples+5];
  //histograms.push_back(h_recoGamma_deltaphi_pion);
  TH1F *h_recoGamma_deltaR_pion[nSamples+5];
  //histograms.push_back(h_recoGamma_deltaR_pion);
  TH1F *h_recoGammasMinDeltaR[nSamples+5];
  //histograms.push_back(h_recoGammasMinDeltaR);
  TH1F *h_recoPiZeroMinDeltaM[nSamples+5];
  //histograms.push_back(h_recoPiZeroMinDeltaM);
  TH1F *h_recoPiZeroDeltaM[nSamples+5];
  //histograms.push_back(h_recoPiZeroDeltaM);
  TH1F *h_NrecoPiZero[nSamples+5];
  //histograms.push_back(h_NrecoPiZero);
  TH1F *h_recoPiZero_pt[nSamples+5];
  //histograms.push_back(h_recoPiZero_pt);
  TH1F *h_recoPiZero_eta[nSamples+5];
  //histograms.push_back(h_recoPiZero_eta);
  TH1F *h_recoPiZero_phi[nSamples+5];
  //histograms.push_back(h_recoPiZero_phi);
  TH1F *h_recoPiZero_deltaphi_muon[nSamples+5];
  //histograms.push_back(h_recoPiZero_deltaphi_muon);
  TH1F *h_recoPiZero_deltaphi_pion[nSamples+5];
  //histograms.push_back(h_recoPiZero_deltaphi_pion);
  TH2F *h_reco_pion_energy_HCAL_ECAL[nSamples];
  
  TH1F *A_highNch_highHF[nSamples+2+nNchCategories];
  TH1F *B_lowNch_highHF[nSamples+2+nNchCategories];
  TH1F *C_highNch_lowHF[nSamples+2+nNchCategories];
  TH1F *D_lowNch_lowHF[nSamples+2+nNchCategories];
  
  // ABCD validation
  TH1F *ABCD_validation[6];
  ABCD_validation[0] = new TH1F("ABCD_validation_highNch_highHF","ABCD validation highNch highHF;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[1] = new TH1F("ABCD_validation_lowNch_highHF","ABCD validation lowNch highHF;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[2] = new TH1F("ABCD_validation_highNch_lowHF","ABCD validation highNch lowHF;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[3] = new TH1F("ABCD_validation_lowNch_lowHF","ABCD validation lowNch lowHF;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[4] = new TH1F("ABCD_validation_highNch_HF_ratio","ABCD validation highNch HF ratio;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  ABCD_validation[5] = new TH1F("ABCD_validation_lowNch_HF_ratio","ABCD validation lowNch HF ratio;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  
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
  TH1F *h_MET[nSamples+5];
  histograms.push_back(h_MET);
  plotNames.push_back("MET");
  
  TH1F *h_sys_muon_SF_Down = new TH1F("muon-SFDown","muon SF uncertainty down;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  TH1F *h_sys_muon_SF_Up = new TH1F("muon-SFUp","muon SF uncertainty up;#Delta#phi(#tau_{#mu}, #tau_{1prong})", deltaphi_bins, 0, TMath::Pi());
  h_sys_muon_SF_Down->Sumw2(); h_sys_muon_SF_Up->Sumw2();
  
  
  for (int i = 0; i < nSamples; i++){
    
    h_cutflow[i] = new TH1F(("h_cutflow_" + tag[i]).c_str(),("Analysis cutflow - " + tag[i]).c_str(),8, 0, 8);
    std::string cutflow_bins_string[] = {"input from Ntuplizer", "#tau #mu", "nCh", "#tau_{1prong}", "HF", "HF & #tau_{1prong}", "HF & #tau_{1prong} & nch", "..."};
    for(size_t j=0; j< 8; j++){
      h_cutflow[0]->GetXaxis()->SetBinLabel(j+1, (cutflow_bins_string[j]).c_str());
    }
    h_cutflow[i]->Sumw2();
    /*if (i != 0) {*/h_cutflow[i]->SetLineColor(colors[i]); h_cutflow[i]->SetMarkerStyle(styles[i]);
  
    A_highNch_highHF[i] = new TH1F(("A_highNch_highHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - high HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    A_highNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - high HF"); A_highNch_highHF[i]->Sumw2();
    A_highNch_highHF[i]->SetLineColor(colors[i]); A_highNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    B_lowNch_highHF[i] = new TH1F(("B_lowNch_highHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - high HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    B_lowNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - high HF"); B_lowNch_highHF[i]->Sumw2();
    B_lowNch_highHF[i]->SetLineColor(colors[i]); B_lowNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    C_highNch_lowHF[i] = new TH1F(("C_highNch_lowHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - low HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    C_highNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - low HF"); C_highNch_lowHF[i]->Sumw2();
    C_highNch_lowHF[i]->SetLineColor(colors[i]); C_highNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
    D_lowNch_lowHF[i] = new TH1F(("D_lowNch_lowHF_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - low HF " + tag[i]).c_str(), deltaphi_bins, 0, TMath::Pi());
    D_lowNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - low HF"); D_lowNch_lowHF[i]->Sumw2();
    D_lowNch_lowHF[i]->SetLineColor(colors[i]); D_lowNch_lowHF[i]->SetMarkerStyle(styles[i]);
    
    /*h_pion_leading_pt[i] = new TH1F(("h_pion_leading_pt_" + tag[i]).c_str(),("leading pion pt " + tag[i]).c_str(),10, 0, 10*pionLeadingCut);
    h_pion_leading_pt[i]->SetXTitle("leading pion pt [GeV]"); h_pion_leading_pt[i]->Sumw2();
    h_pion_leading_pt[i]->SetLineColor(colors[i]); h_pion_leading_pt[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_pt[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_leading_eta[i] = new TH1F(("h_pion_leading_eta_" + tag[i]).c_str(),("leading pion eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_pion_leading_eta[i]->SetXTitle("leading pion eta"); h_pion_leading_eta[i]->Sumw2();
    h_pion_leading_eta[i]->SetLineColor(colors[i]); h_pion_leading_eta[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_eta[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_pion_leading_phi[i] = new TH1F(("h_pion_leading_phi_" + tag[i]).c_str(),("leading pion phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_pion_leading_phi[i]->SetXTitle("leading pion phi"); h_pion_leading_phi[i]->Sumw2();
    h_pion_leading_phi[i]->SetLineColor(colors[i]); h_pion_leading_phi[i]->SetMarkerStyle(styles[i]);
    h_pion_leading_phi[i]->SetBinErrorOption(TH1::kPoisson);*/
    
    h_residual_leadingGenGamma_ditau_pt[i] = new TH1F(("h_residual_leadingGenGamma_ditau_pt_" + tag[i]).c_str(),("leading gen gamma - ditau p_{T} " + tag[i]).c_str(),49, -6, 6);
    h_residual_leadingGenGamma_ditau_pt[i]->SetXTitle("leading gen #gamma p_{T} - visible ditau p_{T} [GeV]"); h_residual_leadingGenGamma_ditau_pt[i]->Sumw2();
    /*if (i != 0) {*/h_residual_leadingGenGamma_ditau_pt[i]->SetLineColor(colors[i]); h_residual_leadingGenGamma_ditau_pt[i]->SetMarkerStyle(styles[i]);
    
    h_leadingGenGamma_ditau_deltaPhi[i] = new TH1F(("h_leadingGenGamma_ditau_deltaPhi_" + tag[i]).c_str(),("#Delta#Phi(leading gen gamma, ditau)  " + tag[i]).c_str(),120, 0, TMath::Pi());
    h_leadingGenGamma_ditau_deltaPhi[i]->SetXTitle("#Delta#Phi(leading gen #gamma , visible ditau)"); h_leadingGenGamma_ditau_deltaPhi[i]->Sumw2();
    /*if (i != 0) {*/h_leadingGenGamma_ditau_deltaPhi[i]->SetLineColor(colors[i]); h_leadingGenGamma_ditau_deltaPhi[i]->SetMarkerStyle(styles[i]);
    
    h_leadingGenGamma_ditau_deltaR[i] = new TH1F(("h_leadingGenGamma_ditau_deltaR_" + tag[i]).c_str(),("#DeltaR(leading gen gamma, ditau)  " + tag[i]).c_str(),48, -6, 6);
    h_leadingGenGamma_ditau_deltaR[i]->SetXTitle("#DeltaR(leading gen #gamma , visible ditau)"); h_leadingGenGamma_ditau_deltaR[i]->Sumw2();
    /*if (i != 0) {*/h_leadingGenGamma_ditau_deltaR[i]->SetLineColor(colors[i]); h_leadingGenGamma_ditau_deltaR[i]->SetMarkerStyle(styles[i]);
    
    h_residual_leadingRecoGamma_ditau_pt[i] = new TH1F(("h_residual_leadingRecoGamma_ditau_pt_" + tag[i]).c_str(),("leading reco gamma - ditau p_{T} " + tag[i]).c_str(),90, -6, 3);
    h_residual_leadingRecoGamma_ditau_pt[i]->SetXTitle("leading reco #gamma p_{T} - visible ditau p_{T} [GeV]"); h_residual_leadingRecoGamma_ditau_pt[i]->Sumw2();
    /*if (i != 0) {*/h_residual_leadingRecoGamma_ditau_pt[i]->SetLineColor(colors[i]); h_residual_leadingRecoGamma_ditau_pt[i]->SetMarkerStyle(styles[i]);
    
    h_leadingRecoGamma_ditau_deltaPhi[i] = new TH1F(("h_leadingRecoGamma_ditau_deltaPhi_" + tag[i]).c_str(),("#Delta#Phi(leading reco gamma, ditau)  " + tag[i]).c_str(),120, 0, TMath::Pi());
    h_leadingRecoGamma_ditau_deltaPhi[i]->SetXTitle("#Delta#Phi(leading reco #gamma , visible ditau)"); h_leadingRecoGamma_ditau_deltaPhi[i]->Sumw2();
    /*if (i != 0) {*/h_leadingRecoGamma_ditau_deltaPhi[i]->SetLineColor(colors[i]); h_leadingRecoGamma_ditau_deltaPhi[i]->SetMarkerStyle(styles[i]);
    
    h_leadingRecoGamma_ditau_deltaR[i] = new TH1F(("h_leadingRecoGamma_ditau_deltaR_" + tag[i]).c_str(),("#DeltaR(leading reco gamma, ditau)  " + tag[i]).c_str(),48, -6, 6);
    h_leadingRecoGamma_ditau_deltaR[i]->SetXTitle("#DeltaR(leading reco #gamma , visible ditau)"); h_leadingRecoGamma_ditau_deltaR[i]->Sumw2();
    /*if (i != 0) {*/h_leadingRecoGamma_ditau_deltaR[i]->SetLineColor(colors[i]); h_leadingRecoGamma_ditau_deltaR[i]->SetMarkerStyle(styles[i]);
    
    h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma[i] = new TH2F(("h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma_" + tag[i]).c_str(),"leading reco gamma - ditau p_{T} vs leading reco gamma p_{T};leading reco gamma p_{T} [GeV];leading reco #gamma p_{T} - visible ditau p_{T} [GeV]",24,0,6,120, -6, 6);
    
    h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma[i] = new TH2F(("h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma_" + tag[i]).c_str(),"leading reco gamma - ditau p_{T} vs leading reco gamma energy;leading reco gamma energy [GeV];leading reco #gamma energy - visible ditau p_{T} [GeV]",32,0,8,120, -6, 6);
    
    h_residual_leadingRecoGamma_ditau_pt_aco[i] = new TH2F(("h_residual_leadingRecoGamma_ditau_pt_aco_" + tag[i]).c_str(),"leading reco gamma - ditau p_{T} vs #alpha(leading reco #gamma,ditau);#alpha(leading reco #gamma,ditau);leading reco #gamma p_{T} - visible ditau p_{T} [GeV]",50,0,1,100, -5, 5);
    
    h_leadingGammaPt_reco_gen[i] = new TH2F(("h_leadingGammaPt_reco_gen_" + tag[i]).c_str(),"leading #gamma p_{T}: reco vs gen;leading gen #gamma p_{T} [GeV];leading reco #gamma p_{T} [GeV]",40,0,8,40, 0, 8);
    
    h_ditau_pt_aco[i] = new TH2F(("h_ditau_pt_aco_" + tag[i]).c_str(),"ditau p_{T} vs #alpha(#mu,#pi^{#pm});#alpha(#mu,#pi^{#pm});visible ditau p_{T} [GeV]",40,0,1,28, 0, 7);
    
    h_NrecoGamma[i] = new TH1F(("h_NrecoGamma_" + tag[i]).c_str(),("total number of reco gammas " + tag[i]).c_str(),30, -0.5, 29.5);
    h_NrecoGamma[i]->SetXTitle("number of reco gammas"); h_NrecoGamma[i]->Sumw2();
    /*if (i != 0) {*/h_NrecoGamma[i]->SetLineColor(colors[i]); h_NrecoGamma[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_pt[i] = new TH1F(("h_recoGamma_pt_" + tag[i]).c_str(),("gamma reco p_{T} " + tag[i]).c_str(),40, 0, 5);
    h_recoGamma_pt[i]->SetXTitle("gamma reco p_{T} [GeV]"); h_recoGamma_pt[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_pt[i]->SetLineColor(colors[i]); h_recoGamma_pt[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_eta[i] = new TH1F(("h_recoGamma_eta_" + tag[i]).c_str(),("gamma reco eta " + tag[i]).c_str(),29, -3, 3);
    h_recoGamma_eta[i]->SetXTitle("gamma reco eta"); h_recoGamma_eta[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_eta[i]->SetLineColor(colors[i]); h_recoGamma_eta[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_phi[i] = new TH1F(("h_recoGamma_phi_" + tag[i]).c_str(),("gamma reco phi " + tag[i]).c_str(),20, -TMath::Pi(), TMath::Pi());
    h_recoGamma_phi[i]->SetXTitle("gamma reco phi"); h_recoGamma_phi[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_phi[i]->SetLineColor(colors[i]); h_recoGamma_phi[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_deltaphi_muon[i] = new TH1F(("h_recoGamma_deltaphi_muon_" + tag[i]).c_str(),("#Delta#phi(reco #mu, reco gamma) " + tag[i]).c_str(),40, 0, TMath::Pi());
    h_recoGamma_deltaphi_muon[i]->SetXTitle("#Delta#phi(reco #mu, reco gamma)"); h_recoGamma_deltaphi_muon[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_deltaphi_muon[i]->SetLineColor(colors[i]); h_recoGamma_deltaphi_muon[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_deltaR_muon[i] = new TH1F(("h_recoGamma_deltaR_muon_" + tag[i]).c_str(),("#DeltaR(reco #mu, reco gamma) " + tag[i]).c_str(),40, 0, 4);
    h_recoGamma_deltaR_muon[i]->SetXTitle("#DeltaR(reco #mu, reco gamma)"); h_recoGamma_deltaR_muon[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_deltaR_muon[i]->SetLineColor(colors[i]); h_recoGamma_deltaR_muon[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_deltaphi_pion[i] = new TH1F(("h_recoGamma_deltaphi_pion_" + tag[i]).c_str(),("#Delta#phi(reco charged pion, reco gamma) " + tag[i]).c_str(),40, 0, TMath::Pi());
    h_recoGamma_deltaphi_pion[i]->SetXTitle("#Delta#phi(reco charged pion, reco gamma)"); h_recoGamma_deltaphi_pion[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_deltaphi_pion[i]->SetLineColor(colors[i]); h_recoGamma_deltaphi_pion[i]->SetMarkerStyle(styles[i]);
    
    h_recoGamma_deltaR_pion[i] = new TH1F(("h_recoGamma_deltaR_pion_" + tag[i]).c_str(),("#DeltaR(reco charged pion, reco gamma) " + tag[i]).c_str(),40, 0, 4);
    h_recoGamma_deltaR_pion[i]->SetXTitle("#DeltaR(reco charged pion, reco gamma)"); h_recoGamma_deltaR_pion[i]->Sumw2();
    /*if (i != 0) {*/h_recoGamma_deltaR_pion[i]->SetLineColor(colors[i]); h_recoGamma_deltaR_pion[i]->SetMarkerStyle(styles[i]);
    
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
    
    h_reco_pion_energy_HCAL_ECAL[i] = new TH2F(("h_reco_pion_energy_HCAL_ECAL_" + tag[i]).c_str(),("reco charged pion HCAL vs ECAL energy" + tag[i]).c_str(),40, 0, 20,40, 0, 20);
    h_reco_pion_energy_HCAL_ECAL[i]->SetXTitle("ECAL energy"); h_reco_pion_energy_HCAL_ECAL[i]->SetYTitle("HCAL energy");
    
    h_NgenGamma[i] = new TH1F(("h_NgenGamma_" + tag[i]).c_str(),("total number of gen gammas " + tag[i]).c_str(),30, -0.5, 29.5);
    h_NgenGamma[i]->SetXTitle("number of gen gammas"); h_NgenGamma[i]->Sumw2();
    /*if (i != 0) {*/h_NgenGamma[i]->SetLineColor(colors[i]); h_NgenGamma[i]->SetMarkerStyle(styles[i]);
    
    h_genGamma_pt[i] = new TH1F(("h_genGamma_pt_" + tag[i]).c_str(),("gamma gen p_{T} " + tag[i]).c_str(),40, 0, 5);
    h_genGamma_pt[i]->SetXTitle("gamma gen p_{T} [GeV]"); h_genGamma_pt[i]->Sumw2();
    /*if (i != 0) {*/h_genGamma_pt[i]->SetLineColor(colors[i]); h_genGamma_pt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_FSR_pt[i] = new TH1F(("h_eff_FSR_pt_" + tag[i]).c_str(),("FSR efficiency " + tag[i]).c_str(),40, 0, 5);
    h_eff_FSR_pt[i]->SetXTitle("FSR gen p_{T} [GeV]"); h_eff_FSR_pt[i]->SetYTitle("efficiency"); h_eff_FSR_pt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_FSR_pt[i]->SetLineColor(colors[i]); h_eff_FSR_pt[i]->SetMarkerStyle(styles[i]);
    
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
    
    h_N_matched_genPionPt[i] = new TH1F(("h_N_matched_genPionPt_" + tag[i]).c_str(),("Number of matched gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_matched_genPionPt[i]->SetXTitle("matched pion gen p_{T} [GeV]"); h_N_matched_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_matched_genPionPt[i]->SetLineColor(colors[i]); h_N_matched_genPionPt[i]->SetMarkerStyle(styles[i]);
    
    h_N_matched_leading_genPionPt[i] = new TH1F(("h_N_matched_leading_genPionPt_" + tag[i]).c_str(),("Number of matched leading gen pions - " + tag[i]).c_str(),40, -0.1, 7.9);
    h_N_matched_leading_genPionPt[i]->SetXTitle("leading matched pion gen p_{T} [GeV]"); h_N_matched_leading_genPionPt[i]->Sumw2();
    /*if (i != 0) {*/h_N_matched_leading_genPionPt[i]->SetLineColor(colors[i]); h_N_matched_leading_genPionPt[i]->SetMarkerStyle(styles[i]);
    
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
    h_N_genTauHadpt[i]->SetXTitle("#tau_{1prong} gen p_{T} [GeV]"); h_N_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_N_genTauHadpt[i]->SetLineColor(colors[i]); h_N_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_3prongReco_genTauHadpt[i] = new TH1F(("h_eff_3prongReco_genTauHadpt_" + tag[i]).c_str(),("Efficiency of reconstructing 3 pions - " + tag[i]).c_str(),20, 0, 20);
    h_eff_3prongReco_genTauHadpt[i]->SetXTitle("#tau_{1prong} gen p_{T} [GeV]");
    h_eff_3prongReco_genTauHadpt[i]->SetYTitle("Efficiency of reconstructing 3 pions");
    h_eff_3prongReco_genTauHadpt[i]->Sumw2();
    /*if (i != 0) {*/h_eff_3prongReco_genTauHadpt[i]->SetLineColor(colors[i]); h_eff_3prongReco_genTauHadpt[i]->SetMarkerStyle(styles[i]);
    
    h_eff_tauReco_genTauHadpt[i] = new TH1F(("h_eff_tauReco_genTauHadpt_" + tag[i]).c_str(),("Efficiency of reconstructing a hadronic tau - " + tag[i]).c_str(),20, 0, 20);
    h_eff_tauReco_genTauHadpt[i]->SetXTitle("#tau_{1prong} gen p_{T} [GeV]");
    h_eff_tauReco_genTauHadpt[i]->SetYTitle("#tau_{1prong} reconstruction efficiency");
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
        
    h_tau_mu_eta[i] = new TH1F(("h_tau_mu_eta_" + tag[i]).c_str(),("visible #tau_{#mu} #eta " + tag[i]).c_str(),16, -2.4, 2.4);
    h_tau_mu_eta[i]->SetXTitle("visible #tau_{#mu} #eta"); h_tau_mu_eta[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_eta[i]->SetLineColor(colors[i]); h_tau_mu_eta[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_eta[i]->SetBinErrorOption(TH1::kPoisson);
        
    h_tau_mu_theta[i] = new TH1F(("h_tau_mu_theta_" + tag[i]).c_str(),("visible #tau_{#mu} #theta " + tag[i]).c_str(),16, 0, TMath::Pi());
    h_tau_mu_theta[i]->SetXTitle("visible #tau_{#mu} #theta"); h_tau_mu_theta[i]->Sumw2();
    h_tau_mu_theta[i]->SetLineColor(colors[i]); h_tau_mu_theta[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_theta[i]->SetBinErrorOption(TH1::kPoisson);
        
    h_tau_mu_phi[i] = new TH1F(("h_tau_mu_phi_" + tag[i]).c_str(),("visible #tau_{#mu} #phi " + tag[i]).c_str(),tau_mu_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_mu_phi[i]->SetXTitle("visible #tau_{#mu} #phi"); h_tau_mu_phi[i]->Sumw2();
    /*if (i != 0) {*/h_tau_mu_phi[i]->SetLineColor(colors[i]); h_tau_mu_phi[i]->SetMarkerStyle(styles[i]);
    h_tau_mu_phi[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_tau_mu_tau_hadron_phi[i] = new TH2F(("h_tau_mu_tau_hadron_phi_" + tag[i]).c_str(),("visible #tau_{#mu}-#tau_{1prong} #phi " + tag[i] +";visible #tau_{#mu} #phi;visible #tau_{1prong} #phi").c_str(),tau_mu_phi_bins, -TMath::Pi(), TMath::Pi(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    
    // hadron related
    
    h_tau_hadron_p[i] = new TH1F(("h_tau_hadron_p_" + tag[i]).c_str(),("#tau_{1prong} vector sum p " + tag[i]).c_str(),10, 0, 20);
    h_tau_hadron_p[i]->SetXTitle("#tau_{1prong} vector sum p [GeV]"); h_tau_hadron_p[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_p[i]->SetLineColor(colors[i]); h_tau_hadron_p[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_pz[i] = new TH1F(("h_tau_hadron_pz_" + tag[i]).c_str(),("#tau_{1prong} p_{z} " + tag[i]).c_str(),20, -10, 20);
    h_tau_hadron_pz[i]->SetXTitle("#tau_{1prong} p_{z} [GeV]"); h_tau_hadron_pz[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_pz[i]->SetLineColor(colors[i]); h_tau_hadron_pz[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_ptScalar[i] = new TH1F(("h_tau_hadron_ptScalar_" + tag[i]).c_str(),("3prong scalar sum p_{T} " + tag[i]).c_str(),10, 0, 20);
    h_tau_hadron_ptScalar[i]->SetXTitle("3prong scalar sum p_{T} [GeV]"); h_tau_hadron_ptScalar[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_ptScalar[i]->SetLineColor(colors[i]); h_tau_hadron_ptScalar[i]->SetMarkerStyle(styles[i]);
    
    h_tau_hadron_pt[i] = new TH1F(("h_tau_hadron_pt_" + tag[i]).c_str(),("#tau_{1prong} p_{T} " + tag[i]).c_str(),21, 0, 10.5);
    h_tau_hadron_pt[i]->SetXTitle("#tau_{1prong} p_{T} [GeV]"); h_tau_hadron_pt[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_pt[i]->SetLineColor(colors[i]); h_tau_hadron_pt[i]->SetMarkerStyle(styles[i]);
    
    h_gen_tau_hadron_visible_pt[i] = new TH1F(("h_gen_tau_hadron_visible_pt_" + tag[i]).c_str(),("gen #tau_{1prong} visible p_{T} " + tag[i]).c_str(),20, 0, 20);
    h_gen_tau_hadron_visible_pt[i]->SetXTitle("gen #tau_{1prong} p_{T} [GeV]"); h_gen_tau_hadron_visible_pt[i]->Sumw2();
    /*if (i != 0) {*/h_gen_tau_hadron_visible_pt[i]->SetLineColor(colors[i]); h_gen_tau_hadron_visible_pt[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_eta[i] = new TH1F(("h_tau_hadron_eta_" + tag[i]).c_str(),("#tau_{1prong} #eta " + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_tau_hadron_eta[i]->SetXTitle("#tau_{1prong} #eta"); h_tau_hadron_eta[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_eta[i]->SetLineColor(colors[i]); h_tau_hadron_eta[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_theta[i] = new TH1F(("h_tau_hadron_theta_" + tag[i]).c_str(),("#tau_{1prong} #theta " + tag[i]).c_str(),16, 0, TMath::Pi());
    h_tau_hadron_theta[i]->SetXTitle("#tau_{1prong} #theta"); h_tau_hadron_theta[i]->Sumw2();
    h_tau_hadron_theta[i]->SetLineColor(colors[i]); h_tau_hadron_theta[i]->SetMarkerStyle(styles[i]);
        
    h_tau_hadron_phi[i] = new TH1F(("h_tau_hadron_phi_" + tag[i]).c_str(),("#tau_{1prong} #phi " + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_hadron_phi[i]->SetXTitle("#tau_{1prong} #phi"); h_tau_hadron_phi[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_phi[i]->SetLineColor(colors[i]); h_tau_hadron_phi[i]->SetMarkerStyle(styles[i]);
        
    for (int j=0; j < 2;j++){
      h_tau_hadron_rhomass[i][j] = new TH1F(("h_tau_hadron_rhomass" + to_string(j) + "_" + tag[i]).c_str(),("#tau_{1prong} #rho_{"+to_string(j+1)+"} mass [GeV] " + tag[i]).c_str(),tau_hadron_rhomass_bins, 0.2, 1.4);
      h_tau_hadron_rhomass[i][j]->SetXTitle(("#rho_{"+to_string(j+1)+"} mass [GeV]").c_str()); h_tau_hadron_rhomass[i][j]->Sumw2();
      /*if (i != 0) {*/h_tau_hadron_rhomass[i][j]->SetLineColor(colors[i]); h_tau_hadron_rhomass[i][j]->SetMarkerStyle(styles[i]);
    }
    
    h_tau_hadron_rhomass2D[i] = new TH2F(("h_tau_hadron_rhomass2D_" + tag[i]).c_str(),("#tau_{1prong} #rho mass [GeV] " + tag[i]).c_str(),2*tau_hadron_rhomass_bins, 0.2, 1.4, 2*tau_hadron_rhomass_bins, 0.2, 1.4);
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
    
    h_tau_hadron_mass[i] = new TH1F(("h_tau_hadron_mass_" + tag[i]).c_str(),("#tau_{1prong} visible mass " + tag[i]).c_str(), 16, 0.1, 1.7);
    h_tau_hadron_mass[i]->SetXTitle("visible #tau_{1prong} mass [GeV]"); h_tau_hadron_mass[i]->Sumw2();
    /*if (i != 0) {*/h_tau_hadron_mass[i]->SetLineColor(colors[i]); h_tau_hadron_mass[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_mass[i] = new TH1F(("h_ditau_mass_" + tag[i]).c_str(),("#tau_{#mu}#tau_{1prong} visible invariant mass " + tag[i]).c_str(), 36, 2, 14);
    h_ditau_mass[i]->SetXTitle("visible #tau_{#mu}#tau_{1prong} invariant mass [GeV]"); h_ditau_mass[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_mass[i]->SetLineColor(colors[i]); h_ditau_mass[i]->SetMarkerStyle(styles[i]);
    
    h_muon_cosThetaStar[i] = new TH1F(("h_muon_cosThetaStar_" + tag[i]).c_str(),("muon cos(#theta^{*}) " + tag[i]).c_str(), 80, -1, 1);
    h_muon_cosThetaStar[i]->SetXTitle("muon cos(#theta^{*})"); h_muon_cosThetaStar[i]->Sumw2();
    h_muon_cosThetaStar[i]->SetLineColor(colors[i]); h_muon_cosThetaStar[i]->SetMarkerStyle(styles[i]);
    
    h_pion_cosThetaStar[i] = new TH1F(("h_pion_cosThetaStar_" + tag[i]).c_str(),("pion cos(#theta^{*}) " + tag[i]).c_str(), 80, -1, 1);
    h_pion_cosThetaStar[i]->SetXTitle("pion cos(#theta^{*})"); h_pion_cosThetaStar[i]->Sumw2();
    h_pion_cosThetaStar[i]->SetLineColor(colors[i]); h_pion_cosThetaStar[i]->SetMarkerStyle(styles[i]);
    
    h_muon_RF_cosDeltaPhi[i] = new TH1F(("h_muon_RF_cosDeltaPhi_" + tag[i]).c_str(),("cos(#Delta#phi(#mu^{RF},RF)) " + tag[i]).c_str(), 80, -1, 1);
    h_muon_RF_cosDeltaPhi[i]->SetXTitle("cos(#Delta#phi(#mu^{RF},RF))"); h_muon_RF_cosDeltaPhi[i]->Sumw2();
    h_muon_RF_cosDeltaPhi[i]->SetLineColor(colors[i]); h_muon_RF_cosDeltaPhi[i]->SetMarkerStyle(styles[i]);
    
    h_pion_RF_cosDeltaPhi[i] = new TH1F(("h_pion_RF_cosDeltaPhi_" + tag[i]).c_str(),("cos(#Delta#phi(#pi^{RF},RF)) " + tag[i]).c_str(), 80, -1, 1);
    h_pion_RF_cosDeltaPhi[i]->SetXTitle("cos(#Delta#phi(#pi^{RF},RF))"); h_pion_RF_cosDeltaPhi[i]->Sumw2();
    h_pion_RF_cosDeltaPhi[i]->SetLineColor(colors[i]); h_pion_RF_cosDeltaPhi[i]->SetMarkerStyle(styles[i]);
    
    h_muon_RF_pz_RF_pz[i] = new TH2F(("h_muon_RF_pz_RF_pz_" + tag[i]).c_str(),("#mu^{RF} p_{z} vs RF p_{z} " + tag[i]).c_str(),40, -10, 10, 40, -10, 10);
    h_muon_RF_pz_RF_pz[i]->SetXTitle("RF p_{z} [GeV]"); h_muon_RF_pz_RF_pz[i]->SetYTitle("#mu^{RF} p_{z} [GeV]");
    
    h_muon_pz_pion_pz[i] = new TH2F(("h_muon_pz_pion_pz_" + tag[i]).c_str(),("p^{#mu}_{z} vs p^{#pi}_{z} " + tag[i]).c_str(),52, -13, 13, 52, -13, 13);
    h_muon_pz_pion_pz[i]->SetXTitle("p^{#pi}_{z} [GeV]"); h_muon_pz_pion_pz[i]->SetYTitle("p^{#mu}_{z} [GeV]");
    
    h_ditau_p[i] = new TH1F(("h_ditau_p_" + tag[i]).c_str(),("#tau_{#mu}#tau_{1prong} visible total p " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_p[i]->SetXTitle("#tau_{#mu}#tau_{1prong} total p [GeV]"); h_ditau_p[i]->Sumw2();
    h_ditau_p[i]->SetLineColor(colors[i]); h_ditau_p[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_pz[i] = new TH1F(("h_ditau_pz_" + tag[i]).c_str(),("#tau_{#mu}#tau_{1prong} visible p_{z} " + tag[i]).c_str(), 40, -40, 40);
    h_ditau_pz[i]->SetXTitle("visible #tau_{#mu}#tau_{1prong} p_{z} [GeV]"); h_ditau_pz[i]->Sumw2();
    h_ditau_pz[i]->SetLineColor(colors[i]); h_ditau_pz[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_pt[i] = new TH1F(("h_ditau_pt_" + tag[i]).c_str(),("#tau_{#mu}#tau_{1prong} visible p_{T} " + tag[i]).c_str(), 42, 0., 10.5);
    h_ditau_pt[i]->SetXTitle("visible #tau_{#mu}#tau_{1prong} p_{T} [GeV]"); h_ditau_pt[i]->Sumw2();
    h_ditau_pt[i]->SetLineColor(colors[i]); h_ditau_pt[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_ptScalar[i] = new TH1F(("h_ditau_ptScalar_" + tag[i]).c_str(),("#mu + 1prong scalar sum p_{T} " + tag[i]).c_str(), 20, 0, 10);
    h_ditau_ptScalar[i]->SetXTitle("#mu + 1prong scalar sum p_{T} [GeV]"); h_ditau_ptScalar[i]->Sumw2();
    //h_ditau_ptScalar[i]->SetXTitle("scalar p_{T}^{#mu} - p_{T}^{1prong} [GeV]"); h_ditau_ptScalar[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_ptScalar[i]->SetLineColor(colors[i]); h_ditau_ptScalar[i]->SetMarkerStyle(styles[i]);
    
    h_ditau_HF_deltaphi[i] = new TH1F(("h_ditau_HF_deltaphi_" + tag[i]).c_str(),("#Delta#phi(leading HF, visible ditau) " + tag[i]).c_str(), 100, 0.0*TMath::Pi(), TMath::Pi());
    h_ditau_HF_deltaphi[i]->SetXTitle("#Delta#phi(leading HF, visible ditau)"); h_ditau_HF_deltaphi[i]->Sumw2();
    /*if (i != 0) {*/h_ditau_HF_deltaphi[i]->SetLineColor(colors[i]); h_ditau_HF_deltaphi[i]->SetMarkerStyle(styles[i]);
    
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
        
    //h_deltaphi_tau_mu_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) " + tag[i]).c_str(),deltaphi_bins/4, 0.75*TMath::Pi(), TMath::Pi());
    h_acoplanarity_tau_mu_tau_hadron[i] = new TH1F(("h_acoplanarity_tau_mu_tau_hadron_" + tag[i]).c_str(),("#alpha(#tau_{#mu}, #tau_{1prong}) " + tag[i]).c_str(),deltaphi_bins, 0.0, 1);
    h_acoplanarity_tau_mu_tau_hadron[i]->SetXTitle("#alpha(#tau_{#mu}, #tau_{1prong})"); h_acoplanarity_tau_mu_tau_hadron[i]->Sumw2();
    h_acoplanarity_tau_mu_tau_hadron[i]->SetLineColor(colors[i]); h_acoplanarity_tau_mu_tau_hadron[i]->SetMarkerStyle(styles[i]);
    h_acoplanarity_tau_mu_tau_hadron[i]->SetBinErrorOption(TH1::kPoisson);
    
    h_deltaphi_tau_mu_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) " + tag[i]).c_str(),deltaphi_bins, 0.0*TMath::Pi(), TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong})"); h_deltaphi_tau_mu_tau_hadron[i]->Sumw2();
    h_deltaphi_tau_mu_tau_hadron[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_tau_hadron[i]->SetMarkerStyle(styles[i]);
    h_deltaphi_tau_mu_tau_hadron[i]->SetBinErrorOption(TH1::kPoisson);
    /*h_deltaphi_tau_mu_tau_hadron[i]->SetFillColor(i);*/
        
    h_deltaphi_tau_mu_tau_hadron_zoomed[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_zoomed_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) " + tag[i]).c_str(),10, 70*TMath::Pi()/80, TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong})"); h_deltaphi_tau_mu_tau_hadron_zoomed[i]->Sumw2();
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetMarkerStyle(styles[i]);
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetBinErrorOption(TH1::kPoisson);
    /*h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetFillColor(i);*/
        
    h_deltaphi_tau_mu_full_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_full_tau_hadron_" + tag[i]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}+#pi^{0}(s)) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}+#pi^{0}(s))"); h_deltaphi_tau_mu_full_tau_hadron[i]->Sumw2();
    h_deltaphi_tau_mu_full_tau_hadron[i]->SetLineColor(colors[i]); h_deltaphi_tau_mu_full_tau_hadron[i]->SetMarkerStyle(styles[i]);
    /*h_deltaphi_tau_mu_full_tau_hadron[i]->SetFillColor(i);*/
        
        
    h_PionMuDeltaPhi[i] = new TH1F(("h_PionMuDeltaPhi_" + tag[i]).c_str(),("#Delta#phi(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_PionMuDeltaPhi[i]->SetXTitle("#Delta#phi(#mu,#pi^{#pm})"); h_PionMuDeltaPhi[i]->Sumw2();
    h_PionMuDeltaPhi[i]->SetLineColor(colors[i]); h_PionMuDeltaPhi[i]->SetMarkerStyle(styles[i]);
        
    h_PionMuDeltaEta[i] = new TH1F(("h_PionMuDeltaEta_" + tag[i]).c_str(),("#Delta#eta(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaEta_bins, 0, 5);
    h_PionMuDeltaEta[i]->SetXTitle("#Delta#eta(#mu,#pi^{#pm})"); h_PionMuDeltaEta[i]->Sumw2();
    h_PionMuDeltaEta[i]->SetLineColor(colors[i]); h_PionMuDeltaEta[i]->SetMarkerStyle(styles[i]);
        
    h_PionMuDeltaR[i] = new TH1F(("h_PionMuDeltaR_" + tag[i]).c_str(),("#DeltaR(#mu,#pi^{#pm}) " + tag[i]).c_str(),deltaR_bins, 0, 5);
    h_PionMuDeltaR[i]->SetXTitle("#DeltaR(#mu,#pi^{#pm})"); h_PionMuDeltaR[i]->Sumw2();
    h_PionMuDeltaR[i]->SetLineColor(colors[i]); h_PionMuDeltaR[i]->SetMarkerStyle(styles[i]);
    
    h_deltaphi_tau_mu_tau_hadron_mueta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_mueta_" + tag[i]).c_str(),("#eta_{#mu} vs #Delta#phi(#tau_{#mu}, #tau_{1prong}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong})"); h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetYTitle("#eta_{#mu}");
    
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_deltaeta_" + tag[i]).c_str(),("#Delta(abs(#eta))_{#mu_#tau} vs #Delta#phi(#tau_{#mu}, #tau_{1prong}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong})"); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetYTitle("abs(#tau_{1prong} #eta) - abs(#mu #eta)");
    
    h_mueta_taueta[i] = new TH2F(("h_mueta_taueta_" + tag[i]).c_str(),("#tau_{1prong} #eta vs #tau_{#mu} #eta " + tag[i]).c_str(),18,-3,3, 18,-3,3);
    h_mueta_taueta[i]->SetXTitle("#tau_{1prong} #eta"); h_mueta_taueta[i]->SetYTitle("#tau_{#mu} #eta");
    
    h_PV_N[i] = new TH1F(("h_PV_N_" + tag[i]).c_str(),("h_PV_N_" + tag[i]).c_str(),5, -.5, 4.5);
    h_PV_N[i]->SetXTitle("N_{PV}"); h_PV_N[i]->Sumw2();
    /*if (i != 0) {*/h_PV_N[i]->SetLineColor(colors[i]); h_PV_N[i]->SetMarkerStyle(styles[i]);
    
    h_AP[i] = new TH2F(("h_AP_" + tag[i]).c_str(),("Armenteros-Podolansky " + tag[i]).c_str(),40, -1, 1, 40,0,0.8);
    h_AP[i]->SetXTitle("#alpha"); h_AP[i]->SetYTitle("#epsilon");
    
    // ZDC
    
    h_sumZDCplus[i] = new TH1F(("h_sumZDCplus_" + tag[i]).c_str(),("sum ZDC+ " + tag[i]).c_str(), 80, 0, 40);
    h_sumZDCplus[i]->SetXTitle("sum ZDC+ / 1000"); h_sumZDCplus[i]->Sumw2();
    /*if (i != 0) {*/h_sumZDCplus[i]->SetLineColor(colors[i]); //h_sumZDCplus[i]->SetMarkerStyle(styles[i]);
    
    h_sumZDCminus[i] = new TH1F(("h_sumZDCminus_" + tag[i]).c_str(),("sum ZDC- " + tag[i]).c_str(), 80, 0, 40);
    h_sumZDCminus[i]->SetXTitle("sum ZDC- / 1000"); h_sumZDCminus[i]->Sumw2();
    /*if (i != 0) {*/h_sumZDCminus[i]->SetLineColor(colors[i]); //h_sumZDCminus[i]->SetMarkerStyle(styles[i]);
    
    h_sumZDC_pm[i] = new TH2F(("h_sumZDC_pm_" + tag[i]).c_str(),("sum ZDC+ vs ZDC- " + tag[i]).c_str(),40, 0, 20, 40, 0, 20);
    h_sumZDC_pm[i]->SetXTitle("sum ZDC- / 1000"); h_sumZDC_pm[i]->SetYTitle("sum ZDC+ / 1000");
    
    h_averageZDCside_averageHFeta[i] = new TH2F(("h_averageZDCside_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs average ZDC side " + tag[i]).c_str(),21, -1, 1, 20, -5, 5);
    h_averageZDCside_averageHFeta[i]->SetXTitle("average ZDC side"); h_averageZDCside_averageHFeta[i]->SetYTitle("average HF #eta");
    
    h_muEta_averageZDCside[i] = new TH2F(("h_muEta_averageZDCside_" + tag[i]).c_str(),("average ZDC side vs muon #eta " + tag[i]).c_str(),19, -2.5, 2.5, 21, -1, 1);
    h_muEta_averageZDCside[i]->SetXTitle("muon #eta"); h_muEta_averageZDCside[i]->SetYTitle("average ZDC side");
    
    h_tauEta_averageZDCside[i] = new TH2F(("h_tauEta_averageZDCside_" + tag[i]).c_str(),("average ZDC side vs #tau_{1prong} #eta " + tag[i]).c_str(),19, -2.5, 2.5, 21, -1, 1);
    h_tauEta_averageZDCside[i]->SetXTitle("#tau_{1prong} #eta"); h_tauEta_averageZDCside[i]->SetYTitle("average ZDC side");
    
    //HF
    
    h_calo_energyHFp[i] = new TH1F(("h_calo_energyHFp_" + tag[i]).c_str(),("energy HF+ " + tag[i]).c_str(),40, -0.5, 19.5);
    h_calo_energyHFp[i]->SetXTitle("energy HF+ [GeV]"); h_calo_energyHFp[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFp[i]->SetLineColor(colors[i]); h_calo_energyHFp[i]->SetMarkerStyle(styles[i]);
    
    h_calo_energyHFm[i] = new TH1F(("h_calo_energyHFm_" + tag[i]).c_str(),("energy HF- " + tag[i]).c_str(),40, -0.5, 19.5);
    h_calo_energyHFm[i]->SetXTitle("energy HF- [GeV]"); h_calo_energyHFm[i]->Sumw2();
    /*if (i != 0) {*/h_calo_energyHFm[i]->SetLineColor(colors[i]); h_calo_energyHFm[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFp[i] = new TH1F(("h_calo_leadingHFp_" + tag[i]).c_str(),("energy leading tower HF+ " + tag[i]).c_str(),48, 0, 16);
    h_calo_leadingHFp[i]->SetXTitle("leading tower energy HF+ [GeV]"); h_calo_leadingHFp[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFp[i]->SetLineColor(colors[i]); h_calo_leadingHFp[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHE[i] = new TH1F(("h_calo_leadingHE_" + tag[i]).c_str(),("energy leading tower HE " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHE[i]->SetXTitle("leading tower energy HE [GeV]"); h_calo_leadingHE[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHE[i]->SetLineColor(colors[i]); h_calo_leadingHE[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHCAL[i] = new TH1F(("h_calo_leadingHCAL_" + tag[i]).c_str(),("energy leading tower HE " + tag[i]).c_str(),30, 0, 9);
    h_calo_leadingHCAL[i]->SetXTitle("leading HCAL tower energy [GeV]"); h_calo_leadingHCAL[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHCAL[i]->SetLineColor(colors[i]); h_calo_leadingHCAL[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingEE[i] = new TH1F(("h_calo_leadingEE_" + tag[i]).c_str(),("energy leading tower HE " + tag[i]).c_str(),60, 0, 18);
    h_calo_leadingEE[i]->SetXTitle("leading tower energy EE [GeV]"); h_calo_leadingEE[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingEE[i]->SetLineColor(colors[i]); h_calo_leadingEE[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingECAL[i] = new TH1F(("h_calo_leadingECAL_" + tag[i]).c_str(),("energy leading tower HE " + tag[i]).c_str(),60, 0, 18);
    h_calo_leadingECAL[i]->SetXTitle("leading ECAL tower energy [GeV]"); h_calo_leadingECAL[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingECAL[i]->SetLineColor(colors[i]); h_calo_leadingECAL[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFp_highNch[i] = new TH1F(("h_calo_leadingHFp_highNch_" + tag[i]).c_str(),("energy leading tower HF+ " + tag[i]).c_str(),48, 0, 16);
    h_calo_leadingHFp_highNch[i]->SetXTitle("leading tower energy HF+ [GeV]"); h_calo_leadingHFp_highNch[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFp_highNch[i]->SetLineColor(colors[i]); h_calo_leadingHFp_highNch[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFm[i] = new TH1F(("h_calo_leadingHFm_" + tag[i]).c_str(),("energy leading tower HF- " + tag[i]).c_str(),48, 0, 16);
    h_calo_leadingHFm[i]->SetXTitle("leading tower energy HF- [GeV]"); h_calo_leadingHFm[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFm[i]->SetLineColor(colors[i]); h_calo_leadingHFm[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingHFm_highNch[i] = new TH1F(("h_calo_leadingHFm_highNch_" + tag[i]).c_str(),("energy leading tower HF- " + tag[i]).c_str(),48, 0, 16);
    h_calo_leadingHFm_highNch[i]->SetXTitle("leading tower energy HF- [GeV]"); h_calo_leadingHFm_highNch[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingHFm_highNch[i]->SetLineColor(colors[i]); h_calo_leadingHFm_highNch[i]->SetMarkerStyle(styles[i]);
    
    h_calo_Et[i] = new TH1F(("h_calo_Et_" + tag[i]).c_str(),("calo E_{T} " + tag[i]).c_str(),20, 0, 20);
    h_calo_Et[i]->SetXTitle("calo E_{T} [GeV]"); h_calo_Et[i]->Sumw2();
    /*if (i != 0) {*/h_calo_Et[i]->SetLineColor(colors[i]); h_calo_Et[i]->SetMarkerStyle(styles[i]);
    
    h_calo_Et_eta[i] = new TH2F(("h_calo_Et_eta_" + tag[i]).c_str(),("calo E_{T} vs #eta " + tag[i] +";#eta;calo E_{T} [GeV]").c_str(),50, -5, 5, 36, -1, 71);
    
    h_calo_leadingEt_eta[i] = new TH2F(("h_calo_leadingEt_eta_" + tag[i]).c_str(),("calo leading E_{T} vs #eta " + tag[i] +";#eta of leading calo tower;E_{T} of leading calo tower [GeV]").c_str(),50, -5, 5, 26, -1, 51);
    
    h_calo_sumEt[i] = new TH1F(("h_calo_sumEt_" + tag[i]).c_str(),("calo E_{T} sum " + tag[i]).c_str(),40, 0, 80);
    h_calo_sumEt[i]->SetXTitle("calo E_{T} sum [GeV]"); h_calo_sumEt[i]->Sumw2();
    /*if (i != 0) {*/h_calo_sumEt[i]->SetLineColor(colors[i]); h_calo_sumEt[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingEt[i] = new TH1F(("h_calo_leadingEt_" + tag[i]).c_str(),("leading calo E_{T} " + tag[i]).c_str(),10, 0, 10);
    h_calo_leadingEt[i]->SetXTitle("leading calo E_{T} [GeV]"); h_calo_leadingEt[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingEt[i]->SetLineColor(colors[i]); h_calo_leadingEt[i]->SetMarkerStyle(styles[i]);
    
    h_calo_E[i] = new TH1F(("h_calo_E_" + tag[i]).c_str(),("calo energy " + tag[i]).c_str(),50, -1, 99);
    h_calo_E[i]->SetXTitle("calo energy [GeV]"); h_calo_E[i]->Sumw2();
    /*if (i != 0) {*/h_calo_E[i]->SetLineColor(colors[i]); h_calo_E[i]->SetMarkerStyle(styles[i]);
    
    h_calo_E_eta[i] = new TH2F(("h_calo_E_eta_" + tag[i]).c_str(),("calo energy vs #eta " + tag[i] +";#eta;calo energy [GeV]").c_str(),50, -5, 5, 36, -1, 71);
    
    h_calo_leadingE_eta[i] = new TH2F(("h_calo_leadingE_eta_" + tag[i]).c_str(),("calo leading energy vs #eta " + tag[i] +";#eta of leading calo tower;energy of leading calo tower [GeV]").c_str(),50, -5, 5, 26, -1, 51);
    
    
    h_calo_energy_muon_deltaR[i] = new TH2F(("h_calo_energy_muon_deltaR_" + tag[i]).c_str(),("calo energy vs #DeltaR(#mu,calo) " + tag[i] +";#DeltaR(#mu,calo);energy of calo tower [GeV]").c_str(),50, 0, 5, 40, 0, 10);
    
    
    h_calo_energy_pion_deltaR[i] = new TH2F(("h_calo_energy_pion_deltaR_" + tag[i]).c_str(),("calo energy vs #DeltaR(#pi,calo) " + tag[i] +";#DeltaR(#pi,calo);energy of calo tower [GeV]").c_str(),50, 0, 5, 40, 0, 10);
    
    
    h_calo_energy_muon_deltaEta[i] = new TH2F(("h_calo_energy_muon_deltaEta_" + tag[i]).c_str(),("calo energy vs #Delta#eta(#mu,calo) " + tag[i] +";#Delta#eta(#mu,calo);energy of calo tower [GeV]").c_str(),50, 0, 5, 40, 0, 10);
    
    
    h_calo_energy_pion_deltaEta[i] = new TH2F(("h_calo_energy_pion_deltaEta_" + tag[i]).c_str(),("calo energy vs #Delta#eta(#pi,calo) " + tag[i] +";#Delta#eta(#pi,calo);energy of calo tower [GeV]").c_str(),50, 0, 5, 40, 0, 10);
    
    h_calo_sumE[i] = new TH1F(("h_calo_sumE_" + tag[i]).c_str(),("calo energy sum " + tag[i]).c_str(),60, 0, 600);
    h_calo_sumE[i]->SetXTitle("calo energy sum [GeV]"); h_calo_sumE[i]->Sumw2();
    /*if (i != 0) {*/h_calo_sumE[i]->SetLineColor(colors[i]); h_calo_sumE[i]->SetMarkerStyle(styles[i]);
    
    h_calo_leadingE[i] = new TH1F(("h_calo_leadingE_" + tag[i]).c_str(),("leading calo energy " + tag[i]).c_str(),15, -2, 58);
    h_calo_leadingE[i]->SetXTitle("leading calo energy [GeV]"); h_calo_leadingE[i]->Sumw2();
    /*if (i != 0) {*/h_calo_leadingE[i]->SetLineColor(colors[i]); h_calo_leadingE[i]->SetMarkerStyle(styles[i]);
    
    h_nCaloTowers[i] = new TH1F(("h_nCaloTowers_" + tag[i]).c_str(),("number of calo towers " + tag[i]).c_str(),60, -5, 595);
    h_nCaloTowers[i]->SetXTitle("number of calo towers"); h_nCaloTowers[i]->Sumw2();
    /*if (i != 0) {*/h_nCaloTowers[i]->SetLineColor(colors[i]); h_nCaloTowers[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersMuon[i] = new TH1F(("h_nTowersMuon_" + tag[i]).c_str(),("number of towers matching #mu " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersMuon[i]->SetXTitle("number of towers mathing #mu"); h_nTowersMuon[i]->Sumw2();
    /*if (i != 0) {*/h_nTowersMuon[i]->SetLineColor(colors[i]); h_nTowersMuon[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersPion[i] = new TH1F(("h_nTowersPion_" + tag[i]).c_str(),("number of towers matching #pi^{#pm} " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersPion[i]->SetXTitle("number of towers mathing #pi^{#pm}"); h_nTowersPion[i]->Sumw2();
    /*if (i != 0) {*/h_nTowersPion[i]->SetLineColor(colors[i]); h_nTowersPion[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersECALMuon[i] = new TH1F(("h_nTowersECALMuon_" + tag[i]).c_str(),("number of ECAL towers matching #mu " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersECALMuon[i]->SetXTitle("number of ECAL towers mathing #mu"); h_nTowersECALMuon[i]->Sumw2();
    h_nTowersECALMuon[i]->SetLineColor(colors[i]); h_nTowersECALMuon[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersECALPion[i] = new TH1F(("h_nTowersECALPion_" + tag[i]).c_str(),("number of ECAL towers matching #pi^{#pm} " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersECALPion[i]->SetXTitle("number of ECAL towers mathing #pi^{#pm}"); h_nTowersECALPion[i]->Sumw2();
    h_nTowersECALPion[i]->SetLineColor(colors[i]); h_nTowersECALPion[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersECALFSR[i] = new TH1F(("h_nTowersECALFSR_" + tag[i]).c_str(),("# ECAL towers in FSR region " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersECALFSR[i]->SetXTitle("# ECAL towers in FSR region"); h_nTowersECALFSR[i]->Sumw2();
    h_nTowersECALFSR[i]->SetLineColor(colors[i]); h_nTowersECALFSR[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersHCALMuon[i] = new TH1F(("h_nTowersHCALMuon_" + tag[i]).c_str(),("number of HCAL towers matching #mu " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersHCALMuon[i]->SetXTitle("number of HCAL towers mathing #mu"); h_nTowersHCALMuon[i]->Sumw2();
    h_nTowersHCALMuon[i]->SetLineColor(colors[i]); h_nTowersHCALMuon[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersHCALPion[i] = new TH1F(("h_nTowersHCALPion_" + tag[i]).c_str(),("number of HCAL towers matching #pi^{#pm} " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersHCALPion[i]->SetXTitle("number of HCAL towers mathing #pi^{#pm}"); h_nTowersHCALPion[i]->Sumw2();
    h_nTowersHCALPion[i]->SetLineColor(colors[i]); h_nTowersHCALPion[i]->SetMarkerStyle(styles[i]);
    
    h_nTowersHCALFSR[i] = new TH1F(("h_nTowersHCALFSR_" + tag[i]).c_str(),("# HCAL towers in FSR region " + tag[i]).c_str(),10, -0.5, 9.5);
    h_nTowersHCALFSR[i]->SetXTitle("# HCAL towers in FSR region"); h_nTowersHCALFSR[i]->Sumw2();
    h_nTowersHCALFSR[i]->SetLineColor(colors[i]); h_nTowersHCALFSR[i]->SetMarkerStyle(styles[i]);
    
    h_ratio_ECAL_HCAL_FSR[i] = new TH1F(("h_ratio_ECAL_HCAL_FSR_" + tag[i]).c_str(),("ECAL/(ECAL+HCAL) in FSR region " + tag[i]).c_str(),20, 0.0, 1);
    h_ratio_ECAL_HCAL_FSR[i]->SetXTitle("ECAL/(ECAL+HCAL) in FSR region"); h_ratio_ECAL_HCAL_FSR[i]->Sumw2();
    h_ratio_ECAL_HCAL_FSR[i]->SetLineColor(colors[i]); h_ratio_ECAL_HCAL_FSR[i]->SetMarkerStyle(styles[i]);
    
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
    
    h_calo_tauEta_leadingHFeta[i] = new TH2F(("h_calo_tauEta_leadingHFeta_" + tag[i]).c_str(),("leading HF #eta vs #tau_{1prong} #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_tauEta_leadingHFeta[i]->SetXTitle("#tau_{1prong} #eta"); h_calo_tauEta_leadingHFeta[i]->SetYTitle("leading HF #eta");
    
    h_calo_muEta_averageHFeta[i] = new TH2F(("h_calo_muEta_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs muon #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_muEta_averageHFeta[i]->SetXTitle("muon #eta"); h_calo_muEta_averageHFeta[i]->SetYTitle("average HF #eta");
    
    h_calo_tauEta_averageHFeta[i] = new TH2F(("h_calo_tauEta_averageHFeta_" + tag[i]).c_str(),("average HF #eta vs #tau_{1prong} #eta " + tag[i]).c_str(),19, -3, 3, 20, -5, 5);
    h_calo_tauEta_averageHFeta[i]->SetXTitle("#tau_{1prong} #eta"); h_calo_tauEta_averageHFeta[i]->SetYTitle("average HF #eta");
    
    // MET
    
    h_MET[i] = new TH1F(("h_MET_" + tag[i]).c_str(),("MET " + tag[i]).c_str(),35, -0.5, 34.5);
    h_MET[i]->SetXTitle("MET"); h_MET[i]->Sumw2();
    /*if (i != 0) {*/h_MET[i]->SetLineColor(colors[i]); h_MET[i]->SetMarkerStyle(styles[i]);
    
  } //loop on nSamples
  
  for (int i = nSamples; i < nSamples+2+nNchCategories; i++){
    
    A_highNch_highHF[i] = new TH1F(("A_highNch_highHF_" + to_string(i)).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - high HF " + to_string(i)).c_str(), deltaphi_bins, 0, TMath::Pi());
    A_highNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - high HF"); A_highNch_highHF[i]->Sumw2();
    A_highNch_highHF[i]->SetLineColor(colors[i]); A_highNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    B_lowNch_highHF[i] = new TH1F((ABCDsysNames[i-nSamples]).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) - " + ABCDsysNames[i-nSamples]).c_str(), deltaphi_bins, 0, TMath::Pi());
    B_lowNch_highHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - high HF"); B_lowNch_highHF[i]->Sumw2();
    B_lowNch_highHF[i]->SetLineColor(colors[i]); B_lowNch_highHF[i]->SetMarkerStyle(styles[i]);
  
    C_highNch_lowHF[i] = new TH1F(("C_highNch_lowHF_" + to_string(i)).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - low HF " + to_string(i)).c_str(), deltaphi_bins, 0, TMath::Pi());
    C_highNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) high Nch - low HF"); C_highNch_lowHF[i]->Sumw2();
    C_highNch_lowHF[i]->SetLineColor(colors[i]); C_highNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
    D_lowNch_lowHF[i] = new TH1F(("D_lowNch_lowHF_" + to_string(i)).c_str(),("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - low HF " + to_string(i)).c_str(), deltaphi_bins, 0, TMath::Pi());
    D_lowNch_lowHF[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong}) low Nch - low HF"); D_lowNch_lowHF[i]->Sumw2();
    D_lowNch_lowHF[i]->SetLineColor(colors[i]); D_lowNch_lowHF[i]->SetMarkerStyle(styles[i]);
  
  }
  
  for (int h = 0; h < histograms.size(); h++){
    string histogramsName = histograms.at(h)[0]->GetName();
    histograms.at(h)[nSamples] = (TH1F*)histograms.at(h)[0]->Clone(("A - "+histogramsName).c_str());
    histograms.at(h)[nSamples+1] = (TH1F*)histograms.at(h)[0]->Clone(("B - "+histogramsName).c_str());
    histograms.at(h)[nSamples+2] = (TH1F*)histograms.at(h)[0]->Clone(("C - "+histogramsName).c_str());
    histograms.at(h)[nSamples+4] = (TH1F*)histograms.at(h)[0]->Clone(("same sign - "+histogramsName).c_str());
  }
  

  TH2F *h_calo_energyHFp_nch = new TH2F("h_calo_energyHFp_nch","leading tower HF+ vs nch",9 , 2.5, 11.5, 12, 1, 5);
  h_calo_energyHFp_nch->SetXTitle("nch"); h_calo_energyHFp_nch->SetYTitle("leading tower HF+ [GeV]");
  
  TH2F *h_calo_energyHFm_nch = new TH2F("h_calo_energyHFm_nch","leading tower HF- vs nch",9 , 2.5, 11.5, 12, 1, 5);
  h_calo_energyHFm_nch->SetXTitle("nch"); h_calo_energyHFm_nch->SetYTitle("leading tower HF- [GeV]");

    
  TH1F *h_SB_deltaphi = new TH1F("h_SB_deltaphi","S/sqrt(S+B) for #Delta#phi(#tau_{#mu}, #tau_{1prong})",deltaphi_bins, 0, TMath::Pi());
  h_SB_deltaphi->SetXTitle("S/sqrt(S+B) for #Delta#phi(#tau_{#mu}, #tau_{1prong})"); h_SB_deltaphi->Sumw2();
  h_SB_deltaphi->SetLineColor(colors[0]); h_SB_deltaphi->SetMarkerStyle(styles[1]);
  h_SB_deltaphi->SetBinErrorOption(TH1::kPoisson);
    
  // other
  TH1F *h_track_activity_pt     = new TH1F("h_track_activity_pt",     "h_track_activity_pt",     25, 0, 10); h_track_activity_pt->Sumw2(); h_track_activity_pt->SetXTitle("track p_{T} [GeV]");
  TH2F *h_track_activity_pt_eta = new TH2F("h_track_activity_pt_eta", "h_track_activity_pt_eta", 25, 0, 10, 10, -2.5, 2.5);
  h_track_activity_pt_eta->SetXTitle("track p_{T} [GeV]"); h_track_activity_pt_eta->SetYTitle("track #eta");
  
  TH2F *h_deltaphi_tau_mu_tau_hadron_nch = new TH2F("h_deltaphi_tau_mu_tau_hadron_nch", "nch vs #Delta#phi(#tau_{#mu}, #tau_{1prong})", 8, 0, TMath::Pi(), 5 , 0.5, 5.5);
  h_deltaphi_tau_mu_tau_hadron_nch->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{1prong})"); h_deltaphi_tau_mu_tau_hadron_nch->SetYTitle("nch");


// end of histogram declaration
// *******************************************

//  double mass1, mass2, mass3;
//  TFile *ntuple = new TFile("ntuple.root", "RECREATE");
//  TTree *aux;
//  aux = new TTree("tree", "tree");
//  aux->Branch("mass1", &mass1);
//  aux->Branch("mass2", &mass2);
//  aux->Branch("mass3", &mass3);

  //double muon_mass = 105.6583 / 1000.; //GeV
  //double pionMass = 139.570 / 1000.; //GeV
  double muon_mass = 105.6583 / 1000.; //GeV
  double pionMass = 139.570 / 1000.; //GeV
  int entries = (TREE->fChain)->GetEntries();
  
  TLorentzVector tau_muon, tau_hadron, full_tau_hadron, gamma, tempGamma, pion[3];
  int charge_counter[nSamples][3];
  for (int s = 0; s < nSamples; s++){
    for (int i = 0; i < 3; i++) charge_counter[s][i] = 0;
  }
  
  int skimmedEvent;
  int nSkimmed;
  int category = -1;

#ifdef writeSkimData
  TFile *skimmedEventsFile = new TFile((tagSkimData+".root").c_str(),"recreate");
  TTree *skimmedEventTree = new TTree ("skimmedEventTree","skimmedEventTree");
  skimmedEventTree->Branch("skimmedEvent",&skimmedEvent,"skimmedEvent/I");
  skimmedEventTree->Branch("category",&category,"category/I");
#endif

#ifdef readSkimData
  TFile *skimmedEventsFile = new TFile((tagSkimData+".root").c_str(),"READ");
  TTree *skimmedEventTree = (TTree*) skimmedEventsFile->Get("skimmedEventTree");
  nSkimmed  = skimmedEventTree->GetEntries();
  cout << "number of skimmed events: " << nSkimmed << endl;
  cout << "dataset tag:  " << tagSkimData << endl;
  skimmedEventTree->SetBranchAddress("skimmedEvent", &skimmedEvent);
  skimmedEventTree->SetBranchAddress("category", &category);
#endif  
  
  
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
  float pionLeadingPt = -1;
  //float muonSF = tnp_weight_trk_pbpb(0);
  //float muonSFerror = TMath::Max(TMath::Abs(tnp_weight_trk_pbpb(-1)-tnp_weight_trk_pbpb(0)),TMath::Abs(tnp_weight_trk_pbpb(0)-tnp_weight_trk_pbpb(-2)));
  //cout << "muon inner tracker SF: " << muonSF << endl;
  //cout << "uncertainty on muon inner tracker SF: " << muonSFerror << endl;
  int muon_charge = 0;
  int tauh_charge = 0;
  float muon_weight = 1;
  float muon_weight_error = 0;
  bool passedNch = false;
  int candCounter = 0;
  float tau_hadron_ptScalar = 0;
  float ditau_ptScalar = 0;
  float delta_phi = 0;
  float full_delta_phi = 0;
  float maxSB = 0;
  int maxSBbin = 1;
  int nBatches = 20;
#ifdef readSkimData
  int printing_every = nSkimmed/nBatches;
  int lastPrintedEntry = -1*printing_every;
#else
  int printing_every = entries/nBatches;
  int lastPrintedEntry = -1*printing_every;
#endif
  int temp_counter = 0;

  cout << "Running on Data ..." << endl;
#ifdef readSkimData
  for(int iSkimmedEntry=0; iSkimmedEntry<nSkimmed; iSkimmedEntry++){
    skimmedEventTree->GetEntry(iSkimmedEntry);
    (TREE->fChain)->GetEntry(skimmedEvent);
#else
  for(int iEntry=0; iEntry<entries; iEntry++){
    (TREE->fChain)->GetEntry(iEntry);
#endif
    
    //if (TREE->BsTauTau_mu1_pt->at(0) < 2.5) continue; // fix me
    //if (TREE->BsTauTau_nch->at(0) > 8) continue; // fix me
    
    TLorentzVector tau_hadron_pi1, tau_hadron_pi2, tau_hadron_pi3;
    
    bool triggered = TREE->triggered->at(0);
    if (!triggered) continue; // fix me
    int temp_nch = TREE->BsTauTau_nch->at(0); 
    //if (temp_nch>(firstNchCategory+nNchCategories-1)) continue; // fix me if ABCD categorization changed
    int nPions=TREE->BsTauTau_nPions->at(0);
    
    
    double maxHFp = TREE->calo_leading_energy_HFp->at(0);
    double maxHFm = TREE->calo_leading_energy_HFm->at(0);
    if (maxHFp > 15 || maxHFm > 15) continue; // fix me
    double sumHFp = TREE->calo_sum_energy_HFp->at(0);
    double sumHFm = TREE->calo_sum_energy_HFm->at(0);
    int sizeHFp = TREE->calo_nTowerHFp->at(0);
    int sizeHFm = TREE->calo_nTowerHFm->at(0);
    
    //if (!(iEntry%printing_every)) cout << "\r" << int((100.0/nBatches)*(iEntry/printing_every)) << "%" << flush;
    
    pionLeadingIndex = -1;
    pionLeadingPt = -1;
    bool passedmu  = false;
    bool passedtau = false;
    bool passedcalo = true;
    bool passedcaloDown = true;
    bool passedcaloUp = true;
    bool passedMET = true;
    bool passedGamma = true;
    bool sameSign = false;
    passedNch = false;
    category = -1;
    tau_hadron_ptScalar = 0;
    ditau_ptScalar = 0;
    double temp_tau_pt_comparison = 0.;
    muon_charge = 0;
    /*double temp_npv = 0.;
    double temp_pvz = 0.;
    double temp_pv_trk_1 = 0.;
    double temp_pv_trk_2 = 0.;
    double temp_pv_trk_3 = 0.;*/
    
    float sumZDCp = -1;
    float sumZDCm = -1;
    bool passedZDC = true;
    if (hasZDCinfo){
      sumZDCp = TREE->BsTauTau_calo_zdcSumPlus->at(0);
      sumZDCm = TREE->BsTauTau_calo_zdcSumMinus->at(0);
      if (sumZDCp > maxZDCp || sumZDCm > maxZDCm || sumZDCp < minZDCp || sumZDCm < minZDCm) passedZDC = false;
      //passedZDC = !passedZDC;
    }
    if (!passedZDC) continue;
    
    if(TREE->unfiltered_muon_pt->size()>1) continue;
    
    for (int i=0; i<(int)TREE->BsTauTau_mu1_pt->size(); i++) {
if( (TREE->BsTauTau_mu1_pt->at(i) < barrel_mu_pt_cut && TMath::Abs(TREE->BsTauTau_mu1_eta->at(i)) < 1.2) || (TREE->BsTauTau_mu1_pt->at(i) < endcap_mu_pt_cut && TMath::Abs(TREE->BsTauTau_mu1_eta->at(i)) > 1.2))  continue;
      if (tau_muon_isSoft){
        if (!TREE->BsTauTau_mu1_isSoft->at(i)) continue;
      }
      if (tau_muon_isGlobal){
        if (!TREE->BsTauTau_mu1_isGlobal->at(i)) continue;
      }
      if (tau_muon_isTracker){
        if (!TREE->BsTauTau_mu1_isTracker->at(i)) continue;
      }
      //if (TREE->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }
      if (TMath::Abs(TREE->BsTauTau_mu1_eta->at(i)) < 2.4) passedmu = true;
      if (TREE->BsTauTau_mu1_pt->at(i) > temp_tau_pt_comparison){
        temp_tau_pt_comparison = TREE->BsTauTau_mu1_pt->at(i);
        tau_muon.SetPtEtaPhiM (TREE->BsTauTau_mu1_pt->at(i), TREE->BsTauTau_mu1_eta->at(i), TREE->BsTauTau_mu1_phi->at(i), muon_mass);
        muon_charge = TREE->BsTauTau_mu1_q->at(i);
      }
    } // loop over the size of the muons
    
    //if(TREE->BsTauTau_mu1_pt->size()>1) passedmu = false;
    if(TREE->unfiltered_muon_pt->size()>1) passedmu = false;
    if (!passedmu) continue;

    
    
    //if (temp_nch != 3) continue;   
    
    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    candCounter = 0;
    int chosenTauIndex = -1;
    float vprob = -1;
    float temp_pt = 0;
    double temp_tau_rho1 = 0.; double temp_tau_rho2 = 0.;
    
    int nAccPion = 0;
    
    for (int i=0; i<(int)TREE->BsTauTau_tau_pt->size(); i++) {
      //if (muon_charge*TREE->BsTauTau_tau_q->at(i) != MuTauCharge) continue;
      if (muon_charge*TREE->BsTauTau_tau_q->at(i) != MuTauCharge) sameSign = true;
      if (TREE->BsTauTau_tau_pt->at(i) > pionExtraCut) nAccPion++;
      if (TREE->BsTauTau_tau_pt->at(i) < pionLeadingCut) continue;
      if (TREE->BsTauTau_tau_pt->at(i) > temp_pt) {
        vprob = 100*TREE->BsTauTau_B_vprob->at(i);
        if (temp_nch == 1) vprob = 100*TREE->BsTauTau_B_vprob->at(i);
        temp_pt = TREE->BsTauTau_tau_pt->at(i);
      }
      if (TREE->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      passedtau = true;
      //if (temp_nch == 1 && TREE->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      
      if (TREE->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREE->BsTauTau_tau_pt->at(i);
        tauh_charge = TREE->BsTauTau_tau_q->at(i);
        tau_hadron.SetPtEtaPhiM (TREE->BsTauTau_tau_pt->at(i), TREE->BsTauTau_tau_eta->at(i), TREE->BsTauTau_tau_phi->at(i), pionMass);
      }
    } // loop over the size of the tau candidates
    //if (candCounter != 1) passedNch = false;
    
    
    //if (temp_nch>=min_nch && temp_nch<=max_nch && nPions==temp_nch) passedNch = true; 
    if (nAccPion>=min_nch && nAccPion<=max_nch) passedNch = true;
    //else passedtau = false;
    
    //if (nAccPion>(firstNchCategory+nNchCategories-1)) continue; // fix me if ABCD categorization changed
    
    ditau_ptScalar = tau_muon.Pt() + tau_hadron.Pt();
    bool passedDitauMass = true;
    bool passedDitauPt = true;
    TLorentzVector ditau = tau_muon+tau_hadron;
    TLorentzVector MET;
    MET.SetPtEtaPhiM(ditau.Pt(),-ditau.Eta(),TMath::Pi()-ditau.Phi(),0);
    //TLorentzVector unboosted_muon = tau_muon + 0.5*MET;
    //TLorentzVector unboosted_pion = tau_hadron + 0.5*MET;
    
    delta_phi = TMath::Abs(tau_muon.DeltaPhi(tau_hadron));
    float acoplanarity = 1-delta_phi/TMath::Pi();
    full_delta_phi = TMath::Abs(tau_muon.DeltaPhi(full_tau_hadron));
    //float unboosted_delta_phi = TMath::Abs(unboosted_muon.DeltaPhi(unboosted_pion));
    
    if (ditau.M() < ditauMassCut) passedDitauMass = false;
    //if (ditau.M() < 3.0 || ditau.M() > 3.3) {passedDitauMass = false; continue;}
    
    //if (ditau.Pt() < ditau_pt_cut) passedDitauPt = false;
    if (ditau.Pt() < ditau_pt_cut && acoplanarity < aco_cut_for_low_ditau_pt) passedDitauPt = false;
    //if (ditau.Pt() > 0.3){ passedDitauPt = false; continue;};
    //if (!passedDitauPt) continue;
    
    TVector3 muon_v3 = tau_muon.Vect();
    TVector3 pion_v3 = tau_hadron.Vect();
    TVector3 restFrame = 0.5*muon_v3 + 0.5*pion_v3;
    TVector3 muon_RF = muon_v3 - restFrame;
    TVector3 pion_RF = pion_v3 - restFrame;
    
    float muon_cosThetaStar = cos(muon_RF.Angle(restFrame));
    float pion_cosThetaStar = cos(pion_RF.Angle(restFrame));
    
    float muon_RF_cosDeltaPhi = cos(muon_RF.DeltaPhi(restFrame));
    float pion_RF_cosDeltaPhi = cos(pion_RF.DeltaPhi(restFrame));
    
    float muon_AP_phi = muon_v3.Angle(restFrame);
    float pion_AP_phi = pion_v3.Angle(restFrame);
    float phiAP = muon_AP_phi+pion_AP_phi;
    float alphaAP = sin(pion_AP_phi-muon_AP_phi) / sin(phiAP);
    float epsilonAP = 2*sin(pion_AP_phi)*sin(muon_AP_phi)/sin(phiAP);
    
    //cout << "muon angle: "<< 180*muon_v3.DeltaPhi(restFrame)/TMath::Pi() << "  muon RF angle: "<< 180*muon_RF.DeltaPhi(restFrame)/TMath::Pi() << endl;
    //cout << "pion angle: "<< 180*pion_v3.DeltaPhi(restFrame)/TMath::Pi() << "  pion RF angle: "<< 180*pion_RF.DeltaPhi(restFrame)/TMath::Pi() << endl;
    
    
    
    //if (ditau.Pt() > MET_cut) passedMET = false;

    bool found_mu_track = false;
    bool found_pi1_track = false;
    bool found_pi2_track = false;
    bool found_pi3_track = false;
    bool high_activity = false;

    //temp_pvz = TREE->BsTauTau_bbPV_vz->at(0);
    //temp_npv = TREE->PV_N;
    double maxHF = 0;
    if (maxHFp > HFpLeading_high || maxHFm > HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcalo = false;
    if (maxHFp > 0.9*HFpLeading_high || maxHFm > 0.9*HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcaloDown = false;
    if (maxHFp > 1.1*HFpLeading_high || maxHFm > 1.1*HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcaloUp = false;
    double leadingHFeta = 0;
    double leadingHFphi = 0;
    double averageHFeta = 0;
    double sumEt = TREE->calo_sum_eT->at(0);
    double sumE = TREE->calo_sum_energy->at(0);
    double leadingE = TREE->calo_leading_energy->at(0);
    double leadingE_eta = TREE->calo_leading_energy_eta->at(0);
    double leadingE_nonHF = 0;
    double leadingE_eta_nonHF = 1000;
    double leadingE_phi_nonHF = 1000;
    double leadingEt = TREE->calo_leading_eT->at(0);
    double leadingEt_Eta = TREE->calo_leading_eT_eta->at(0);
    int nTowersMuon = 0;
    int nTowersPion = 0;
    int nTowersECALMuon = 0;
    int nTowersECALPion = 0;
    int nTowersECALFSR = 0;
    int nTowersHCALMuon = 0;
    int nTowersHCALPion = 0;
    int nTowersHCALFSR = 0;
    double ECAL_energy_FSR = 0;
    double HCAL_energy_FSR = 0;
    
    bool SR = false;
    if (passedmu && triggered && passedtau /*&& muon_charge*tauh_charge == -1*/ && passedDitauPt && passedDitauMass) SR = true;
    if (SR){
      for (int i=0; i<(int)TREE->EB_eta->size(); i++) {
        double eta = TREE->EB_eta->at(i);
        double energy = TREE->EB_energy->at(i);
        if (abs(eta) > 1.442) continue;
        //if (energy < 0.7) continue;
        //double deltaR_muon = sqrt(pow(deltaphi(TREE->EB_phi->at(i),tau_muon.Phi()),2)+pow(TREE->EB_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREE->EB_phi->at(i),tau_hadron.Phi()),2)+pow(TREE->EB_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREE->EB_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREE->EB_eta->at(i)-tau_hadron.Eta());
        //double deltaR_MET = sqrt(pow(deltaphi(TREE->EB_phi->at(i),MET.Phi()),2)+pow(TREE->EB_eta->at(i)-MET.Eta(),2));
        if (deltaEta_muon < muonEB_DEta) nTowersECALMuon += 1;
        if (deltaEta_pion < pionEB_DEta) nTowersECALPion += 1;
        if (1-deltaphi(TREE->EB_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersECALFSR += 1;
          ECAL_energy_FSR += energy;
        }
        //if (deltaR_MET < MET_EB_DEta) nTowersECALMET += 1;
      }
      for (int i=0; i<(int)TREE->EE_eta->size(); i++) {
        double eta = TREE->EE_eta->at(i);
        double energy = TREE->EE_energy->at(i);
        if (abs(eta) < 1.566 || abs(eta) > 2.6) continue;
        //if (energy < 3) continue;
        //double deltaR_muon = sqrt(pow(deltaphi(TREE->EE_phi->at(i),tau_muon.Phi()),2)+pow(TREE->EE_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREE->EE_phi->at(i),tau_hadron.Phi()),2)+pow(TREE->EE_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREE->EE_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREE->EE_eta->at(i)-tau_hadron.Eta());
        if (deltaEta_muon < muonEE_DEta) nTowersECALMuon += 1;
        if (deltaEta_pion < pionEE_DEta) nTowersECALPion += 1;
        if (1-deltaphi(TREE->EE_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersECALFSR += 1;
          ECAL_energy_FSR += energy;
        }
      }
      for (int i=0; i<(int)TREE->HB_eta->size(); i++) {
        double eta = TREE->HB_eta->at(i);
        double energy = TREE->HB_energy->at(i);
        if (abs(eta) > 1.305) continue;
        //if (energy < 2.8) continue; // true threshold
        //double deltaR_muon = sqrt(pow(deltaphi(TREE->HB_phi->at(i),tau_muon.Phi()),2)+pow(TREE->HB_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREE->HB_phi->at(i),tau_hadron.Phi()),2)+pow(TREE->HB_eta->at(i)-tau_hadron.Eta(),2));
        //double deltaR_pion = sqrt(pow(TREE->HB_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREE->HB_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREE->HB_eta->at(i)-tau_hadron.Eta());
        if (deltaEta_muon < muonHB_DEta) nTowersHCALMuon += 1;
        if (deltaEta_pion < pionHB_DEta) nTowersHCALPion += 1;
        if (1-deltaphi(TREE->HB_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersHCALFSR += 1;
          HCAL_energy_FSR += energy;
        }
      }
      for (int i=0; i<(int)TREE->HE_eta->size(); i++) {
        double eta = TREE->HE_eta->at(i);
        double energy = TREE->HE_energy->at(i);
        if (abs(eta) < 1.41 || abs(eta) > 3) continue;
        //if (energy < 1) continue; // true threshold
        //double deltaR_muon = sqrt(pow(deltaphi(TREE->HE_phi->at(i),tau_muon.Phi()),2)+pow(TREE->HE_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREE->HE_phi->at(i),tau_hadron.Phi()),2)+pow(TREE->HE_eta->at(i)-tau_hadron.Eta(),2));
        //double deltaR_pion = sqrt(pow(TREE->HE_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREE->HE_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREE->HE_eta->at(i)-tau_hadron.Eta());
        if (deltaEta_muon < muonHE_DEta) nTowersHCALMuon += 1;
        if (deltaEta_pion < pionHE_DEta) nTowersHCALPion += 1;
        if (1-deltaphi(TREE->HE_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersHCALFSR += 1;
          HCAL_energy_FSR += energy;
        }
      }
    } // if(SR)
    
    
#ifdef writeSkimData
    if (nAccPion >= firstNchCategory && nAccPion <= firstNchCategory+nNchCategories-1 && !passedcalo) category = nSamples;
    if (passedNch && !passedcalo) category = nSamples+1;
    if (nAccPion >= firstNchCategory && nAccPion <= firstNchCategory+nNchCategories-1 && passedcalo) category = nSamples+2;
    if (passedNch && passedcalo){
      if (muon_charge*tauh_charge == MuTauCharge) category = 0;
      if (muon_charge*tauh_charge == -1*MuTauCharge) category = nSamples+4;
    }
#endif
    
    
    TLorentzVector leadingGamma, leadingFSR, matchingFSR;
    leadingGamma.SetPtEtaPhiM (0, 100, 0, 0);
    leadingFSR.SetPtEtaPhiM (0, 100, 0, 0);
    matchingFSR.SetPtEtaPhiM (1000, 100, 0, 0);
    bool passedGammaDitau = true;

    bool counted = false;
    
    SR = false;
    if (passedmu /*&& passedNch && passedcalo*/ && triggered && passedtau && muon_charge*tauh_charge == -1 && passedDitauPt && passedDitauMass && passedZDC) SR = true;
#ifdef writeSkimData
    if(SR)
    {
    TLorentzVector recoPiZero;
    int Ngamma = 0;
    float minDeltaR = 1000;
    float minDeltaM = 100000; //MeV
    int NrecoPiZero = 0;
    if (!threeProng[0]){
    Ngamma = TREE->BsTauTau_nGammas->at(0);
    full_tau_hadron = tau_hadron;
    for (int i=0; i<TREE->BsTauTau_nGammas->at(0); i++){
      gamma.SetPtEtaPhiM (TREE->reco_gamma_pt->at(i),TREE->reco_gamma_eta->at(i),TREE->reco_gamma_phi->at(i),0);
      if (gamma.Pt()<gammaPtCut || TMath::Abs(gamma.Eta()) > gammaEtaCut){
        Ngamma -= 1;
        continue;
      }
      if (gamma.Pt() > leadingGamma.Pt()) leadingGamma = gamma;
      if (gamma.Pt() > leadingFSR.Pt() && acop(gamma,ditau) < FSR_ditau_acoCut) leadingFSR = gamma;
      if ((ditau+gamma).Pt() < (ditau+matchingFSR).Pt() && acop(gamma,ditau) < FSR_ditau_acoCut) matchingFSR = gamma;
      
      if (TMath::Abs(gamma.Pt()-ditau.Pt()) < 5 && !category){
        //cout << temp_counter << " gamma pt: " << gamma.Pt() << " eta: " << gamma.Eta() << " phi: " << gamma.Phi() << " ditau aco:         " << int(1000*acop(gamma,ditau)) << " res pt: " << gamma.Pt()-ditau.Pt() << endl;
        counted = true;
      }
      /*if (TMath::Abs(gamma.DeltaR(tau_hadron)) > piZeroDeltaRCut){
      //if (TMath::Abs(gamma.DeltaR(tau_muon)) > piZeroDeltaRCut){
        Ngamma -= 1;
        continue;
      }*/
      //if (i == 0) recoNeutralPion = gamma;
      //else recoNeutralPion += gamma;
      for (int j=i+1; j<TREE->BsTauTau_nGammas->at(0); j++){
        if (TREE->reco_gamma_pt->at(j)<gammaPtCut || TMath::Abs(TREE->reco_gamma_eta->at(j)) > gammaEtaCut) continue;
        tempGamma.SetPtEtaPhiM (TREE->reco_gamma_pt->at(j),TREE->reco_gamma_eta->at(j),TREE->reco_gamma_phi->at(j),0);
        //if (TMath::Abs(tempGamma.DeltaPhi(tau_hadron)) > piZeroDeltaPhiCut) continue;
        //if (TMath::Abs(tempGamma.DeltaPhi(tau_muon)) > piZeroDeltaPhiCut) continue;
        //if (TMath::Abs(tempGamma.DeltaR(tau_hadron)) > piZeroDeltaRCut) continue;
        //if (TMath::Abs(tempGamma.DeltaR(tau_muon)) > piZeroDeltaRCut) continue;
        float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
        if (deltaR < minDeltaR) minDeltaR = deltaR;
        if (deltaR > gammaDeltaRCut) continue;
        recoPiZero = gamma+tempGamma;
        if (recoPiZero.Pt() < piZeroPtCut) continue;
        double recoPiZeroMass = 1000*recoPiZero.M(); // MeV
        if (TMath::Abs(recoPiZeroMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = recoPiZeroMass-PiZeroMass;
        if (SR && category == 0){
          h_recoPiZeroDeltaM[0]->Fill(recoPiZeroMass-PiZeroMass);
          h_recoPiZero_pt[0]->Fill(recoPiZero.Pt());
          h_recoPiZero_eta[0]->Fill(recoPiZero.Eta());
          h_recoPiZero_phi[0]->Fill(recoPiZero.Phi());
          h_recoPiZero_deltaphi_muon[0]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_muon)));
          h_recoPiZero_deltaphi_pion[0]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)));
        }
        NrecoPiZero += 1;
        if (gamma.Pt()>gammaPtCut && tempGamma.Pt()>gammaPtCut) full_tau_hadron += recoPiZero;
      } // loop on secondary gammas
      
    } // loop on reco gammas
    
    
    //leadingFSR=matchingFSR; // fix me
    
    float acoGammaDitau = acop(leadingFSR,ditau);
    float resGammaDitauPt = leadingFSR.Pt()-ditau.Pt();
    if (leadingFSR.Pt()==0){
      acoGammaDitau = -1;
      resGammaDitauPt = -1000;
    }
    /*for (int i=0; i<2; i++){
      if(acoGammaDitau<acoGammaDitau_cut[i]){
        if(resGammaDitauPt<resDitauGammaPt_cut_up[i] && resGammaDitauPt>resDitauGammaPt_cut_down[i]) passedGammaDitau = false;
      }
    }*/
    if (counted) temp_counter++;
    //if (temp_counter == 20) break;
    
    
    if(leadingGamma.E() > gammaECut && resGammaDitauPt<resDitauGammaPt_cut_up && resGammaDitauPt>resDitauGammaPt_cut_down) passedGammaDitau = false;
    
    
    //if (!passedGammaDitau) SR = false;

    if (SR && category != -1) h_residual_leadingRecoGamma_ditau_pt[category]->Fill(resGammaDitauPt);
    
    if (SR && passedNch && passedcalo){
      h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma[0]->Fill(leadingFSR.Pt(),resGammaDitauPt);
      h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma[0]->Fill(leadingFSR.E(),resGammaDitauPt);
      h_residual_leadingRecoGamma_ditau_pt_aco[0]->Fill(acoGammaDitau,resGammaDitauPt);
      
      if (passedGammaDitau){
        if (leadingGamma.Pt() != 0){
          h_leadingRecoGamma_ditau_deltaPhi[0]->Fill(TMath::Abs(ditau.DeltaPhi(leadingGamma)));
          h_leadingRecoGamma_ditau_deltaR[0]->Fill(MET.DeltaR(leadingGamma));
          h_recoGamma_pt[0]->Fill(leadingGamma.Pt());
          h_recoGamma_eta[0]->Fill(leadingGamma.Eta());
          h_recoGamma_phi[0]->Fill(leadingGamma.Phi());
          h_recoGamma_deltaphi_muon[0]->Fill(TMath::Abs(leadingGamma.DeltaPhi(tau_muon)));
          h_recoGamma_deltaphi_pion[0]->Fill(TMath::Abs(leadingGamma.DeltaPhi(tau_hadron)));
          h_recoGamma_deltaR_muon[0]->Fill(TMath::Abs(leadingGamma.DeltaR(tau_muon)));
          h_recoGamma_deltaR_pion[0]->Fill(TMath::Abs(leadingGamma.DeltaR(tau_hadron)));
        }
        h_NrecoPiZero[0]->Fill(NrecoPiZero);
        h_recoGammasMinDeltaR[0]->Fill(minDeltaR);
        h_recoPiZeroMinDeltaM[0]->Fill(minDeltaM);
        h_NrecoGamma[0]->Fill(Ngamma);
      }
    } // if in the signal region
    } // if not three prong
    } // if true
#endif
    
    
    
    
    
    SR = false;
    if (passedmu && triggered && passedtau && muon_charge*tauh_charge == -1 && passedDitauPt && passedDitauMass && passedGammaDitau && passedZDC) SR = true;
    if(SR)
    {
      for (int i=0; i<(int)TREE->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREE->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREE->BsTauTau_calo_energyHFm->at(i);
        if (passedNch && SR) h_calo_E[0]->Fill(TREE->BsTauTau_calo_energy->at(i));
        if (passedNch && SR) h_calo_E_eta[0]->Fill(TREE->BsTauTau_calo_eta->at(i),TREE->BsTauTau_calo_energy->at(i));
        if (passedNch && SR) h_calo_Et[0]->Fill(TREE->BsTauTau_calo_eT->at(i));
        if (passedNch && SR) h_calo_Et_eta[0]->Fill(TREE->BsTauTau_calo_eta->at(i),TREE->BsTauTau_calo_eT->at(i));
        double etaCal = TREE->BsTauTau_calo_eta->at(i);
        double phiCal = TREE->BsTauTau_calo_phi->at(i);
        
        if (TREE->BsTauTau_calo_energy->at(i) > 0){
          double deltaR_muon = sqrt(pow(deltaphi(phiCal,tau_muon.Phi()),2)+pow(etaCal-tau_muon.Eta(),2));
          double deltaR_pion = sqrt(pow(deltaphi(phiCal,tau_hadron.Phi()),2)+pow(etaCal-tau_hadron.Eta(),2));
          double deltaEta_muon = abs(etaCal-tau_muon.Eta());
          double deltaEta_pion = abs(etaCal-tau_hadron.Eta());
          if (deltaEta_muon < muon_DEta) nTowersMuon += 1;
          if (deltaEta_pion < pion_DEta) nTowersPion += 1;
          if (passedNch && passedcalo && abs(etaCal) < 5){
            h_calo_energy_muon_deltaR[0]->Fill(deltaR_muon,TREE->BsTauTau_calo_energy->at(i));
            h_calo_energy_pion_deltaR[0]->Fill(deltaR_pion,TREE->BsTauTau_calo_energy->at(i));
            h_calo_energy_muon_deltaEta[0]->Fill(deltaEta_muon,TREE->BsTauTau_calo_energy->at(i));
            h_calo_energy_pion_deltaEta[0]->Fill(deltaEta_pion,TREE->BsTauTau_calo_energy->at(i));
          }
        }
        
        if (abs(etaCal) < 2.5 && leadingE_nonHF < TREE->BsTauTau_calo_energy->at(i)){
          leadingE_nonHF = TREE->BsTauTau_calo_energy->at(i);
          leadingE_eta_nonHF = etaCal;
          leadingE_phi_nonHF = phiCal;
        }
        
        if (eHFp != -1){
          if (eHFp > maxHFp) maxHFp = eHFp;
          if (maxHFp > maxHF) {maxHF = maxHFp; leadingHFeta = etaCal; leadingHFphi = phiCal;};
          averageHFeta += eHFp * etaCal;
          if (passedNch && SR) h_calo_energyHFp[0]->Fill(eHFp);
          if (passedNch && SR) h_calo_HF_eta[0]->Fill(etaCal);
          if (passedNch && SR) h_calo_HF_energy_eta[0]->Fill(etaCal,eHFp);
        }
        if (eHFm != -1){
          if (eHFm > maxHFm) maxHFm = eHFm;
          if (maxHFm > maxHF) {maxHF = maxHFm; leadingHFeta = etaCal; leadingHFphi = phiCal;};
          averageHFeta += eHFm * etaCal;
          if (passedNch && SR) h_calo_energyHFm[0]->Fill(eHFm);
          if (passedNch && SR) h_calo_HF_eta[0]->Fill(etaCal);
          if (passedNch && SR) h_calo_HF_energy_eta[0]->Fill(etaCal,eHFm);
        }
      }
      averageHFeta /= (sumHFp+sumHFm);
      if (passedcalo && passedNch && SR) h_calo_sumE[0]->Fill(sumE);
      if (passedcalo && passedNch && SR) h_calo_leadingE[0]->Fill(leadingE);
      if (passedNch && SR) h_calo_leadingE_eta[0]->Fill(leadingE_eta_nonHF,leadingE_nonHF);
      if (passedcalo && passedNch && SR) h_calo_sumEt[0]->Fill(sumEt);
      if (passedcalo && passedNch && SR) h_calo_leadingEt[0]->Fill(leadingEt);
      if (passedNch) h_calo_leadingEt_eta[0]->Fill(leadingEt_Eta,leadingEt);
      if (passedcalo && passedNch && SR) h_nCaloTowers[0]->Fill(TREE->BsTauTau_calo_eta->size());
      if (passedcalo && passedNch && SR) h_nTowersMuon[0]->Fill(nTowersMuon);
      if (passedcalo && passedNch && SR) h_nTowersPion[0]->Fill(nTowersPion);
      if (passedcalo && passedNch && SR) h_nTowersECALMuon[0]->Fill(nTowersECALMuon);
      if (passedcalo && passedNch && SR) h_nTowersECALPion[0]->Fill(nTowersECALPion);
      if (passedcalo && passedNch && SR && leadingFSR.Pt()!=0) h_nTowersECALFSR[0]->Fill(nTowersECALFSR);
      if (passedcalo && passedNch && SR) h_nTowersHCALMuon[0]->Fill(nTowersHCALMuon);
      if (passedcalo && passedNch && SR) h_nTowersHCALPion[0]->Fill(nTowersHCALPion);
      if (passedcalo && passedNch && SR && leadingFSR.Pt()!=0) h_nTowersHCALFSR[0]->Fill(nTowersHCALFSR);
      if (passedcalo && passedNch && SR && (ECAL_energy_FSR+HCAL_energy_FSR) != 0 && leadingFSR.Pt()!=0) h_ratio_ECAL_HCAL_FSR[0]->Fill(ECAL_energy_FSR/(ECAL_energy_FSR+HCAL_energy_FSR));
      if (passedcalo && passedNch && SR) h_calo_energyHFp_sum[0]->Fill(sumHFp);
      if (passedcalo && passedNch && SR) h_calo_energyHFm_sum[0]->Fill(sumHFm);
      if (passedcalo && passedNch && SR) h_calo_energyHF_pm[0]->Fill(sumHFm,sumHFp);
      if (passedcalo && passedNch && SR) h_calo_energyHFp_size[0]->Fill(sizeHFp);
      if (passedcalo && passedNch && SR) h_calo_energyHFm_size[0]->Fill(sizeHFm);
      if (passedNch && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low && SR) h_calo_leadingHFp[0]->Fill(maxHFp);
      if (passedNch && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low && SR) h_calo_leadingHFm[0]->Fill(maxHFm);
      if (passedNch && passedcalo && SR) h_calo_leadingHE[0]->Fill(TREE->calo_leading_energy_HE->at(0));
      if (passedNch && passedcalo && SR) h_calo_leadingHCAL[0]->Fill(TMath::Max(TREE->calo_leading_energy_HB->at(0),TREE->calo_leading_energy_HE->at(0)));
      if (passedNch && passedcalo && SR) h_calo_leadingEE[0]->Fill(TREE->calo_leading_energy_EE->at(0));
      if (passedNch && passedcalo && SR) h_calo_leadingECAL[0]->Fill(TMath::Max(TREE->calo_leading_energy_EB->at(0),TREE->calo_leading_energy_EE->at(0)));
      if (nAccPion>=firstNchCategory && nAccPion<=firstNchCategory+nNchCategories-1 && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low && SR) h_calo_leadingHFp_highNch[0]->Fill(maxHFp);
      if (nAccPion>=firstNchCategory && nAccPion<=firstNchCategory+nNchCategories-1 && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low && SR) h_calo_leadingHFm_highNch[0]->Fill(maxHFm);
      if (passedNch && SR) h_calo_leadingHF_pm[0]->Fill(maxHFm,maxHFp);
    } // if(SR)
    
    

    
    
    h_cutflow[0]->Fill(0); //input from Ntuplizer
    if (passedmu) h_cutflow[0]->Fill(1); //tau mu
    if (passedNch) h_cutflow[0]->Fill(2); //nCh
    if (passedtau) h_cutflow[0]->Fill(3); //tau hadron
    if (passedcalo) h_cutflow[0]->Fill(4); //HF
    if (passedcalo && passedtau) h_cutflow[0]->Fill(5); //HF and tau hadron
    if (passedcalo && passedtau && passedNch) h_cutflow[0]->Fill(6); //HF and tau hadron and nch
    //if () h_cutflow[0]->Fill(4); //vertex prob

    
    bool keepEvent = false;
    
    if (passedmu && passedMET && passedGamma && triggered && passedGammaDitau)
    {
      keepEvent = true;
    }
    if (keepEvent && passedtau && /*muon_charge*tauh_charge != 0 &&*/ passedZDC && delta_phi >= deltaPhi_cut && category != -1){
      if (passedDitauPt) h_ditau_mass[category]->Fill(ditau.M());
      if (passedDitauMass){
        h_ditau_pt[category]->Fill(ditau.Pt());
        if (!category) h_ditau_pt_aco[0]->Fill(acoplanarity,ditau.Pt());
      }
    }
    if (!passedDitauPt || !passedDitauMass) keepEvent = false;
    if (keepEvent && category != -1 && /*muon_charge*tauh_charge != 0 &&*/ passedZDC && delta_phi >= deltaPhi_cut){
      h_tau_hadron_vprob[category]->Fill(vprob);
    }
    if (keepEvent && passedtau && (category == 0 || category == nSamples+4) && passedZDC && delta_phi >= deltaPhi_cut){
      charge_counter[0][muon_charge*tauh_charge + 1] += 1;
    }
    if (!passedtau) keepEvent = false;
    //if (muon_charge*tauh_charge == 0) keepEvent = false;
    if (keepEvent && category != -1 && delta_phi >= deltaPhi_cut){
      h_sumZDCplus[category]->Fill(sumZDCp/1000);
      h_sumZDCminus[category]->Fill(sumZDCm/1000);
      h_sumZDC_pm[0]->Fill(sumZDCm/1000,sumZDCp/1000);
      float averageZDC = (sumZDCp-ratioZDCpm*sumZDCm)/(sumZDCp+ratioZDCpm*sumZDCm);
      h_muEta_averageZDCside[0]->Fill(tau_muon.Eta(),averageZDC);
      h_tauEta_averageZDCside[0]->Fill(tau_hadron.Eta(),averageZDC);
      h_averageZDCside_averageHFeta[0]->Fill(averageZDC,averageHFeta);
    }
    if(!passedZDC) keepEvent = false;
    if (keepEvent && category != -1){
      h_deltaphi_tau_mu_tau_hadron[category]->Fill(delta_phi);
      h_acoplanarity_tau_mu_tau_hadron[category]->Fill(acoplanarity);
      //h_deltaphi_tau_mu_tau_hadron[category]->Fill(delta_phi);
      h_deltaphi_tau_mu_tau_hadron_zoomed[category]->Fill(delta_phi);
      if (!category) h_deltaphi_tau_mu_full_tau_hadron[0]->Fill(full_delta_phi);
      h_deltaphi_tau_mu_tau_hadron_nch->Fill(delta_phi,temp_nch);
      if (!category) h_deltaphi_tau_mu_tau_hadron_mueta[0]->Fill(delta_phi,tau_muon.Eta());
      if (!category) h_deltaphi_tau_mu_tau_hadron_deltaeta[0]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()));
    }
    if(delta_phi < deltaPhi_cut) keepEvent = false;
    if (keepEvent){
      //if (passedcalo) h_tau_hadron_nch[0]->Fill(temp_nch);
      if (category==0 || category==nSamples+2) h_tau_hadron_nch[0]->Fill(nAccPion);
      if (!passedcalo) h_tau_hadron_nch_highHF[0]->Fill(temp_nch);
      if (category==0 || category==nSamples+2) h_tau_hadron_ncand_final[0]->Fill(candCounter);
      h_calo_energyHFp_nch->Fill(temp_nch,maxHFp); h_calo_energyHFm_nch->Fill(temp_nch,maxHFm);
    }
    if (category == -1) keepEvent = false;
    if (keepEvent){
#ifdef readSkimData
      if (iSkimmedEntry>=(lastPrintedEntry+printing_every)){
        cout << "\r" << int(iSkimmedEntry*100.0/entries) << "%" << flush;
        lastPrintedEntry = iSkimmedEntry;
      }
#else
      if (iEntry>=(lastPrintedEntry+printing_every)){
        cout << "\r" << int(iEntry*100.0/entries) << "%" << flush;
        lastPrintedEntry = iEntry;
      }
#endif
            
      if (!category) h_tau_hadron_nPions[0]->Fill(nPions);
      h_PionMuDeltaPhi[category]->Fill(TMath::Abs(tau_muon.DeltaPhi(tau_hadron)));
      h_PionMuDeltaEta[category]->Fill(TMath::Abs(tau_muon.Eta()-tau_hadron.Eta()));
      h_PionMuDeltaR[category]->Fill(tau_muon.DeltaR(tau_hadron));
      for (int i = 0; i < (int)TREE->reco_pion_ecalEnergy->size(); i++){
        //if (TREE->reco_pion_ecalEnergy->at(i) != 0 || TREE->reco_pion_hcalEnergy->at(i) != 0)
        h_reco_pion_energy_HCAL_ECAL[0]->Fill(TREE->reco_pion_ecalEnergy->at(i),TREE->reco_pion_hcalEnergy->at(i));
      }
      h_tau_mu_p[category]->Fill(tau_muon.P());
      h_tau_mu_pz[category]->Fill(tau_muon.Pz());
      h_tau_mu_dz[category]->Fill(10*(TREE->BsTauTau_mu1_vz->at(0)-TREE->BsTauTau_PV_vz->at(0)));
      h_tau_mu_pt[category]->Fill(tau_muon.Pt());
      h_tau_mu_eta[category]->Fill(tau_muon.Eta());
      h_tau_mu_theta[category]->Fill(tau_muon.Theta());
      h_tau_mu_phi[category]->Fill(tau_muon.Phi());
      h_tau_hadron_p[category]->Fill(tau_hadron.P());
      h_tau_hadron_pz[category]->Fill(tau_hadron.Pz());
      h_tau_hadron_ptScalar[category]->Fill(tau_hadron_ptScalar);
      h_tau_hadron_pt[category]->Fill(tau_hadron.Pt());
      h_tau_hadron_eta[category]->Fill(tau_hadron.Eta());
      h_tau_hadron_theta[category]->Fill(tau_hadron.Theta());
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
        h_mueta_taueta[0]->Fill(tau_hadron.Eta(),tau_muon.Eta());
        //h_tau_hadron_track_pvz[0][0]->Fill(temp_pv_trk_1 - temp_pvz);
        //h_PV_N[0]->Fill(temp_npv);
        //h_AP[0]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2);
        h_ditau_HF_deltaphi[0]->Fill(deltaphi(ditau.Phi(),leadingHFphi));
        h_AP[0]->Fill(alphaAP,epsilonAP);
        h_muon_RF_pz_RF_pz[0]->Fill(restFrame.Pz(),muon_RF.Pz());
        h_muon_pz_pion_pz[0]->Fill(tau_hadron.Pz(),tau_muon.Pz());
      }
      
      h_MET[category]->Fill(MET.Pt());
      h_muon_cosThetaStar[category]->Fill(muon_cosThetaStar);
      h_pion_cosThetaStar[category]->Fill(pion_cosThetaStar);
      h_muon_RF_cosDeltaPhi[category]->Fill(muon_RF_cosDeltaPhi);
      h_pion_RF_cosDeltaPhi[category]->Fill(pion_RF_cosDeltaPhi);
      h_ditau_p[category]->Fill(ditau.P());
      h_ditau_pz[category]->Fill(ditau.Pz());
      h_ditau_ptScalar[category]->Fill(ditau_ptScalar);
      //cout << TREE->BsTauTau_calo_zdcSumPlus->at(0) << "," << TREE->BsTauTau_calo_zdcSumMinus->at(0) << endl;
      
      //cout << (tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()) << " ---- " << tau_muon.Pt() << endl;
      //h_MET[0]->Fill(TREE->MET_sumEt->at(0));
      

#ifdef writeSkimData
      skimmedEvent = iEntry;
      skimmedEventTree->Fill();
#endif
      
      if (!category && delta_phi > 2.78816){
        event = TREE->EVENT_event;
        run = TREE->EVENT_run;
        lumiBlock = TREE->EVENT_lumiBlock;
        deltaPhi = delta_phi;
        eventID->Fill();
      }
    } //if (keepEvent)
    
    if (passedmu && triggered && passedtau /*&& muon_charge*tauh_charge == -1*/ && passedDitauPt && passedDitauMass && passedGammaDitau && passedZDC /*&& delta_phi >= deltaPhi_cut*/){
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && !passedcalo) A_highNch_highHF[0]->Fill(delta_phi);
      if (nAccPion == 1 && passedNch && !passedcalo) B_lowNch_highHF[0]->Fill(delta_phi);
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && passedcalo) C_highNch_lowHF[0]->Fill(delta_phi);
      if (nAccPion == 1 && passedNch && passedcalo) D_lowNch_lowHF[0]->Fill(delta_phi);
      
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && !passedcaloDown) A_highNch_highHF[nSamples]->Fill(delta_phi);
      if (nAccPion == 1 && passedNch && !passedcaloDown) B_lowNch_highHF[nSamples]->Fill(delta_phi);
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && passedcaloDown) C_highNch_lowHF[nSamples]->Fill(delta_phi);
      if (nAccPion == 1 && passedNch && passedcaloDown) D_lowNch_lowHF[nSamples]->Fill(delta_phi);
      
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && !passedcaloUp) A_highNch_highHF[nSamples+1]->Fill(delta_phi);
      if (nAccPion == 1 && passedNch && !passedcaloUp) B_lowNch_highHF[nSamples+1]->Fill(delta_phi);
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && passedcaloUp) C_highNch_lowHF[nSamples+1]->Fill(delta_phi);
      if (nAccPion == 1 && passedNch && passedcaloUp) D_lowNch_lowHF[nSamples+1]->Fill(delta_phi);
      
      for (int cat = 0; cat < nNchCategories; cat++){      
        if (nAccPion == firstNchCategory+cat && !passedcalo) A_highNch_highHF[nSamples+2+cat]->Fill(delta_phi);
        if (nAccPion == 1 && passedNch && !passedcalo) B_lowNch_highHF[nSamples+2+cat]->Fill(delta_phi);
        if (nAccPion == firstNchCategory+cat && passedcalo) C_highNch_lowHF[nSamples+2+cat]->Fill(delta_phi);
        if (nAccPion == 1 && passedNch && passedcalo) D_lowNch_lowHF[nSamples+2+cat]->Fill(delta_phi);
      }
      
      
      // ABCD validation
      if (category == nSamples){
        bool ultraHighHF = (maxHFp > 20 || maxHFm > 20);
        if (temp_nch > firstNchCategory+1 && ultraHighHF) ABCD_validation[0]->Fill(delta_phi);
        if (temp_nch <= firstNchCategory+1 && ultraHighHF) ABCD_validation[1]->Fill(delta_phi);
        if (temp_nch > firstNchCategory+1 && !ultraHighHF) ABCD_validation[2]->Fill(delta_phi);
        if (temp_nch > firstNchCategory+1 && !ultraHighHF) ABCD_validation[4]->Fill(delta_phi);
        if (temp_nch <= firstNchCategory+1 && !ultraHighHF) ABCD_validation[3]->Fill(delta_phi);
        if (temp_nch <= firstNchCategory+1 && !ultraHighHF) ABCD_validation[5]->Fill(delta_phi);
      }
    }

//      aux->Fill();

  } // loop over the entries
  
  
  

  eventFile->Write();
  eventFile->Close();
  
#ifdef writeSkimData
  skimmedEventsFile->Write();
  skimmedEventsFile->Close();
#endif
  
  
  TCanvas *c_ABCD_subcategories = new TCanvas("c_ABCD_subcategories", "c_ABCD_subcategories", 800, 800);
  ABCD_validation[0]->SetMinimum(0.9);
  for (int cat = 0; cat < 4; cat++) ABCD_validation[cat]->Draw("he0e1x0same");
  CMS_lumi( c_ABCD_subcategories, iPeriod, iPos );
  c_ABCD_subcategories->SetLogy(1);
  
  TLegend *legend_ABCD_subcategories = new TLegend(0.4,0.75,0.55,0.9);
  legend_ABCD_subcategories->SetFillStyle(0); legend_ABCD_subcategories->SetFillColor(0);
  legend_ABCD_subcategories->SetLineColor(0); legend_ABCD_subcategories->SetShadowColor(0); legend_ABCD_subcategories->SetTextSize(0.035);
  legend_ABCD_subcategories->AddEntry(ABCD_validation[0],    "N_{ch} = 6,7 & HF > 20", "l");
  legend_ABCD_subcategories->AddEntry(ABCD_validation[1],    "N_{ch} = 4,5 & HF > 20", "l");
  legend_ABCD_subcategories->AddEntry(ABCD_validation[2],    "N_{ch} = 6,7 & 6 #leq HF < 20", "l");
  legend_ABCD_subcategories->AddEntry(ABCD_validation[3],    "N_{ch} = 4,5 & 6 #leq HF < 20", "l");
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
  legend_ABCD_validation->AddEntry((TObject*)0, "#frac{6 #leq HF < 20}{HF > 20}", "");
  legend_ABCD_validation->AddEntry((TObject*)0, "", "");
  legend_ABCD_validation->AddEntry(ABCD_validation[4],    "N_{ch} = 6,7", "l");
  legend_ABCD_validation->AddEntry(ABCD_validation[5],    "N_{ch} = 4,5", "l");
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
  
  /*TCanvas *c_pionSF = new TCanvas("c_pionSF", "c_pionSF", 800, 800);
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
  c_pionSF->SaveAs((basePlotDir+"/singlePlot/pion_SF.pdf").c_str());*/
  
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
  float gen_muon_pt_cut = 0;
  double temp_pv_trk_1;
  
  for (int s = 0; s < nSamples-1; s++){
  cout << "Starting " << tag[s+1] << endl;
  entriesMC = (TREEMCs[s]->fChain)->GetEntries();
  nMu3prong = cutflow[1]->GetBinContent(8);
  acceptedReco = 0;
  acceptedRecoRaw = 0;
  acceptedRecoSquare = 0;
  
  temp_counter = 0;
  for(int iEntry=0; iEntry<entriesMC; iEntry++) {
    (TREEMCs[s]->fChain)->GetEntry(iEntry);
    
    if (!(iEntry%(entriesMC/10))) cout << "\r" << int(100*iEntry/entriesMC) << "%" << flush;
    
    //if (TREEMCs[s]->gen_tautau_to_mu1prong->at(0) == 0) continue;
    bool gen_tautau_to_mu1prong = TREEMCs[s]->gen_tautau_to_mu1prong->at(0);
    
    //if (s==3 && !TREEMCs[s]->gen_tautau_to_mu3prong->at(0)) continue;
    //if (s==4 && !TREEMCs[s]->gen_tautau_to_mumu->at(0)) continue;
    if (s==1 && TREEMCs[s]->gen_tautau_to_mu1prong->at(0)) continue; 
    
    if (!threeProng[s+1]){
      gen_tautau_to_mu1prong = TREEMCs[s]->gen_tautau_to_mu1prong->at(0);
      bool gen_tautau_to_mu1prong_0npi = TREEMCs[s]->gen_tautau_to_mu1prong_0npi->at(0);
      bool gen_tautau_to_mu1prong_1npi = TREEMCs[s]->gen_tautau_to_mu1prong_1npi->at(0);
      bool gen_tautau_to_mu1prong_npis = TREEMCs[s]->gen_tautau_to_mu1prong_npis->at(0);
      bool passGen[] = {1,gen_tautau_to_mu1prong,gen_tautau_to_mu1prong_0npi,gen_tautau_to_mu1prong_1npi,gen_tautau_to_mu1prong_npis,gen_tautau_to_mu1prong};
      //if (!passGen[s]) continue;
    }
    bool triggered = TREEMCs[s]->triggered->at(0);
    
    passedmu  = false;
    passedtau = false;
    passedcalo = true;
    bool passedMET = true;
    passedNch = false;
    temp_tau_pt_comparison = 0.;
    muon_charge = 0;
    tau_hadron_ptScalar = 0;
    ditau_ptScalar = 0;
    temp_nch = 0.;
    temp_npv = 0.;
    temp_pvz = 0.;
    matchedpT_rank = 0;
    matchedTauHad = false;
    pionLeadingIndex = -1;
    pionLeadingPt = -1;
    
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
      if( (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < barrel_mu_pt_cut && TMath::Abs(TREEMCs[s]->BsTauTau_mu1_eta->at(i)) < 1.2) || (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < endcap_mu_pt_cut && TMath::Abs(TREEMCs[s]->BsTauTau_mu1_eta->at(i)) > 1.2))  continue;
      //if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }
      if(TMath::Abs(TREEMCs[s]->BsTauTau_mu1_eta->at(i)) < 2.4) passedmu = true;
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_mu1_pt->at(i);
        tau_muon.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_mu1_pt->at(i), TREEMCs[s]->BsTauTau_mu1_eta->at(i), TREEMCs[s]->BsTauTau_mu1_phi->at(i), muon_mass);
        muon_charge = TREEMCs[s]->BsTauTau_mu1_q->at(i);
      }
    }    
    
    //if(TREEMCs[s]->BsTauTau_mu1_pt->size()>1) passedmu = false;
    if(TREEMCs[s]->unfiltered_muon_pt->size()>1) passedmu = false;
    
    
    if (!TREEMCs[s]->BsTauTau_nch->empty()) temp_nch = TREEMCs[s]->BsTauTau_nch->at(0);
    int nPions=TREEMCs[s]->BsTauTau_nPions->at(0);
    
    if (!TREEMCs[s]->BsTauTau_bbPV_vz->empty()) temp_pvz = TREEMCs[s]->BsTauTau_bbPV_vz->at(0);
    //temp_npv = TREEMCs[s]->PV_N;
    
    
    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    float muon_weight_trg = tnp_weight_trg_pbpb_upc(tau_muon.Pt(),tau_muon.Eta(),0);
    float muon_weight_trk = tnp_weight_trk_pbpb(0);
    //muon_weight = muon_weight_trg*muon_weight_trk;
    muon_weight = 1; // temporary - to be replaced with above
    muon_weight_error = 0;
    candCounter = 0;
    int chosenTauIndex = -1;
    temp_tau_rho1 = 0.; temp_tau_rho2 = 0.;
    
    //if(!TREEMCs[s]->BsTauTau_tau_pt) continue;
    float vprob = -1;
    float temp_pt = 0;
    
    int nAccPion = 0;
    
    for (int i=0; i<(int)TREEMCs[s]->BsTauTau_tau_pt->size(); i++) {
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > pionExtraCut) nAccPion++;
      //if (muon_charge*TREEMCs[s]->BsTauTau_tau_q->at(i) != MuTauCharge) continue;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) < pionLeadingCut) continue;
      //if (temp_nch != 1 && TREEMCs[s]->BsTauTau_tau_pt->at(i) < 2.0) continue;
      //if(TMath::Abs(TREEMCs[s]->BsTauTau_tau_eta->at(i)) > 2) continue;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_pt) {
        vprob = 100*TREEMCs[s]->BsTauTau_B_vprob->at(i);
        if (temp_nch == 1) vprob = 100*TREEMCs[s]->BsTauTau_B_vprob->at(i);
        temp_pt = TREEMCs[s]->BsTauTau_tau_pt->at(i);
      }
      if (TREEMCs[s]->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      passedtau = true;
      //if (temp_nch == 1 && TREEMCs[s]->BsTauTau_B_vprob->at(i) < tau_hadron_vertexprob) continue;
      
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) / TREEMCs[s]->BsTauTau_tau_matched_gentaupt->at(i) > 0.85) matchedpT_rank += 1;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        tauh_charge = TREEMCs[s]->BsTauTau_tau_q->at(i);
        // Mixing tau and muon SF
        //tau_weight *= muonSF;
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        tau_hadron.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_tau_pt->at(i), TREEMCs[s]->BsTauTau_tau_eta->at(i), TREEMCs[s]->BsTauTau_tau_phi->at(i), pionMass);
        if (threeProng[s+1]){
          temp_tau_rho1 = TREEMCs[s]->BsTauTau_tau_rhomass1->at(i);
          if (temp_nch>=3) temp_tau_rho2 = TREEMCs[s]->BsTauTau_tau_rhomass2->at(i);
        }
      }
      temp_pv_trk_1 = TREEMCs[s]->BsTauTau_tau_pi1_z->at(i);
      
    } // loop over the size of the tau candidates
    //if (candCounter != 1) passedNch = false;
    
    
    //if (temp_nch>=min_nch && temp_nch<=max_nch && nPions==temp_nch) passedNch = true;
    if (nAccPion>=min_nch && nAccPion<=max_nch) passedNch = true;
    
    
    //if (triggered&&passedmu&&passedNch&&passedtau&&s) cout << "triggered + passedmu + passedNch + passedtau" << endl;
    bool passedDitauMass = true;
    bool passedDitauPt = true;
    TLorentzVector ditau = tau_muon+tau_hadron;
    TLorentzVector MET;
    MET.SetPtEtaPhiM(ditau.Pt(),-ditau.Eta(),TMath::Pi()-ditau.Phi(),0);
    TLorentzVector unboosted_muon = tau_muon + 0.5*MET;
    TLorentzVector unboosted_pion = tau_hadron + 0.5*MET;
    
    
    delta_phi = TMath::Abs(tau_muon.DeltaPhi(tau_hadron));
    float acoplanarity = 1-delta_phi/TMath::Pi();
    full_delta_phi = TMath::Abs(tau_muon.DeltaPhi(full_tau_hadron));
    float unboosted_delta_phi = TMath::Abs(unboosted_muon.DeltaPhi(unboosted_pion));
    
    
    if (ditau.M() < ditauMassCut) passedDitauMass = false;
    //if (ditau.M() < 3.0 || ditau.M() > 3.3) passedDitauMass = false;
    //if (ditau.Pt() < ditau_pt_cut) passedDitauPt = false;
    if (ditau.Pt() < ditau_pt_cut && acoplanarity < aco_cut_for_low_ditau_pt) passedDitauPt = false;
    
    
    TVector3 muon_v3 = tau_muon.Vect();
    TVector3 pion_v3 = tau_hadron.Vect();
    TVector3 restFrame = 0.5*muon_v3 + 0.5*pion_v3;
    TVector3 muon_RF = muon_v3 - restFrame;
    TVector3 pion_RF = pion_v3 - restFrame;
    
    float muon_cosThetaStar = cos(muon_RF.Angle(restFrame));
    float pion_cosThetaStar = cos(pion_RF.Angle(restFrame));
    
    float muon_RF_cosDeltaPhi = cos(muon_RF.DeltaPhi(restFrame));
    float pion_RF_cosDeltaPhi = cos(pion_RF.DeltaPhi(restFrame));
    
    float muon_AP_phi = muon_v3.Angle(restFrame);
    float pion_AP_phi = pion_v3.Angle(restFrame);
    float phiAP = muon_AP_phi+pion_AP_phi;
    float alphaAP = sin(pion_AP_phi-muon_AP_phi) / sin(phiAP);
    float epsilonAP = 2*sin(pion_AP_phi)*sin(muon_AP_phi)/sin(phiAP);
    
    //cout << "muon transverse angle: "<< 180*muon_v3.DeltaPhi(restFrame)/TMath::Pi() << "  muon RF angle: "<< 180*muon_RF.DeltaPhi(restFrame)/TMath::Pi() << endl;
    //cout << "pion transverse angle: "<< 180*pion_v3.DeltaPhi(restFrame)/TMath::Pi() << "  pion RF angle: "<< 180*pion_RF.DeltaPhi(restFrame)/TMath::Pi() << endl;
    //cout << "rest frame pz: " << restFrame.Pz() << endl;
    //cout << "muon pz: " << muon_v3.Pz() << "  muon RF pz: " << muon_RF.Pz() << endl;
    //cout << "pion pz: " << pion_v3.Pz() << "  pion RF pz: " << pion_RF.Pz() << endl;
    
    /*if(muon_RF_cosDeltaPhi < -0.95){
      cout << "muon_RF_cosDeltaPhi: " << muon_RF_cosDeltaPhi << endl;
      cout << "\tRF pt: " << restFrame.Pt() << endl;
      cout << "\tmuon RF pt: "  << muon_RF.Pt() << endl;
      cout << "\tmuon lab pt: " << muon_v3.Pt() << endl;
      cout << "\tpion lab pt: " << pion_v3.Pt() << endl;
    }*/
    
    
    
    TLorentzVector leadingGamma, leadingGenGamma, leadingFSR, matchingFSR;
    leadingGamma.SetPtEtaPhiM (0, 100, 0, 0);
    leadingGenGamma.SetPtEtaPhiM (0, 100, 0, 0);
    leadingFSR.SetPtEtaPhiM (0, 100, 0, 0);
    matchingFSR.SetPtEtaPhiM (1000, 100, 0, 0);
    
    // Applying pion SF and up/down variations
    float tauSFup[3][14];
    float tauSFdown[3][14];
    float tauSF = 1;
    float tauSFmean = 1; //should be equal to tauSF. This is a crosscheck
    
    ditau_ptScalar = tau_muon.Pt() + tau_hadron.Pt();
    
    //if (ditau.Pt() > MET_cut) passedMET = false;
    
    
    double maxHFp = TREEMCs[s]->calo_leading_energy_HFp->at(0);
    double maxHFm = TREEMCs[s]->calo_leading_energy_HFm->at(0);
    double sumHFp = TREEMCs[s]->calo_sum_energy_HFp->at(0);
    double sumHFm = TREEMCs[s]->calo_sum_energy_HFm->at(0);
    int sizeHFp = TREEMCs[s]->calo_nTowerHFp->at(0);
    int sizeHFm = TREEMCs[s]->calo_nTowerHFm->at(0);
    double maxHF = 0;
    double leadingHFeta = 0;
    double leadingHFphi = 0;
    double averageHFeta = 0;
    double sumE = TREEMCs[s]->calo_sum_energy->at(0);
    double leadingE = TREEMCs[s]->calo_leading_energy->at(0);
    double leadingE_eta = TREEMCs[s]->calo_leading_energy_eta->at(0);
    double leadingE_nonHF = 0;
    double leadingE_eta_nonHF = 1000;
    double leadingE_phi_nonHF = 1000;
    double sumEt = TREEMCs[s]->calo_sum_eT->at(0);
    double leadingEt = TREEMCs[s]->calo_leading_eT->at(0);
    double leadingEt_Eta = TREEMCs[s]->calo_leading_eT_eta->at(0);
    int nTowersMuon = 0;
    int nTowersPion = 0;
    int nTowersECALMuon = 0;
    int nTowersECALPion = 0;
    int nTowersECALFSR = 0;
    int nTowersHCALMuon = 0;
    int nTowersHCALPion = 0;
    int nTowersHCALFSR = 0;
    double ECAL_energy_FSR = 0;
    double HCAL_energy_FSR = 0;
    
    bool SR = false;
    if (passedmu && triggered && passedtau /*&& muon_charge*tauh_charge == -1*/ && passedDitauPt && passedDitauMass) SR = true;
    if (SR){
      for (int i=0; i<(int)TREEMCs[s]->EB_eta->size(); i++) {
        double eta = TREEMCs[s]->EB_eta->at(i);
        double energy = TREEMCs[s]->EB_energy->at(i);
        if (abs(eta) > 1.442) continue;
        //if (energy < 0.7) continue;
        //double deltaR_muon = sqrt(pow(deltaphi(TREEMCs[s]->EB_phi->at(i),tau_muon.Phi()),2)+pow(TREEMCs[s]->EB_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREEMCs[s]->EB_phi->at(i),tau_hadron.Phi()),2)+pow(TREEMCs[s]->EB_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREEMCs[s]->EB_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREEMCs[s]->EB_eta->at(i)-tau_hadron.Eta());
        //double deltaR_MET = sqrt(pow(deltaphi(TREEMCs[s]->EB_phi->at(i),MET.Phi()),2)+pow(TREEMCs[s]->EB_eta->at(i)-MET.Eta(),2));
        if (deltaEta_muon < muonEB_DEta) nTowersECALMuon += 1;
        if (deltaEta_pion < pionEB_DEta) nTowersECALPion += 1;
        if (1-deltaphi(TREEMCs[s]->EB_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersECALFSR += 1;
          ECAL_energy_FSR += energy;
        }
        //if (deltaR_MET < MET_EB_DEta) nTowersECALMET += 1;
      }
      for (int i=0; i<(int)TREEMCs[s]->EE_eta->size(); i++) {
        double eta = TREEMCs[s]->EE_eta->at(i);
        double energy = TREEMCs[s]->EE_energy->at(i);
        if (abs(eta) < 1.566 || abs(eta) > 2.6) continue;
        //if (energy < 3) continue;
        //double deltaR_muon = sqrt(pow(deltaphi(TREEMCs[s]->EE_phi->at(i),tau_muon.Phi()),2)+pow(TREEMCs[s]->EE_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREEMCs[s]->EE_phi->at(i),tau_hadron.Phi()),2)+pow(TREEMCs[s]->EE_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREEMCs[s]->EE_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREEMCs[s]->EE_eta->at(i)-tau_hadron.Eta());
        if (deltaEta_muon < muonEE_DEta) nTowersECALMuon += 1;
        if (deltaEta_pion < pionEE_DEta) nTowersECALPion += 1;
        if (1-deltaphi(TREEMCs[s]->EE_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersECALFSR += 1;
          ECAL_energy_FSR += energy;
        }
      }
      for (int i=0; i<(int)TREEMCs[s]->HB_eta->size(); i++) {
        double eta = TREEMCs[s]->HB_eta->at(i);
        double energy = TREEMCs[s]->HB_energy->at(i);
        if (abs(eta) > 1.305) continue;
        //if (energy < 2.8) continue; // true threshold
        //double deltaR_muon = sqrt(pow(deltaphi(TREEMCs[s]->HB_phi->at(i),tau_muon.Phi()),2)+pow(TREEMCs[s]->HB_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREEMCs[s]->HB_phi->at(i),tau_hadron.Phi()),2)+pow(TREEMCs[s]->HB_eta->at(i)-tau_hadron.Eta(),2));
        //double deltaR_pion = sqrt(pow(TREEMCs[s]->HB_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREEMCs[s]->HB_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREEMCs[s]->HB_eta->at(i)-tau_hadron.Eta());
        if (deltaEta_muon < muonHB_DEta) nTowersHCALMuon += 1;
        if (deltaEta_pion < pionHB_DEta) nTowersHCALPion += 1;
        if (1-deltaphi(TREEMCs[s]->HB_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersHCALFSR += 1;
          HCAL_energy_FSR += energy;
        }
      }
      for (int i=0; i<(int)TREEMCs[s]->HE_eta->size(); i++) {
        double eta = TREEMCs[s]->HE_eta->at(i);
        double energy = TREEMCs[s]->HE_energy->at(i);
        if (abs(eta) < 1.41 || abs(eta) > 3) continue;
        //if (energy < 1) continue; // true threshold
        //double deltaR_muon = sqrt(pow(deltaphi(TREEMCs[s]->HE_phi->at(i),tau_muon.Phi()),2)+pow(TREEMCs[s]->HE_eta->at(i)-tau_muon.Eta(),2));
        //double deltaR_pion = sqrt(pow(deltaphi(TREEMCs[s]->HE_phi->at(i),tau_hadron.Phi()),2)+pow(TREEMCs[s]->HE_eta->at(i)-tau_hadron.Eta(),2));
        //double deltaR_pion = sqrt(pow(TREEMCs[s]->HE_eta->at(i)-tau_hadron.Eta(),2));
        double deltaEta_muon = abs(TREEMCs[s]->HE_eta->at(i)-tau_muon.Eta());
        double deltaEta_pion = abs(TREEMCs[s]->HE_eta->at(i)-tau_hadron.Eta());
        if (deltaEta_muon < muonHE_DEta) nTowersHCALMuon += 1;
        if (deltaEta_pion < pionHE_DEta) nTowersHCALPion += 1;
        if (1-deltaphi(TREEMCs[s]->HE_phi->at(i),ditau.Phi())/TMath::Pi() < FSR_ditau_acoCut){
          nTowersHCALFSR += 1;
          HCAL_energy_FSR += energy;
        }
      }
    } // if(SR)
    
    
    float gen_tau_hadron_pt = -1000;
    float gen_tau_hadron_visible_pt = 0;
    TLorentzVector gen_tau_hadron_visible;
    int nAccGenPions = 0;
    int nLeadingPion = 0;
    int nMatchedPions = 0;
    float gen_muon_pt = -1;
    if (gen_tautau_to_mu1prong){
      float leadingGenPt = 0;
      for (int i=0; i<(int)TREEMCs[s]->gen_tau_daughter_pdgId->size(); i++) { // loop on gen daughters
        if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==13){
          if((TREEMCs[s]->gen_tau_daughter_pt->at(i) < barrel_mu_pt_cut && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i)) < 1.2) || (TREEMCs[s]->gen_tau_daughter_pt->at(i) < endcap_mu_pt_cut && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i)) >= 1.2)) continue;
          if(TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))>2.4) continue;
          gen_muon_pt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
        }
        if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<2.5){
        //if (TMath::Abs(TREEMCs[s]->gen_tau_daughter_pdgId->at(i))==211 && TMath::Abs(TREEMCs[s]->gen_tau_daughter_eta->at(i))<1){
          if (TREEMCs[s]->gen_tau_daughter_pt->at(i) > pionLeadingCut) nLeadingPion++;
          nAccGenPions++;
          if (leadingGenPt < TREEMCs[s]->gen_tau_daughter_pt->at(i)) leadingGenPt = TREEMCs[s]->gen_tau_daughter_pt->at(i);
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
            //if (TMath::Abs(temp_reco.Pt()/temp_gen.Pt()-1) < 0.07 && temp_gen.DeltaR(temp_reco) < 0.015){
            if (TMath::Abs(temp_reco.Pt()-temp_gen.Pt()) < 0.15 && temp_gen.DeltaR(temp_reco) < 0.015){
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
      if (gen_tautau_to_mu1prong && nLeadingPion >= 1 && gen_tau_hadron_visible.Pt() >= 2 && gen_muon_pt > gen_muon_pt_cut/*2.5*/){
        h_N_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
        acceptedGen++;
      }
      if(gen_tautau_to_mu1prong && nLeadingPion >= 1){
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
        h_N_matched_leading_genPionPt[s+1]->Fill(leadingGenPt);
      } // if passedtau and matchedTauHad and passedmu
    }
    
    TLorentzVector recoPiZero, genPiZero;
    int Ngamma = 0;
    int NgenPiZero = 0;
    float minDeltaR = 1000;
    float minDeltaM = 100000; //MeV
    SR = false;
    if (!threeProng[s+1] && passedmu && triggered && passedtau && muon_charge*tauh_charge == -1 && passedDitauPt && passedDitauMass && passedcalo && passedNch) SR = true;
    if (true){
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
                if (SR){
                  h_genPiZero_pt[s+1]->Fill(genPiZero.Pt());
                  h_genPiZero_eta[s+1]->Fill(genPiZero.Eta());
                  h_genPiZero_phi[s+1]->Fill(genPiZero.Phi());
                  h_genPiZero_deltaphi_muon[s+1]->Fill(TMath::Abs(genPiZero.DeltaPhi(tau_muon)));
                  h_genPiZero_deltaphi_pion[s+1]->Fill(TMath::Abs(genPiZero.DeltaPhi(tau_hadron)));
                }
                NgenPiZero += 1;
              }
            }
          }
          Ngamma += 1;
          //h_genGamma_pt[s+1]->Fill(TREEMCs[s]->gen_neutral_pion_daughter_pt->at(i));
          //h_genGamma_eta[s+1]->Fill(TREEMCs[s]->gen_neutral_pion_daughter_eta->at(i));
          //h_genGamma_phi[s+1]->Fill(TREEMCs[s]->gen_neutral_pion_daughter_phi->at(i));
          //h_genGamma_deltaphi_muon[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_muon)));
          //h_genGamma_deltaphi_pion[s+1]->Fill(TMath::Abs(gamma.DeltaPhi(tau_hadron)));
        }
      } // loop on gen neutral pion daughters
      //h_NgenGamma[s+1]->Fill(Ngamma);
      if (SR){
        h_NgenPiZero[s+1]->Fill(NgenPiZero);
        h_genGammasMinDeltaR[s+1]->Fill(minDeltaR);
        h_genPiZeroMinDeltaM[s+1]->Fill(minDeltaM);
      }
    }
    
    SR = false;
    if (passedmu && triggered && passedtau && muon_charge*tauh_charge == -1 && passedDitauPt && passedDitauMass && passedcalo && passedNch && !TREEMCs[s]->gen_gamma_pt->empty()) SR = true;
    if (SR){
      for (int i=0; i<(int)TREEMCs[s]->gen_gamma_pt->size(); i++) {
        if (TREEMCs[s]->gen_gamma_pt->at(i)<gammaPtCut || TMath::Abs(TREEMCs[s]->gen_gamma_eta->at(i)) > gammaEtaCut) continue;
        gamma.SetPtEtaPhiM (TREEMCs[s]->gen_gamma_pt->at(i),TREEMCs[s]->gen_gamma_eta->at(i),TREEMCs[s]->gen_gamma_phi->at(i),0);
        if (gamma.Pt() > leadingGenGamma.Pt())leadingGenGamma = gamma;
      } // loop on gen gammas, if present
      if (SR){
        h_NgenGamma[s+1]->Fill(TREEMCs[s]->gen_gamma_pt->size());
        h_residual_leadingGenGamma_ditau_pt[s+1]->Fill(leadingGenGamma.Pt()-ditau.Pt());
        h_leadingGenGamma_ditau_deltaPhi[s+1]->Fill(TMath::Abs(ditau.DeltaPhi(leadingGenGamma)));
        h_leadingGenGamma_ditau_deltaR[s+1]->Fill(MET.DeltaR(leadingGenGamma));
        h_genGamma_pt[s+1]->Fill(leadingGenGamma.Pt());
        h_genGamma_eta[s+1]->Fill(leadingGenGamma.Eta());
        h_genGamma_phi[s+1]->Fill(leadingGenGamma.Phi());
        h_genGamma_deltaphi_muon[s+1]->Fill(TMath::Abs(leadingGenGamma.DeltaPhi(tau_muon)));
        h_genGamma_deltaphi_pion[s+1]->Fill(TMath::Abs(leadingGenGamma.DeltaPhi(tau_hadron)));
      }
    } // event selection
    
    int NrecoPiZero = 0;
    minDeltaR = 1000;
    minDeltaM = 100000; //MeV
    
    bool counted = false;
    
    bool passedGammaDitau = true;
    //SR = false;
    //if (!threeProng[s+1] && passedmu && triggered && passedtau && muon_charge*tauh_charge == -1 && passedDitauPt && passedDitauMass && passedcalo && passedNch) SR = true;
    if (SR){
    Ngamma = TREEMCs[s]->BsTauTau_nGammas->at(0);
    full_tau_hadron = tau_hadron;
    for (int i=0; i<TREEMCs[s]->BsTauTau_nGammas->at(0); i++){
      gamma.SetPtEtaPhiM (TREEMCs[s]->reco_gamma_pt->at(i),TREEMCs[s]->reco_gamma_eta->at(i),TREEMCs[s]->reco_gamma_phi->at(i),0);
      if (gamma.Pt()<gammaPtCut || TMath::Abs(gamma.Eta()) > gammaEtaCut){
        Ngamma -= 1;
        continue;
      }
      if (gamma.Pt() > leadingGamma.Pt()) leadingGamma = gamma;
      if (gamma.Pt() > leadingFSR.Pt() && acop(gamma,ditau) < FSR_ditau_acoCut) leadingFSR = gamma;
      if ((ditau+gamma).Pt() < (ditau+matchingFSR).Pt() && acop(gamma,ditau) < FSR_ditau_acoCut) matchingFSR = gamma;
      
      if (TMath::Abs(gamma.Pt()-ditau.Pt()) < 5){
        //cout << temp_counter << " gamma pt: " << gamma.Pt() << " eta: " << gamma.Eta() << " phi: " << gamma.Phi() << " ditau aco:         " << int(1000*acop(gamma,ditau)) << " res pt: " << gamma.Pt()-ditau.Pt() << endl;
        counted = true;
      }
      /*if (TMath::Abs(gamma.DeltaR(tau_hadron)) > piZeroDeltaRCut){
      //if (TMath::Abs(gamma.DeltaR(tau_muon)) > piZeroDeltaRCut){
        Ngamma -= 1;
        continue;
      }*/
      //if (i == 0) recoNeutralPion = gamma;
      //else recoNeutralPion += gamma;
      for (int j=i+1; j<TREEMCs[s]->BsTauTau_nGammas->at(0); j++){
        if (TREEMCs[s]->reco_gamma_pt->at(j)<gammaPtCut || TMath::Abs(TREEMCs[s]->reco_gamma_eta->at(j)) > gammaEtaCut) continue;
        tempGamma.SetPtEtaPhiM (TREEMCs[s]->reco_gamma_pt->at(j),TREEMCs[s]->reco_gamma_eta->at(j),TREEMCs[s]->reco_gamma_phi->at(j),0);
        //if (TMath::Abs(tempGamma.DeltaPhi(tau_hadron)) > piZeroDeltaPhiCut) continue;
        //if (TMath::Abs(tempGamma.DeltaPhi(tau_muon)) > piZeroDeltaPhiCut) continue;
        //if (TMath::Abs(tempGamma.DeltaR(tau_hadron)) < piZeroDeltaRCut) continue;
        //if (TMath::Abs(tempGamma.DeltaR(tau_muon)) < piZeroDeltaRCut) continue;
        recoPiZero = gamma+tempGamma;
        if (recoPiZero.Pt() < piZeroPtCut) continue;
        //if (TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)) > piZeroDeltaPhiCut) continue;
        float deltaR = TMath::Abs(gamma.DeltaR(tempGamma));
        if (deltaR < minDeltaR) minDeltaR = deltaR;
        if (deltaR > gammaDeltaRCut) continue;
        double recoPiZeroMass = 1000*recoPiZero.M(); // MeV
        if (TMath::Abs(recoPiZeroMass-PiZeroMass) < TMath::Abs(minDeltaM)) minDeltaM = recoPiZeroMass-PiZeroMass;
        if (SR && passedcalo && passedNch){
          h_recoPiZeroDeltaM[s+1]->Fill(recoPiZeroMass-PiZeroMass);
          h_recoPiZero_pt[s+1]->Fill(recoPiZero.Pt());
          h_recoPiZero_eta[s+1]->Fill(recoPiZero.Eta());
          h_recoPiZero_phi[s+1]->Fill(recoPiZero.Phi());
          h_recoPiZero_deltaphi_muon[s+1]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_muon)));
          h_recoPiZero_deltaphi_pion[s+1]->Fill(TMath::Abs(recoPiZero.DeltaPhi(tau_hadron)));
        }
        NrecoPiZero += 1;
        if (gamma.Pt()>gammaPtCut && tempGamma.Pt()>gammaPtCut) full_tau_hadron += recoPiZero;
      } // loop on secondary gammas
    } // loop on reco gammas
    
    bool genMatchedFSR = false;
    if (TMath::Abs(leadingGenGamma.DeltaR(leadingFSR)) < 0.1 && TMath::Abs(leadingGenGamma.Pt()-leadingFSR.Pt())<0.2) genMatchedFSR = true;
    //leadingFSR=matchingFSR; // fix me
    
    float acoGammaDitau = acop(leadingFSR,ditau);
    float resGammaDitauPt = leadingFSR.Pt()-ditau.Pt();
    if (leadingFSR.Pt()==0 || !genMatchedFSR){
      acoGammaDitau = -1;
      resGammaDitauPt = -1000;
    }
    /*for (int i=0; i<2; i++){
      if(acoGammaDitau<acoGammaDitau_cut[i]){
        if(resGammaDitauPt<resDitauGammaPt_cut_up[i] && resGammaDitauPt>resDitauGammaPt_cut_down[i]) passedGammaDitau = false;
      }
    }*/
    
    if (counted) temp_counter++;
    //if (temp_counter == 20) break;
    
    if(leadingGamma.E() > gammaECut && resGammaDitauPt<resDitauGammaPt_cut_up && resGammaDitauPt>resDitauGammaPt_cut_down) passedGammaDitau = false;
    
    
    if (SR){
      h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma[s+1]->Fill(leadingFSR.Pt(),resGammaDitauPt);
      h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma[s+1]->Fill(leadingFSR.E(),resGammaDitauPt);
      h_residual_leadingRecoGamma_ditau_pt_aco[s+1]->Fill(acoGammaDitau,resGammaDitauPt);
      h_residual_leadingRecoGamma_ditau_pt[s+1]->Fill(resGammaDitauPt);
      if (leadingGenGamma.Pt() != 0) h_leadingGammaPt_reco_gen[s+1]->Fill(leadingGenGamma.Pt(),leadingFSR.Pt());
      if (TMath::Abs(leadingGenGamma.DeltaR(leadingFSR)) < 0.1 && TMath::Abs(leadingGenGamma.Pt()-leadingFSR.Pt())<0.2){ // if gen and reco FSR match
        h_eff_FSR_pt[s+1]->Fill(leadingGenGamma.Pt());
      }
      
      if (passedGammaDitau){
        if (leadingGamma.Pt() != 0){
          h_leadingRecoGamma_ditau_deltaPhi[s+1]->Fill(TMath::Abs(ditau.DeltaPhi(leadingGamma)));
          h_leadingRecoGamma_ditau_deltaR[s+1]->Fill(MET.DeltaR(leadingGamma));
          h_recoGamma_pt[s+1]->Fill(leadingGamma.Pt());
          h_recoGamma_eta[s+1]->Fill(leadingGamma.Eta());
          h_recoGamma_phi[s+1]->Fill(leadingGamma.Phi());
          h_recoGamma_deltaphi_muon[s+1]->Fill(TMath::Abs(leadingGamma.DeltaPhi(tau_muon)));
          h_recoGamma_deltaphi_pion[s+1]->Fill(TMath::Abs(leadingGamma.DeltaPhi(tau_hadron)));
          h_recoGamma_deltaR_muon[s+1]->Fill(TMath::Abs(leadingGamma.DeltaR(tau_muon)));
          h_recoGamma_deltaR_pion[s+1]->Fill(TMath::Abs(leadingGamma.DeltaR(tau_hadron)));
        }
        h_NrecoPiZero[s+1]->Fill(NrecoPiZero);
        h_recoGammasMinDeltaR[s+1]->Fill(minDeltaR);
        h_recoPiZeroMinDeltaM[s+1]->Fill(minDeltaM);
        h_NrecoGamma[s+1]->Fill(Ngamma);
      }
    } // if (SR)
    }
    
    
    
    SR = false;
    if (passedmu && triggered && passedtau /*&& muon_charge*tauh_charge == -1*/ && passedDitauPt && passedDitauMass && passedGammaDitau) SR = true;
    if (SR)
    {
      for (int i=0; i<(int)TREEMCs[s]->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREEMCs[s]->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREEMCs[s]->BsTauTau_calo_energyHFm->at(i);
        if (passedNch && SR) h_calo_E[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_energy->at(i));
        if (passedNch && SR) h_calo_E_eta[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eta->at(i),TREEMCs[s]->BsTauTau_calo_energy->at(i));
        if (passedNch && SR) h_calo_Et[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eT->at(i));
        if (passedNch && SR) h_calo_Et_eta[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eta->at(i),TREEMCs[s]->BsTauTau_calo_eT->at(i));
        double etaCal = TREEMCs[s]->BsTauTau_calo_eta->at(i);
        double phiCal = TREEMCs[s]->BsTauTau_calo_phi->at(i);
        
        if (TREEMCs[s]->BsTauTau_calo_energy->at(i) > 0){
          double deltaR_muon = sqrt(pow(deltaphi(phiCal,tau_muon.Phi()),2)+pow(etaCal-tau_muon.Eta(),2));
          double deltaR_pion = sqrt(pow(deltaphi(phiCal,tau_hadron.Phi()),2)+pow(etaCal-tau_hadron.Eta(),2));
          double deltaEta_muon = abs(etaCal-tau_muon.Eta());
          double deltaEta_pion = abs(etaCal-tau_hadron.Eta());
          if (deltaEta_muon < muon_DEta) nTowersMuon += 1;
          if (deltaEta_pion < pion_DEta) nTowersPion += 1;
          if (passedNch && passedcalo && abs(etaCal) < 5){
            h_calo_energy_muon_deltaR[s+1]->Fill(deltaR_muon,TREEMCs[s]->BsTauTau_calo_energy->at(i));
            h_calo_energy_pion_deltaR[s+1]->Fill(deltaR_pion,TREEMCs[s]->BsTauTau_calo_energy->at(i));
            h_calo_energy_muon_deltaEta[s+1]->Fill(deltaEta_muon,TREEMCs[s]->BsTauTau_calo_energy->at(i));
            h_calo_energy_pion_deltaEta[s+1]->Fill(deltaEta_pion,TREEMCs[s]->BsTauTau_calo_energy->at(i));
          }
        }
        
        if (abs(etaCal) < 2.5 && leadingE_nonHF < TREEMCs[s]->BsTauTau_calo_energy->at(i)){
          leadingE_nonHF = TREEMCs[s]->BsTauTau_calo_energy->at(i);
          leadingE_eta_nonHF = etaCal;
          leadingE_phi_nonHF = phiCal;
        }
        
        if (eHFp != -1){
          if (eHFp > maxHFp) maxHFp = eHFp;
          if (maxHFp > maxHF) {maxHF = maxHFp; leadingHFeta = etaCal; leadingHFphi = phiCal;};
          averageHFeta += eHFp * etaCal;
          if (passedNch && SR) h_calo_energyHFp[s+1]->Fill(eHFp);
          if (passedNch && SR) h_calo_HF_eta[s+1]->Fill(etaCal);
          if (passedNch && SR) h_calo_HF_energy_eta[s+1]->Fill(etaCal,eHFp);        }
        if (eHFm != -1){
          if (eHFm > maxHFm) maxHFm = eHFm;
          if (maxHFm > maxHF) {maxHF = maxHFm; leadingHFeta = etaCal; leadingHFphi = phiCal;};
          averageHFeta += eHFm * etaCal;
          if (passedNch && SR) h_calo_energyHFm[s+1]->Fill(eHFm);
          if (passedNch && SR) h_calo_HF_eta[s+1]->Fill(etaCal);
          if (passedNch && SR) h_calo_HF_energy_eta[s+1]->Fill(etaCal,eHFm);        }
      }
      averageHFeta /= (sumHFp+sumHFm);
      if (maxHFp > HFpLeading_high || maxHFm > HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcalo = false;
      if (SR){
        if (passedcalo && passedNch) h_calo_sumE[s+1]->Fill(sumE);
        if (passedcalo && passedNch) h_calo_leadingE[s+1]->Fill(leadingE);
        if (passedNch) h_calo_leadingE_eta[s+1]->Fill(leadingE_eta_nonHF,leadingE_nonHF);
        if (passedcalo && passedNch) h_calo_sumEt[s+1]->Fill(sumEt);
        if (passedcalo && passedNch) h_calo_leadingEt[s+1]->Fill(leadingEt);
        if (passedNch) h_calo_leadingEt_eta[s+1]->Fill(leadingEt_Eta,leadingEt);
        if (passedcalo && passedNch) h_nCaloTowers[s+1]->Fill(TREEMCs[s]->BsTauTau_calo_eta->size());
        if (passedcalo && passedNch) h_nTowersMuon[s+1]->Fill(nTowersMuon);
        if (passedcalo && passedNch) h_nTowersPion[s+1]->Fill(nTowersPion);
        if (passedcalo && passedNch) h_nTowersECALMuon[s+1]->Fill(nTowersECALMuon);
        if (passedcalo && passedNch) h_nTowersECALPion[s+1]->Fill(nTowersECALPion);
        if (passedcalo && passedNch && leadingFSR.Pt()!=0) h_nTowersECALFSR[s+1]->Fill(nTowersECALFSR);
        if (passedcalo && passedNch) h_nTowersHCALMuon[s+1]->Fill(nTowersHCALMuon);
        if (passedcalo && passedNch) h_nTowersHCALPion[s+1]->Fill(nTowersHCALPion);
        if (passedcalo && passedNch && leadingFSR.Pt()!=0) h_nTowersHCALFSR[s+1]->Fill(nTowersHCALFSR);
        if (passedcalo && passedNch && SR && (ECAL_energy_FSR+HCAL_energy_FSR) != 0 && leadingFSR.Pt()!=0) h_ratio_ECAL_HCAL_FSR[s+1]->Fill(ECAL_energy_FSR/(ECAL_energy_FSR+HCAL_energy_FSR));
        if (passedcalo && passedNch) h_calo_energyHFp_sum[s+1]->Fill(sumHFp,tauSF*muon_weight);
        if (passedcalo && passedNch) h_calo_energyHFm_sum[s+1]->Fill(sumHFm,tauSF*muon_weight);
        if (passedcalo && passedNch) h_calo_energyHF_pm[s+1]->Fill(sumHFm,sumHFp,tauSF*muon_weight);
        //h_calo_energyHFp_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),maxHFp); h_calo_energyHFm_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),maxHFm,tauSF*muon_weight);
        if (passedcalo && passedNch) h_calo_energyHFp_size[s+1]->Fill(sizeHFp,tauSF*muon_weight);
        if (passedcalo && passedNch) h_calo_energyHFm_size[s+1]->Fill(sizeHFm,tauSF*muon_weight);
        if (passedNch && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFp[s+1]->Fill(maxHFp,tauSF*muon_weight);
        if (passedNch && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low) h_calo_leadingHFm[s+1]->Fill(maxHFm,tauSF*muon_weight);
        if (passedNch && passedcalo) h_calo_leadingHE[s+1]->Fill(TREEMCs[s]->calo_leading_energy_HE->at(0),tauSF*muon_weight);
      if (passedNch && passedcalo && SR) h_calo_leadingHCAL[s+1]->Fill(TMath::Max(TREEMCs[s]->calo_leading_energy_HB->at(0),TREEMCs[s]->calo_leading_energy_HE->at(0)),tauSF*muon_weight);
      if (passedNch && passedcalo && SR) h_calo_leadingEE[s+1]->Fill(TREEMCs[s]->calo_leading_energy_EE->at(0),tauSF*muon_weight);
      if (passedNch && passedcalo && SR) h_calo_leadingECAL[s+1]->Fill(TMath::Max(TREEMCs[s]->calo_leading_energy_EB->at(0),TREEMCs[s]->calo_leading_energy_EE->at(0)),tauSF*muon_weight);
        if (nAccPion>=firstNchCategory && nAccPion<=firstNchCategory+nNchCategories-1 && maxHFm < HFmLeading_high && maxHFm > HFmLeading_low) h_calo_leadingHFp_highNch[s+1]->Fill(maxHFp,tauSF*muon_weight);
        if (nAccPion>=firstNchCategory && nAccPion<=firstNchCategory+nNchCategories-1 && maxHFp < HFpLeading_high && maxHFp > HFpLeading_low) h_calo_leadingHFm_highNch[s+1]->Fill(maxHFm,tauSF*muon_weight);
        if (passedNch) h_calo_leadingHF_pm[s+1]->Fill(maxHFm,maxHFp,tauSF*muon_weight);
      } // if (SR)
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    
    } // if SR
    
    
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
    
    
    if (passedmu && passedcalo && passedGamma && gen_tautau_to_mu1prong && passedtau && muon_charge*tauh_charge != 0 && delta_phi >= deltaPhi_cut && passedNch && nLeadingPion >= 1) h_N_recoMuonPt[s+1]->Fill(tau_muon.Pt());
    
    //if (passedmu && passedcalo && passedGamma && gen_tautau_to_mu1prong && passedtau && muon_charge*tauh_charge != 0 && delta_phi >= deltaPhi_cut && passedNch && nLeadingPion >= 1) h_N_recoMuonPt[s+1]->Fill(gen_muon_pt);
    
    
    
    if (passedmu && triggered && passedtau /*&& muon_charge*tauh_charge == -1*/ && passedDitauPt && passedDitauMass && passedGammaDitau /*&& passedGamma*/ && ((gen_tautau_to_mu1prong && matchedTauHad && TREEMCs[s]->gen_tau_pt->at(0) > 0 && TREEMCs[s]->gen_tau_pt->at(1) > 0 && gen_muon_pt > gen_muon_pt_cut) || (s > 0))){
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && !passedcalo) A_highNch_highHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
      if (nAccPion == 1 && passedNch  && !passedcalo) B_lowNch_highHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
      if (nAccPion >= firstNchCategory && nAccPion < firstNchCategory+nNchCategories && passedcalo) C_highNch_lowHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
      if (nAccPion == 1 && passedNch  && passedcalo) D_lowNch_lowHF[s+1]->Fill(delta_phi,tauSF*muon_weight);
    }    
    
    
    
    
    if (passedmu /*&& passedGamma*/ && triggered && passedGammaDitau && ((gen_tautau_to_mu1prong && matchedTauHad && TREEMCs[s]->gen_tau_pt->at(0) > 0 && TREEMCs[s]->gen_tau_pt->at(1) > 0 && gen_muon_pt > gen_muon_pt_cut) || (s > 0)))
    {
      if (passedtau && muon_charge*tauh_charge != 0 && delta_phi >= deltaPhi_cut && passedcalo && passedNch){
        if (passedDitauPt) h_ditau_mass[s+1]->Fill(ditau.M(),tauSF*muon_weight);
        if (passedDitauMass){
          h_ditau_pt[s+1]->Fill(ditau.Pt(),tauSF*muon_weight);
          h_ditau_pt_aco[s+1]->Fill(acoplanarity,ditau.Pt());
        }
      }
      if (!passedDitauPt || !passedDitauMass) continue;
      //if (!TREEMCs[s]->BsTauTau_tau_isRight->at(0)) continue;
      if (passedNch && passedcalo && muon_charge*tauh_charge == MuTauCharge && delta_phi >= deltaPhi_cut) h_tau_hadron_vprob[s+1]->Fill(vprob,tauSF*muon_weight);
      if (passedNch && passedcalo && delta_phi >= deltaPhi_cut) charge_counter[s+1][muon_charge*tauh_charge + 1] += 1;
      if (!passedtau && passedcalo) continue;
      if (muon_charge*tauh_charge != MuTauCharge && passedcalo && passedNch) continue;
      if(passedNch && gen_tautau_to_mu1prong && nLeadingPion >= 1 && passedcalo) h_eff_tauReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
      //if(gen_tautau_to_mu1prong && nAccGenPions == 3 && TREEMCs[s]->BsTauTau_tau_isRight->at(0)) h_eff_tauReco_genTauHadpt[s+1]->Fill(gen_tau_hadron_pt);
      if (passedNch && passedcalo){
        h_deltaphi_tau_mu_tau_hadron_zoomed[s+1]->Fill(delta_phi,tauSF*muon_weight);
        h_deltaphi_tau_mu_full_tau_hadron[s+1]->Fill(full_delta_phi,tauSF*muon_weight);
        h_deltaphi_tau_mu_tau_hadron_mueta[s+1]->Fill(delta_phi,tau_muon.Eta(),tauSF*muon_weight);
        h_deltaphi_tau_mu_tau_hadron_deltaeta[s+1]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()),tauSF*muon_weight);
      }
      
      /*
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
        muon_weight_error = sqrt(pow(muon_weight_sys_trigger,2)+pow(muon_weight_sys_tracking,2)+pow(muon_weight_sys_muid,2)+pow(muon_weight_sys_binned,2)+muon_weight_stat_trigger+muon_weight_stat_muid);
        
        
        //float muon_weight_stat = sqrt(muon_weight_stat_trigger+muon_weight_stat_muid+muon_weight_stat_sta);
        //float muon_weight_sys = sqrt(pow(muon_weight_sys_trigger,2)+pow(muon_weight_sys_tracking,2)+pow(muon_weight_sys_muid,2)+pow(muon_weight_sys_sta,2)+pow(muon_weight_sys_binned,2));
        //cout << "muon SF: " << muon_weight << " +- " << muon_weight_stat << "(stat) +- " << muon_weight_sys << "(sys)" << endl;
        
        h_sys_muon_SF_Down->Fill(delta_phi,tauSF*(muon_weight-muon_weight_error));
        h_sys_muon_SF_Up->Fill(delta_phi,tauSF*(muon_weight+muon_weight_error));
        
        // #### 3prong tau SF systematic uncertainties ####
        for (int ebin = 0; ebin < 3; ebin++){
          for (int pbin = 0; pbin < 14; pbin++){
            h_tauSFup[ebin][pbin]->Fill(delta_phi,tauSFup[ebin][pbin]*muon_weight);
            h_tauSFdown[ebin][pbin]->Fill(delta_phi,tauSFdown[ebin][pbin]*muon_weight);
          }
        }
        h_tauSFmean->Fill(delta_phi,tauSFmean*muon_weight);
        
      } // if(s==0) for muon and pion SF uncertainties */
      
      if (passedNch && passedcalo){
        h_deltaphi_tau_mu_tau_hadron[s+1]->Fill(delta_phi,tauSF*muon_weight);
        h_acoplanarity_tau_mu_tau_hadron[s+1]->Fill(acoplanarity,tauSF*muon_weight);
      }
      
      
      if(delta_phi < deltaPhi_cut && passedcalo) continue;
      if (!passedcalo) h_tau_hadron_nch_highHF[s+1]->Fill(temp_nch,tauSF*muon_weight);
      if (!passedcalo) continue;
      //if (!TREEMCs[s]->BsTauTau_nPions->empty()) temp_nch = TREEMCs[s]->BsTauTau_nPions->at(0);
      //h_tau_hadron_nch[s+1]->Fill(temp_nch,tauSF*muon_weight);
      h_tau_hadron_nch[s+1]->Fill(nAccPion,tauSF*muon_weight);
      h_tau_hadron_ncand_final[s+1]->Fill(candCounter,tauSF*muon_weight);
    
      if (!passedNch) continue;
      h_tau_hadron_nPions[s+1]->Fill(nPions);
      h_PionMuDeltaPhi[s+1]->Fill(TMath::Abs(tau_muon.DeltaPhi(tau_hadron)),tauSF*muon_weight);
      h_PionMuDeltaEta[s+1]->Fill(TMath::Abs(tau_muon.Eta()-tau_hadron.Eta()),tauSF*muon_weight);
      h_PionMuDeltaR[s+1]->Fill(tau_muon.DeltaR(tau_hadron),tauSF*muon_weight);
      //if(nLeadingPion >= 1) h_eff_trigger[s+1]->Fill(gen_muon_pt);
      if(nLeadingPion >= 1) h_eff_trigger[s+1]->Fill(tau_muon.Pt());
      for (int i = 0; i < (int)TREEMCs[s]->reco_pion_ecalEnergy->size(); i++){
        //if (TREEMCs[s]->reco_pion_ecalEnergy->at(i) != 0 || TREEMCs[s]->reco_pion_hcalEnergy->at(i) != 0)
        h_reco_pion_energy_HCAL_ECAL[s+1]->Fill(TREEMCs[s]->reco_pion_ecalEnergy->at(i),TREEMCs[s]->reco_pion_hcalEnergy->at(i));
      }
      
      h_tau_hadron_matched_pt_index[s+1]->Fill(matchedpT_rank,tauSF*muon_weight);
      h_tau_mu_p[s+1]->Fill(tau_muon.P(),tauSF*muon_weight);
      h_tau_mu_pz[s+1]->Fill(tau_muon.Pz(),tauSF*muon_weight);
      h_tau_mu_dz[s+1]->Fill(10*(TREEMCs[s]->BsTauTau_mu1_vz->at(0)-TREEMCs[s]->BsTauTau_PV_vz->at(0)),tauSF*muon_weight);
      h_tau_mu_pt[s+1]->Fill(tau_muon.Pt(),tauSF*muon_weight);
      if (s==0) h_muon_pt_before_SF->Fill(tau_muon.Pt(),tauSF);
      if (s==0) h_muon_pt_after_SF->Fill(tau_muon.Pt(),tauSF*muon_weight);
      h_tau_mu_eta[s+1]->Fill(tau_muon.Eta(),tauSF*muon_weight);
      h_tau_mu_theta[s+1]->Fill(tau_muon.Theta(),tauSF*muon_weight);
      h_tau_mu_phi[s+1]->Fill(tau_muon.Phi(),tauSF*muon_weight);
      h_tau_mu_tau_hadron_phi[s+1]->Fill(tau_muon.Phi(),tau_hadron.Phi());
      h_tau_hadron_p[s+1]->Fill(tau_hadron.P(),tauSF*muon_weight);
      h_tau_hadron_pz[s+1]->Fill(tau_hadron.Pz(),tauSF*muon_weight);
      h_tau_hadron_ptScalar[s+1]->Fill(tau_hadron_ptScalar,tauSF*muon_weight);
      h_tau_hadron_pt[s+1]->Fill(tau_hadron.Pt(),tauSF*muon_weight);
      h_gen_tau_hadron_visible_pt[s+1]->Fill(gen_tau_hadron_visible.Pt(),tauSF*muon_weight);
      //if (gen_tautau_to_mu1prong/* && TREEMCs[s]->BsTauTau_tau_isRight->at(0)*/) h_resVisTauPt[s+1]->Fill(tau_hadron.Pt() - gen_tau_hadron_visible_pt);
      h_tau_hadron_eta[s+1]->Fill(tau_hadron.Eta(),tauSF*muon_weight);
      h_tau_hadron_theta[s+1]->Fill(tau_hadron.Theta(),tauSF*muon_weight);
      h_tau_hadron_phi[s+1]->Fill(tau_hadron.Phi(),tauSF*muon_weight);
      if (gen_tautau_to_mu1prong/* && !TREEMCs[s]->BsTauTau_tau_isRight->empty()*/){
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
      h_tau_hadron_track_pvz[s+1][0]->Fill(temp_pv_trk_1 - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][1]->Fill(tau_track2.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][2]->Fill(tau_track3.Z() - temp_pvz);
      h_calo_muEta_leadingHFeta[s+1]->Fill(tau_muon.Eta(),leadingHFeta);
      h_calo_tauEta_leadingHFeta[s+1]->Fill(tau_hadron.Eta(),leadingHFeta);
      h_calo_muEta_averageHFeta[s+1]->Fill(tau_muon.Eta(),averageHFeta);
      h_calo_tauEta_averageHFeta[s+1]->Fill(tau_hadron.Eta(),averageHFeta);
      h_mueta_taueta[s+1]->Fill(tau_hadron.Eta(),tau_muon.Eta(),tauSF*muon_weight);
      h_PV_N[s+1]->Fill(temp_npv,tauSF*muon_weight);
      TLorentzVector MET;
      MET.SetPtEtaPhiM(ditau.Pt(),ditau.Eta(),-ditau.Phi(),0);
      h_MET[s+1]->Fill(MET.Pt(),tauSF*muon_weight);
      //h_ditau_mass[s+1]->Fill((ditau+MET).M(),tauSF*muon_weight);
      h_muon_RF_pz_RF_pz[s+1]->Fill(restFrame.Pz(),muon_RF.Pz());
      h_muon_pz_pion_pz[s+1]->Fill(tau_hadron.Pz(),tau_muon.Pz());
      h_muon_cosThetaStar[s+1]->Fill(muon_cosThetaStar,tauSF*muon_weight);
      h_pion_cosThetaStar[s+1]->Fill(pion_cosThetaStar,tauSF*muon_weight);
      h_muon_RF_cosDeltaPhi[s+1]->Fill(muon_RF_cosDeltaPhi,tauSF*muon_weight);
      h_pion_RF_cosDeltaPhi[s+1]->Fill(pion_RF_cosDeltaPhi,tauSF*muon_weight);
      h_ditau_p[s+1]->Fill(ditau.P(),tauSF*muon_weight);
      h_ditau_pz[s+1]->Fill(ditau.Pz(),tauSF*muon_weight);
      h_ditau_ptScalar[s+1]->Fill(ditau_ptScalar,tauSF*muon_weight);
      h_ditau_HF_deltaphi[s+1]->Fill(deltaphi(ditau.Phi(),leadingHFphi),tauSF*muon_weight);
      
      
      h_AP[s+1]->Fill(alphaAP,epsilonAP);
      //h_AP[s+1]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2,tauSF*muon_weight);
      //h_MET[s+1]->Fill(TREEMCs[s]->MET_sumEt->at(0));
      //if (matchedTauHad){
      if (true){ // Wrong statement: For efficiency, we make sure gen is solid but reco can be a fake just as in data efficiency.
        acceptedReco += tauSF*muon_weight;
        acceptedRecoSquare += pow(tauSF*muon_weight,2);
        acceptedRecoRaw++;
      }
    }
    
    
    
  } // loop on entries
  
  cout << "\nMuon charge times hadronic tau charge for " << tag[s+1] << ":\n -1: " << charge_counter[s+1][0] << " -> " << charge_counter[s+1][0]*SF[s+1] << "\n  0: " << charge_counter[s+1][1] << "\n +1: " << charge_counter[s+1][2] << endl;
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
    h_tau_mu_p[i]->Scale(SF[i]); h_tau_mu_pz[i]->Scale(SF[i]); h_tau_mu_pt[i]->Scale(SF[i]);
    h_tau_mu_eta[i]->Scale(SF[i]); h_tau_mu_theta[i]->Scale(SF[i]); h_tau_mu_phi[i]->Scale(SF[i]);
    h_tau_mu_dz[i]->Scale(SF[i]);
    h_tau_hadron_p[i]->Scale(SF[i]); h_tau_hadron_pz[i]->Scale(SF[i]); h_tau_hadron_pt[i]->Scale(SF[i]);
    //h_gen_tau_hadron_visible_pt[i]->Scale(SF[i]);
    h_gen_tau_hadron_visible_pt[i]->Scale(1./h_gen_tau_hadron_visible_pt[i]->GetMaximum());
    h_tau_hadron_eta[i]->Scale(SF[i]); h_tau_hadron_theta[i]->Scale(SF[i]); h_tau_hadron_phi[i]->Scale(SF[i]);
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
    h_acoplanarity_tau_mu_tau_hadron[i]->Scale(SF[i]);
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->Scale(SF[i]);
    h_deltaphi_tau_mu_full_tau_hadron[i]->Scale(SF[i]);
    h_PionMuDeltaPhi[i]->Scale(SF[i]);
    h_PionMuDeltaEta[i]->Scale(SF[i]);
    h_PionMuDeltaR[i]->Scale(SF[i]);
    
    h_residual_leadingGenGamma_ditau_pt[i]->Scale(SF[i]);
    h_leadingGenGamma_ditau_deltaPhi[i]->Scale(SF[i]);
    h_leadingGenGamma_ditau_deltaR[i]->Scale(SF[i]);
    h_NgenGamma[i]->Scale(SF[i]);
    h_genGamma_pt[i]->Scale(SF[i]);
    h_eff_FSR_pt[i]->Scale(SF[i]);
    h_genGamma_eta[i]->Scale(SF[i]);
    h_genGamma_phi[i]->Scale(SF[i]);
    h_genGamma_deltaphi_muon[i]->Scale(SF[i]);
    h_genGamma_deltaphi_pion[i]->Scale(SF[i]);
    h_residual_leadingRecoGamma_ditau_pt[i]->Scale(SF[i]);
    h_leadingRecoGamma_ditau_deltaPhi[i]->Scale(SF[i]);
    h_leadingRecoGamma_ditau_deltaR[i]->Scale(SF[i]);
    h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma[i]->Scale(SF[i]);
    h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma[i]->Scale(SF[i]);
    h_residual_leadingRecoGamma_ditau_pt_aco[i]->Scale(SF[i]);
    h_leadingGammaPt_reco_gen[i]->Scale(SF[i]);
    h_ditau_pt_aco[i]->Scale(SF[i]);
    h_recoGamma_pt[i]->Scale(SF[i]);
    h_recoGamma_eta[i]->Scale(SF[i]);
    h_recoGamma_phi[i]->Scale(SF[i]);
    h_recoGamma_deltaphi_muon[i]->Scale(SF[i]);
    h_recoGamma_deltaR_muon[i]->Scale(SF[i]);
    h_recoGamma_deltaphi_pion[i]->Scale(SF[i]);
    h_recoGamma_deltaR_pion[i]->Scale(SF[i]);
    
    h_reco_pion_energy_HCAL_ECAL[i]->Scale(SF[i]);
    
    //if (i == 0) h_deltaphi_tau_mu_tau_hadron[i]->Scale(0.1);
    
    h_PV_N[i]->Scale(SF[i]);
    h_calo_energyHFp[i]->Scale(SF[i]);
    //h_calo_energyHFp_sum[i]->Scale(1./h_calo_energyHFp_sum[i]->GetMaximum());
    h_calo_E[i]->Scale(SF[i]);
    h_calo_sumE[i]->Scale(SF[i]);
    h_calo_leadingE[i]->Scale(SF[i]);
    h_calo_Et[i]->Scale(SF[i]);
    h_calo_sumEt[i]->Scale(SF[i]);
    h_calo_leadingEt[i]->Scale(SF[i]);
    if (i>=0){
      h_nCaloTowers[i]->Scale(1./h_nCaloTowers[i]->GetMaximum());
      h_nTowersMuon[i]->Scale(1./h_nTowersMuon[i]->GetMaximum());
      h_nTowersPion[i]->Scale(1./h_nTowersPion[i]->GetMaximum());
      h_nTowersECALMuon[i]->Scale(1./h_nTowersECALMuon[i]->GetMaximum());
      h_nTowersECALPion[i]->Scale(1./h_nTowersECALPion[i]->GetMaximum());
      h_nTowersHCALMuon[i]->Scale(1./h_nTowersHCALMuon[i]->GetMaximum());
      h_nTowersHCALPion[i]->Scale(1./h_nTowersHCALPion[i]->GetMaximum());
      h_ratio_ECAL_HCAL_FSR[i]->Scale(1./h_ratio_ECAL_HCAL_FSR[i]->GetMaximum());
    }
    h_nTowersECALFSR[i]->Scale(SF[i]);
    h_nTowersHCALFSR[i]->Scale(SF[i]);
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
    h_calo_leadingHE[i]->Scale(SF[i]);
    h_calo_leadingHCAL[i]->Scale(SF[i]);
    h_calo_leadingEE[i]->Scale(SF[i]);
    h_calo_leadingECAL[i]->Scale(SF[i]);
    h_calo_leadingHFp_highNch[i]->Scale(SF[i]);
    h_calo_leadingHFm_highNch[i]->Scale(SF[i]);
    h_MET[i]->Scale(SF[i]);
    h_muon_cosThetaStar[i]->Scale(SF[i]);
    h_pion_cosThetaStar[i]->Scale(SF[i]);
    h_muon_RF_cosDeltaPhi[i]->Scale(SF[i]);
    h_pion_RF_cosDeltaPhi[i]->Scale(SF[i]);
    h_ditau_p[i]->Scale(SF[i]);
    h_ditau_pz[i]->Scale(SF[i]);
    h_ditau_pt[i]->Scale(SF[i]);
    //h_ditau_mass[i]->Scale(1./h_ditau_mass[i]->GetMaximum());
    h_ditau_mass[i]->Scale(SF[i]);
    //for (int bin=1; bin <=4; bin++) cutflow[i]->SetBinContent(bin,SF[i]*cutflow[i]->GetBinContent(bin));
    h_resVisTauPt[i]->Scale(SF[i]);
    h_tau_hadron_ptScalar[i]->Scale(SF[i]);
    h_ditau_ptScalar[i]->Scale(SF[i]);
    
    /*h_pion_leading_pt[i]->Scale(SF[i]);
    h_pion_leading_eta[i]->Scale(SF[i]);
    h_pion_leading_phi[i]->Scale(SF[i]);*/
    
    A_highNch_highHF[i]->Scale(SF[i]);
    B_lowNch_highHF[i]->Scale(SF[i]);
    C_highNch_lowHF[i]->Scale(SF[i]);
    D_lowNch_lowHF[i]->Scale(SF[i]);
    h_NrecoGamma[i]->Scale(SF[i]);
    h_recoPiZeroMinDeltaM[i]->Scale(SF[i]);
    h_recoGammasMinDeltaR[i]->Scale(SF[i]);
    h_NrecoPiZero[i]->Scale(SF[i]);
    h_recoPiZeroDeltaM[i]->Scale(SF[i]);
    h_recoPiZero_pt[i]->Scale(SF[i]);
    h_recoPiZero_eta[i]->Scale(SF[i]);
    h_recoPiZero_phi[i]->Scale(SF[i]);
    h_recoPiZero_deltaphi_muon[i]->Scale(SF[i]);
    h_recoPiZero_deltaphi_pion[i]->Scale(SF[i]);
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
  TH1F *post_background = new TH1F("post_background","#Delta#phi(#tau_{#mu}, #tau_{1prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  TH1F *post_signal = new TH1F("post_signal","#Delta#phi(#tau_{#mu}, #tau_{1prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  
  if ((1+deltaphi_bins-maxSBbin) != temp_post_signal->GetNbinsX()) shift = maxSBbin-1;
  for (int i = 1; i <= 1+deltaphi_bins-maxSBbin; i++){
    post_background->SetBinContent(i,temp_post_background->GetBinContent(i+shift));
    post_signal->SetBinContent(i,temp_post_signal->GetBinContent(i+shift));
  }
  nPostSignal = post_signal->Integral();
  nPostBackground = post_background->Integral();
#endif
  
  // subtracting the MC samples from control regions of deltaphi
  for (int s = 1; s < nSamples; s++){
    A_highNch_highHF[0]->Add(A_highNch_highHF[s],-1);
    B_lowNch_highHF[0]->Add(B_lowNch_highHF[s],-1);
    C_highNch_lowHF[0]->Add(C_highNch_lowHF[s],-1);
  }
  
  TH1F *background = (TH1F*)B_lowNch_highHF[0]->Clone("background");
  background->SetTitle("Estimated background from ABCD;#Delta#phi(#tau_{#mu}, #tau_{1prong});background in the signal region");
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
  
  
  TH1F *pre_background = new TH1F("pre_background","#Delta#phi(#tau_{#mu}, #tau_{1prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  TH1F *pre_signal = new TH1F("pre_signal","#Delta#phi(#tau_{#mu}, #tau_{1prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
  TH1F *pre_data = new TH1F("pre_data","#Delta#phi(#tau_{#mu}, #tau_{1prong})",1+deltaphi_bins-maxSBbin, (maxSBbin-1)*TMath::Pi()/deltaphi_bins, TMath::Pi());
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
  
  
  /*cout << "\nABCD sys prefit nch variations: ";
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
  //background_sysUp->Divide(A_highNch_highHF[nSamples+1]);*/


#ifdef plotting
  
  TCanvas *c_tau_had_pv = new TCanvas("c_tau_had_pv", "c_tau_had_pv", 1500, 500); c_tau_had_pv->Divide(3,1);
  for (int j = 0; j < 3; j++){
    for (int i = 1; i < nSamples; i++) hs_tau_hadron_track_pvz[j]->Add(h_tau_hadron_track_pvz[i][j]);
    c_tau_had_pv->cd(j+1);
    //hs_tau_hadron_track_pvz[j]->Draw("he");
    h_tau_hadron_track_pvz[0][j]->Draw("e0e1x0");
    //c_tau_had_pv->cd(j+1); h_tau_hadron_track_pvz[0][j]->Draw("e0e1x0"); h_tau_hadron_track_pvz[0][j]->GetYaxis()->SetRangeUser(0., 0.6);
    for (int i = 1; i < nSamples; i++){
      h_tau_hadron_track_pvz[i][j]->Draw("hesame");
      //h_tau_hadron_track_pvz[i][j]->GetYaxis()->SetRangeUser(0., 0.1);
    }
  }

  c_tau_had_pv->SaveAs((basePlotDir+"/collective/tau_had_pv."+plotFormat).c_str());
    
  TCanvas *c_delta_phi = new TCanvas("c_delta_phi", "c_delta_phi", 800, 800);  
  //for (int i = 1; i < nSamples; i++) hs_deltaphi_tau_mu_tau_hadron->Add(h_deltaphi_tau_mu_tau_hadron[i]);
  //if (hs_deltaphi_tau_mu_tau_hadron->GetMaximum() < h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum()) hs_deltaphi_tau_mu_tau_hadron->SetMaximum(h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum());
  //c_delta_phi->cd(1); hs_deltaphi_tau_mu_tau_hadron->Draw("he"); h_deltaphi_tau_mu_tau_hadron[0]->Draw("e1same");
  
  hs_deltaphi_tau_mu_tau_hadron->Add(background); 
  for (int sample = nSamples-1; sample > 0; sample--) hs_deltaphi_tau_mu_tau_hadron->Add(h_deltaphi_tau_mu_tau_hadron[sample]); 
  if (hs_deltaphi_tau_mu_tau_hadron->GetMaximum() > h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum()) h_deltaphi_tau_mu_tau_hadron[0]->SetMaximum(1.2*hs_deltaphi_tau_mu_tau_hadron->GetMaximum());
  /*for (int i = 1; i < nSamples; i++){
    //hs_deltaphi_tau_mu_tau_hadron->Add(h_deltaphi_tau_mu_tau_hadron[i]);
    if (h_deltaphi_tau_mu_tau_hadron[i]->GetMaximum() > h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum()) h_deltaphi_tau_mu_tau_hadron[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron[i]->GetMaximum());
    if (h_deltaphi_tau_mu_tau_hadron[i]->GetMinimum() < h_deltaphi_tau_mu_tau_hadron[0]->GetMinimum() && h_deltaphi_tau_mu_tau_hadron[i]->GetMinimum()  > 0) h_deltaphi_tau_mu_tau_hadron[0]->SetMinimum(0.8*h_deltaphi_tau_mu_tau_hadron[i]->GetMinimum());
    h_deltaphi_tau_mu_tau_hadron[i]->SetLineWidth(3);
  }*/
  
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
  
  //for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_tau_hadron[i]->Draw("hesame"); 
  hs_deltaphi_tau_mu_tau_hadron->Draw("noclearhesame");  
  
  h_deltaphi_tau_mu_tau_hadron[0]->Draw("e0e1x0same"); 
  
  CMS_lumi( c_delta_phi, iPeriod, iPos );
  
  //background->Draw("hesame");
  
  TLegend *tempGend = new TLegend(0.65,0.75,0.9,0.93);
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
  TLegend *legend = new TLegend(0.65,0.65,0.9,0.93);
  legend->SetFillStyle(0);
  legend->SetFillColor(0); legend->SetLineColor(0); legend->SetShadowColor(0); legend->SetTextSize(0.035);
  //legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    (tag[0]).c_str(), "ep");
  legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data", "ep");
  for (int i = 1; i < nSamples; i++) legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tags[i]).c_str(), "l");
  legend->AddEntry(background,    "ABCD background", "l");
  //legend->Draw();
  //gStyle->SetErrorX(0.);
  
  TLegend *nobackgend = new TLegend(0.65,0.75,0.9,0.93);
  nobackgend->SetFillStyle(0);
  nobackgend->SetFillColor(0); nobackgend->SetLineColor(0); nobackgend->SetShadowColor(0); nobackgend->SetTextSize(0.03);
  nobackgend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    "Data", "ep");
  for (int i = 1; i < nSamples; i++) nobackgend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tags[i]).c_str(), "l");
  
  TCanvas *tempCanvas;
  
  // loop on histograms
  for (int h = 0; h < histograms.size(); h++){
    /*
    // subtracting the MC samples from ABCD regions --> This is wrong since we still don't save CR distributions of the MC
    for (int region = nSamples; region <= nSamples+2; region++){
      for (int s = 1; s < nSamples; s++) histograms.at(h)[region]->Add(histograms.at(h)[s],-1);
    }
    */
    
    string histogramsName = histograms.at(h)[0]->GetName();
    histograms.at(h)[nSamples+3] = (TH1F*)histograms.at(h)[nSamples+1]->Clone(("background - "+histogramsName).c_str());
    histograms.at(h)[nSamples+3]->Multiply(histograms.at(h)[nSamples+2]);
    histograms.at(h)[nSamples+3]->Divide(histograms.at(h)[nSamples]);
    histograms.at(h)[nSamples+3]->SetLineColor(colors[nSamples]);
    histograms.at(h)[nSamples+3]->SetMarkerStyle(styles[nSamples]);
    histograms.at(h)[nSamples+4]->SetLineColor(colors[nSamples+1]);
    histograms.at(h)[nSamples+4]->SetMarkerStyle(styles[nSamples]);
    histograms.at(h)[nSamples+1]->SetBinContent(1,h); //hack to save the index of histograms. Needed later to plot the stacked histogram.
    histograms.at(h)[0]->SetBinErrorOption(TH1::kPoisson);
#ifdef post
    histograms.at(h)[nSamples+3]->Scale(nPostBackground/histograms.at(h)[nSamples+3]->Integral());
    if (nSamples > 1) histograms.at(h)[1]->Scale(nPostSignal/histograms.at(h)[1]->Integral());
#endif
    THStack *tempStack = new THStack(("stacked - "+histogramsName).c_str(),histogramsName.c_str());
    //tempStack->Add(histograms.at(h)[nSamples+4]);
    tempStack->Add(histograms.at(h)[nSamples+3]);
    for (int sample = nSamples-1; sample>0; sample--) tempStack->Add(histograms.at(h)[sample]);
    //if (nSamples > 1) tempStack->Add(histograms.at(h)[1]);
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
      //if (histogramsName == "h_deltaphi_tau_mu_tau_hadron_zoomed_"){
        //cout << val_n << " " << err_n << " " << val_d << " " << err_d << endl;
      //}
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
    
    TLegend *loopgend = new TLegend(0.65,0.75,0.9,0.93);
    loopgend->SetFillStyle(0);
    loopgend->SetFillColor(0); loopgend->SetLineColor(0); loopgend->SetShadowColor(0); loopgend->SetTextSize(0.035);
    loopgend->AddEntry(histograms.at(h)[0],    "Data", "ep");
    for (int i = 1; i < nSamples; i++) loopgend->AddEntry(histograms.at(h)[i], (tags[i]).c_str(), "l");
    loopgend->AddEntry(background,    "ABCD background", "l");
    //loopgend->AddEntry(histograms.at(h)[nSamples+4],    "same sign", "l");
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
  
  
  
  cout << endl << endl << "delta phi event count:" << endl;
  for (int i = 0; i < nSamples+5; i++) cout << h_deltaphi_tau_mu_tau_hadron[i]->Integral() << endl;
  
  
  cout << endl << "counting number of CR and SR events:" << endl;
  for (int i = 0; i < nSamples+4; i++){
    if (i >= nSamples && i <= nSamples+2) continue;
    float count_CR = 0;
    float count_SR = 0;
    for (int b = 1; b <= deltaphi_bins; b++){
      if (h_deltaphi_tau_mu_tau_hadron[i]->GetXaxis()->GetBinCenter(b) < deltaPhi_cut) count_CR += h_deltaphi_tau_mu_tau_hadron[i]->GetBinContent(b);
      else count_SR += h_deltaphi_tau_mu_tau_hadron[i]->GetBinContent(b);
    }
    if (i == 0) cout << "data: CR:  " << count_CR << "  SR:  " << count_SR << endl;
    else if (i == nSamples+3) cout << "ABCD: CR:  " << count_CR << "  SR:  " << count_SR << endl;
    else cout << tag[i] << " CR:  " << count_CR << "  SR:  " << count_SR << endl;
  }
  
  //legend->AddEntry(histograms.at(0)[nSamples+4],    "same sign", "l");
  
  
  TCanvas *c_acoplanarity = new TCanvas("c_acoplanarity", "c_acoplanarity", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_acoplanarity_tau_mu_tau_hadron[i]->GetMaximum() > h_acoplanarity_tau_mu_tau_hadron[0]->GetMaximum()) h_acoplanarity_tau_mu_tau_hadron[0]->SetMaximum(1.2*h_acoplanarity_tau_mu_tau_hadron[i]->GetMaximum());
    if (h_acoplanarity_tau_mu_tau_hadron[i]->GetMinimum() < h_acoplanarity_tau_mu_tau_hadron[0]->GetMinimum() && h_acoplanarity_tau_mu_tau_hadron[i]->GetMinimum()  > 0) h_acoplanarity_tau_mu_tau_hadron[0]->SetMinimum(0.8*h_acoplanarity_tau_mu_tau_hadron[i]->GetMinimum());
    h_acoplanarity_tau_mu_tau_hadron[i]->SetLineWidth(3);
  }
  
  
  //h_acoplanarity_tau_mu_tau_hadron[0]->SetMaximum(1.2*h_acoplanarity_tau_mu_tau_hadron[0]->GetMaximum());
  
  //h_acoplanarity_tau_mu_tau_hadron[0]->GetYaxis()->SetTitle((t+binSize+")").c_str());
  h_acoplanarity_tau_mu_tau_hadron[0]->Draw("e0e1x0");
  //c_acoplanarity->SetGrid();
  c_acoplanarity->SetLogy(1);
  
  hs_deltaphi_tau_mu_tau_hadron->Add(background); 
  for (int sample = nSamples-1; sample > 0; sample--) hs_acoplanarity_tau_mu_tau_hadron->Add(h_acoplanarity_tau_mu_tau_hadron[sample]); 
  hs_acoplanarity_tau_mu_tau_hadron->Draw("noclearhesame");  
  
  h_acoplanarity_tau_mu_tau_hadron[0]->Draw("e0e1x0same"); 
  
  //for (int i = 1; i < nSamples; i++) h_acoplanarity_tau_mu_tau_hadron[i]->Draw("hesame"); 
  
  CMS_lumi( c_acoplanarity, iPeriod, iPos );
  
  tempGend->Draw();
  
  c_acoplanarity->SaveAs((basePlotDir+"/singlePlot/acoplanarity.png").c_str());
  c_acoplanarity->SaveAs((basePlotDir+"/singlePlot/acoplanarity.pdf").c_str());
  
  c_acoplanarity->SaveAs((basePlotDir+"/collective/acoplanarity."+plotFormat).c_str());



  
  TCanvas *c_delta_phi_zoomed = new TCanvas("c_delta_phi_zoomed", "c_delta_phi_zoomed", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMaximum() > h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetMaximum()) h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMaximum());
    if (h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMinimum() < h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetMinimum() && h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMinimum()  > 0) h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMinimum(0.8*h_deltaphi_tau_mu_tau_hadron_zoomed[i]->GetMinimum());
    h_deltaphi_tau_mu_tau_hadron_zoomed[i]->SetLineWidth(3);
  }
  h_deltaphi_tau_mu_tau_hadron_zoomed[1]->SetMaximum(140);
  //h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMaximum(1.2*h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetMaximum());
  //h_deltaphi_tau_mu_tau_hadron_zoomed[0]->SetMinimum(0.005);
  
  h_deltaphi_tau_mu_tau_hadron_zoomed[0]->GetYaxis()->SetTitle((t+binSize+")").c_str());
  //h_deltaphi_tau_mu_tau_hadron_zoomed[0]->Draw("e0e1x0");
  //c_delta_phi_zoomed->SetGrid();
  //c_delta_phi_zoomed->SetLogy(1);
  for (int i = 1; i < nSamples; i++) h_deltaphi_tau_mu_tau_hadron_zoomed[i]->Draw("hesame"); 
  
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
  
  
  TCanvas *c_PionMuDeltaPhi = new TCanvas("c_PionMuDeltaPhi", "c_PionMuDeltaPhi", 800, 800);
  h_PionMuDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_PionMuDeltaPhi[i]->Draw("hesame");
  stacks.at(h_PionMuDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_PionMuDeltaPhi, iPeriod, iPos );
  c_PionMuDeltaPhi->SaveAs((basePlotDir+"/singlePlot/PionMuDeltaPhi.png").c_str());
  c_PionMuDeltaPhi->SaveAs((basePlotDir+"/singlePlot/PionMuDeltaPhi.pdf").c_str());
  
  TCanvas *c_PionMuDeltaEta = new TCanvas("c_PionMuDeltaEta", "c_PionMuDeltaEta", 800, 800);
  h_PionMuDeltaEta[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_PionMuDeltaEta[i]->Draw("hesame");
  stacks.at(h_PionMuDeltaEta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_PionMuDeltaEta, iPeriod, iPos );
  c_PionMuDeltaEta->SaveAs((basePlotDir+"/singlePlot/PionMuDeltaEta.png").c_str());
  c_PionMuDeltaEta->SaveAs((basePlotDir+"/singlePlot/PionMuDeltaEta.pdf").c_str());
  
  TCanvas *c_PionMuDeltaR = new TCanvas("c_PionMuDeltaR", "c_PionMuDeltaR", 800, 800);
  h_PionMuDeltaR[0]->Draw("e0e1x0"); legend->Draw();
  //for (int i = 1; i < nSamples; i++) h_PionMuDeltaR[i]->Draw("hesame");
  stacks.at(h_PionMuDeltaR[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_PionMuDeltaR, iPeriod, iPos );
  c_PionMuDeltaR->SaveAs((basePlotDir+"/singlePlot/PionMuDeltaR.png").c_str());
  c_PionMuDeltaR->SaveAs((basePlotDir+"/singlePlot/PionMuDeltaR.pdf").c_str());
  
  
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

  TCanvas *c_AP = new TCanvas("c_AP", "c_AP", 800*nSamples, 800); c_AP->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_AP->cd(i+1); h_AP[i]->Draw("colz");}
  c_AP->SaveAs((basePlotDir+"/collective/Armenteros_Podolansky."+plotFormat).c_str());
  
  
  TCanvas *c_ABCD_ditau_mass = new TCanvas("c_ABCD_ditau_mass", "c_ABCD_ditau_mass", 1600, 2400); c_ABCD_ditau_mass->Divide(2,3);
  
  c_ABCD_ditau_mass->cd(1); h_ditau_mass[nSamples]->Draw("e0e1x0");
  
  c_ABCD_ditau_mass->cd(2); h_ditau_mass[nSamples+1]->Draw("e0e1x0");
  
  c_ABCD_ditau_mass->cd(3); h_ditau_mass[nSamples+2]->Draw("e0e1x0");
  
  c_ABCD_ditau_mass->cd(4); h_ditau_mass[0]->Draw("e0e1x0");
  
  c_ABCD_ditau_mass->cd(5); h_ditau_mass[nSamples+3]->Draw("e0e1x0");
  
  c_ABCD_ditau_mass->SaveAs((basePlotDir+"/collective/ABCD_ditau_mass."+plotFormat).c_str());
  
  
  TCanvas *c_ABCD_ditau_pt = new TCanvas("c_ABCD_ditau_pt", "c_ABCD_ditau_pt", 1600, 2400); c_ABCD_ditau_pt->Divide(2,3);
  
  c_ABCD_ditau_pt->cd(1); h_ditau_pt[nSamples]->Draw("e0e1x0");
  
  c_ABCD_ditau_pt->cd(2); h_ditau_pt[nSamples+1]->Draw("e0e1x0");
  
  c_ABCD_ditau_pt->cd(3); h_ditau_pt[nSamples+2]->Draw("e0e1x0");
  
  c_ABCD_ditau_pt->cd(4); h_ditau_pt[0]->Draw("e0e1x0");
  
  c_ABCD_ditau_pt->cd(5); h_ditau_pt[nSamples+3]->Draw("e0e1x0");
  
  c_ABCD_ditau_pt->SaveAs((basePlotDir+"/collective/ABCD_ditau_pt."+plotFormat).c_str());
  
  
  TCanvas *c_ditau_HF_deltaphi = new TCanvas("c_ditau_HF_deltaphi", "c_ditau_HF_deltaphi", 800, 800);
  h_ditau_HF_deltaphi[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++) h_ditau_HF_deltaphi[i]->Draw("hesame");
  //c_ditau_HF_deltaphi->SetLogy(1);
  CMS_lumi( c_ditau_HF_deltaphi, iPeriod, iPos );
  c_ditau_HF_deltaphi->SaveAs((basePlotDir+"/singlePlot/ditau_HF_deltaphi.png").c_str());
  c_ditau_HF_deltaphi->SaveAs((basePlotDir+"/singlePlot/ditau_HF_deltaphi.pdf").c_str());

  TCanvas *c_muon_RF_pz_RF_pz = new TCanvas("c_muon_RF_pz_RF_pz", "c_muon_RF_pz_RF_pz", 800*nSamples, 800); c_muon_RF_pz_RF_pz->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_muon_RF_pz_RF_pz->cd(i+1); h_muon_RF_pz_RF_pz[i]->Draw("colz");}
  c_muon_RF_pz_RF_pz->SaveAs((basePlotDir+"/collective/muon_RF_pz_RF_pz."+plotFormat).c_str());

  TCanvas *c_muon_pz_pion_pz = new TCanvas("c_muon_pz_pion_pz", "c_muon_pz_pion_pz", 800*nSamples, 800); c_muon_pz_pion_pz->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_muon_pz_pion_pz->cd(i+1); h_muon_pz_pion_pz[i]->Draw("colz");}
  c_muon_pz_pion_pz->SaveAs((basePlotDir+"/collective/muon_pz_pion_pz."+plotFormat).c_str());
  
  TCanvas *c_muon_cosThetaStar = new TCanvas("c_muon_cosThetaStar", "c_muon_cosThetaStar", 800, 800);
  h_muon_cosThetaStar[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_muon_cosThetaStar[i]->Draw("hesame");
  //stacks.at(h_muon_cosThetaStar[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_muon_cosThetaStar, iPeriod, iPos );
  c_muon_cosThetaStar->SaveAs((basePlotDir+"/singlePlot/muon_cosThetaStar.png").c_str());
  c_muon_cosThetaStar->SaveAs((basePlotDir+"/singlePlot/muon_cosThetaStar.pdf").c_str());
  
  TCanvas *c_pion_cosThetaStar = new TCanvas("c_pion_cosThetaStar", "c_pion_cosThetaStar", 800, 800);
  h_pion_cosThetaStar[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_cosThetaStar[i]->Draw("hesame");
  //stacks.at(h_pion_cosThetaStar[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_cosThetaStar, iPeriod, iPos );
  c_pion_cosThetaStar->SaveAs((basePlotDir+"/singlePlot/pion_cosThetaStar.png").c_str());
  c_pion_cosThetaStar->SaveAs((basePlotDir+"/singlePlot/pion_cosThetaStar.pdf").c_str());
  
  TCanvas *c_muon_RF_cosDeltaPhi = new TCanvas("c_muon_RF_cosDeltaPhi", "c_muon_RF_cosDeltaPhi", 800, 800);
  h_muon_RF_cosDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_muon_RF_cosDeltaPhi[i]->Draw("hesame");
  //stacks.at(h_muon_RF_cosDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_muon_RF_cosDeltaPhi, iPeriod, iPos );
  c_muon_RF_cosDeltaPhi->SaveAs((basePlotDir+"/singlePlot/muon_RF_cosDeltaPhi.png").c_str());
  c_muon_RF_cosDeltaPhi->SaveAs((basePlotDir+"/singlePlot/muon_RF_cosDeltaPhi.pdf").c_str());
  
  TCanvas *c_pion_RF_cosDeltaPhi = new TCanvas("c_pion_RF_cosDeltaPhi", "c_pion_RF_cosDeltaPhi", 800, 800);
  h_pion_RF_cosDeltaPhi[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_pion_RF_cosDeltaPhi[i]->Draw("hesame");
  //stacks.at(h_pion_RF_cosDeltaPhi[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_pion_RF_cosDeltaPhi, iPeriod, iPos );
  c_pion_RF_cosDeltaPhi->SaveAs((basePlotDir+"/singlePlot/pion_RF_cosDeltaPhi.png").c_str());
  c_pion_RF_cosDeltaPhi->SaveAs((basePlotDir+"/singlePlot/pion_RF_cosDeltaPhi.pdf").c_str());
  
  
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
  

  TCanvas *c_residual_leadingGenGamma_ditau_pt = new TCanvas("c_residual_leadingGenGamma_ditau_pt", "c_residual_leadingGenGamma_ditau_pt", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_residual_leadingGenGamma_ditau_pt[i]->GetMaximum() > h_residual_leadingGenGamma_ditau_pt[1]->GetMaximum()) h_residual_leadingGenGamma_ditau_pt[1]->SetMaximum(1.2*h_residual_leadingGenGamma_ditau_pt[i]->GetMaximum());
    h_residual_leadingGenGamma_ditau_pt[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_residual_leadingGenGamma_ditau_pt, iPeriod, iPos );
  c_residual_leadingGenGamma_ditau_pt->SaveAs((basePlotDir+"/singlePlot/residual_leadingGenGamma_ditau_pt.png").c_str());
  c_residual_leadingGenGamma_ditau_pt->SaveAs((basePlotDir+"/singlePlot/residual_leadingGenGamma_ditau_pt.pdf").c_str());
  

  TCanvas *c_leadingGenGamma_ditau_deltaPhi = new TCanvas("c_leadingGenGamma_ditau_deltaPhi", "c_leadingGenGamma_ditau_deltaPhi", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_leadingGenGamma_ditau_deltaPhi[i]->GetMaximum() > h_leadingGenGamma_ditau_deltaPhi[1]->GetMaximum()) h_leadingGenGamma_ditau_deltaPhi[1]->SetMaximum(1.2*h_leadingGenGamma_ditau_deltaPhi[i]->GetMaximum());
    h_leadingGenGamma_ditau_deltaPhi[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_leadingGenGamma_ditau_deltaPhi, iPeriod, iPos );
  c_leadingGenGamma_ditau_deltaPhi->SetLogy(1);
  c_leadingGenGamma_ditau_deltaPhi->SaveAs((basePlotDir+"/singlePlot/leadingGenGamma_ditau_deltaPhi.png").c_str());
  c_leadingGenGamma_ditau_deltaPhi->SaveAs((basePlotDir+"/singlePlot/leadingGenGamma_ditau_deltaPhi.pdf").c_str());
  

  TCanvas *c_leadingGenGamma_ditau_deltaR = new TCanvas("c_leadingGenGamma_ditau_deltaR", "c_leadingGenGamma_ditau_deltaR", 800, 800);
  for (int i = 1; i < nSamples; i++){
    if (h_leadingGenGamma_ditau_deltaR[i]->GetMaximum() > h_leadingGenGamma_ditau_deltaR[1]->GetMaximum()) h_leadingGenGamma_ditau_deltaR[1]->SetMaximum(1.2*h_leadingGenGamma_ditau_deltaR[i]->GetMaximum());
    h_leadingGenGamma_ditau_deltaR[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_leadingGenGamma_ditau_deltaR, iPeriod, iPos );
  c_leadingGenGamma_ditau_deltaR->SaveAs((basePlotDir+"/singlePlot/leadingGenGamma_ditau_deltaR.png").c_str());
  c_leadingGenGamma_ditau_deltaR->SaveAs((basePlotDir+"/singlePlot/leadingGenGamma_ditau_deltaR.pdf").c_str());
  
  TCanvas *c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma = new TCanvas("c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma", "c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma", nSamples*800, 800); c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma->cd(i+1);
    h_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma[i]->Draw("colz");
  }
  c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma->SaveAs((basePlotDir+"/collective/residual_leadingRecoGamma_ditau_pt_leadingRecoGamma.png").c_str());
  c_residual_leadingRecoGamma_ditau_pt_leadingRecoGamma->SaveAs((basePlotDir+"/collective/residual_leadingRecoGamma_ditau_pt_leadingRecoGamma.pdf").c_str());
  
  TCanvas *c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma = new TCanvas("c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma", "c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma", nSamples*800, 800); c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma->cd(i+1);
    h_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma[i]->Draw("colz");
  }
  c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma->SaveAs((basePlotDir+"/collective/residual_leadingRecoGamma_ditau_energy_leadingRecoGamma.png").c_str());
  c_residual_leadingRecoGamma_ditau_energy_leadingRecoGamma->SaveAs((basePlotDir+"/collective/residual_leadingRecoGamma_ditau_energy_leadingRecoGamma.pdf").c_str());
  
  TCanvas *c_residual_leadingRecoGamma_ditau_pt_aco = new TCanvas("c_residual_leadingRecoGamma_ditau_pt_aco", "c_residual_leadingRecoGamma_ditau_pt_aco", nSamples*800, 800); c_residual_leadingRecoGamma_ditau_pt_aco->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_residual_leadingRecoGamma_ditau_pt_aco->cd(i+1);
    h_residual_leadingRecoGamma_ditau_pt_aco[i]->Draw("colz");
  }
  c_residual_leadingRecoGamma_ditau_pt_aco->SaveAs((basePlotDir+"/collective/residual_leadingRecoGamma_ditau_pt_aco.png").c_str());
  c_residual_leadingRecoGamma_ditau_pt_aco->SaveAs((basePlotDir+"/collective/residual_leadingRecoGamma_ditau_pt_aco.pdf").c_str());
  
  TCanvas *c_leadingGammaPt_reco_gen = new TCanvas("c_leadingGammaPt_reco_gen", "c_leadingGammaPt_reco_gen", nSamples*800, 800); c_leadingGammaPt_reco_gen->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_leadingGammaPt_reco_gen->cd(i+1);
    h_leadingGammaPt_reco_gen[i]->Draw("colz");
  }
  c_leadingGammaPt_reco_gen->SaveAs((basePlotDir+"/collective/leadingGammaPt_reco_gen.png").c_str());
  c_leadingGammaPt_reco_gen->SaveAs((basePlotDir+"/collective/leadingGammaPt_reco_gen.pdf").c_str());
  
  TCanvas *c_eff_FSR_pt = new TCanvas("c_eff_FSR_pt", "c_eff_FSR_pt", (nSamples-1)*800, 800); c_eff_FSR_pt->Divide(nSamples-1,1);
  for (int i = 1; i < nSamples; i++){
    c_eff_FSR_pt->cd(i);
    h_eff_FSR_pt[i]->Divide(h_genGamma_pt[i]);
    h_eff_FSR_pt[i]->Draw("hesame");
  }
  c_eff_FSR_pt->SaveAs((basePlotDir+"/collective/eff_FSR_pt.png").c_str());
  c_eff_FSR_pt->SaveAs((basePlotDir+"/collective/eff_FSR_pt.pdf").c_str());
  
  
  TCanvas *c_ditau_pt_aco = new TCanvas("c_ditau_pt_aco", "c_ditau_pt_aco", nSamples*800, 800); c_ditau_pt_aco->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_ditau_pt_aco->cd(i+1);
    c_ditau_pt_aco->cd(i+1)->SetLogz(1);
    h_ditau_pt_aco[i]->Draw("colz");
  }
  c_ditau_pt_aco->SaveAs((basePlotDir+"/collective/ditau_pt_aco.png").c_str());
  c_ditau_pt_aco->SaveAs((basePlotDir+"/collective/ditau_pt_aco.pdf").c_str());

  TCanvas *c_residual_leadingRecoGamma_ditau_pt = new TCanvas("c_residual_leadingRecoGamma_ditau_pt", "c_residual_leadingRecoGamma_ditau_pt", 800, 800);
  h_residual_leadingRecoGamma_ditau_pt[0]->Draw("e0e1x0");
  /*for (int i = 1; i < nSamples; i++){
    if (h_residual_leadingRecoGamma_ditau_pt[i]->GetMaximum() > h_residual_leadingRecoGamma_ditau_pt[0]->GetMaximum()) h_residual_leadingRecoGamma_ditau_pt[0]->SetMaximum(1.2*h_residual_leadingRecoGamma_ditau_pt[i]->GetMaximum());
    h_residual_leadingRecoGamma_ditau_pt[i]->Draw("hesame");
  }
  nobackgend->Draw();*/
  stacks.at(h_residual_leadingRecoGamma_ditau_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  legend->Draw();
  CMS_lumi( c_residual_leadingRecoGamma_ditau_pt, iPeriod, iPos );
  c_residual_leadingRecoGamma_ditau_pt->SaveAs((basePlotDir+"/singlePlot/residual_leadingRecoGamma_ditau_pt.png").c_str());
  c_residual_leadingRecoGamma_ditau_pt->SaveAs((basePlotDir+"/singlePlot/residual_leadingRecoGamma_ditau_pt.pdf").c_str());
  

  TCanvas *c_leadingRecoGamma_ditau_deltaPhi = new TCanvas("c_leadingRecoGamma_ditau_deltaPhi", "c_leadingRecoGamma_ditau_deltaPhi", 800, 800);
  h_leadingRecoGamma_ditau_deltaPhi[0]->Draw("e0e1x0");
  for (int i = 1; i < nSamples; i++){
    if (h_leadingRecoGamma_ditau_deltaPhi[i]->GetMaximum() > h_leadingRecoGamma_ditau_deltaPhi[0]->GetMaximum()) h_leadingRecoGamma_ditau_deltaPhi[0]->SetMaximum(1.2*h_leadingRecoGamma_ditau_deltaPhi[i]->GetMaximum());
    h_leadingRecoGamma_ditau_deltaPhi[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_leadingRecoGamma_ditau_deltaPhi, iPeriod, iPos );
  c_leadingRecoGamma_ditau_deltaPhi->SetLogy(1);
  c_leadingRecoGamma_ditau_deltaPhi->SaveAs((basePlotDir+"/singlePlot/leadingRecoGamma_ditau_deltaPhi.png").c_str());
  c_leadingRecoGamma_ditau_deltaPhi->SaveAs((basePlotDir+"/singlePlot/leadingRecoGamma_ditau_deltaPhi.pdf").c_str());
  

  TCanvas *c_leadingRecoGamma_ditau_deltaR = new TCanvas("c_leadingRecoGamma_ditau_deltaR", "c_leadingRecoGamma_ditau_deltaR", 800, 800);
  h_leadingRecoGamma_ditau_deltaR[0]->Draw("e0e1x0");
  for (int i = 1; i < nSamples; i++){
    if (h_leadingRecoGamma_ditau_deltaR[i]->GetMaximum() > h_leadingRecoGamma_ditau_deltaR[0]->GetMaximum()) h_leadingRecoGamma_ditau_deltaR[0]->SetMaximum(1.2*h_leadingRecoGamma_ditau_deltaR[i]->GetMaximum());
    h_leadingRecoGamma_ditau_deltaR[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_leadingRecoGamma_ditau_deltaR, iPeriod, iPos );
  c_leadingRecoGamma_ditau_deltaR->SaveAs((basePlotDir+"/singlePlot/leadingRecoGamma_ditau_deltaR.png").c_str());
  c_leadingRecoGamma_ditau_deltaR->SaveAs((basePlotDir+"/singlePlot/leadingRecoGamma_ditau_deltaR.pdf").c_str());

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
  /*
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
  
  } // if nSamples == 1
  */
  
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
  //h_tau_mu_pt[1]->SetMaximum(230);
  //for (int i = 1; i < nSamples; i++) h_tau_mu_pt[i]->Draw("hesame");
  cout << endl << endl << "tau muon event count:" << endl;
  for (int i = 0; i < nSamples+5; i++) cout << h_tau_mu_pt[i]->Integral() << endl;
  //tempGend->Draw();
  //h_tau_mu_pt[nSamples+3]->Draw("hesame");
  stacks.at(h_tau_mu_pt[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_pt, iPeriod, iPos );
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt.png").c_str());
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt.pdf").c_str());
  c_tau_muon_pt->SaveAs((basePlotDir+"/singlePlot/tau_muon_pt.C").c_str());
  
  
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
    if (h_tau_mu_eta[i]->GetMinimum() < h_tau_mu_eta[0]->GetMinimum()) h_tau_mu_eta[0]->SetMinimum(0.8*h_tau_mu_eta[i]->GetMinimum());
  }
  //h_tau_mu_eta[0]->SetMaximum(1.2*h_tau_mu_eta[0]->GetMaximum());
  h_tau_mu_eta[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_eta[i]->Draw("hesame");
  
  
  TCanvas *c_tau_muon_eta = new TCanvas("c_tau_muon_eta", "c_tau_muon_eta", 800, 800);
  h_tau_mu_eta[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_eta[i]->Draw("hesame");
  //stacks.at(h_tau_mu_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_eta, iPeriod, iPos );
  c_tau_muon_eta->SaveAs((basePlotDir+"/singlePlot/tau_muon_eta.png").c_str());
  c_tau_muon_eta->SaveAs((basePlotDir+"/singlePlot/tau_muon_eta.pdf").c_str());
  
  TCanvas *c_tau_muon_theta = new TCanvas("c_tau_muon_theta", "c_tau_muon_theta", 800, 800);
  h_tau_mu_theta[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_mu_theta[i]->Draw("hesame");
  //stacks.at(h_tau_mu_theta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_muon_theta, iPeriod, iPos );
  c_tau_muon_theta->SaveAs((basePlotDir+"/singlePlot/tau_muon_theta.png").c_str());
  c_tau_muon_theta->SaveAs((basePlotDir+"/singlePlot/tau_muon_theta.pdf").c_str());
  
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
  
  

  /*TCanvas *c_tau_muon_tau_hadron_phi_data = new TCanvas("c_tau_muon_tau_hadron_phi_data", "c_tau_muon_tau_hadron_phi_data", 800, 800);
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
  c_tau_muon_tau_hadron_phi_data->SaveAs((basePlotDir+"/singlePlot/tau_muon_tau_hadron_phi_data.pdf").c_str());*/
  
  /*TCanvas *c_tau_muon_tau_hadron_phi_signal = new TCanvas("c_tau_muon_tau_hadron_phi_signal", "c_tau_muon_tau_hadron_phi_signal", 800, 800);
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
  c_tau_muon_tau_hadron_phi_signal->SaveAs((basePlotDir+"/singlePlot/tau_muon_tau_hadron_phi_signal.pdf").c_str());*/

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
  for (int i = 1; i < nSamples; i++) h_tau_hadron_eta[i]->Draw("hesame");
  //stacks.at(h_tau_hadron_eta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_eta, iPeriod, iPos );
  c_tau_hadron_eta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_eta.png").c_str());
  c_tau_hadron_eta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_eta.pdf").c_str());
  
  TCanvas *c_tau_hadron_theta = new TCanvas("c_tau_hadron_theta", "c_tau_hadron_theta", 800, 800);
  h_tau_hadron_theta[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_tau_hadron_theta[i]->Draw("hesame");
  //stacks.at(h_tau_hadron_theta[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_tau_hadron_theta, iPeriod, iPos );
  c_tau_hadron_theta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_theta.png").c_str());
  c_tau_hadron_theta->SaveAs((basePlotDir+"/singlePlot/tau_hadron_theta.pdf").c_str());
  
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
  //h_tau_hadron_mass[0]->SetTitle("#tau_{1prong} visible mass");
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
  
  cout << endl << endl << "ditau mass event count:" << endl;
  for (int i = 0; i < nSamples+5; i++) cout << h_ditau_mass[i]->Integral(6,36) << endl;
  
  stacks.at(h_ditau_mass[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_ditau_mass, iPeriod, iPos );
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass.png").c_str());
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass.pdf").c_str());
  c_ditau_mass->SaveAs((basePlotDir+"/singlePlot/ditau_mass.C").c_str());
  
  TCanvas *c_ditau_ptScalar = new TCanvas("c_ditau_ptScalar", "c_ditau_ptScalar", 800, 800);
  h_ditau_ptScalar[0]->Draw("e0e1x0"); legend->Draw();
  for (int i = 1; i < nSamples; i++) h_ditau_ptScalar[i]->Draw("hesame");
  //stacks.at(h_ditau_ptScalar[nSamples+1]->GetBinContent(1))->Draw("noclearhesame");
  CMS_lumi( c_ditau_ptScalar, iPeriod, iPos );
  c_ditau_ptScalar->SaveAs((basePlotDir+"/singlePlot/ditau_scalar_pt.png").c_str());
  c_ditau_ptScalar->SaveAs((basePlotDir+"/singlePlot/ditau_scalar_pt.pdf").c_str());
  
  //TH1F *stackedDitauScalarPt = new TH1F(*(TH1F*)((stacks.at(h_ditau_ptScalar[nSamples+1]->GetBinContent(1))->GetStack()->Last())));
  //cout << " **** **** **** ****" << endl << "KS test for ditau scalar pt: " << h_ditau_ptScalar[0]->KolmogorovTest(stackedDitauScalarPt) << endl << endl;
  
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
  c_ditau_pt->SetLogy(1);
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
  
  TCanvas *c_reco_pion_energy_HCAL_ECAL = new TCanvas("c_reco_pion_energy_HCAL_ECAL", "c_reco_pion_energy_HCAL_ECAL", nSamples*800, 800); c_reco_pion_energy_HCAL_ECAL->Divide(nSamples,1);
  for (int s = 0; s < nSamples; s++){
    c_reco_pion_energy_HCAL_ECAL->cd(s+1);
    c_reco_pion_energy_HCAL_ECAL->cd(s+1)->SetLogz(1);
    h_reco_pion_energy_HCAL_ECAL[s]->Draw("colz");
  }
  c_reco_pion_energy_HCAL_ECAL->SaveAs((basePlotDir+"/collective/reco_pion_energy_HCAL_ECAL.png").c_str());
  c_reco_pion_energy_HCAL_ECAL->SaveAs((basePlotDir+"/collective/reco_pion_energy_HCAL_ECAL.pdf").c_str());
  
  /*TCanvas *c_recoPion = new TCanvas("c_recoPion", "c_recoPion", 3*800, 3*800); c_recoPion->Divide(3,3);
  
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
  
  c_recoPion->SaveAs((basePlotDir+"/collective/reco_pion."+plotFormat).c_str());*/

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
  
  TCanvas *c_calo_leadingHE = new TCanvas("c_calo_leadingHE", "c_calo_leadingHE", 800, 800);
  h_calo_leadingHE[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    if (h_calo_leadingHE[i]->GetMaximum() > h_calo_leadingHE[0]->GetMaximum()) h_calo_leadingHE[0]->SetMaximum(1.2*h_calo_leadingHE[i]->GetMaximum());
    h_calo_leadingHE[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_leadingHE, iPeriod, iPos );
  c_calo_leadingHE->cd()->SetLogy(1);
  c_calo_leadingHE->SaveAs((basePlotDir+"/singlePlot/calo_leadingHE.png").c_str());
  c_calo_leadingHE->SaveAs((basePlotDir+"/singlePlot/calo_leadingHE.pdf").c_str());
  
  TCanvas *c_calo_leadingHCAL = new TCanvas("c_calo_leadingHCAL", "c_calo_leadingHCAL", 800, 800);
  h_calo_leadingHCAL[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    if (h_calo_leadingHCAL[i]->GetMaximum() > h_calo_leadingHCAL[0]->GetMaximum()) h_calo_leadingHCAL[0]->SetMaximum(1.2*h_calo_leadingHCAL[i]->GetMaximum());
    h_calo_leadingHCAL[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_leadingHCAL, iPeriod, iPos );
  c_calo_leadingHCAL->cd()->SetLogy(1);
  c_calo_leadingHCAL->SaveAs((basePlotDir+"/singlePlot/calo_leadingHCAL.png").c_str());
  c_calo_leadingHCAL->SaveAs((basePlotDir+"/singlePlot/calo_leadingHCAL.pdf").c_str());
  
  TCanvas *c_calo_leadingEE = new TCanvas("c_calo_leadingEE", "c_calo_leadingEE", 800, 800);
  h_calo_leadingEE[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    if (h_calo_leadingEE[i]->GetMaximum() > h_calo_leadingEE[0]->GetMaximum()) h_calo_leadingEE[0]->SetMaximum(1.2*h_calo_leadingEE[i]->GetMaximum());
    h_calo_leadingEE[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_leadingEE, iPeriod, iPos );
  c_calo_leadingEE->cd()->SetLogy(1);
  c_calo_leadingEE->SaveAs((basePlotDir+"/singlePlot/calo_leadingEE.png").c_str());
  c_calo_leadingEE->SaveAs((basePlotDir+"/singlePlot/calo_leadingEE.pdf").c_str());
  
  TCanvas *c_calo_leadingECAL = new TCanvas("c_calo_leadingECAL", "c_calo_leadingECAL", 800, 800);
  h_calo_leadingECAL[0]->Draw("e0e1x0"); nobackgend->Draw();
  for (int i = 1; i < nSamples; i++){
    if (h_calo_leadingECAL[i]->GetMaximum() > h_calo_leadingECAL[0]->GetMaximum()) h_calo_leadingECAL[0]->SetMaximum(1.2*h_calo_leadingECAL[i]->GetMaximum());
    h_calo_leadingECAL[i]->Draw("hesame");
  }
  CMS_lumi( c_calo_leadingECAL, iPeriod, iPos );
  c_calo_leadingECAL->cd()->SetLogy(1);
  c_calo_leadingECAL->SaveAs((basePlotDir+"/singlePlot/calo_leadingECAL.png").c_str());
  c_calo_leadingECAL->SaveAs((basePlotDir+"/singlePlot/calo_leadingECAL.pdf").c_str());
  
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
  
  TCanvas *c_calo_energy_muon_deltaR = new TCanvas("c_calo_energy_muon_deltaR", "c_calo_energy_muon_deltaR", nSamples*800, 800); c_calo_energy_muon_deltaR->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_energy_muon_deltaR->cd(i+1);
    h_calo_energy_muon_deltaR[i]->Draw("COLZ");
    c_calo_energy_muon_deltaR->cd(i+1)->SetLogz(1);
  }
  c_calo_energy_muon_deltaR->SaveAs((basePlotDir+"/collective/calo_energy_muon_deltaR.png").c_str());
  c_calo_energy_muon_deltaR->SaveAs((basePlotDir+"/collective/calo_energy_muon_deltaR.pdf").c_str());
  
  TCanvas *c_calo_energy_pion_deltaR = new TCanvas("c_calo_energy_pion_deltaR", "c_calo_energy_pion_deltaR", nSamples*800, 800); c_calo_energy_pion_deltaR->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_energy_pion_deltaR->cd(i+1);
    h_calo_energy_pion_deltaR[i]->Draw("COLZ");
    c_calo_energy_pion_deltaR->cd(i+1)->SetLogz(1);
  }
  c_calo_energy_pion_deltaR->SaveAs((basePlotDir+"/collective/calo_energy_pion_deltaR.png").c_str());
  c_calo_energy_pion_deltaR->SaveAs((basePlotDir+"/collective/calo_energy_pion_deltaR.pdf").c_str());
  
  TCanvas *c_calo_energy_muon_deltaEta = new TCanvas("c_calo_energy_muon_deltaEta", "c_calo_energy_muon_deltaEta", nSamples*800, 800); c_calo_energy_muon_deltaEta->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_energy_muon_deltaEta->cd(i+1);
    h_calo_energy_muon_deltaEta[i]->Draw("COLZ");
    c_calo_energy_muon_deltaEta->cd(i+1)->SetLogz(1);
  }
  c_calo_energy_muon_deltaEta->SaveAs((basePlotDir+"/collective/calo_energy_muon_deltaEta.png").c_str());
  c_calo_energy_muon_deltaEta->SaveAs((basePlotDir+"/collective/calo_energy_muon_deltaEta.pdf").c_str());
  
  TCanvas *c_calo_energy_pion_deltaEta = new TCanvas("c_calo_energy_pion_deltaEta", "c_calo_energy_pion_deltaEta", nSamples*800, 800); c_calo_energy_pion_deltaEta->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++){
    c_calo_energy_pion_deltaEta->cd(i+1);
    h_calo_energy_pion_deltaEta[i]->Draw("COLZ");
    c_calo_energy_pion_deltaEta->cd(i+1)->SetLogz(1);
  }
  c_calo_energy_pion_deltaEta->SaveAs((basePlotDir+"/collective/calo_energy_pion_deltaEta.png").c_str());
  c_calo_energy_pion_deltaEta->SaveAs((basePlotDir+"/collective/calo_energy_pion_deltaEta.pdf").c_str());
  
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
  h_nCaloTowers[0]->Draw("e0e1x0"); 
  for (int i = 1; i < nSamples; i++){
    h_nCaloTowers[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_nCaloTowers, iPeriod, iPos );
  //c_nCaloTowers->cd()->SetLogy(1);
  c_nCaloTowers->SaveAs((basePlotDir+"/singlePlot/nCaloTowers.png").c_str());
  c_nCaloTowers->SaveAs((basePlotDir+"/singlePlot/nCaloTowers.pdf").c_str());
  
  TCanvas *c_nTowersMuon = new TCanvas("c_nTowersMuon", "c_nTowersMuon", 800, 800);
  h_nTowersMuon[0]->Draw("e0e1x0"); 
  for (int i = 1; i < nSamples; i++){
    h_nTowersMuon[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_nTowersMuon, iPeriod, iPos );
  //c_nTowersMuon->cd()->SetLogy(1);
  c_nTowersMuon->SaveAs((basePlotDir+"/singlePlot/nTowersMuon.png").c_str());
  c_nTowersMuon->SaveAs((basePlotDir+"/singlePlot/nTowersMuon.pdf").c_str());
  
  TCanvas *c_nTowersPion = new TCanvas("c_nTowersPion", "c_nTowersPion", 800, 800);
  h_nTowersPion[0]->Draw("e0e1x0"); 
  for (int i = 1; i < nSamples; i++){
    h_nTowersPion[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_nTowersPion, iPeriod, iPos );
  //c_nTowersPion->cd()->SetLogy(1);
  c_nTowersPion->SaveAs((basePlotDir+"/singlePlot/nTowersPion.png").c_str());
  c_nTowersPion->SaveAs((basePlotDir+"/singlePlot/nTowersPion.pdf").c_str());
  
  TCanvas *c_nTowersECALMuon = new TCanvas("c_nTowersECALMuon", "c_nTowersECALMuon", 800, 800);
  h_nTowersECALMuon[0]->Draw("e0e1x0"); 
  for (int i = 1; i < nSamples; i++){
    h_nTowersECALMuon[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_nTowersECALMuon, iPeriod, iPos );
  //c_nTowersECALMuon->cd()->SetLogy(1);
  c_nTowersECALMuon->SaveAs((basePlotDir+"/singlePlot/nTowersECALMuon.png").c_str());
  c_nTowersECALMuon->SaveAs((basePlotDir+"/singlePlot/nTowersECALMuon.pdf").c_str());
  
  TCanvas *c_nTowersECALPion = new TCanvas("c_nTowersECALPion", "c_nTowersECALPion", 800, 800);
  h_nTowersECALPion[0]->Draw("e0e1x0"); 
  for (int i = 1; i < nSamples; i++){
    h_nTowersECALPion[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_nTowersECALPion, iPeriod, iPos );
  //c_nTowersECALPion->cd()->SetLogy(1);
  c_nTowersECALPion->SaveAs((basePlotDir+"/singlePlot/nTowersECALPion.png").c_str());
  c_nTowersECALPion->SaveAs((basePlotDir+"/singlePlot/nTowersECALPion.pdf").c_str());
  
  TCanvas *c_nTowersECALFSR = new TCanvas("c_nTowersECALFSR", "c_nTowersECALFSR", 800, 800);
  h_nTowersECALFSR[0]->Draw("e0e1x0");
  //h_nTowersECALFSR[0]->GetYaxis()->SetRangeUser(0,400);
  for (int i = 1; i < nSamples; i++) h_nTowersECALFSR[i]->Draw("hesame"); 
  nobackgend->Draw();
  CMS_lumi( c_nTowersECALFSR, iPeriod, iPos );
  c_nTowersECALFSR->SaveAs((basePlotDir+"/singlePlot/nTowersECALFSR.png").c_str());
  c_nTowersECALFSR->SaveAs((basePlotDir+"/singlePlot/nTowersECALFSR.pdf").c_str());
  
  TCanvas *c_nTowersHCALMuon = new TCanvas("c_nTowersHCALMuon", "c_nTowersHCALMuon", 800, 800);
  h_nTowersHCALMuon[0]->Draw("e0e1x0"); 
  for (int i = 1; i < nSamples; i++){
    h_nTowersHCALMuon[i]->Draw("hesame");
  }
  nobackgend->Draw();
  CMS_lumi( c_nTowersHCALMuon, iPeriod, iPos );
  //c_nTowersHCALMuon->cd()->SetLogy(1);
  c_nTowersHCALMuon->SaveAs((basePlotDir+"/singlePlot/nTowersHCALMuon.png").c_str());
  c_nTowersHCALMuon->SaveAs((basePlotDir+"/singlePlot/nTowersHCALMuon.pdf").c_str());
  
  TCanvas *c_nTowersHCALPion = new TCanvas("c_nTowersHCALPion", "c_nTowersHCALPion", 800, 800);
  h_nTowersHCALPion[0]->Draw("e0e1x0");
  h_nTowersHCALPion[0]->GetYaxis()->SetRangeUser(0,1.2);
  for (int i = 1; i < nSamples; i++){
    h_nTowersHCALPion[i]->Draw("hesame");
  } 
  nobackgend->Draw();
  CMS_lumi( c_nTowersHCALPion, iPeriod, iPos );
  //c_nTowersHCALPion->cd()->SetLogy(1);
  c_nTowersHCALPion->SaveAs((basePlotDir+"/singlePlot/nTowersHCALPion.png").c_str());
  c_nTowersHCALPion->SaveAs((basePlotDir+"/singlePlot/nTowersHCALPion.pdf").c_str());
  
  TCanvas *c_nTowersHCALFSR = new TCanvas("c_nTowersHCALFSR", "c_nTowersHCALFSR", 800, 800);
  h_nTowersHCALFSR[0]->Draw("e0e1x0");
  //h_nTowersHCALFSR[0]->GetYaxis()->SetRangeUser(0,400);
  for (int i = 1; i < nSamples; i++) h_nTowersHCALFSR[i]->Draw("hesame"); 
  nobackgend->Draw();
  CMS_lumi( c_nTowersHCALFSR, iPeriod, iPos );
  c_nTowersHCALFSR->SaveAs((basePlotDir+"/singlePlot/nTowersHCALFSR.png").c_str());
  c_nTowersHCALFSR->SaveAs((basePlotDir+"/singlePlot/nTowersHCALFSR.pdf").c_str());
  
  TCanvas *c_ratio_ECAL_HCAL_FSR = new TCanvas("c_ratio_ECAL_HCAL_FSR", "c_ratio_ECAL_HCAL_FSR", 800, 800);
  h_ratio_ECAL_HCAL_FSR[0]->Draw("e0e1x0");
  h_ratio_ECAL_HCAL_FSR[0]->GetYaxis()->SetRangeUser(0,1.2);
  for (int i = 1; i < nSamples; i++) h_ratio_ECAL_HCAL_FSR[i]->Draw("hesame"); 
  nobackgend->Draw();
  CMS_lumi( c_ratio_ECAL_HCAL_FSR, iPeriod, iPos );
  c_ratio_ECAL_HCAL_FSR->SaveAs((basePlotDir+"/singlePlot/ratio_ECAL_HCAL_FSR.png").c_str());
  c_ratio_ECAL_HCAL_FSR->SaveAs((basePlotDir+"/singlePlot/ratio_ECAL_HCAL_FSR.pdf").c_str());

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

#endif


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
  
  //background_sysNchDown->SetDirectory(gDirectory); background_sysNchDown->Write();
  //background_sysNchUp->SetDirectory(gDirectory); background_sysNchUp->Write();
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