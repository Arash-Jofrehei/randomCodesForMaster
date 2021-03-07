#include <vector>
#include <iostream>
#include <string>
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
//#include "/Applications/root/macros/AtlasStyle.C"

//#define nch_cut

// **************************
// cut definition
  bool tau_muon_isSoft = true;
  double tau_hadron_vertexprob = 0.1;
  float HFpLeading_low = 0;
  float HFmLeading_low = 0;
  float HFpLeading_high = 200;
  float HFmLeading_high = 200;
  int max_nch = 40;
  float MET_cut = 200;
// end of cut definition
// **************************

double my_pi=TMath::Pi();
double deltaphi(double phi1, double phi2)
{
 double dphi=fabs(phi1-phi2);
 return (dphi<=my_pi)? dphi : 2.*my_pi-dphi;
}

//void run(string fileData, string fileMCS){//, string fileMCB1, string fileMCB2){
  //string files[2] = {fileData,fileMCS};
  //const int nSamples = 2;
  
//string inputFiles[] = {"../ntuples/flatTuple_data_2015.root", "/eos/user/a/ajofrehe/temporary_barn/flatTuple_mcggTauTau_gen_2015.root"};

//string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root", "../ntuples/flatTuple_mcggTauTau_fixedAOD_2015.root"};

string inputFiles[] = {"../ntuples/flatTuple_fixedAOD_data_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2015.root", "../ntuples/flatTuple_mcggTauTau_AOD_MadGraph_2018.root", "../ntuples/flatTuple_mcggTauTau_AOD_SuperChic_2018.root", "/eos/user/a/ajofrehe/gtau/ntuples/flatTuple_mcggBBbar_4f_AOD_MadGraph_2015.root"};
  
//string inputFiles[] = {"../ntuples/flatTuple_AOD_data_2015.root", "../ntuples/flatTuple_fixedAOD_data_2015.root","../ntuples/flatTuple_data_2015.root", "../ntuples/flatTuple_mcggTauTau_2015.root","../ntuples/flatTuple_mcggTauTau_fixedAOD_2015.root","../ntuples/flatTuple_mcggTauTau_2018.root"};

//string inputFiles[] = {"../ntuples/flatTuple_data_2015.root","../ntuples/flatTuple_mcggTauTau_2015.root","../ntuples/flatTuple_mcggBBbar_2015.root","../ntuples/flatTuple_mcggCCbar_2015.root","../ntuples/flatTuple_mcggTauTau_HydjetDrumMB_2018.root"};

//string inputFiles[] = {"../ntuples/flatTuple_data_2015.root","../ntuples/flatTuple_mcggTauTau_2018.root","../ntuples/flatTuple_mcggBBbar_2018.root","../ntuples/flatTuple_mcggCCbar_2018.root","../ntuples/flatTuple_mcggTauTau_HydjetDrumMB_2018.root"};

void run(const int nSamples, string files[] = inputFiles){
//void run(vector<string> files){
  //const int nSamples = files.size();
  TFile *root_file[nSamples];
  TREE_DATA *TREE;
  TREE_MC *TREEMCs[nSamples-1];
  TH1F *cutflow[nSamples];
  //SetAtlasStyle();

  Double_t Red[2]   = { 1.00, .00};
  Double_t Green[2] = { 1.00, .00};
  Double_t Blue[2]  = { 1.00, .00};
  Double_t Stops[2] = { .00, 1.00};
                                                                               
  Int_t nb=5;
  TColor::CreateGradientColorTable(2,Stops,Red,Green,Blue,nb);

  root_file[0] = new TFile(files[0].c_str(),"READ");
  TREE = new TREE_DATA((TTree*)(root_file[0])->Get("tree"),root_file[0]);
  cutflow[0] = (TH1F*) root_file[0]->Get("ntuplizer/cutflow_perevt");
  for (int i = 0; i < nSamples-1; i++){
    root_file[i+1] = new TFile(files[i+1].c_str(),"READ");
    TREEMCs[i] = new TREE_MC((TTree*)(root_file[i+1])->Get("tree"),root_file[i+1]);
    cutflow[i+1] = (TH1F*) root_file[i+1]->Get("ntuplizer/cutflow_perevt");
    cutflow[i+1]->SetLineColor(i+1); cutflow[i+1]->SetMarkerSize(0);
  }
  if (nSamples>5) cutflow[5]->SetLineColor(6);
  
// *******************************************
// histogram declaration
  
  float dataLumi = 0.55;
  //string tag[5] = {"Data 2015","MC-Signal 2015","MC-BBbar 2015","MC-CCbar 2015","MC-MB 2015"};
  string tag[5] = {"Data 2015","MC Signal 2015 MadGraph","MC Signal 2018 MadGraph","MC Signal 2018 SuperChic", "MC BBbar 2015 MadGraph"};
  //string tag[6] = {"Data AOD 2015","Data modified AOD 2015","Data miniAOD 2015","MC-Signal AOD 2015","MC-Signal modified AOD 2015","MC-Signal miniAOD 2018"};
  float nEvents[5] = {1,1184600,2000000,623000,30000000};
  //for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1);
  for (int s = 0; s < nSamples; s++) nEvents[s] = cutflow[s]->GetBinContent(1)/200;
  //float crossSectionMC[4] = {570000,1500,300000,570000};
  float crossSectionMC[4] = {570000,570000,570000,1500};
  float SF[5] = {1,0,0,0,0};
  for (int s = 1; s < nSamples; s++) if(nEvents[s] != 0) SF[s] = dataLumi * crossSectionMC[s-1] / nEvents[s];
  float SF2[5] = {1,0,0,0,0};
  for (int s = 1; s < nSamples; s++) if(cutflow[s]->GetBinContent(1) != 0) SF2[s] = dataLumi * crossSectionMC[s-1] / cutflow[s]->GetBinContent(1);
  //SF[1] = nEvents[0]/nEvents[1];
  //SF[2] = nEvents[0]/nEvents[2];

  int tau_mu_pt_bins = 40;
  int tau_mu_eta_bins = 20;
  int tau_mu_phi_bins = 20;
  int tau_hadron_pt_bins = 40;
  int tau_hadron_eta_bins = 20;
  int tau_hadron_phi_bins = 20;
  int tau_hadron_rhomass_bins = 12;
  int tau_hadron_nch_bins = 20;
  int tau_hadron_pvz_bins = 50;
  int calo_energy_bins = 125;
  int deltaphi_bins = 8;
  
  THStack *hs_tau_mu_pt = new THStack("hs_tau_mu_pt","#tau_{#mu} p_{T} [GeV]");
  TH1F *h_tau_mu_pt[nSamples];
  THStack *hs_tau_mu_eta = new THStack("hs_tau_mu_eta","#tau_{#mu} #eta");
  TH1F *h_tau_mu_eta[nSamples];
  THStack *hs_tau_mu_phi = new THStack("hs_tau_mu_phi","#tau_{#mu} #phi");
  TH1F *h_tau_mu_phi[nSamples];
  THStack *hs_tau_hadron_pt = new THStack("hs_tau_hadron_pt","#tau_{hadron} p_{T} [GeV]");
  TH1F *h_tau_hadron_pt[nSamples];
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
  THStack *hs_tau_hadron_ncand_final = new THStack("hs_tau_hadron_ncand_final","#tau_{hadron} cands final");
  TH1F *h_tau_hadron_ncand_final[nSamples];
  THStack *hs_tau_hadron_vprob = new THStack("hs_tau_hadron_vprob","#tau_{hadron} vprob (%)");
  TH1F *h_tau_hadron_vprob[nSamples];
  TH1F *h_tau_hadron_matched_pt_index[nSamples];
  THStack *hs_tau_hadron_mass = new THStack("hs_tau_hadron_mass","#tau_{hadron} mass [GeV]");
  TH1F *h_tau_hadron_mass[nSamples];
  THStack *hs_ditau_mass = new THStack("hs_tau_hadron_mass","#tau#tau mass [GeV]");
  TH1F *h_ditau_mass[nSamples];
  THStack *hs_tau_hadron_track_pvz[3];
  for (int j=0; j < 3;j++) hs_tau_hadron_track_pvz[j] = new THStack(("h_tau_hadron_track" + to_string(j) + "_pvz").c_str(),("#tau_{hadron} track" + to_string(j) + " pvz").c_str());
  TH1F *h_tau_hadron_track_pvz[nSamples][3];
  THStack *hs_deltaphi_tau_mu_tau_hadron = new THStack("hs_deltaphi_tau_mu_tau_hadron","#Delta#phi(#tau_{#mu}, #tau_{hadron})");
  TH1F *h_deltaphi_tau_mu_tau_hadron[nSamples];
  THStack *hs_PV_N = new THStack("hs_PV_N","number of PV");
  TH1F *h_PV_N[nSamples];
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
  THStack *hs_calo_HF_eta = new THStack("hs_calo_HF_eta","HF #eta");
  TH1F *h_calo_HF_eta[nSamples];
  TH2F *h_calo_HF_energy_eta[nSamples];
  TH2F *h_deltaphi_tau_mu_tau_hadron_mueta[nSamples];
  TH2F *h_deltaphi_tau_mu_tau_hadron_deltaeta[nSamples];
  TH2F *h_mueta_taueta[nSamples];
  TH2F *h_AP[nSamples];
  
  // MET
  TH1F *h_MET[nSamples];
  
  for (int i = 0; i < nSamples; i++){
    h_tau_mu_pt[i] = new TH1F(("h_tau_mu_pt_" + tag[i]).c_str(),("h_tau_mu_pt_" + tag[i]).c_str(),tau_mu_pt_bins, 0, 20);
    h_tau_mu_pt[i]->SetXTitle("#tau_{#mu} p_{T} [GeV]"); h_tau_mu_pt[i]->Sumw2();
    if (i != 0) {h_tau_mu_pt[i]->SetLineColor(i); h_tau_mu_pt[i]->SetMarkerSize(0);}
        
    h_tau_mu_eta[i] = new TH1F(("h_tau_mu_eta_" + tag[i]).c_str(),("h_tau_mu_eta_" + tag[i]).c_str(),tau_mu_eta_bins, -2.5, 2.5);
    h_tau_mu_eta[i]->SetXTitle("#tau_{#mu} #eta"); h_tau_mu_eta[i]->Sumw2();
    if (i != 0) {h_tau_mu_eta[i]->SetLineColor(i); h_tau_mu_eta[i]->SetMarkerSize(0);}
        
    h_tau_mu_phi[i] = new TH1F(("h_tau_mu_phi_" + tag[i]).c_str(),("h_tau_mu_phi_" + tag[i]).c_str(),tau_mu_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_mu_phi[i]->SetXTitle("#tau_{#mu} #phi"); h_tau_mu_phi[i]->Sumw2();
    if (i != 0) {h_tau_mu_phi[i]->SetLineColor(i); h_tau_mu_phi[i]->SetMarkerSize(0);}
    
    // hadron related
    
    h_tau_hadron_pt[i] = new TH1F(("h_tau_hadron_pt_" + tag[i]).c_str(),("h_tau_hadron_pt_" + tag[i]).c_str(),tau_hadron_pt_bins, 0, 20);
    h_tau_hadron_pt[i]->SetXTitle("#tau_{hadron} p_{T} [GeV]"); h_tau_hadron_pt[i]->Sumw2();
    if (i != 0) {h_tau_hadron_pt[i]->SetLineColor(i); h_tau_hadron_pt[i]->SetMarkerSize(0);}
        
    h_tau_hadron_eta[i] = new TH1F(("h_tau_hadron_eta_" + tag[i]).c_str(),("h_tau_hadron_eta_" + tag[i]).c_str(),tau_hadron_eta_bins, -2.5, 2.5);
    h_tau_hadron_eta[i]->SetXTitle("#tau_{hadron} #eta"); h_tau_hadron_eta[i]->Sumw2();
    if (i != 0) {h_tau_hadron_eta[i]->SetLineColor(i); h_tau_hadron_eta[i]->SetMarkerSize(0);}
        
    h_tau_hadron_phi[i] = new TH1F(("h_tau_hadron_phi_" + tag[i]).c_str(),("h_tau_hadron_phi_" + tag[i]).c_str(),tau_hadron_phi_bins, -TMath::Pi(), TMath::Pi());
    h_tau_hadron_phi[i]->SetXTitle("#tau_{hadron} #phi"); h_tau_hadron_phi[i]->Sumw2();
    if (i != 0) {h_tau_hadron_phi[i]->SetLineColor(i); h_tau_hadron_phi[i]->SetMarkerSize(0);}
        
    for (int j=0; j < 2;j++){
      h_tau_hadron_rhomass[i][j] = new TH1F(("h_tau_hadron_rhomass" + to_string(j) + "_" + tag[i]).c_str(),("#tau_{hadron} #rho_{"+to_string(j+1)+"} mass [GeV] " + tag[i]).c_str(),tau_hadron_rhomass_bins, 0.2, 1.4);
      h_tau_hadron_rhomass[i][j]->SetXTitle(("#rho_{"+to_string(j+1)+"} mass [GeV]").c_str()); h_tau_hadron_rhomass[i][j]->Sumw2();
      if (i != 0) {h_tau_hadron_rhomass[i][j]->SetLineColor(i); h_tau_hadron_rhomass[i][j]->SetMarkerSize(0);}
    }
    
    h_tau_hadron_rhomass2D[i] = new TH2F(("h_tau_hadron_rhomass2D_" + tag[i]).c_str(),("#tau_{hadron} #rho mass [GeV] " + tag[i]).c_str(),2*tau_hadron_rhomass_bins, 0.2, 1.4, 2*tau_hadron_rhomass_bins, 0.2, 1.4);
    h_tau_hadron_rhomass2D[i]->SetXTitle("#rho_{1} mass [GeV]"); h_tau_hadron_rhomass2D[i]->SetYTitle("#rho_{2} mass [GeV]");
    
    h_tau_hadron_nch[i] = new TH1F(("h_tau_hadron_nch_" + tag[i]).c_str(),("h_tau_hadron_nch_" + tag[i]).c_str(),tau_hadron_nch_bins, 0.5, 20.5);
    h_tau_hadron_nch[i]->SetXTitle("nch"); h_tau_hadron_nch[i]->Sumw2();
    if (i != 0) {h_tau_hadron_nch[i]->SetLineColor(i); h_tau_hadron_nch[i]->SetMarkerSize(0);}
    
    h_tau_hadron_ncand_final[i] = new TH1F(("h_tau_hadron_ncand_final_" + tag[i]).c_str(),("h_tau_hadron_ncand_final_" + tag[i]).c_str(),20, 0.5, 20.5);
    h_tau_hadron_ncand_final[i]->SetXTitle("ncand final"); h_tau_hadron_ncand_final[i]->Sumw2();
    if (i != 0) {h_tau_hadron_ncand_final[i]->SetLineColor(i); h_tau_hadron_ncand_final[i]->SetMarkerSize(0);}
    
    h_tau_hadron_vprob[i] = new TH1F(("h_tau_hadron_vprob_" + tag[i]).c_str(),("h_tau_hadron_vprob_" + tag[i]).c_str(),20, 0, 100);
    h_tau_hadron_vprob[i]->SetXTitle("vprob (%)"); h_tau_hadron_vprob[i]->Sumw2();
    if (i != 0) {h_tau_hadron_vprob[i]->SetLineColor(i); h_tau_hadron_vprob[i]->SetMarkerSize(0);}
    
    h_tau_hadron_matched_pt_index[i] = new TH1F(("h_tau_hadron_matched_pt_index_" + tag[i]).c_str(),("h_tau_hadron_matched_pt_index_" + tag[i]).c_str(),10, 0.5, 10.5);
    h_tau_hadron_matched_pt_index[i]->SetXTitle("p_T rank of the matched tau"); h_tau_hadron_matched_pt_index[i]->Sumw2();
    if (i != 0) {h_tau_hadron_matched_pt_index[i]->SetLineColor(i); h_tau_hadron_matched_pt_index[i]->SetMarkerSize(0);}
    
    h_tau_hadron_mass[i] = new TH1F(("h_tau_hadron_mass_" + tag[i]).c_str(),("#tau_{hadron} mass [GeV] " + tag[i]).c_str(), 13, 0.4, 1.7);
    h_tau_hadron_mass[i]->SetXTitle("#tau_{hadron} mass [GeV]"); h_tau_hadron_mass[i]->Sumw2();
    if (i != 0) {h_tau_hadron_mass[i]->SetLineColor(i); h_tau_hadron_mass[i]->SetMarkerSize(0);}
    
    h_ditau_mass[i] = new TH1F(("h_ditau_mass_" + tag[i]).c_str(),("#tau#tau invariant mass [GeV] " + tag[i]).c_str(), 20, 0, 40);
    h_ditau_mass[i]->SetXTitle("#tau#tau invariant mass [GeV]"); h_ditau_mass[i]->Sumw2();
    if (i != 0) {h_ditau_mass[i]->SetLineColor(i); h_ditau_mass[i]->SetMarkerSize(0);}
        
    for (int j=0; j < 3;j++){
      h_tau_hadron_track_pvz[i][j] = new TH1F(("h_tau_hadron_track" + to_string(j) + "_pvz_" + tag[i]).c_str(),("h_tau_hadron_track" + to_string(j) + "_pvz_" + tag[i]).c_str(),tau_hadron_pvz_bins, -0.5, 0.5);
      h_tau_hadron_track_pvz[i][j]->SetXTitle(("PV_{z} - trk"+to_string(j)+"_{z}").c_str()); h_tau_hadron_track_pvz[i][j]->Sumw2();
      if (i != 0) {h_tau_hadron_track_pvz[i][j]->SetLineColor(i); h_tau_hadron_track_pvz[i][j]->SetMarkerSize(0);}
    }
    
    // other
        
    h_deltaphi_tau_mu_tau_hadron[i] = new TH1F(("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),("h_deltaphi_tau_mu_tau_hadron_" + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi());
    h_deltaphi_tau_mu_tau_hadron[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron[i]->Sumw2();
    if (i != 0) {h_deltaphi_tau_mu_tau_hadron[i]->SetLineColor(i); h_deltaphi_tau_mu_tau_hadron[i]->SetMarkerSize(0); /*h_deltaphi_tau_mu_tau_hadron[i]->SetFillColor(i);*/}
    if (i == 5) h_deltaphi_tau_mu_tau_hadron[i]->SetLineColor(6);
    
    h_deltaphi_tau_mu_tau_hadron_mueta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_mueta_" + tag[i]).c_str(),("#eta_{#mu} vs #Delta#phi(#tau_{#mu}, #tau_{hadron}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron_mueta[i]->SetYTitle("#eta_{#mu}");
    
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i] = new TH2F(("h_deltaphi_tau_mu_tau_hadron_deltaeta_" + tag[i]).c_str(),("#Delta(abs(#eta))_{#mu_#tau} vs #Delta#phi(#tau_{#mu}, #tau_{hadron}) " + tag[i]).c_str(),deltaphi_bins, 0, TMath::Pi(), 18,-3,3);
    h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetXTitle("#Delta#phi(#tau_{#mu}, #tau_{hadron})"); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->SetYTitle("#Delta(abs(#eta))_{#mu_#tau}");
    
    h_mueta_taueta[i] = new TH2F(("h_mueta_taueta_" + tag[i]).c_str(),("#tau_{hadron} #eta vs #tau_{#mu} #eta " + tag[i]).c_str(),18,-3,3, 18,-3,3);
    h_mueta_taueta[i]->SetXTitle("#tau_{hadron} #eta"); h_mueta_taueta[i]->SetYTitle("#tau_{#mu} #eta");
    
    h_PV_N[i] = new TH1F(("h_PV_N_" + tag[i]).c_str(),("h_PV_N_" + tag[i]).c_str(),5, -.5, 4.5);
    h_PV_N[i]->SetXTitle("N_{PV}"); h_PV_N[i]->Sumw2();
    if (i != 0) {h_PV_N[i]->SetLineColor(i); h_PV_N[i]->SetMarkerSize(0);}
    
    h_AP[i] = new TH2F(("h_AP_" + tag[i]).c_str(),("AP " + tag[i]).c_str(),40, -1, 1, 45,0,15);
    h_AP[i]->SetXTitle("p_{z}^{+} - p_{z}^{-} / p_{z}^{+} + p_{z}^{-}"); h_AP[i]->SetYTitle("#tau_{#mu} p_{T} [GeV]");
    
    // calo
    
    h_calo_energyHFp[i] = new TH1F(("h_calo_energyHFp_" + tag[i]).c_str(),("energy HF+ " + tag[i]).c_str(),40, -0.5, 19.5);
    h_calo_energyHFp[i]->SetXTitle("energy HF+ [GeV]"); h_calo_energyHFp[i]->Sumw2();
    if (i != 0) {h_calo_energyHFp[i]->SetLineColor(i); h_calo_energyHFp[i]->SetMarkerSize(0);}
    
    h_calo_energyHFm[i] = new TH1F(("h_calo_energyHFm_" + tag[i]).c_str(),("energy HF- " + tag[i]).c_str(),40, -0.5, 19.5);
    h_calo_energyHFm[i]->SetXTitle("energy HF- [GeV]"); h_calo_energyHFm[i]->Sumw2();
    if (i != 0) {h_calo_energyHFm[i]->SetLineColor(i); h_calo_energyHFm[i]->SetMarkerSize(0);}
    
    h_calo_leadingHFp[i] = new TH1F(("h_calo_leadingHFp_" + tag[i]).c_str(),("energy leading tower HF+ " + tag[i]).c_str(),80, -0.5, 39.5);
    h_calo_leadingHFp[i]->SetXTitle("energy HF+ [GeV]"); h_calo_leadingHFp[i]->Sumw2();
    if (i != 0) {h_calo_leadingHFp[i]->SetLineColor(i); h_calo_leadingHFp[i]->SetMarkerSize(0);}
    
    h_calo_leadingHFm[i] = new TH1F(("h_calo_leadingHFm_" + tag[i]).c_str(),("energy leading tower HF- " + tag[i]).c_str(),80, -0.5, 39.5);
    h_calo_leadingHFm[i]->SetXTitle("energy HF- [GeV]"); h_calo_leadingHFm[i]->Sumw2();
    if (i != 0) {h_calo_leadingHFm[i]->SetLineColor(i); h_calo_leadingHFm[i]->SetMarkerSize(0);}
    
    h_calo_energyHFp_sum[i] = new TH1F(("h_calo_energyHFp_sum_" + tag[i]).c_str(),("energy sum HF+ " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFp_sum[i]->SetXTitle("energy HF+ [GeV]"); h_calo_energyHFp_sum[i]->Sumw2();
    if (i != 0) {h_calo_energyHFp_sum[i]->SetLineColor(i); h_calo_energyHFp_sum[i]->SetMarkerSize(0);}
    
    h_calo_energyHFm_sum[i] = new TH1F(("h_calo_energyHFm_sum_" + tag[i]).c_str(),("energy sum HF- " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFm_sum[i]->SetXTitle("energy HF- [GeV]"); h_calo_energyHFm_sum[i]->Sumw2();
    if (i != 0) {h_calo_energyHFm_sum[i]->SetLineColor(i); h_calo_energyHFm_sum[i]->SetMarkerSize(0);}
    
    h_calo_energyHFp_size[i] = new TH1F(("h_calo_energyHFp_size_" + tag[i]).c_str(),("size HF+ " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFp_size[i]->SetXTitle("size HF+"); h_calo_energyHFp_size[i]->Sumw2();
    if (i != 0) {h_calo_energyHFp_size[i]->SetLineColor(i); h_calo_energyHFp_size[i]->SetMarkerSize(0);}
    
    h_calo_energyHFm_size[i] = new TH1F(("h_calo_energyHFm_size_" + tag[i]).c_str(),("size HF- " + tag[i]).c_str(),calo_energy_bins, -0.5, 249.5);
    h_calo_energyHFm_size[i]->SetXTitle("size HF-"); h_calo_energyHFm_size[i]->Sumw2();
    if (i != 0) {h_calo_energyHFm_size[i]->SetLineColor(i); h_calo_energyHFm_size[i]->SetMarkerSize(0);}
    
    h_calo_energyHF_pm[i] = new TH2F(("h_calo_energyHF_pm_" + tag[i]).c_str(),("total energy HF+ vs HF- " + tag[i]).c_str(),50, -0.5, 249.5, 50, -0.5, 249.5);
    h_calo_energyHF_pm[i]->SetXTitle("total energy HF-"); h_calo_energyHF_pm[i]->SetYTitle("total energy HF+");
    
    h_calo_leadingHF_pm[i] = new TH2F(("h_calo_leadingHF_pm_" + tag[i]).c_str(),("leading tower HF+ vs HF- " + tag[i]).c_str(),20, -0.5, 9.5, 20, -0.5, 9.5);
    h_calo_leadingHF_pm[i]->SetXTitle("leading tower HF-"); h_calo_leadingHF_pm[i]->SetYTitle("leading tower HF+");
    
    h_calo_HF_eta[i] = new TH1F(("h_calo_HF_eta_" + tag[i]).c_str(),("HF eta " + tag[i]).c_str(), 20, -5, 5);
    h_calo_HF_eta[i]->SetXTitle("#eta"); h_calo_HF_eta[i]->Sumw2();
    if (i != 0) {h_calo_HF_eta[i]->SetLineColor(i); h_calo_HF_eta[i]->SetMarkerSize(0);}
    
    h_calo_HF_energy_eta[i] = new TH2F(("h_calo_HF_energy_eta_" + tag[i]).c_str(),("energy HF vs #eta " + tag[i]).c_str(),20, -5, 5, 25, -0.5, 24.5);
    h_calo_HF_energy_eta[i]->SetXTitle("#eta"); h_calo_HF_energy_eta[i]->SetYTitle("energy HF");
    
    // MET
    
    h_MET[i] = new TH1F(("h_MET_" + tag[i]).c_str(),("MET " + tag[i]).c_str(),15, -0.5, 14.5);
    h_MET[i]->SetXTitle("MET"); h_MET[i]->Sumw2();
    if (i != 0) {h_MET[i]->SetLineColor(i); h_MET[i]->SetMarkerSize(0);}
    
  }

  TH2F *h_calo_energyHFp_nch = new TH2F("h_calo_energyHFp_nch","energy HF+ vs nch",4 , 2.5, 6.5, 30, 100, 160);
  h_calo_energyHFp_nch->SetXTitle("nch"); h_calo_energyHFp_nch->SetYTitle("energy HF+ [GeV]");
  
  TH2F *h_calo_energyHFm_nch = new TH2F("h_calo_energyHFm_nch","energy HF- vs nch",4 , 2.5, 6.5, 30, 100, 160);
  h_calo_energyHFm_nch->SetXTitle("nch"); h_calo_energyHFm_nch->SetYTitle("energy HF- [GeV]");

  // other
  TH1F *h_track_activity_pt     = new TH1F("h_track_activity_pt",     "h_track_activity_pt",     25, 0, 10); h_track_activity_pt->Sumw2(); h_track_activity_pt->SetXTitle("track p_{T} [GeV]");
  TH2F *h_track_activity_pt_eta = new TH2F("h_track_activity_pt_eta", "h_track_activity_pt_eta", 25, 0, 10, 10, -2.5, 2.5);
  h_track_activity_pt_eta->SetXTitle("track p_{T} [GeV]"); h_track_activity_pt_eta->SetYTitle("track #eta");
  
  TH2F *h_deltaphi_tau_mu_tau_hadron_nch = new TH2F("h_deltaphi_tau_mu_tau_hadron_nch", "nch vs #Delta#phi(#tau_{#mu}, #tau_{hadron})", 8, 0, TMath::Pi(), 4 , 2.5, 6.5);
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
  
  TLorentzVector tau_muon, tau_hadron;
  int charge_counter[nSamples][3];
  for (int s = 0; s < nSamples; s++){
    for (int i = 0; i < 3; i++) charge_counter[s][i] = 0;
  }
  int muon_charge = 0;
  int tauh_charge = 0;
  int candCounter = 0;
  for(int iEntry=0; iEntry<entries; iEntry++) {
    (TREE->fChain)->GetEntry(iEntry);

    TLorentzVector tau_hadron_pi1, tau_hadron_pi2, tau_hadron_pi3;

    bool passedmu  = false;
    bool passedtau = false;
    bool passedcalo = true;
    bool passedMET = true;
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

    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    candCounter = 0;
    double temp_tau_rho1 = 0.; double temp_tau_rho2 = 0.;
    for (int i=0; i<(int)TREE->BsTauTau_tau_pt->size(); i++) {
      h_tau_hadron_vprob[0]->Fill(100*TREE->BsTauTau_tau_vprob->at(i));
      if (TREE->BsTauTau_tau_vprob->at(i) < tau_hadron_vertexprob) { continue; }
      if (TREE->BsTauTau_tau_pt->at(i) < 3.0) { continue; }
      if (muon_charge*TREE->BsTauTau_tau_q->at(i) == 1) candCounter += 1;
      passedtau = true;
      if (TREE->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison && muon_charge*TREE->BsTauTau_tau_q->at(i) != 1) {
        temp_tau_pt_comparison = TREE->BsTauTau_tau_pt->at(i);
        tauh_charge = TREE->BsTauTau_tau_q->at(i);
        //tau_hadron.SetPtEtaPhiM (TREE->BsTauTau_tau_pt->at(i), TREE->BsTauTau_tau_eta->at(i), TREE->BsTauTau_tau_phi->at(i), muon_mass);
        tau_hadron.SetPtEtaPhiM (TREE->BsTauTau_tau_pt->at(i), TREE->BsTauTau_tau_eta->at(i), TREE->BsTauTau_tau_phi->at(i), TREE->BsTauTau_tau_mass->at(i));
        tau_hadron_pi1.SetPtEtaPhiM (TREE->BsTauTau_tau_pi1_pt->at(i), TREE->BsTauTau_tau_pi1_eta->at(i), TREE->BsTauTau_tau_pi1_phi->at(i), .13957018);
        tau_hadron_pi2.SetPtEtaPhiM (TREE->BsTauTau_tau_pi2_pt->at(i), TREE->BsTauTau_tau_pi2_eta->at(i), TREE->BsTauTau_tau_pi2_phi->at(i), .13957018);
        tau_hadron_pi3.SetPtEtaPhiM (TREE->BsTauTau_tau_pi3_pt->at(i), TREE->BsTauTau_tau_pi3_eta->at(i), TREE->BsTauTau_tau_pi3_phi->at(i), .13957018);
        temp_tau_rho1 = TREE->BsTauTau_tau_rhomass1->at(i); temp_tau_rho2 = TREE->BsTauTau_tau_rhomass2->at(i);
        temp_pv_trk_1 = TREE->BsTauTau_tau_pi1_z->at(i); temp_pv_trk_2 = TREE->BsTauTau_tau_pi2_z->at(i); temp_pv_trk_3 = TREE->BsTauTau_tau_pi3_z->at(i);
      }
    } // loop over the size of the tau candidates
    
    h_tau_hadron_ncand_final[0]->Fill(candCounter);
    
    //if ((tau_muon+tau_hadron).Pt() > MET_cut) passedMET = false;

    bool found_mu_track = false;
    bool found_pi1_track = false;
    bool found_pi2_track = false;
    bool found_pi3_track = false;
    bool high_activity = false;
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
    } // check activity

    temp_nch = TREE->BsTauTau_nch->at(0);
    temp_pvz = TREE->BsTauTau_bbPV_vz->at(0);
    temp_npv = TREE->PV_N;
    double sumHFp = 0;
    double maxHFp = 0;
    int sizeHFp = 0;
    double sumHFm = 0;
    double maxHFm = 0;
    int sizeHFm = 0;
    
#ifdef nch_cut
    if (passedmu && passedtau && temp_nch<=max_nch)
#else
    if (passedmu && passedtau)
#endif
    {
      for (int i=0; i<(int)TREE->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREE->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREE->BsTauTau_calo_energyHFm->at(i);
        double etaCal = TREE->BsTauTau_calo_eta->at(i);
        if (eHFp != -1)
          {h_calo_energyHFp[0]->Fill(eHFp); sizeHFp++; sumHFp += eHFp; h_calo_HF_eta[0]->Fill(etaCal); h_calo_HF_energy_eta[0]->Fill(etaCal,eHFp);}
        if (eHFm != -1)
          {h_calo_energyHFm[0]->Fill(eHFm); sizeHFm++; sumHFm += eHFm; h_calo_HF_eta[0]->Fill(etaCal); h_calo_HF_energy_eta[0]->Fill(etaCal,eHFm);}
        if (eHFp > maxHFp) maxHFp = eHFp;
        if (eHFm > maxHFm) maxHFm = eHFm;
      }
      if (maxHFp > HFpLeading_high || maxHFm > HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcalo = false;
      if (passedcalo) h_calo_energyHFp_sum[0]->Fill(sumHFp);
      if (passedcalo) h_calo_energyHFm_sum[0]->Fill(sumHFm);
      if (passedcalo) h_calo_energyHF_pm[0]->Fill(sumHFm,sumHFp);
      if (passedcalo) h_calo_energyHFp_nch->Fill(TREE->BsTauTau_nch->at(0),sumHFp); h_calo_energyHFm_nch->Fill(TREE->BsTauTau_nch->at(0),sumHFm);
      if (passedcalo) h_calo_energyHFp_size[0]->Fill(sizeHFp);
      if (passedcalo) h_calo_energyHFm_size[0]->Fill(sizeHFm);
      h_calo_leadingHFp[0]->Fill(maxHFp);
      h_calo_leadingHFm[0]->Fill(maxHFm);
      if (passedcalo) h_calo_leadingHF_pm[0]->Fill(maxHFm,maxHFp);
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    }

    
#ifdef nch_cut
    if (passedmu && passedtau && temp_nch<=max_nch && passedcalo /*&& !high_activity*/ && passedMET)
#else
    if (passedmu && passedtau && passedcalo /*&& !high_activity*/ && passedMET)
#endif
    {
      charge_counter[0][muon_charge*tauh_charge + 1] += 1;
      if (muon_charge*tauh_charge == 0) continue;
      h_tau_mu_pt[0] ->Fill(tau_muon.Pt());
      h_tau_mu_eta[0]->Fill(tau_muon.Eta());
      h_tau_mu_phi[0]->Fill(tau_muon.Phi());
      h_tau_hadron_pt[0]->Fill(tau_hadron.Pt());
      h_tau_hadron_eta[0]->Fill(tau_hadron.Eta());
      h_tau_hadron_phi[0]->Fill(tau_hadron.Phi());
      h_tau_hadron_rhomass[0][0]->Fill(temp_tau_rho1);
      h_tau_hadron_rhomass[0][1]->Fill(temp_tau_rho2);
      h_tau_hadron_rhomass2D[0]->Fill(temp_tau_rho1,temp_tau_rho2);
      h_tau_hadron_nch[0]->Fill(temp_nch);
      h_tau_hadron_mass[0]->Fill(tau_hadron.M());
      float delta_phi = TMath::Abs(tau_muon.DeltaPhi(tau_hadron));
      h_deltaphi_tau_mu_tau_hadron_nch->Fill(delta_phi,temp_nch);
      h_deltaphi_tau_mu_tau_hadron_mueta[0]->Fill(delta_phi,tau_muon.Eta());
      h_deltaphi_tau_mu_tau_hadron_deltaeta[0]->Fill(delta_phi,TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()));
      h_mueta_taueta[0]->Fill(tau_muon.Eta(),tau_hadron.Eta());
      h_tau_hadron_track_pvz[0][0]->Fill(temp_pvz - temp_pv_trk_1);
      h_tau_hadron_track_pvz[0][1]->Fill(temp_pvz - temp_pv_trk_2);
      h_tau_hadron_track_pvz[0][2]->Fill(temp_pvz - temp_pv_trk_3);

      h_deltaphi_tau_mu_tau_hadron[0]->Fill(TMath::Abs(tau_muon.DeltaPhi(tau_hadron)));
      h_PV_N[0]->Fill(temp_npv);
      TLorentzVector MET;
      MET.SetPtEtaPhiM((tau_muon+tau_hadron).Pt(),(tau_muon+tau_hadron).Eta(),-(tau_muon+tau_hadron).Phi(),0);
      h_MET[0]->Fill(MET.Pt());
      h_ditau_mass[0]->Fill((tau_muon+tau_hadron+MET).M());
      h_AP[0]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2);
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
    }
    

//      aux->Fill();

  } // loop over the entries
  cout << "Muon charge times hadronic tau charge for Data:\n -1: " << charge_counter[0][0] << "\n 0: " << charge_counter[0][1] << "\n +1: " << charge_counter[0][2] << endl;

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
  double temp_tau_pt_comparison = 0.;
  double temp_nch = 0.;
  double temp_npv = 0.;
  double temp_pvz = 0.;
  int matchedpT_rank = 0;
  
  for (int s = 0; s < nSamples-1; s++){
  entriesMC = (TREEMCs[s]->fChain)->GetEntries();
  for(int iEntry=0; iEntry<entriesMC; iEntry++) {
    (TREEMCs[s]->fChain)->GetEntry(iEntry);
    
    //cout << "number of taus: " << (int)TREEMCs[s]->gen_tau_pt->size() << endl;
    if (TREEMCs[s]->gen_tautau_to_mu3prong->at(0) == 0) continue;
    if (TREEMCs[s]->triggered->at(0) == 0) continue;

    passedmu  = false;
    passedtau = false;
    passedcalo = true;
    bool passedMET = true;
    temp_tau_pt_comparison = 0.;
    muon_charge = 0;
    temp_nch = 0.;
    temp_npv = 0.;
    temp_pvz = 0.;
    matchedpT_rank = 0;
    
    for (int i=0; i<(int)TREEMCs[s]->BsTauTau_mu1_pt->size(); i++) {
      if (TREEMCs[s]->BsTauTau_mu1_isSoft->at(i) != tau_muon_isSoft) { continue; }
      //if (TREEMCs[s]->BsTauTau_mu1_isTracker->at(i) != tau_muon_isSoft) { continue; }
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) < 3.5) { continue; }
      passedmu = true;
      if (TREEMCs[s]->BsTauTau_mu1_pt->at(i) > temp_tau_pt_comparison) {
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_mu1_pt->at(i);
        tau_muon.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_mu1_pt->at(i), TREEMCs[s]->BsTauTau_mu1_eta->at(i), TREEMCs[s]->BsTauTau_mu1_phi->at(i), muon_mass);
        muon_charge = TREEMCs[s]->BsTauTau_mu1_q->at(i);
      }
    }
    
    temp_tau_pt_comparison = 0.;
    tauh_charge = 0;
    candCounter = 0;
    temp_tau_rho1 = 0.; temp_tau_rho2 = 0.;
    
    //if(!TREEMCs[s]->BsTauTau_tau_pt) continue;
    float vprob = -1;
    float temp_pt = -10;
    for (int i=0; i<(int)TREEMCs[s]->BsTauTau_tau_pt->size(); i++) {
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_pt){
        vprob = 100*TREEMCs[s]->BsTauTau_tau_vprob->at(i);
        temp_pt = TREEMCs[s]->BsTauTau_tau_pt->at(i);
      }
      if (TREEMCs[s]->BsTauTau_tau_vprob->at(i) < tau_hadron_vertexprob) { continue; }
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) < 3.0) { continue; }
      if (muon_charge*TREEMCs[s]->BsTauTau_tau_q->at(i) == 1) candCounter += 1;
      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) / TREEMCs[s]->BsTauTau_tau_matched_gentaupt->at(i) > 0.85) matchedpT_rank += 1;
      passedtau = true;

      if (TREEMCs[s]->BsTauTau_tau_pt->at(i) > temp_tau_pt_comparison && muon_charge*TREEMCs[s]->BsTauTau_tau_q->at(i) != 1) {
        tauh_charge = TREEMCs[s]->BsTauTau_tau_q->at(i);
        temp_tau_pt_comparison = TREEMCs[s]->BsTauTau_tau_pt->at(i);
        //tau_hadron.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_tau_pt->at(i), TREEMCs[s]->BsTauTau_tau_eta->at(i), TREEMCs[s]->BsTauTau_tau_phi->at(i), muon_mass);
        tau_hadron.SetPtEtaPhiM (TREEMCs[s]->BsTauTau_tau_pt->at(i), TREEMCs[s]->BsTauTau_tau_eta->at(i), TREEMCs[s]->BsTauTau_tau_phi->at(i), TREEMCs[s]->BsTauTau_tau_mass->at(i));
        temp_tau_rho1 = TREEMCs[s]->BsTauTau_tau_rhomass1->at(i); temp_tau_rho2 = TREEMCs[s]->BsTauTau_tau_rhomass2->at(i);
      }
    } // loop over the size of the reco taus
    h_tau_hadron_vprob[s+1]->Fill(vprob);
    
    h_tau_hadron_ncand_final[s+1]->Fill(candCounter);
    h_tau_hadron_matched_pt_index[s+1]->Fill(matchedpT_rank);
    
    //if ((tau_muon+tau_hadron).Pt() > MET_cut) passedMET = false;
    
    if (!TREEMCs[s]->BsTauTau_nch->empty()) temp_nch = TREEMCs[s]->BsTauTau_nch->at(0);
    if (!TREEMCs[s]->BsTauTau_bbPV_vz->empty()) temp_pvz = TREEMCs[s]->BsTauTau_bbPV_vz->at(0);
    temp_npv = TREEMCs[s]->PV_N;

#ifdef nch_cut
    if (passedmu && passedtau && temp_nch<=max_nch)
#else
    if (passedmu && passedtau)
#endif
    {
      double sumHFp = 0;
      double maxHFp = 0;
      int sizeHFp = 0;
      double sumHFm = 0;
      double maxHFm = 0;
      int sizeHFm = 0;
      for (int i=0; i<(int)TREEMCs[s]->BsTauTau_calo_eta->size(); i++) {
        double eHFp = TREEMCs[s]->BsTauTau_calo_energyHFp->at(i);
        double eHFm = TREEMCs[s]->BsTauTau_calo_energyHFm->at(i);
        double etaCal = TREEMCs[s]->BsTauTau_calo_eta->at(i);
        if (eHFp != -1){
          h_calo_energyHFp[s+1]->Fill(eHFp); sizeHFp++; sumHFp += eHFp; h_calo_HF_eta[s+1]->Fill(etaCal); h_calo_HF_energy_eta[s+1]->Fill(etaCal,eHFp);}
        if (eHFm != -1){
          h_calo_energyHFm[s+1]->Fill(eHFm); sizeHFm++; sumHFm += eHFm; h_calo_HF_eta[s+1]->Fill(etaCal); h_calo_HF_energy_eta[s+1]->Fill(etaCal,eHFm);}
        if (eHFp > maxHFp) maxHFp = eHFp;
        if (eHFm > maxHFm) maxHFm = eHFm;
      }
      if (maxHFp > HFpLeading_high || maxHFm > HFmLeading_high || maxHFp < HFpLeading_low || maxHFm < HFmLeading_low) passedcalo = false;
      if (passedcalo) h_calo_energyHFp_sum[s+1]->Fill(sumHFp);
      if (passedcalo) h_calo_energyHFm_sum[s+1]->Fill(sumHFm);
      if (passedcalo) h_calo_energyHF_pm[s+1]->Fill(sumHFm,sumHFp);
      if (passedcalo) h_calo_energyHFp_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),sumHFp); h_calo_energyHFm_nch->Fill(TREEMCs[s]->BsTauTau_nch->at(0),sumHFm);
      if (passedcalo) h_calo_energyHFp_size[s+1]->Fill(sizeHFp);
      if (passedcalo) h_calo_energyHFm_size[s+1]->Fill(sizeHFm);
      h_calo_leadingHFp[s+1]->Fill(maxHFp);
      h_calo_leadingHFm[s+1]->Fill(maxHFm);
      if (passedcalo) h_calo_leadingHF_pm[s+1]->Fill(maxHFm,maxHFp);
      //if (sumHFp < 95 || sumHFp > 160 || sumHFm < 95 || sumHFm > 160) passedcalo = false;
    }

#ifdef nch_cut
    if (passedmu && passedtau && temp_nch<=max_nch && passedcalo && passedMET)
#else
    if (passedmu && passedtau && passedcalo && passedMET)
#endif
    {
      if (!TREEMCs[s]->BsTauTau_tau_isRight->at(0)) continue;
      if (temp_nch == 14) cout << "event: " << TREEMCs[s]->EVENT_event << "  Run: " << TREEMCs[s]->EVENT_run << "  lumi block: " << TREEMCs[s]->EVENT_lumiBlock << endl;
      charge_counter[s+1][muon_charge*tauh_charge + 1] += 1;
      h_tau_mu_pt[s+1]->Fill(tau_muon.Pt());
      h_tau_mu_eta[s+1]->Fill(tau_muon.Eta());
      h_tau_mu_phi[s+1]->Fill(tau_muon.Phi());
      h_tau_hadron_pt[s+1]->Fill(tau_hadron.Pt());
      h_tau_hadron_eta[s+1]->Fill(tau_hadron.Eta());
      h_tau_hadron_phi[s+1]->Fill(tau_hadron.Phi());
      h_tau_hadron_rhomass[s+1][0]->Fill(temp_tau_rho1);
      h_tau_hadron_rhomass[s+1][1]->Fill(temp_tau_rho2);
      h_tau_hadron_rhomass2D[s+1]->Fill(temp_tau_rho1,temp_tau_rho2);
      h_tau_hadron_nch[s+1]->Fill(temp_nch);
      h_tau_hadron_mass[s+1]->Fill(tau_hadron.M());
//      h_tau_hadron_track_pvz[s+1][0]->Fill(tau_track1.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][1]->Fill(tau_track2.Z() - temp_pvz);
//      h_tau_hadron_track_pvz[s+1][2]->Fill(tau_track3.Z() - temp_pvz);
      h_deltaphi_tau_mu_tau_hadron[s+1]->Fill(TMath::Abs(tau_muon.DeltaPhi(tau_hadron)));
      h_deltaphi_tau_mu_tau_hadron_mueta[s+1]->Fill(TMath::Abs(tau_muon.DeltaPhi(tau_hadron)),tau_muon.Eta());
      h_deltaphi_tau_mu_tau_hadron_deltaeta[s+1]->Fill(TMath::Abs(tau_muon.DeltaPhi(tau_hadron)),TMath::Abs(tau_hadron.Eta())-TMath::Abs(tau_muon.Eta()));
      h_mueta_taueta[s+1]->Fill(tau_muon.Eta(),tau_hadron.Eta());
      h_PV_N[s+1]->Fill(temp_npv);
      TLorentzVector MET;
      MET.SetPtEtaPhiM((tau_muon+tau_hadron).Pt(),(tau_muon+tau_hadron).Eta(),-(tau_muon+tau_hadron).Phi(),0);
      h_MET[s+1]->Fill(MET.Pt());
      h_ditau_mass[s+1]->Fill((tau_muon+tau_hadron+MET).M());
      h_AP[s+1]->Fill((tau_muon.Pz()-tau_hadron.Pz()) / (tau_muon.Pz()+tau_hadron.Pz()),(tau_muon.Pt()+tau_hadron.Pt())/2);
      //h_MET[s+1]->Fill(TREEMCs[s]->MET_sumEt->at(0));
    }

  } // loop on entries
  cout << "Muon charge times hadronic tau charge for " << tag[s+1] << ":\n -1: " << charge_counter[s+1][0] << "\n 0: " << charge_counter[s+1][1] << "\n +1: " << charge_counter[s+1][2] << endl;
  } // loop on samples

  // *******************************************
  // *******************************************
  // Plotting starts here
  // *******************************************
  // *******************************************
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0,1);
  
  TLegend *legend = new TLegend(0.25,0.6,0.6,0.9);
  
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
    cout << SF2[i] << " / " << SF2[i]*int(h_tau_mu_pt[i]->GetEntries());
    if (i != nSamples-1) cout << ", ";
  }
  cout << endl;

  for (int i = 0; i < nSamples; i++){
    if (i != 0 && h_tau_hadron_nch[i]->GetMaximum()!=0) SF[i] = 0.8 * h_tau_hadron_nch[0]->GetMaximum()/h_tau_hadron_nch[i]->GetMaximum();
    cout << "Temporary SF: " << SF[i] << endl;
    h_tau_mu_pt[i]->Scale(SF[i]); h_tau_mu_eta[i]->Scale(SF[i]); h_tau_mu_phi[i]->Scale(SF[i]);
    h_tau_hadron_pt[i]->Scale(SF[i]); h_tau_hadron_eta[i]->Scale(SF[i]); h_tau_hadron_phi[i]->Scale(SF[i]);
    for (int j = 0; j < 2; j++) h_tau_hadron_rhomass[i][j]->Scale(SF[i]);
    h_tau_hadron_nch[i]->Scale(SF[i]); h_tau_hadron_mass[i]->Scale(SF[i]);
    h_tau_hadron_ncand_final[i]->Scale(h_tau_hadron_ncand_final[0]->GetBinContent(1) / h_tau_hadron_ncand_final[i]->GetBinContent(1));
    h_tau_hadron_vprob[i]->Scale(h_tau_hadron_vprob[0]->GetBinContent(1) / h_tau_hadron_vprob[i]->GetBinContent(1));
    for (int j = 0; j < 3; j++) h_tau_hadron_track_pvz[i][j]->Scale(SF[i]);
    h_deltaphi_tau_mu_tau_hadron[i]->Scale(SF[i]);
    h_PV_N[i]->Scale(SF[i]);
    h_calo_energyHFp[i]->Scale(SF[i]);
    //h_calo_energyHFp_sum[i]->Scale(1./h_calo_energyHFp_sum[i]->GetMaximum());
    h_calo_energyHFp_sum[i]->Scale(SF[i]);
    h_calo_energyHFp_size[i]->Scale(1./h_calo_energyHFp_size[i]->GetMaximum());
    h_calo_energyHFm[i]->Scale(SF[i]);
    //h_calo_energyHFm_sum[i]->Scale(1./h_calo_energyHFm_sum[i]->GetMaximum());
    h_calo_energyHFm_sum[i]->Scale(SF[i]);
    h_calo_energyHFm_size[i]->Scale(1./h_calo_energyHFm_size[i]->GetMaximum());
    h_calo_HF_eta[i]->Scale(SF[i]);
    h_calo_leadingHFp[i]->Scale(SF[i]);
    h_calo_leadingHFm[i]->Scale(SF[i]);
    h_MET[i]->Scale(SF[i]);
    h_ditau_mass[i]->Scale(SF[i]);
    cutflow[i]->Scale(SF2[i]);
  }

  TCanvas *c_tau_had_pv = new TCanvas("c_tau_had_pv", "c_tau_had_pv", 1500, 500); c_tau_had_pv->Divide(3,1);
  for (int j = 0; j < 3; j++){
    for (int i = 1; i < nSamples; i++) hs_tau_hadron_track_pvz[j]->Add(h_tau_hadron_track_pvz[i][j]);
    c_tau_had_pv->cd(j+1); hs_tau_hadron_track_pvz[j]->Draw("he"); h_tau_hadron_track_pvz[0][j]->Draw("e1same");
    //c_tau_had_pv->cd(j+1); h_tau_hadron_track_pvz[0][j]->Draw("e1"); h_tau_hadron_track_pvz[0][j]->GetYaxis()->SetRangeUser(0., 0.6);
    //for (int i = 1; i < nSamples; i++) h_tau_hadron_track_pvz[i][j]->Draw("hesame"); h_tau_hadron_track_pvz[i][j]->GetYaxis()->SetRangeUser(0., 0.1);
  }
#ifdef nch_cut
  c_tau_had_pv->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/tau_had_pv.png");
#else
  c_tau_had_pv->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/tau_had_pv.png");
#endif
    
  TCanvas *c_other_plots = new TCanvas("c_other_plots", "c_deltaPhi_PV", 1300, 500); c_other_plots->Divide(3,1);
  for (int i = 1; i < nSamples; i++) hs_deltaphi_tau_mu_tau_hadron->Add(h_deltaphi_tau_mu_tau_hadron[i]);
  if (hs_deltaphi_tau_mu_tau_hadron->GetMaximum() < h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum()) hs_deltaphi_tau_mu_tau_hadron->SetMaximum(h_deltaphi_tau_mu_tau_hadron[0]->GetMaximum());
  c_other_plots->cd(1); hs_deltaphi_tau_mu_tau_hadron->Draw("he"); h_deltaphi_tau_mu_tau_hadron[0]->Draw("e1same"); 
  legend = new TLegend(0.5,0.6,0.75,0.88);
  legend->SetFillColor(0); legend->SetLineColor(0); legend->SetShadowColor(0); legend->SetTextSize(0.035);
  legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[0],    (tag[0]).c_str(), "ep");
  for (int i = 1; i < nSamples; i++) legend->AddEntry(h_deltaphi_tau_mu_tau_hadron[i], (tag[i]).c_str(), "f");
  legend->Draw();
  c_other_plots->cd(2); h_deltaphi_tau_mu_tau_hadron_nch->Draw("COLZ");
  for (int i = 1; i < nSamples; i++) hs_PV_N->Add(h_PV_N[i]);
  c_other_plots->cd(3); h_PV_N[0]->Draw("e1"); hs_PV_N->Draw("hesame");
#ifdef nch_cut
  c_other_plots->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/deltaPhi_PV.png");
#else
  c_other_plots->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/deltaPhi_PV.png");
#endif

  TCanvas *c_tau_muon_kinematics = new TCanvas("c_tau_muon_kinematics", "c_tau_muon_kinematics", 1500, 500); c_tau_muon_kinematics->Divide(3,1);
  c_tau_muon_kinematics->cd(1); for (int i = 1; i < nSamples; i++) hs_tau_mu_pt->Add(h_tau_mu_pt[i]);
  if (hs_tau_mu_pt->GetMaximum() < h_tau_mu_pt[0]->GetMaximum()) hs_tau_mu_pt->SetMaximum(h_tau_mu_pt[0]->GetMaximum());
  hs_tau_mu_pt->Draw("he"); h_tau_mu_pt[0]->Draw("e1same"); legend->Draw();
  c_tau_muon_kinematics->cd(2); for (int i = 1; i < nSamples; i++) hs_tau_mu_eta->Add(h_tau_mu_eta[i]);
  if (hs_tau_mu_eta->GetMaximum() < h_tau_mu_eta[0]->GetMaximum()) hs_tau_mu_eta->SetMaximum(h_tau_mu_eta[0]->GetMaximum());
  hs_tau_mu_eta->Draw("he"); h_tau_mu_eta[0]->Draw("e1same");
  c_tau_muon_kinematics->cd(3); for (int i = 1; i < nSamples; i++) hs_tau_mu_phi->Add(h_tau_mu_phi[i]);
  if (hs_tau_mu_phi->GetMaximum() < h_tau_mu_phi[0]->GetMaximum()) hs_tau_mu_phi->SetMaximum(h_tau_mu_phi[0]->GetMaximum());
  hs_tau_mu_phi->Draw("he"); h_tau_mu_phi[0]->Draw("e1same");
#ifdef nch_cut
  c_tau_muon_kinematics->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/tau_muon_kinematics.png");
#else
  c_tau_muon_kinematics->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/tau_muon_kinematics.png");
#endif

TCanvas *c_tau_had = new TCanvas("c_tau_had", "c_tau_had", 2000, 1500); c_tau_had->Divide(4,3);
  c_tau_had->cd(1); for (int i = 1; i < nSamples; i++) hs_tau_hadron_pt->Add(h_tau_hadron_pt[i]);
  if (hs_tau_hadron_pt->GetMaximum() < h_tau_hadron_pt[0]->GetMaximum()) hs_tau_hadron_pt->SetMaximum(h_tau_hadron_pt[0]->GetMaximum());
  hs_tau_hadron_pt->Draw("he"); h_tau_hadron_pt[0]->Draw("e1same");
  c_tau_had->cd(2); for (int i = 1; i < nSamples; i++) hs_tau_hadron_eta->Add(h_tau_hadron_eta[i]);
  if (hs_tau_hadron_eta->GetMaximum() < h_tau_hadron_eta[0]->GetMaximum()) hs_tau_hadron_eta->SetMaximum(h_tau_hadron_eta[0]->GetMaximum());
  hs_tau_hadron_eta->Draw("he"); h_tau_hadron_eta[0]->Draw("e1same");
  c_tau_had->cd(3); for (int i = 1; i < nSamples; i++) hs_tau_hadron_phi->Add(h_tau_hadron_phi[i]);
  if (hs_tau_hadron_phi->GetMaximum() < h_tau_hadron_phi[0]->GetMaximum()) hs_tau_hadron_phi->SetMaximum(h_tau_hadron_phi[0]->GetMaximum());
  hs_tau_hadron_phi->Draw("he"); h_tau_hadron_phi[0]->Draw("e1same");
  for (int j = 0; j < 2; j++){
    c_tau_had->cd(j+4);
    for (int i = 1; i < nSamples; i++) hs_tau_hadron_rhomass[j]->Add(h_tau_hadron_rhomass[i][j]);
    if (hs_tau_hadron_rhomass[j]->GetMaximum() < h_tau_hadron_rhomass[0][j]->GetMaximum()) hs_tau_hadron_rhomass[j]->SetMaximum(h_tau_hadron_rhomass[0][j]->GetMaximum());
    hs_tau_hadron_rhomass[j]->Draw("he"); h_tau_hadron_rhomass[0][j]->Draw("e1same");
  }
  c_tau_had->cd(6); for (int i = 1; i < nSamples; i++) hs_tau_hadron_nch->Add(h_tau_hadron_nch[i]);
  if (hs_tau_hadron_nch->GetMaximum() < h_tau_hadron_nch[0]->GetMaximum()) hs_tau_hadron_nch->SetMaximum(h_tau_hadron_nch[0]->GetMaximum());
  /*hs_tau_hadron_nch->Draw("he");*/ h_tau_hadron_nch[0]->Draw("e1same"); for (int i = 1; i < nSamples; i++) h_tau_hadron_nch[i]->Draw("hesame"); legend->Draw();
  c_tau_had->cd(7); h_tau_hadron_ncand_final[0]->Draw("e1same"); for (int i = 1; i < nSamples; i++) h_tau_hadron_ncand_final[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(8); h_tau_hadron_vprob[0]->Draw("e1same"); for (int i = 1; i < nSamples; i++) h_tau_hadron_vprob[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(9); for (int i = 1; i < nSamples; i++) h_tau_hadron_matched_pt_index[i]->Draw("hesame"); legend->Draw(); 
  c_tau_had->cd(10); for (int i = 1; i < nSamples; i++) hs_tau_hadron_mass->Add(h_tau_hadron_mass[i]);
  if (hs_tau_hadron_mass->GetMaximum() < h_tau_hadron_mass[0]->GetMaximum()) hs_tau_hadron_mass->SetMaximum(h_tau_hadron_mass[0]->GetMaximum());
  hs_tau_hadron_mass->Draw("he"); h_tau_hadron_mass[0]->Draw("e1same"); //legend->Draw();
  c_tau_had->cd(11); for (int i = 1; i < nSamples; i++) hs_ditau_mass->Add(h_ditau_mass[i]);
  if (hs_ditau_mass->GetMaximum() < h_ditau_mass[0]->GetMaximum()) hs_ditau_mass->SetMaximum(h_ditau_mass[0]->GetMaximum());
  hs_ditau_mass->Draw("he"); h_ditau_mass[0]->Draw("e1same"); //legend->Draw();
#ifdef nch_cut
  c_tau_had->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/tau_had.png");
#else
  c_tau_had->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/tau_had.png");
#endif

  TCanvas *c_track_activity = new TCanvas("c_track_activity", "c_track_activity", 1200, 600); c_track_activity->Divide(2,1);
//  c_track_activity->cd(1); h_track_activity_pt->Draw("e1");
  c_track_activity->cd(1); h_track_activity_pt->DrawNormalized("e1");
//  c_track_activity->cd(2); h_track_activity_pt_eta->Draw("COLZ");
  c_track_activity->cd(2); h_track_activity_pt_eta->DrawNormalized("COLZ");
#ifdef nch_cut
  c_track_activity->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/track_activity.png");
#else
  c_track_activity->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/track_activity.png");
#endif

  TCanvas *c_calo = new TCanvas("c_calo", "c_calo", 1600, 5*800); c_calo->Divide(2,5);
  for (int i = 1; i < nSamples; i++) hs_calo_energyHFp->Add(h_calo_energyHFp[i]);
  c_calo->cd(1); h_calo_energyHFp[0]->Draw("e1"); h_calo_energyHFp[0]->GetYaxis()->SetRangeUser(1,2*h_calo_energyHFp[0]->GetMaximum());
  c_calo->cd(1)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_energyHFp[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++) hs_calo_energyHFm->Add(h_calo_energyHFm[i]);
  c_calo->cd(2); h_calo_energyHFm[0]->Draw("e1"); h_calo_energyHFm[0]->GetYaxis()->SetRangeUser(1,2*h_calo_energyHFm[0]->GetMaximum());
  c_calo->cd(2)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_energyHFm[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++) hs_calo_leadingHFp->Add(h_calo_leadingHFp[i]);
  c_calo->cd(3); h_calo_leadingHFp[0]->Draw("e1"); h_calo_leadingHFp[0]->GetYaxis()->SetRangeUser(1,2*h_calo_leadingHFp[0]->GetMaximum());
  c_calo->cd(3)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_leadingHFp[i]->Draw("hesame"); legend->Draw();
  for (int i = 1; i < nSamples; i++) hs_calo_leadingHFm->Add(h_calo_leadingHFm[i]);
  c_calo->cd(4); h_calo_leadingHFm[0]->Draw("e1"); h_calo_leadingHFm[0]->GetYaxis()->SetRangeUser(1,2*h_calo_leadingHFm[0]->GetMaximum());
  c_calo->cd(4)->SetLogy(1); for (int i = 1; i < nSamples; i++) h_calo_leadingHFm[i]->Draw("hesame"); legend->Draw();
  c_calo->cd(5); h_calo_energyHFp_nch->Draw("COLZ");
  c_calo->cd(6); h_calo_energyHFm_nch->Draw("COLZ");
  c_calo->cd(7); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFp_sum[i]->Draw("e1same");
  c_calo->cd(8); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFm_sum[i]->Draw("e1same");
  c_calo->cd(9); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFp_size[i]->Draw("e1same");
  c_calo->cd(10); for (int i = 0; i < nSamples; i++) if (i<4) h_calo_energyHFm_size[i]->Draw("e1same");
#ifdef nch_cut
  c_calo->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/calo.png");
#else
  c_calo->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/calo.png");
#endif

  TCanvas *c_HF_eta = new TCanvas("c_HF_eta", "c_HF_eta", 500*(nSamples+1), 500); c_HF_eta->Divide(nSamples+1,1);
  for (int i = 1; i < nSamples; i++) hs_calo_HF_eta->Add(h_calo_HF_eta[i]);
  c_HF_eta->cd(1); h_calo_HF_eta[0]->Draw("e1"); for (int i = 1; i < nSamples; i++) h_calo_HF_eta[i]->Draw("hesame"); legend->Draw();
  for (int i = 0; i < nSamples; i++){ c_HF_eta->cd(i+2); h_calo_HF_energy_eta[i]->Draw("colz"); c_HF_eta->cd(i+2)->SetLogz(1);}
#ifdef nch_cut
  c_HF_eta->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/HF_eta.png");
#else
  c_HF_eta->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/HF_eta.png");
#endif

  TCanvas *c_HFpm = new TCanvas("c_HFpm", "c_HFpm", 500*nSamples, 1000); c_HFpm->Divide(nSamples,2);
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(i+1); h_calo_energyHF_pm[i]->Draw("colz"); }
  for (int i = 0; i < nSamples; i++){ c_HFpm->cd(nSamples+i+1); h_calo_leadingHF_pm[i]->Draw("colz"); }
#ifdef nch_cut
  c_HFpm->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/HFpm.png");
#else
  c_HFpm->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/HFpm.png");
#endif

  TCanvas *c_delta_phi_eta = new TCanvas("c_delta_phi_eta", "c_delta_phi_eta", 500*nSamples, 1500); c_delta_phi_eta->Divide(nSamples,3);
  for (int i = 0; i < nSamples; i++){
    c_delta_phi_eta->cd(i+1); h_deltaphi_tau_mu_tau_hadron_mueta[i]->Draw("colz");
    c_delta_phi_eta->cd(nSamples+i+1); h_deltaphi_tau_mu_tau_hadron_deltaeta[i]->Draw("colz");
    c_delta_phi_eta->cd(2*nSamples+i+1); h_mueta_taueta[i]->Draw("colz");
  }
#ifdef nch_cut
  c_delta_phi_eta->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/delta_phi_eta.png");
#else
  c_delta_phi_eta->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/delta_phi_eta.png");
#endif

  TCanvas *c_AP = new TCanvas("c_AP", "c_AP", 500*nSamples, 500); c_AP->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_AP->cd(i+1); h_AP[i]->Draw("colz");}
#ifdef nch_cut
  c_AP->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/AP.png");
#else
  c_AP->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/AP.png");
#endif

  TCanvas *c_rho = new TCanvas("c_rho", "c_rho", 500*nSamples, 500); c_rho->Divide(nSamples,1);
  for (int i = 0; i < nSamples; i++) {c_rho->cd(i+1); h_tau_hadron_rhomass2D[i]->Draw("colz");}
#ifdef nch_cut
  c_rho->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/rho.png");
#else
  c_rho->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/rho.png");
#endif

  TCanvas *c_MET = new TCanvas("c_MET", "c_MET", 500, 500);
  c_MET->cd(1); h_MET[0]->Draw("e1"); for (int i = 1; i < nSamples; i++) h_MET[i]->Draw("hesame");
#ifdef nch_cut
  c_MET->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/MET.png");
#else
  c_MET->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/MET.png");
#endif

  TCanvas *c_cutflow = new TCanvas("c_cutflow", "c_cutflow", 500, 500);
  cutflow[0]->Draw("e1");
  cutflow[0]->GetYaxis()->SetRangeUser(1,2*cutflow[0]->GetMaximum());
  for (int i = 1; i < nSamples; i++) cutflow[i]->Draw("hesame");
  legend->Draw();
  c_cutflow->cd(1)->SetLogy(1);
#ifdef nch_cut
  c_cutflow->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/nch_cut/cutflow.png");
#else
  c_cutflow->SaveAs("/eos/user/a/ajofrehe/www/gtau/plots/cutflow.png");
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
  TCanvas *c_output = new TCanvas("c_output", "c_output", 400*nSamples, 500); c_output->Divide(nSamples,1);
  TH1F *h_deltaphi_tau_mu_tau_hadron_new = (TH1F*)h_deltaphi_tau_mu_tau_hadron[0]->Clone("data_obs");
  c_output->cd(1); h_deltaphi_tau_mu_tau_hadron_new->DrawCopy("e1");
  h_deltaphi_tau_mu_tau_hadron_new->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron_new->Write();
  for (int i = 1; i < nSamples; i++){
    c_output->cd(i+1); h_deltaphi_tau_mu_tau_hadron[i]->DrawCopy("e1");
    h_deltaphi_tau_mu_tau_hadron[i]->SetDirectory(gDirectory);h_deltaphi_tau_mu_tau_hadron[i]->Write();
  }
  //outputFile->Write();
  outputFile->Close();

//  ntuple->Write();
//  ntuple->Close();

}
