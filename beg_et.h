//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 17 01:38:00 2021 by ROOT version 6.22/02
// from TTree event_tree/event_tree
// found on file: prod6/spbeagpb6/out_xaa.root
//////////////////////////////////////////////////////////

#ifndef beg_et_h
#define beg_et_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class beg_et {
public :


  TH1F * h_mass;
  TH1F * h_mass_bc;

  //  TH1F * h_mass;
   TH1F * h_px0;
   TH1F * h_pxb;
   TH1F * h_ttrue;
   TH1F * h_ttrue_bc;
   TH1F * h_t_trans_rot;
   TH1F * h_t_trans;
   TH1F * h_t_trans2;
   TH2F * h2px;
   TH1F * h_ttrueff;
   TH1F * h_t_trans_rotff;
   TH1F * h_t_transff;


   TF1 * f_egausres1;
   //   TF1 * f_egausres2;



   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         cross_section;
   Float_t         event_weight;
   Int_t           n_generator_accepted;
   Int_t           nHits;
   Int_t           hits_layerID[135];   //[nHits]
   Int_t           hits_trueID[135];   //[nHits]
   Float_t         hits_x[135];   //[nHits]
   Float_t         hits_y[135];   //[nHits]
   Float_t         hits_z[135];   //[nHits]
   Float_t         hits_x2[135];   //[nHits]
   Float_t         hits_y2[135];   //[nHits]
   Float_t         hits_z2[135];   //[nHits]
   Float_t         hits_t[135];   //[nHits]
   Float_t         hits_edep[135];   //[nHits]
   Int_t           nTracks;
   Float_t         tracks_ID[29];   //[nTracks]
   Short_t         tracks_charge[29];   //[nTracks]
   Float_t         tracks_px[29];   //[nTracks]
   Float_t         tracks_py[29];   //[nTracks]
   Float_t         tracks_pz[29];   //[nTracks]
   Float_t         tracks_dca[29];   //[nTracks]
   Float_t         tracks_dca_2d[29];   //[nTracks]
   Float_t         tracks_trueID[29];   //[nTracks]
   UShort_t        tracks_source[29];   //[nTracks]
   Float_t         track_pion_LL[29];   //[nTracks]
   Float_t         track_kaon_LL[29];   //[nTracks]
   Float_t         track_proton_LL[29];   //[nTracks]
   Int_t           nProjections;
   Float_t         track_ProjTrackID[49];   //[nProjections]
   Int_t           track_ProjLayer[49];   //[nProjections]
   Float_t         track_TLP_x[49];   //[nProjections]
   Float_t         track_TLP_y[49];   //[nProjections]
   Float_t         track_TLP_z[49];   //[nProjections]
   Float_t         track_TLP_t[49];   //[nProjections]
   Float_t         track_TLP_true_x[49];   //[nProjections]
   Float_t         track_TLP_true_y[49];   //[nProjections]
   Float_t         track_TLP_true_z[49];   //[nProjections]
   Float_t         track_TLP_true_t[49];   //[nProjections]
   Int_t           tower_FHCAL_N;
   Float_t         tower_FHCAL_E[1];   //[tower_FHCAL_N]
   Int_t           tower_FHCAL_iEta[1];   //[tower_FHCAL_N]
   Int_t           tower_FHCAL_iPhi[1];   //[tower_FHCAL_N]
   Int_t           tower_FHCAL_trueID[1];   //[tower_FHCAL_N]
   Int_t           cluster_FHCAL_N;
   Float_t         cluster_FHCAL_E[1];   //[cluster_FHCAL_N]
   Float_t         cluster_FHCAL_Eta[1];   //[cluster_FHCAL_N]
   Float_t         cluster_FHCAL_Phi[1];   //[cluster_FHCAL_N]
   Int_t           cluster_FHCAL_NTower[1];   //[cluster_FHCAL_N]
   Int_t           cluster_FHCAL_trueID[1];   //[cluster_FHCAL_N]
   Int_t           tower_BECAL_N;
   Float_t         tower_BECAL_E[354];   //[tower_BECAL_N]
   Int_t           tower_BECAL_iEta[354];   //[tower_BECAL_N]
   Int_t           tower_BECAL_iPhi[354];   //[tower_BECAL_N]
   Int_t           tower_BECAL_trueID[354];   //[tower_BECAL_N]
   Int_t           tower_HCALIN_N;
   Float_t         tower_HCALIN_E[65];   //[tower_HCALIN_N]
   Int_t           tower_HCALIN_iEta[65];   //[tower_HCALIN_N]
   Int_t           tower_HCALIN_iPhi[65];   //[tower_HCALIN_N]
   Int_t           tower_HCALIN_trueID[65];   //[tower_HCALIN_N]
   Int_t           cluster_HCALIN_N;
   Float_t         cluster_HCALIN_E[1];   //[cluster_HCALIN_N]
   Float_t         cluster_HCALIN_Eta[1];   //[cluster_HCALIN_N]
   Float_t         cluster_HCALIN_Phi[1];   //[cluster_HCALIN_N]
   Int_t           cluster_HCALIN_NTower[1];   //[cluster_HCALIN_N]
   Int_t           cluster_HCALIN_trueID[1];   //[cluster_HCALIN_N]
   Int_t           tower_HCALOUT_N;
   Float_t         tower_HCALOUT_E[11];   //[tower_HCALOUT_N]
   Int_t           tower_HCALOUT_iEta[11];   //[tower_HCALOUT_N]
   Int_t           tower_HCALOUT_iPhi[11];   //[tower_HCALOUT_N]
   Int_t           tower_HCALOUT_trueID[11];   //[tower_HCALOUT_N]
   Int_t           cluster_HCALOUT_N;
   Float_t         cluster_HCALOUT_E[3];   //[cluster_HCALOUT_N]
   Float_t         cluster_HCALOUT_Eta[3];   //[cluster_HCALOUT_N]
   Float_t         cluster_HCALOUT_Phi[3];   //[cluster_HCALOUT_N]
   Int_t           cluster_HCALOUT_NTower[3];   //[cluster_HCALOUT_N]
   Int_t           cluster_HCALOUT_trueID[3];   //[cluster_HCALOUT_N]
   Int_t           tower_FEMC_N;
   Float_t         tower_FEMC_E[623];   //[tower_FEMC_N]
   Int_t           tower_FEMC_iEta[623];   //[tower_FEMC_N]
   Int_t           tower_FEMC_iPhi[623];   //[tower_FEMC_N]
   Int_t           tower_FEMC_trueID[623];   //[tower_FEMC_N]
   Int_t           cluster_FEMC_N;
   Float_t         cluster_FEMC_E[29];   //[cluster_FEMC_N]
   Float_t         cluster_FEMC_Eta[29];   //[cluster_FEMC_N]
   Float_t         cluster_FEMC_Phi[29];   //[cluster_FEMC_N]
   Int_t           cluster_FEMC_NTower[29];   //[cluster_FEMC_N]
   Int_t           cluster_FEMC_trueID[29];   //[cluster_FEMC_N]
   Int_t           tower_EEMC_N;
   Float_t         tower_EEMC_E[143];   //[tower_EEMC_N]
   Int_t           tower_EEMC_iEta[143];   //[tower_EEMC_N]
   Int_t           tower_EEMC_iPhi[143];   //[tower_EEMC_N]
   Int_t           tower_EEMC_trueID[143];   //[tower_EEMC_N]
   Int_t           cluster_EEMC_N;
   Float_t         cluster_EEMC_E[10];   //[cluster_EEMC_N]
   Float_t         cluster_EEMC_Eta[10];   //[cluster_EEMC_N]
   Float_t         cluster_EEMC_Phi[10];   //[cluster_EEMC_N]
   Int_t           cluster_EEMC_NTower[10];   //[cluster_EEMC_N]
   Int_t           cluster_EEMC_trueID[10];   //[cluster_EEMC_N]
   Float_t         vertex_x;
   Float_t         vertex_y;
   Float_t         vertex_z;
   Int_t           vertex_NCont;
   Float_t         vertex_true_x;
   Float_t         vertex_true_y;
   Float_t         vertex_true_z;
   Int_t           nMCPart;
   Int_t           mcpart_ID[199];   //[nMCPart]
   Int_t           mcpart_ID_parent[199];   //[nMCPart]
   Int_t           mcpart_PDG[199];   //[nMCPart]
   Float_t         mcpart_E[199];   //[nMCPart]
   Float_t         mcpart_px[199];   //[nMCPart]
   Float_t         mcpart_py[199];   //[nMCPart]
   Float_t         mcpart_pz[199];   //[nMCPart]
   Int_t           mcpart_BCID[199];   //[nMCPart]
   Int_t           nHepmcp;
   Int_t           hepmcp_procid;
   Float_t         hepmcp_x1;
   Float_t         hepmcp_x2;
   Float_t         hepmcp_Q2;
   Int_t           hepmcp_status[359];   //[nHepmcp]
   Int_t           hepmcp_PDG[359];   //[nHepmcp]
   Float_t         hepmcp_E[359];   //[nHepmcp]
   Float_t         hepmcp_px[359];   //[nHepmcp]
   Float_t         hepmcp_py[359];   //[nHepmcp]
   Float_t         hepmcp_pz[359];   //[nHepmcp]
   Int_t           hepmcp_BCID[359];   //[nHepmcp]
   Int_t           hepmcp_m1[359];   //[nHepmcp]
   Int_t           hepmcp_m2[359];   //[nHepmcp]

   // List of branches
   TBranch        *b_cross_section;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_n_generator_accepted;   //!
   TBranch        *b_nHits;   //!
   TBranch        *b_hits_layerID;   //!
   TBranch        *b_hits_trueID;   //!
   TBranch        *b_hits_x;   //!
   TBranch        *b_hits_y;   //!
   TBranch        *b_hits_z;   //!
   TBranch        *b_hits_x2;   //!
   TBranch        *b_hits_y2;   //!
   TBranch        *b_hits_z2;   //!
   TBranch        *b_hits_t;   //!
   TBranch        *b_hits_edep;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_tracks_ID;   //!
   TBranch        *b_tracks_charge;   //!
   TBranch        *b_tracks_px;   //!
   TBranch        *b_tracks_py;   //!
   TBranch        *b_tracks_pz;   //!
   TBranch        *b_tracks_dca;   //!
   TBranch        *b_tracks_dca_2d;   //!
   TBranch        *b_tracks_trueID;   //!
   TBranch        *b_tracks_source;   //!
   TBranch        *b_track_pion_LL;   //!
   TBranch        *b_track_kaon_LL;   //!
   TBranch        *b_track_proton_LL;   //!
   TBranch        *b_nProjections;   //!
   TBranch        *b_track_ProjTrackID;   //!
   TBranch        *b_track_ProjLayer;   //!
   TBranch        *b_track_TLP_x;   //!
   TBranch        *b_track_TLP_y;   //!
   TBranch        *b_track_TLP_z;   //!
   TBranch        *b_track_TLP_t;   //!
   TBranch        *b_track_TLP_true_x;   //!
   TBranch        *b_track_TLP_true_y;   //!
   TBranch        *b_track_TLP_true_z;   //!
   TBranch        *b_track_TLP_true_t;   //!
   TBranch        *b_tower_FHCAL_N;   //!
   TBranch        *b_tower_FHCAL_E;   //!
   TBranch        *b_tower_FHCAL_iEta;   //!
   TBranch        *b_tower_FHCAL_iPhi;   //!
   TBranch        *b_tower_FHCAL_trueID;   //!
   TBranch        *b_cluster_FHCAL_N;   //!
   TBranch        *b_cluster_FHCAL_E;   //!
   TBranch        *b_cluster_FHCAL_Eta;   //!
   TBranch        *b_cluster_FHCAL_Phi;   //!
   TBranch        *b_cluster_FHCAL_NTower;   //!
   TBranch        *b_cluster_FHCAL_trueID;   //!
   TBranch        *b_tower_BECAL_N;   //!
   TBranch        *b_tower_BECAL_E;   //!
   TBranch        *b_tower_BECAL_iEta;   //!
   TBranch        *b_tower_BECAL_iPhi;   //!
   TBranch        *b_tower_BECAL_trueID;   //!
   TBranch        *b_tower_HCALIN_N;   //!
   TBranch        *b_tower_HCALIN_E;   //!
   TBranch        *b_tower_HCALIN_iEta;   //!
   TBranch        *b_tower_HCALIN_iPhi;   //!
   TBranch        *b_tower_HCALIN_trueID;   //!
   TBranch        *b_cluster_HCALIN_N;   //!
   TBranch        *b_cluster_HCALIN_E;   //!
   TBranch        *b_cluster_HCALIN_Eta;   //!
   TBranch        *b_cluster_HCALIN_Phi;   //!
   TBranch        *b_cluster_HCALIN_NTower;   //!
   TBranch        *b_cluster_HCALIN_trueID;   //!
   TBranch        *b_tower_HCALOUT_N;   //!
   TBranch        *b_tower_HCALOUT_E;   //!
   TBranch        *b_tower_HCALOUT_iEta;   //!
   TBranch        *b_tower_HCALOUT_iPhi;   //!
   TBranch        *b_tower_HCALOUT_trueID;   //!
   TBranch        *b_cluster_HCALOUT_N;   //!
   TBranch        *b_cluster_HCALOUT_E;   //!
   TBranch        *b_cluster_HCALOUT_Eta;   //!
   TBranch        *b_cluster_HCALOUT_Phi;   //!
   TBranch        *b_cluster_HCALOUT_NTower;   //!
   TBranch        *b_cluster_HCALOUT_trueID;   //!
   TBranch        *b_tower_FEMC_N;   //!
   TBranch        *b_tower_FEMC_E;   //!
   TBranch        *b_tower_FEMC_iEta;   //!
   TBranch        *b_tower_FEMC_iPhi;   //!
   TBranch        *b_tower_FEMC_trueID;   //!
   TBranch        *b_cluster_FEMC_N;   //!
   TBranch        *b_cluster_FEMC_E;   //!
   TBranch        *b_cluster_FEMC_Eta;   //!
   TBranch        *b_cluster_FEMC_Phi;   //!
   TBranch        *b_cluster_FEMC_NTower;   //!
   TBranch        *b_cluster_FEMC_trueID;   //!
   TBranch        *b_tower_EEMC_N;   //!
   TBranch        *b_tower_EEMC_E;   //!
   TBranch        *b_tower_EEMC_iEta;   //!
   TBranch        *b_tower_EEMC_iPhi;   //!
   TBranch        *b_tower_EEMC_trueID;   //!
   TBranch        *b_cluster_EEMC_N;   //!
   TBranch        *b_cluster_EEMC_E;   //!
   TBranch        *b_cluster_EEMC_Eta;   //!
   TBranch        *b_cluster_EEMC_Phi;   //!
   TBranch        *b_cluster_EEMC_NTower;   //!
   TBranch        *b_cluster_EEMC_trueID;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_vertex_NCont;   //!
   TBranch        *b_vertex_true_x;   //!
   TBranch        *b_vertex_true_y;   //!
   TBranch        *b_vertex_true_z;   //!
   TBranch        *b_nMCPart;   //!
   TBranch        *b_mcpart_ID;   //!
   TBranch        *b_mcpart_ID_parent;   //!
   TBranch        *b_mcpart_PDG;   //!
   TBranch        *b_mcpart_E;   //!
   TBranch        *b_mcpart_px;   //!
   TBranch        *b_mcpart_py;   //!
   TBranch        *b_mcpart_pz;   //!
   TBranch        *b_mcpart_BCID;   //!
   TBranch        *b_nHepmcp;   //!
   TBranch        *b_hepmcp_procid;   //!
   TBranch        *b_hepmcp_x1;   //!
   TBranch        *b_hepmcp_x2;   //!
   TBranch        *b_hepmcp_Q2;   //!
   TBranch        *b_hepmcp_status;   //!
   TBranch        *b_hepmcp_PDG;   //!
   TBranch        *b_hepmcp_E;   //!
   TBranch        *b_hepmcp_px;   //!
   TBranch        *b_hepmcp_py;   //!
   TBranch        *b_hepmcp_pz;   //!
   TBranch        *b_hepmcp_BCID;   //!
   TBranch        *b_hepmcp_m1;   //!
   TBranch        *b_hepmcp_m2;   //!

   beg_et(TTree *tree=0);
   virtual ~beg_et();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef beg_et_cxx
beg_et::beg_et(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("prod6/spbeagpb6/out_xaa.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("prod6/spbeagpb6/out_xaa.root");
      }
      f->GetObject("event_tree",tree);

   }

   
   h2px = new TH2F("h2px","",200,-3.5,1.5,1000,0.5,1.5);
   //h_px0 = new TH1F("h_px0","",200,-0.5,0.5);
   h_px0 = new TH1F("h_px0","",2000,0.5,1.5);
   h_pxb = new TH1F("h_pxb","",2000,0.5,1.5);

   f_egausres1 = new TF1("fgau_r1","gaus",-10,10);
   f_egausres1->SetParameters(10e4,0,1);
   cout << f_egausres1->GetRandom() << " rand "<< endl;
   //   f_egausres2 = new TF1("fgau_r2","gaus",-10,10);

   h_mass_bc = new TH1F("h_mass_bc","",90,0.5,1.4);
   h_mass = new TH1F("h_mass","",90,0.5,1.4);


   h_ttrue = new TH1F("h_ttrue","",90,0.0,0.25);
   h_ttrue_bc = new TH1F("h_ttrue_bc","",90,0.0,0.25);
   h_t_trans_rot = new TH1F("h_t_trans_rot","",90,0,0.25);
   h_t_trans = new TH1F("h_t_trans","",90,0,0.25);
   h_t_trans2 = new TH1F("h_t_trans2","",90,0,0.25);

   h_ttrueff = new TH1F("h_ttrueff","",90,0.0,0.25);
   h_t_trans_rotff = new TH1F("h_t_trans_rotff","",90,0,0.25);
   h_t_transff = new TH1F("h_t_transff","",90,0,0.25);

   //   h_px0 = new TH1F("h_px0","",200,-0.5,0.5);
   //h_pxb = new TH1F("h_pxb","",200,-0.5,0.5);


   Init(tree);
}

beg_et::~beg_et()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t beg_et::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t beg_et::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void beg_et::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("cross_section", &cross_section, &b_cross_section);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("n_generator_accepted", &n_generator_accepted, &b_n_generator_accepted);
   fChain->SetBranchAddress("nHits", &nHits, &b_nHits);
   fChain->SetBranchAddress("hits_layerID", hits_layerID, &b_hits_layerID);
   fChain->SetBranchAddress("hits_trueID", hits_trueID, &b_hits_trueID);
   fChain->SetBranchAddress("hits_x", hits_x, &b_hits_x);
   fChain->SetBranchAddress("hits_y", hits_y, &b_hits_y);
   fChain->SetBranchAddress("hits_z", hits_z, &b_hits_z);
   fChain->SetBranchAddress("hits_x2", hits_x2, &b_hits_x2);
   fChain->SetBranchAddress("hits_y2", hits_y2, &b_hits_y2);
   fChain->SetBranchAddress("hits_z2", hits_z2, &b_hits_z2);
   fChain->SetBranchAddress("hits_t", hits_t, &b_hits_t);
   fChain->SetBranchAddress("hits_edep", hits_edep, &b_hits_edep);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("tracks_ID", tracks_ID, &b_tracks_ID);
   fChain->SetBranchAddress("tracks_charge", tracks_charge, &b_tracks_charge);
   fChain->SetBranchAddress("tracks_px", tracks_px, &b_tracks_px);
   fChain->SetBranchAddress("tracks_py", tracks_py, &b_tracks_py);
   fChain->SetBranchAddress("tracks_pz", tracks_pz, &b_tracks_pz);
   fChain->SetBranchAddress("tracks_dca", tracks_dca, &b_tracks_dca);
   fChain->SetBranchAddress("tracks_dca_2d", tracks_dca_2d, &b_tracks_dca_2d);
   fChain->SetBranchAddress("tracks_trueID", tracks_trueID, &b_tracks_trueID);
   fChain->SetBranchAddress("tracks_source", tracks_source, &b_tracks_source);
   fChain->SetBranchAddress("track_pion_LL", track_pion_LL, &b_track_pion_LL);
   fChain->SetBranchAddress("track_kaon_LL", track_kaon_LL, &b_track_kaon_LL);
   fChain->SetBranchAddress("track_proton_LL", track_proton_LL, &b_track_proton_LL);
   fChain->SetBranchAddress("nProjections", &nProjections, &b_nProjections);
   fChain->SetBranchAddress("track_ProjTrackID", track_ProjTrackID, &b_track_ProjTrackID);
   fChain->SetBranchAddress("track_ProjLayer", track_ProjLayer, &b_track_ProjLayer);
   fChain->SetBranchAddress("track_TLP_x", track_TLP_x, &b_track_TLP_x);
   fChain->SetBranchAddress("track_TLP_y", track_TLP_y, &b_track_TLP_y);
   fChain->SetBranchAddress("track_TLP_z", track_TLP_z, &b_track_TLP_z);
   fChain->SetBranchAddress("track_TLP_t", track_TLP_t, &b_track_TLP_t);
   fChain->SetBranchAddress("track_TLP_true_x", track_TLP_true_x, &b_track_TLP_true_x);
   fChain->SetBranchAddress("track_TLP_true_y", track_TLP_true_y, &b_track_TLP_true_y);
   fChain->SetBranchAddress("track_TLP_true_z", track_TLP_true_z, &b_track_TLP_true_z);
   fChain->SetBranchAddress("track_TLP_true_t", track_TLP_true_t, &b_track_TLP_true_t);
   fChain->SetBranchAddress("tower_FHCAL_N", &tower_FHCAL_N, &b_tower_FHCAL_N);
   fChain->SetBranchAddress("tower_FHCAL_E", &tower_FHCAL_E, &b_tower_FHCAL_E);
   fChain->SetBranchAddress("tower_FHCAL_iEta", &tower_FHCAL_iEta, &b_tower_FHCAL_iEta);
   fChain->SetBranchAddress("tower_FHCAL_iPhi", &tower_FHCAL_iPhi, &b_tower_FHCAL_iPhi);
   fChain->SetBranchAddress("tower_FHCAL_trueID", &tower_FHCAL_trueID, &b_tower_FHCAL_trueID);
   fChain->SetBranchAddress("cluster_FHCAL_N", &cluster_FHCAL_N, &b_cluster_FHCAL_N);
   fChain->SetBranchAddress("cluster_FHCAL_E", &cluster_FHCAL_E, &b_cluster_FHCAL_E);
   fChain->SetBranchAddress("cluster_FHCAL_Eta", &cluster_FHCAL_Eta, &b_cluster_FHCAL_Eta);
   fChain->SetBranchAddress("cluster_FHCAL_Phi", &cluster_FHCAL_Phi, &b_cluster_FHCAL_Phi);
   fChain->SetBranchAddress("cluster_FHCAL_NTower", &cluster_FHCAL_NTower, &b_cluster_FHCAL_NTower);
   fChain->SetBranchAddress("cluster_FHCAL_trueID", &cluster_FHCAL_trueID, &b_cluster_FHCAL_trueID);
   fChain->SetBranchAddress("tower_BECAL_N", &tower_BECAL_N, &b_tower_BECAL_N);
   fChain->SetBranchAddress("tower_BECAL_E", tower_BECAL_E, &b_tower_BECAL_E);
   fChain->SetBranchAddress("tower_BECAL_iEta", tower_BECAL_iEta, &b_tower_BECAL_iEta);
   fChain->SetBranchAddress("tower_BECAL_iPhi", tower_BECAL_iPhi, &b_tower_BECAL_iPhi);
   fChain->SetBranchAddress("tower_BECAL_trueID", tower_BECAL_trueID, &b_tower_BECAL_trueID);
   fChain->SetBranchAddress("tower_HCALIN_N", &tower_HCALIN_N, &b_tower_HCALIN_N);
   fChain->SetBranchAddress("tower_HCALIN_E", tower_HCALIN_E, &b_tower_HCALIN_E);
   fChain->SetBranchAddress("tower_HCALIN_iEta", tower_HCALIN_iEta, &b_tower_HCALIN_iEta);
   fChain->SetBranchAddress("tower_HCALIN_iPhi", tower_HCALIN_iPhi, &b_tower_HCALIN_iPhi);
   fChain->SetBranchAddress("tower_HCALIN_trueID", tower_HCALIN_trueID, &b_tower_HCALIN_trueID);
   fChain->SetBranchAddress("cluster_HCALIN_N", &cluster_HCALIN_N, &b_cluster_HCALIN_N);
   fChain->SetBranchAddress("cluster_HCALIN_E", cluster_HCALIN_E, &b_cluster_HCALIN_E);
   fChain->SetBranchAddress("cluster_HCALIN_Eta", cluster_HCALIN_Eta, &b_cluster_HCALIN_Eta);
   fChain->SetBranchAddress("cluster_HCALIN_Phi", cluster_HCALIN_Phi, &b_cluster_HCALIN_Phi);
   fChain->SetBranchAddress("cluster_HCALIN_NTower", cluster_HCALIN_NTower, &b_cluster_HCALIN_NTower);
   fChain->SetBranchAddress("cluster_HCALIN_trueID", cluster_HCALIN_trueID, &b_cluster_HCALIN_trueID);
   fChain->SetBranchAddress("tower_HCALOUT_N", &tower_HCALOUT_N, &b_tower_HCALOUT_N);
   fChain->SetBranchAddress("tower_HCALOUT_E", tower_HCALOUT_E, &b_tower_HCALOUT_E);
   fChain->SetBranchAddress("tower_HCALOUT_iEta", tower_HCALOUT_iEta, &b_tower_HCALOUT_iEta);
   fChain->SetBranchAddress("tower_HCALOUT_iPhi", tower_HCALOUT_iPhi, &b_tower_HCALOUT_iPhi);
   fChain->SetBranchAddress("tower_HCALOUT_trueID", tower_HCALOUT_trueID, &b_tower_HCALOUT_trueID);
   fChain->SetBranchAddress("cluster_HCALOUT_N", &cluster_HCALOUT_N, &b_cluster_HCALOUT_N);
   fChain->SetBranchAddress("cluster_HCALOUT_E", cluster_HCALOUT_E, &b_cluster_HCALOUT_E);
   fChain->SetBranchAddress("cluster_HCALOUT_Eta", cluster_HCALOUT_Eta, &b_cluster_HCALOUT_Eta);
   fChain->SetBranchAddress("cluster_HCALOUT_Phi", cluster_HCALOUT_Phi, &b_cluster_HCALOUT_Phi);
   fChain->SetBranchAddress("cluster_HCALOUT_NTower", cluster_HCALOUT_NTower, &b_cluster_HCALOUT_NTower);
   fChain->SetBranchAddress("cluster_HCALOUT_trueID", cluster_HCALOUT_trueID, &b_cluster_HCALOUT_trueID);
   fChain->SetBranchAddress("tower_FEMC_N", &tower_FEMC_N, &b_tower_FEMC_N);
   fChain->SetBranchAddress("tower_FEMC_E", tower_FEMC_E, &b_tower_FEMC_E);
   fChain->SetBranchAddress("tower_FEMC_iEta", tower_FEMC_iEta, &b_tower_FEMC_iEta);
   fChain->SetBranchAddress("tower_FEMC_iPhi", tower_FEMC_iPhi, &b_tower_FEMC_iPhi);
   fChain->SetBranchAddress("tower_FEMC_trueID", tower_FEMC_trueID, &b_tower_FEMC_trueID);
   fChain->SetBranchAddress("cluster_FEMC_N", &cluster_FEMC_N, &b_cluster_FEMC_N);
   fChain->SetBranchAddress("cluster_FEMC_E", cluster_FEMC_E, &b_cluster_FEMC_E);
   fChain->SetBranchAddress("cluster_FEMC_Eta", cluster_FEMC_Eta, &b_cluster_FEMC_Eta);
   fChain->SetBranchAddress("cluster_FEMC_Phi", cluster_FEMC_Phi, &b_cluster_FEMC_Phi);
   fChain->SetBranchAddress("cluster_FEMC_NTower", cluster_FEMC_NTower, &b_cluster_FEMC_NTower);
   fChain->SetBranchAddress("cluster_FEMC_trueID", cluster_FEMC_trueID, &b_cluster_FEMC_trueID);
   fChain->SetBranchAddress("tower_EEMC_N", &tower_EEMC_N, &b_tower_EEMC_N);
   fChain->SetBranchAddress("tower_EEMC_E", tower_EEMC_E, &b_tower_EEMC_E);
   fChain->SetBranchAddress("tower_EEMC_iEta", tower_EEMC_iEta, &b_tower_EEMC_iEta);
   fChain->SetBranchAddress("tower_EEMC_iPhi", tower_EEMC_iPhi, &b_tower_EEMC_iPhi);
   fChain->SetBranchAddress("tower_EEMC_trueID", tower_EEMC_trueID, &b_tower_EEMC_trueID);
   fChain->SetBranchAddress("cluster_EEMC_N", &cluster_EEMC_N, &b_cluster_EEMC_N);
   fChain->SetBranchAddress("cluster_EEMC_E", cluster_EEMC_E, &b_cluster_EEMC_E);
   fChain->SetBranchAddress("cluster_EEMC_Eta", cluster_EEMC_Eta, &b_cluster_EEMC_Eta);
   fChain->SetBranchAddress("cluster_EEMC_Phi", cluster_EEMC_Phi, &b_cluster_EEMC_Phi);
   fChain->SetBranchAddress("cluster_EEMC_NTower", cluster_EEMC_NTower, &b_cluster_EEMC_NTower);
   fChain->SetBranchAddress("cluster_EEMC_trueID", cluster_EEMC_trueID, &b_cluster_EEMC_trueID);
   fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
   fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
   fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
   fChain->SetBranchAddress("vertex_NCont", &vertex_NCont, &b_vertex_NCont);
   fChain->SetBranchAddress("vertex_true_x", &vertex_true_x, &b_vertex_true_x);
   fChain->SetBranchAddress("vertex_true_y", &vertex_true_y, &b_vertex_true_y);
   fChain->SetBranchAddress("vertex_true_z", &vertex_true_z, &b_vertex_true_z);
   fChain->SetBranchAddress("nMCPart", &nMCPart, &b_nMCPart);
   fChain->SetBranchAddress("mcpart_ID", mcpart_ID, &b_mcpart_ID);
   fChain->SetBranchAddress("mcpart_ID_parent", mcpart_ID_parent, &b_mcpart_ID_parent);
   fChain->SetBranchAddress("mcpart_PDG", mcpart_PDG, &b_mcpart_PDG);
   fChain->SetBranchAddress("mcpart_E", mcpart_E, &b_mcpart_E);
   fChain->SetBranchAddress("mcpart_px", mcpart_px, &b_mcpart_px);
   fChain->SetBranchAddress("mcpart_py", mcpart_py, &b_mcpart_py);
   fChain->SetBranchAddress("mcpart_pz", mcpart_pz, &b_mcpart_pz);
   fChain->SetBranchAddress("mcpart_BCID", mcpart_BCID, &b_mcpart_BCID);
   fChain->SetBranchAddress("nHepmcp", &nHepmcp, &b_nHepmcp);
   fChain->SetBranchAddress("hepmcp_procid", &hepmcp_procid, &b_hepmcp_procid);
   fChain->SetBranchAddress("hepmcp_x1", &hepmcp_x1, &b_hepmcp_x1);
   fChain->SetBranchAddress("hepmcp_x2", &hepmcp_x2, &b_hepmcp_x2);
   fChain->SetBranchAddress("hepmcp_Q2", &hepmcp_Q2, &b_hepmcp_Q2);
   fChain->SetBranchAddress("hepmcp_status", hepmcp_status, &b_hepmcp_status);
   fChain->SetBranchAddress("hepmcp_PDG", hepmcp_PDG, &b_hepmcp_PDG);
   fChain->SetBranchAddress("hepmcp_E", hepmcp_E, &b_hepmcp_E);
   fChain->SetBranchAddress("hepmcp_px", hepmcp_px, &b_hepmcp_px);
   fChain->SetBranchAddress("hepmcp_py", hepmcp_py, &b_hepmcp_py);
   fChain->SetBranchAddress("hepmcp_pz", hepmcp_pz, &b_hepmcp_pz);
   fChain->SetBranchAddress("hepmcp_BCID", hepmcp_BCID, &b_hepmcp_BCID);
   fChain->SetBranchAddress("hepmcp_m1", hepmcp_m1, &b_hepmcp_m1);
   fChain->SetBranchAddress("hepmcp_m2", hepmcp_m2, &b_hepmcp_m2);
   Notify();
}

Bool_t beg_et::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void beg_et::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t beg_et::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef beg_et_cxx
