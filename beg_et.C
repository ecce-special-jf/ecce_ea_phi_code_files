#define beg_et_cxx
#include "beg_et.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void beg_et::Loop()
{
//   In a ROOT session, you can do:
//      root> .L beg_et.C
//      root> beg_et t
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


   TH2D* hrp_0_xy = new TH2D("hrp_0_xy","RP 0 YX",200,-20,20,200,-10,10);
   TH2D* hrp_1_xy = new TH2D("hrp_1_xy","RP 1 YX",200,-20,20,200,-10,10);
   TH2D* hrp2_0_xy = new TH2D("hrp2_0_xy","RP 0 YX",200,-5,5,200,-3,3);
   TH2D* hrp2_1_xy = new TH2D("hrp2_1_xy","RP 1 YX",200,-5,5,200,-3,3);
   
   TH2D* hrp_0_xy_cut = new TH2D("hrp_0_xy_cut","RP 0 YX",200,-20,20,200,-10,10);
   TH2D* hrp_1_xy_cut = new TH2D("hrp_1_xy_cut","RP 1 YX",200,-20,20,200,-10,10);
   TH2D* hrp2_0_xy_cut = new TH2D("hrp2_0_xy_cut","RP 0 YX",200,-5,5,200,-3,3);
   TH2D* hrp2_1_xy_cut = new TH2D("hrp2_1_xy_cut","RP 1 YX",200,-5,5,200,-3,3);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      // find phi to KK truth events/info
      int p1i =-1,p2i=-1, k1i=-1, k2i=-1, e1i =-1, e2i =-1;
      //      cout << "ientry " << ientry  << endl;
      for (int i = 0; i < nHepmcp; i++)
	{
	  //	  cout << "  " << hepmcp_PDG[i];
	  if (hepmcp_PDG[i] == -321)
	    k1i = i;
	  if (hepmcp_PDG[i] == 321)
	    k2i = i;
	  if (hepmcp_PDG[i] == 333)
	    {
	      if (p1i < 0)
		p1i = i;
	      else
		p2i = i;
	    }
	  
	  //cout << "hep pdg: " <<  hepmcp_PDG[i] << " " << e1i << " " << i << endl;
	  if (hepmcp_PDG[i] == 11 && hepmcp_status[i] < 5 && i != 0)
	    {
	      if (e1i < 0)
		e1i = i;
	      else
		e2i = i;
	    }

	  //	  cout << "\t2hep pdg: " <<  hepmcp_PDG[i] <<  " " << e1i << " " << e2i << endl;
		  

	  
	}
      //      cout << endl;
      //      cout << k1i << " " << k2i << " " << p1i << "  " << p2i << endl;
      
      TLorentzVector tk1,tp1,tp2,tk2,te1,te2,te0;
      if (k1i > 0 && k2i > 0 && p1i > 0 && p2i > 0 && e1i > 0 && e2i > 0)
	{
	  //	  h_px0->Fill(hepmcp);
	  h_pxb->Fill(hepmcp_E[e1i]/(hepmcp_E[e2i]));	      
	}

      //      continue;


      //      continue;
      te0.SetPxPyPzE(18,-18,0,0);
      //      te2.SetPxPyPzE(hepmcp_px[e1i],hepmcp_py[e1i],hepmcp_pz[e1i],hepmcp_E[e1i]);
      te1.SetPxPyPzE(hepmcp_px[e1i],hepmcp_py[e1i],hepmcp_pz[e1i],hepmcp_E[e1i]);
      te2.SetPxPyPzE(hepmcp_px[e2i],hepmcp_py[e2i],hepmcp_pz[e2i],hepmcp_E[e2i]);
      //      te1.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      tp2.SetPxPyPzE(hepmcp_px[p2i],hepmcp_py[p2i],hepmcp_pz[p2i],hepmcp_E[p2i]);  

      //      cout << tp2.M() << " " <<te1.M() << " " <<  te2.M() << endl; 

      h_ttrue_bc->Fill(    (tp2.Px()+te2.Px())*(tp2.Px()+te2.Px()) +
			(tp2.Py()+te2.Py())*(tp2.Py()+te2.Py()) );            
      
      //  
      // mcpart_PDG == 11 && mcpart_ID_parent == 0 selects scattered electron
      // with beam effects

      map<int,int> mcpMap;
      map<int,int> mcpMap_i;
      //      map<int,int> mcpMap;
      for (int jj=0; jj < nMCPart; jj++)
	{
	  if (mcpart_ID_parent[jj]==0)
	    {
	      mcpMap[mcpart_ID[jj]]=mcpart_PDG[jj];
	      mcpMap_i[mcpart_ID[jj]]=jj;
	    }
	}


      if (cluster_EEMC_N != 1) 
	continue;
      
      int chtrk0= -1, chtrk1=-1, chtrk2=-1;      
      int tchtrk0= -1, tchtrk1=-1, tchtrk2=-1;
      int recotrk0= -1, recotrk1=-1, recotrk2=-1;
      int pidtrk0= -1, pidtrk1=-1, pidtrk2=-1;
      for (int k =0; k< nTracks; k++)
	{
	  if (tracks_source[k] > 0) continue;
	  if (mcpMap[tracks_trueID[k]] == 321)
	    {
	      if (chtrk0 < 0)
		chtrk0 = k;
	      else 
		chtrk0 = -1*chtrk0+-k*1000;
	      // more than one + kaon or -kaon cuts event	      
	    }

	  if (mcpMap[tracks_trueID[k]] == -321)
	    {
	      if (chtrk1 < 0)
		chtrk1 = k;
	      else 
		chtrk1 = -1*chtrk1+-k*1000;
	      // more than one + kaon or -kaon cuts event	      
	    }

	  if (mcpMap[tracks_trueID[k]] == 11)
	    {
	      if (chtrk2 < 0)
		chtrk2 = k;
	      else 
		cout << "found 2 electrons" << endl;
	      // more than one + kaon or -kaon cuts event	      
	    }

	  // now reco
	  

	  if (tracks_charge[k] == -1 ) // electron or kaon will encounter kaon first
	    {
	      if (tracks_trueID[k] == cluster_EEMC_trueID[0] && mcpMap[cluster_EEMC_trueID[0]] == 11)
		//cheat on matching
		{
		  recotrk2 = k;
		  pidtrk2 = k;
		  tchtrk2 = mcpMap_i[tracks_trueID[k]];
		}
	      else
		{
		  recotrk0 = k;	     						   
		  if (mcpMap[tracks_trueID[k]] == -321)
		    tchtrk0 = mcpMap_i[tracks_trueID[k]];
		  if (track_kaon_LL[k] - track_pion_LL[k]  > 1)
		    {
		      pidtrk0 = k;
		    }
		}
	    }

	  if (tracks_charge[k] == 1 ) // electron or kaon will encounter kaon first
	    {
	      recotrk1 = k;	     						   
	      if (mcpMap[tracks_trueID[k]] == 321)
		tchtrk1 = mcpMap_i[tracks_trueID[k]];
	      
	      if (track_kaon_LL[k] - track_pion_LL[k]  > 1)
		{
		  pidtrk1 = k;
		}
	    }
	  	  
	}
      
      if (recotrk0 < 0 ||  recotrk1 < 0)
	continue;


      //      cout << "\t" << ientry << " " << recotrk0 << " " << recotrk1 << " " << pidtrk0 << " " << pidtrk1 << endl; 

      float all_px = tracks_px[0]+tracks_px[1]+tracks_px[2];

      float  p0mag=sqrt(tracks_px[recotrk0]*tracks_px[recotrk0]+tracks_py[recotrk0]*tracks_py[recotrk0]+tracks_pz[recotrk0]*tracks_pz[recotrk0]);    
      float  E0mag=sqrt(p0mag*p0mag+0.495*0.495);
      float  p1mag =sqrt(tracks_px[recotrk1]*tracks_px[recotrk1]+tracks_py[recotrk1]*tracks_py[recotrk1]+tracks_pz[recotrk1]*tracks_pz[recotrk1]);
      float  E1mag=sqrt(p1mag*p1mag+0.495*0.495);
      float  p2mag=sqrt(tracks_px[recotrk2]*tracks_px[recotrk2]+tracks_py[recotrk2]*tracks_py[recotrk2]+tracks_pz[recotrk2]*tracks_pz[recotrk2]);
      
      float  phi_e=E0mag + E1mag;
      float  phi_px=tracks_px[recotrk0]+tracks_px[recotrk1];
      float  phi_py=tracks_py[recotrk0]+tracks_py[recotrk1];
      float  phi_pz=tracks_pz[recotrk0]+tracks_pz[recotrk1];
      
      TLorentzVector phiVb;
      phiVb.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      h_mass_bc->Fill(phiVb.M()); // shouldn't matter if it's before beam correction
     
      if (pidtrk0 < 0 && pidtrk1 < 0)
	continue;
     
      //      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      eVb.SetPxPyPzE(tracks_px[recotrk2],tracks_py[recotrk2],tracks_pz[recotrk2],p2mag);

      TLorentzVector phiVb0;
      phiVb0.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      TLorentzVector eVb0;
      eVb0.SetPxPyPzE(hepmcp_px[e1i],hepmcp_py[e1i],hepmcp_pz[e1i],hepmcp_E[e1i]);
      //      eVb0.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
      
      TLorentzRotation lrTrans;
      //      lrTrans.RotateY(+12.5e-3).Boost(sin(+12.5e-3),0,0);
      lrTrans.RotateY(-12.5e-3).Boost(sin(+12.5e-3),0,0);

      //      TLorentzVector phiVb;
      //phiVb.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);


      phiVb *= lrTrans;
      eVb *= lrTrans;
      
      if (!(tchtrk0 < 0||  tchtrk1 < 0 || tchtrk2 < 0 ) )
	{
	  TLorentzVector mphiVb;
	  mphiVb.SetPxPyPzE(mcpart_px[tchtrk0]+mcpart_px[tchtrk1],
			    mcpart_py[tchtrk0]+mcpart_py[tchtrk1],
			    mcpart_pz[tchtrk0]+mcpart_pz[tchtrk1],
			    mcpart_E[tchtrk0]+mcpart_E[tchtrk1]);
	  
	  TLorentzVector meVb;
	  meVb.SetPxPyPzE(mcpart_px[tchtrk2],
			  mcpart_py[tchtrk2],
			  mcpart_pz[tchtrk2],
			  mcpart_E[tchtrk2]);
	  
	  mphiVb *= lrTrans;
	  meVb *= lrTrans;
	  
	  h_t_trans2->Fill(    (mphiVb.Px()+meVb.Px())*(mphiVb.Px()+meVb.Px()) +
			       (mphiVb.Py()+meVb.Py())*(mphiVb.Py()+meVb.Py()) );            
	}
      
      
      if (fabs(phiVb.M()-1.02) > 0.04) continue;  

      h_mass->Fill(phiVb.M());
				      
      h_px0->Fill(all_px); 
      h_pxb->Fill(phiVb.Px() + eVb.Px() ); 
      
      //      h_ttrue->Fill(-hepmcp_E[2]*hepmcp_E[2]+hepmcp_px[2]*hepmcp_px[2]+hepmcp_py[2]*hepmcp_py[2]+hepmcp_pz[2]*hepmcp_pz[2]);      
      
      h_ttrue->Fill(    (tp2.Px()+te2.Px())*(tp2.Px()+te2.Px()) +
			(tp2.Py()+te2.Py())*(tp2.Py()+te2.Py()) );            
      
      
      float ttr =   (phiVb.Px()+eVb.Px())*(phiVb.Px()+eVb.Px()) +
	(phiVb.Py()+eVb.Py())*(phiVb.Py()+eVb.Py()) ;            
      
      h_t_trans_rot->Fill(  ttr );
      
      h_t_trans->Fill(    (phiVb.Px()+eVb0.Px())*(phiVb.Px()+eVb0.Px()) +
			  (phiVb.Py()+eVb0.Py())*(phiVb.Py()+eVb0.Py()) );            


      //FF analysis
      
      // simple FF analysis
      // ID 70 ZDC 0 
      // ID 71 RP 1 
      // ID 72 B0 2
      // ID 73 RP2 3
      // ID 74 OMD 4


   // definitions of FF detectors
   int nplanes[5] = {1,2,4,2,2}; // ZDC, RP, B0, RP2, OMD
   TString names[5] = {"ZDC","RP","B0","RP2","OMD"};

   
   float planez[5][4] = {
			 {3350,0,0,0},
			 {2600,2800,0,0},
			 {561,585,609,633},
			 {4300,4450,0,0},
			 {3250,3450,0,0}
   };

   float planex[5][4] = {
			 {120,0,0,0},
			 {75.6,78.15,0,0},
			 {0,0,0,0},
			 {101.94+1,106.94+0.4,0,0},
			 {0,0,0,0}
   };

   float sigmax[5][4] = {
			 {0,0,0,0},
			 {4.86,4.34,0,0},
			 {0,0,0,0},
			 {.393,.288,0,0},
			 {0,0,0,0}
   };
   float sigmay[5][4] = {
			 {0,0,0,0},
			 {.709,.637,0,0},
			 {0,0,0,0},
			 {0.2,0.2,0,0}, // way overestimate
			 {0,0,0,0}
   };

   int nhits[5];
   int nhitsCut[5];
   float ehits[5];
 

      for (int i =0;i<5;i++)
	{
	  ehits[i] = 0;
	  nhits[i] = 0;
	  nhitsCut[i] = 0;
	}

      float ezdc_phot = 0;
      float ezdc_neut = 0;
      float ezdc_bg = 0;
      for (int i = 0;i<nHits;i++)
	{
	  if (hits_layerID[i]<70 || hits_layerID[i]>74) continue; // skip non-FF
	  int layer = hits_layerID[i]-70; // remap 70-74 to 0-4
	  int trueid = -1;
	  for (int imc = 0;imc<nMCPart;imc++)
	    {
	      if (mcpart_ID[imc]==hits_trueID[i] && hits_trueID[i]>0) // filter out BG
		{
		  trueid = mcpart_PDG[imc];
		  break;
		}
	    }

	  switch (trueid)
	    {
	    case 22:
	      ezdc_phot += hits_edep[i];
	      break;
	    case 2112:
	      ezdc_neut += hits_edep[i];
	      break;
	    default:
	      ezdc_bg += hits_edep[i];
	    }
	  
	  float edep = hits_edep[i]; 
	  ehits[layer] += edep;
	  nhits[layer]++;
	  int subdet = -1;
	  for (int j = 0;j<nplanes[layer];j++) // find the subdeector and check if it makes r cut
	    {
	      if ( hits_z[i]>planez[layer][j]-10 && hits_z[i]<planez[layer][j]+10 )
		{
		  if (layer==1 && j==0) hrp_0_xy->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
		  if (layer==1 && j==1) hrp_1_xy->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
		  if (layer==3 && j==0) hrp2_0_xy->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
		  if (layer==3 && j==1) hrp2_1_xy->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
		  
		  float r = 0;
		  if (sigmax[layer][j]!=0)
		    {
		      r = sqrt( pow((hits_x[i]-planex[layer][j])/sigmax[layer][j],2) +
				      pow(hits_y[i]/sigmay[layer][j],2) );
		      if (r>1)
			{
			  nhitsCut[layer]++;
			  if (layer==1 && j==0) hrp_0_xy_cut->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
			  if (layer==1 && j==1) hrp_1_xy_cut->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
			  if (layer==3 && j==0) hrp2_0_xy_cut->Fill(hits_x[i]-planex[layer][j],hits_y[i]);
			  if (layer==3 && j==1) hrp2_1_xy_cut->Fill(hits_x[i]-planex[layer][j],hits_y[i]);			  
			}
		    }

		  //cout << "Found hit in layer " << names[layer] << " subdet " << j << " z=" << hits_z[i] << " x=" << hits_x[i]-planex[layer][j] << " y=" << hits_y[i] << " r= " << r << endl;
		  break;
		}
	    }
	}


      if (nhits[0] > 0 || nhits[1] > 0 || nhits[2] > 0  || nhits[4] > 0)
	continue;


      h_ttrueff->Fill(    (tp2.Px()+te2.Px())*(tp2.Px()+te2.Px()) +
			(tp2.Py()+te2.Py())*(tp2.Py()+te2.Py()) );            
      
      
      h_t_trans_rotff->Fill(  ttr );
      
      h_t_transff->Fill(    (phiVb.Px()+eVb0.Px())*(phiVb.Px()+eVb0.Px()) +
			  (phiVb.Py()+eVb0.Py())*(phiVb.Py()+eVb0.Py()) );            

      
   }



}
