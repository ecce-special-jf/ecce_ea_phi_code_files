#define new_et_cxx
#include "new_et.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>


void new_et::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L new_et.C
//      Root > new_et t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //      cout << tracks_px[0] << endl;
      //      if (fabs(tracks_px[2]-hepmcp_px[1])/hepmcp_px[1] > 0.2)
      //	continue;

      //      if (nTracks != 3) continue;
      if (cluster_EEMC_N != 1) continue;
      if (nTracks < 3) continue;
      else if (tracks_trueID[2] != 4 ) continue;
      if (fabs(vertex_z) >3) continue;
      
      

      float all_px = tracks_px[0]+tracks_px[1]+tracks_px[2];

      float  p0mag=sqrt(tracks_px[0]*tracks_px[0]+tracks_py[0]*tracks_py[0]+tracks_pz[0]*tracks_pz[0]);    
      float  E0mag=sqrt(p0mag*p0mag+0.495*0.495);

      float  p1mag =sqrt(tracks_px[1]*tracks_px[1]+tracks_py[1]*tracks_py[1]+tracks_pz[1]*tracks_pz[1]);
      
      //      if (p0mag < 1  || p1mag < 1)  continue;

      TLorentzVector eVb0;
      eVb0.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);

      TLorentzVector eVb02a;
      eVb02a.SetPxPyPzE(mcpart_px[nMCPart- 1],mcpart_py[nMCPart- 1],mcpart_pz[nMCPart- 1],mcpart_E[nMCPart- 1]);
      //      eVb0.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      float  p2mag=sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]+tracks_pz[2]*tracks_pz[2]);
      
      TLorentzVector eVb2a;
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      eVb2a.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      //smearing truth by expected       
      float true_e = mcpart_E[nMCPart -1];
	//float true_e = hepmcp_E[1];
      float strue = sqrt(true_e);
      // PbWO   2.5%/sqrt(E) + 1%
      // reuse PHENIX   8%/sqrt(E) + 2%
      //      float res1 = 0.025, res2 = 0.01;
      float res1 = 0.08, res2 = 0.02;


      //reset to more reasonable smeared value;
      cluster_EEMC_E[0] = true_e 
	+ strue*res1*f_egausres1->GetRandom() 
	+ res2*true_e*f_egausres1->GetRandom();
      
      if (jentry < 1000)
	cout << "true e" << true_e <<  " " << cluster_EEMC_E[0] <<endl;

      //for debug 
      //cluster_EEMC_Eta[0] = //(eVb2a.Eta() + cluster_EEMC_Eta[2])/2;

	cluster_EEMC_Eta[0] = eVb2a.Eta();
      
      //      float pxcal = 0.887, pycal = 0.887, pzcal = 1.07;
      float pxcal = 1., pycal = 1., pzcal = 1.;
      // tracks_px[2] = hepmcp_px[1];
      // tracks_py[2] = hepmcp_py[1];
      // tracks_pz[2] = hepmcp_pz[1];

      //      tracks_px[2] = hepmcp_E[1]*cos(cluster_EEMC_Phi[0]);
      //      tracks_py[2] = hepmcp_E[1]*sin(cluster_EEMC_Phi[0]);
      //      tracks_pz[2] = cluster_EEMC_[0];

      float cosphi_e = tracks_px[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]);
      float sinphi_e = tracks_py[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]);
      //float cosphi_e = cos(cluster_EEMC_Phi[0]);
      //float sinphi_e = sin(cluster_EEMC_Phi[0]);

      //      float cosphi_e = cos(eVb02a.Phi());
      //float sinphi_e = sin(eVb02a.Phi());


      float emc_epx = pxcal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*cosphi_e;
      float emc_epy = pycal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*sinphi_e;

	/*
      float emc_epx = pxcal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_px[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]));

      float emc_epy = pycal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_py[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]));

	*/

      float emc_epz = pzcal*cluster_EEMC_E[0]*tanh(cluster_EEMC_Eta[0]);

      float  E1mag=sqrt(p1mag*p1mag+0.495*0.495);

      //      if (p2mag > 10) continue;

      float  phi_e=E0mag + E1mag;
      float  phi_px=tracks_px[0]+tracks_px[1];
      float  phi_py=tracks_py[0]+tracks_py[1];
      float  phi_pz=tracks_pz[0]+tracks_pz[1];


      TLorentzVector k0Vb;
      TLorentzVector k0Vbt;
      k0Vb.SetPxPyPzE(tracks_px[0],tracks_py[0], tracks_pz[0],E0mag);
      k0Vbt.SetPxPyPzE(hepmcp_px[7], hepmcp_py[7], hepmcp_pz[7], hepmcp_E[7]);


      //      if (fabs(tracks_px[1]/)

      TLorentzVector k1Vb;
      k1Vb.SetPxPyPzE(tracks_px[1],tracks_py[1], tracks_pz[1],E1mag);

      if (fabs(k1Vb.Eta()) > 0.95 || fabs(k0Vb.Eta()) > 0.95) continue;
      if (fabs(k1Vb.Eta()) < 0.17) continue;
      //      if (fabs(tracks_py[1]/tracks_px[1]) < 1.5 || fabs(tracks_py[0]/tracks_px[0]) < 1.5) continue;


      TLorentzVector phiVb;
      phiVb.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      //      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);      

      //      eVb.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
            
      tracks_px[2] = 1.02*tracks_px[2];
      tracks_py[2] = 1.02*tracks_py[2];
      tracks_pz[2] = 1.02*tracks_pz[2];

      float av_epx = (tracks_px[2]+ emc_epx)/2.0;
      float av_epy = (tracks_py[2]+ emc_epy)/2.0;
      float av_epz = (tracks_pz[2]+ emc_epz)/2.0;
      float av_eE = (p2mag+cluster_EEMC_E[0])/2.0;

      if (true)
      if (fabs(tracks_px[2]/emc_epx -1) > 0.08 || 
	  fabs(tracks_py[2]/emc_epy -1) > 0.08 || 
	  fabs(tracks_pz[2]/emc_epz -1) > 0.08)
	//||
	// fabs(p2mag/cluster_EEMC_E[0]-1) > 0.02)
	continue;



      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      eVb.SetPxPyPzE(av_epx,av_epy,av_epz,av_eE);
      //eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      TLorentzVector eVb2;
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      eVb2.SetPxPyPzE(av_epx,av_epy,av_epz,av_eE);

      TLorentzVector phiVb0;
      phiVb0.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      //      TLorentzVector eVb0;
      //      eVb0.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb0.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      if (eVb.Eta() < -3.2) continue;
      
      TLorentzRotation lrTrans;
      lrTrans.RotateY(+12.5e-3).Boost(sin(+12.5e-3),0,0);
      //lrTrans.RotateY(+12.5e-3);
      
      phiVb *= lrTrans;
      eVb *= lrTrans;

      phiVb0 *=lrTrans;
      //phiVb0 *=lrTrans;
				      
      k0Vb *= lrTrans;
      
      // if (!(fabs(k0Vbt.Eta()) > 1.0 && 
      // 	    fabs(k0Vbt.Eta()) < 1.5) 
      // 	  && !(fabs(k0Vbt.Eta()) < 0.2))
      //h2px->Fill(k0Vb.Eta(), k0Vb.Px()/hepmcp_px[7]);
      h2px->Fill(eVb.Eta(), eVb.Px()/hepmcp_px[1]);

      h_px0->Fill(eVb.Px()); 
      h_pxb->Fill(phiVb.Px() + eVb.Px() ); 
      
      //      h_ttrue->Fill(-hepmcp_E[2]*hepmcp_E[2]+hepmcp_px[2]*hepmcp_px[2]+hepmcp_py[2]*hepmcp_py[2]+hepmcp_pz[2]*hepmcp_pz[2]);      
      h_ttrue->Fill(hepmcp_px[2]*hepmcp_px[2]+hepmcp_py[2]*hepmcp_py[2]);      

      //      float ttr2 =   (phiVb.Px()+eVb2.Px())*(phiVb.Px()+eVb2.Px()) +
      //	(phiVb.Py()+eVb2.Py())*(phiVb.Py()+eVb2.Py()) ;            


      //h_ttrue->Fill(ttr2);


      float ttr =   (phiVb.Px()+eVb.Px())*(phiVb.Px()+eVb.Px()) +
	(phiVb.Py()+eVb.Py())*(phiVb.Py()+eVb.Py()) ;            
     
      h_t_trans_rot->Fill(  ttr );

      h_t_trans->Fill(    (phiVb.Px()+eVb0.Px())*(phiVb.Px()+eVb0.Px()) +
			  (phiVb.Py()+eVb0.Py())*(phiVb.Py()+eVb0.Py()) );            


   }
}
