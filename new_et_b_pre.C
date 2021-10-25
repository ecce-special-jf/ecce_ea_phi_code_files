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
      if (!(tracks_trueID[0]==1 && tracks_trueID[1] == 2 && tracks_trueID[2] == 4)) continue; 


      float all_px = tracks_px[0]+tracks_px[1]+tracks_px[2];

      float  p0mag=sqrt(tracks_px[0]*tracks_px[0]+tracks_py[0]*tracks_py[0]+tracks_pz[0]*tracks_pz[0]);    
      float  E0mag=sqrt(p0mag*p0mag+0.495*0.495);

      float  p1mag =sqrt(tracks_px[1]*tracks_px[1]+tracks_py[1]*tracks_py[1]+tracks_pz[1]*tracks_pz[1]);
      
      //      if (p0mag < 1  || p1mag < 1)  continue;

      float emc_epx = 1.11*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_px[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]));

      float emc_epy = 1.09*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_py[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]));

      float emc_epz = 1.09*cluster_EEMC_E[0]*tanh(cluster_EEMC_Eta[0]);

      float  E1mag=sqrt(p1mag*p1mag+0.495*0.495);
      float  p2mag=sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]+tracks_pz[2]*tracks_pz[2]);

      //      if (p2mag > 10) continue;

      float  phi_e=E0mag + E1mag;
      float  phi_px=tracks_px[0]+tracks_px[1];
      float  phi_py=tracks_py[0]+tracks_py[1];
      float  phi_pz=tracks_pz[0]+tracks_pz[1];


      TLorentzVector k0Vb;
      k0Vb.SetPxPyPzE(tracks_px[0],tracks_py[0], tracks_pz[0],E0mag);

      TLorentzVector k1Vb;
      k1Vb.SetPxPyPzE(tracks_px[1],tracks_py[1], tracks_pz[1],E1mag);

      //      if (fabs(k1Vb.Eta()) < 1. || fabs(k0Vb.Eta()) < 1.) continue;


      TLorentzVector phiVb;
      phiVb.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      //      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);


      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
      float av_epx = (tracks_px[2]+ emc_epx)/2.0;
      float av_epy = (tracks_py[2]+ emc_epy)/2.0;
      float av_epz = (tracks_pz[2]+ emc_epz)/2.0;
      float av_eE = (p2mag+1.078*cluster_EEMC_E[0])/2.0;

      if(false)
      if (fabs(tracks_px[2]/emc_epx -1) > 0.03 || 
      	  fabs(tracks_py[2]/emc_epy -1) > 0.03 ||  
	  fabs(tracks_pz[2]/emc_epz -1) > 0.03 ||
	fabs(p2mag/1.078/cluster_EEMC_E[0]-1) > 0.03)
      	continue;

      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      eVb.SetPxPyPzE(av_epx,av_epy,av_epz,av_eE);

      TLorentzVector phiVb0;
      phiVb0.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      TLorentzVector eVb0;
      eVb0.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb0.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      //      if (eVb0.Eta() > -3) continue;
      
      TLorentzRotation lrTrans;
      lrTrans.RotateY(+12.5e-3).Boost(sin(+12.5e-3),0,0);
      //lrTrans.RotateY(+12.5e-3);
      
      phiVb *= lrTrans;
      eVb *= lrTrans;

      phiVb0 *=lrTrans;
      //phiVb0 *=lrTrans;
				      
      h_px0->Fill(all_px); 
      h_pxb->Fill(phiVb.Px() + eVb.Px() ); 
      
      //      h_ttrue->Fill(-hepmcp_E[2]*hepmcp_E[2]+hepmcp_px[2]*hepmcp_px[2]+hepmcp_py[2]*hepmcp_py[2]+hepmcp_pz[2]*hepmcp_pz[2]);      
      h_ttrue->Fill(hepmcp_px[2]*hepmcp_px[2]+hepmcp_py[2]*hepmcp_py[2]);      
      float ttr =   (phiVb.Px()+eVb.Px())*(phiVb.Px()+eVb.Px()) +
	(phiVb.Py()+eVb.Py())*(phiVb.Py()+eVb.Py()) ;            
     
      h_t_trans_rot->Fill(  ttr );

      h_t_trans->Fill(    (phiVb.Px()+eVb0.Px())*(phiVb.Px()+eVb0.Px()) +
			  (phiVb.Py()+eVb0.Py())*(phiVb.Py()+eVb0.Py()) );            


   }
}
