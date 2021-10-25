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

      //for (int i = 0; 
      

      
      // find phi to KK truth events/info
      int p1i =-1,p2i=-1, k1i=-1, k2i=-1;
      cout << "ientry " << ientry ;
      for (int i = 0; i < nHepmcp; i++)
	{
	  cout << "  " << hepmcp_PDG[i];
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
	  
	}
      cout << endl;
      cout << k1i << " " << k2i << " " << p1i << "  " << p2i << endl;
      
      TLorentzVector tk1,tp1,tp2,tk2;
      if (k1i > 0 && k2i > 0 && p1i > 0 && p2i > 0)
	{
	  h_px0->Fill(hepmcp_E[p1i]/(hepmcp_E[k1i]+hepmcp_E[k2i]));
	  h_pxb->Fill(hepmcp_E[p2i]/(hepmcp_E[k1i]+hepmcp_E[k2i]));	      
	}

      continue;

      float all_px = tracks_px[0]+tracks_px[1]+tracks_px[2];

      float  p0mag=sqrt(tracks_px[0]*tracks_px[0]+tracks_py[0]*tracks_py[0]+tracks_pz[0]*tracks_pz[0]);    
      float  E0mag=sqrt(p0mag*p0mag+0.495*0.495);
      float  p1mag =sqrt(tracks_px[1]*tracks_px[1]+tracks_py[1]*tracks_py[1]+tracks_pz[1]*tracks_pz[1]);
      float  E1mag=sqrt(p1mag*p1mag+0.495*0.495);
      float  p2mag=sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]+tracks_pz[2]*tracks_pz[2]);

      float  phi_e=E0mag + E1mag;
      float  phi_px=tracks_px[0]+tracks_px[1];
      float  phi_py=tracks_py[0]+tracks_py[1];
      float  phi_pz=tracks_pz[0]+tracks_pz[1];


      TLorentzVector phiVb;
      phiVb.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      //      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      TLorentzVector phiVb0;
      phiVb0.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      TLorentzVector eVb0;
      eVb0.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb0.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
      
      TLorentzRotation lrTrans;
      lrTrans.RotateY(+12.5e-3).Boost(sin(+12.5e-3),0,0);
      
      phiVb *= lrTrans;
      eVb *= lrTrans;
				      
      h_px0->Fill(all_px); 
      h_pxb->Fill(phiVb.Px() + eVb.Px() ); 
      
      h_ttrue->Fill(-hepmcp_E[2]*hepmcp_E[2]+hepmcp_px[2]*hepmcp_px[2]+hepmcp_py[2]*hepmcp_py[2]+hepmcp_pz[2]*hepmcp_pz[2]);      
      float ttr =   (phiVb.Px()+eVb.Px())*(phiVb.Px()+eVb.Px()) +
	(phiVb.Py()+eVb.Py())*(phiVb.Py()+eVb.Py()) ;            
     
      h_t_trans_rot->Fill(  ttr );

      h_t_trans->Fill(    (phiVb.Px()+eVb0.Px())*(phiVb.Px()+eVb0.Px()) +
			  (phiVb.Py()+eVb0.Py())*(phiVb.Py()+eVb0.Py()) );            


   }
}
