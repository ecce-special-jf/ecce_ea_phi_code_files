#define thnew_et_cxx
#include "thnew_et.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void thnew_et::Loop()
{
//   In a ROOT session, you can do:
//      root> .L thnew_et.C
//      root> thnew_et t
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


     if (cluster_EEMC_N != 1) continue;
      if (tracks_trueID[0] != 1 || tracks_trueID[1] != 2) continue;
      if (nTracks < 2) continue;
      // else if (nTracks == 3 && tracks_trueID[2] != 4 ) continue;
      // if (fabs(vertex_z) >3) continue;
      
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
      float res1 = 0.02, res2 = 0.01;


      //reset to more reasonable smeared value;
      cluster_EEMC_E[0] = true_e 
	+ strue*res1*f_egausres1->GetRandom() 
	+ res2*true_e*f_egausres1->GetRandom();
      

      //for debug 
      //cluster_EEMC_Eta[0] = //(eVb2a.Eta() + cluster_EEMC_Eta[2])/2;

      cluster_EEMC_Eta[0] = eVb2a.Eta();
      if (nTracks == 2 || tracks_trueID[2] != 4)
	cluster_EEMC_Eta[0] = eVb02a.Eta()*(1+0.004*f_egausres1->GetRandom());
	//tracking eta is very good when you have it

      //
      // figure out why h_t_trans is missing some counts they are in the over flow bin



      
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

      //tracking cosphi peters out around 3.2
      if (eVb02a.Eta() < -3.2)
	{
	  cosphi_e = cos(eVb02a.Phi())*(1.+0.017*f_egausres1->GetRandom());
	  sinphi_e = sin(eVb02a.Phi())*(1.+0.017*f_egausres1->GetRandom());
	}

      //      float cosphi_e = cos(eVb02a.Phi());
      //float sinphi_e = sin(eVb02a.Phi());

      if (jentry < 1000)
	if ((nTracks == 2 || tracks_trueID[2] != 4 ) && eVb02a.Eta() < -3.2)
	  cout << "true e" << true_e <<  " " << cluster_EEMC_E[0]  << "  --- cos " << cos(eVb02a.Phi())<< "   " << cosphi_e << " eta : " << eVb02a.Eta() <<  "  " << cluster_EEMC_Eta[0] <<endl;



      float emc_epx = pxcal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*cosphi_e;
      float emc_epy = pycal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*sinphi_e;

	/*
      float emc_epx = pxcal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_px[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]));

      float emc_epy = pycal*cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_py[2]/sqrt(tracks_px[2]*tracks_px[2]+tracks_py[2]*tracks_py[2]));

	*/

      float emc_epz = pzcal*cluster_EEMC_E[0]*tanh(cluster_EEMC_Eta[0]);

      ////////////////////


      float all_px = tracks_px[0]+tracks_px[1]+tracks_px[2];

      float  p0mag=sqrt(tracks_px[0]*tracks_px[0]+tracks_py[0]*tracks_py[0]+tracks_pz[0]*tracks_pz[0]);    
      float  E0mag=sqrt(p0mag*p0mag+0.495*0.495);

      float  p1mag =sqrt(tracks_px[1]*tracks_px[1]+tracks_py[1]*tracks_py[1]+tracks_pz[1]*tracks_pz[1]);
      
      float  E1mag=sqrt(p1mag*p1mag+0.495*0.495);

      //      if (p2mag > 10) continue;


      TLorentzVector k0Vb;
      TLorentzVector k0Vbt;
      k0Vb.SetPxPyPzE(tracks_px[0],tracks_py[0], tracks_pz[0],E0mag);
      k0Vbt.SetPxPyPzE(hepmcp_px[7], hepmcp_py[7], hepmcp_pz[7], hepmcp_E[7]);


      //      if (fabs(tracks_px[1]/)

      TLorentzVector k1Vb;
      k1Vb.SetPxPyPzE(tracks_px[1],tracks_py[1], tracks_pz[1],E1mag);
      TLorentzVector k1Vbt;
      k1Vbt.SetPxPyPzE(hepmcp_px[6], hepmcp_py[6], hepmcp_pz[6], hepmcp_E[6]);


      //      if (fabs(k1Vb.Eta()) > 0.95 || fabs(k0Vb.Eta()) > 0.95) continue;
      if (fabs(k1Vb.Eta()) > 2.0 || fabs(k0Vb.Eta()) > 2.0) continue;
      // if (fabs(k1Vb.Eta()) < 0.17) continue;
      //      if (fabs(tracks_py[1]/tracks_px[1]) < 1.5 || fabs(tracks_py[0]/tracks_px[0]) < 1.5) continue;
      float scale1 = 0.0;
      if (fabs(k1Vb.Eta()) < 2.0 && fabs(k1Vb.Eta()) > 1.0) 
	{
	  if (fabs(k1Vb.Eta()) > 1.5) scale1 = 0.83; else scale1 = 0.61;   
	  tracks_px[1] = mcpart_px[nMCPart - 3] +  scale1 * (tracks_px[1] - mcpart_px[nMCPart-3]);
	  tracks_py[1] = mcpart_py[nMCPart - 3] +  scale1 * (tracks_py[1] - mcpart_py[nMCPart-3]);
	  tracks_pz[1] = mcpart_pz[nMCPart - 3] +  scale1 * (tracks_pz[1] - mcpart_pz[nMCPart-3]);
	  p1mag =sqrt(tracks_px[1]*tracks_px[1]+tracks_py[1]*tracks_py[1]+tracks_pz[1]*tracks_pz[1]);
	  //	  tracks_px[1] = mcpart_px[nMCPart - 4] +  0.5 * (tracks_px[1] - mcpart_px[nMCPart-4]);
	  E1mag=sqrt(p1mag*p1mag+0.495*0.495);
	}
      
      if (fabs(k0Vb.Eta()) < 2.0 && fabs(k0Vb.Eta()) > 1.0) 
	{
	  if (fabs(k0Vb.Eta()) > 1.5) scale1 = 0.83; else scale1 = 0.61;   
	  tracks_px[0] = mcpart_px[nMCPart - 4] +  scale1 * (tracks_px[0] - mcpart_px[nMCPart-4]);
	  tracks_py[0] = mcpart_py[nMCPart - 4] +  scale1 * (tracks_py[0] - mcpart_py[nMCPart-4]);
	  tracks_pz[0] = mcpart_pz[nMCPart - 4] +  scale1 * (tracks_pz[0] - mcpart_pz[nMCPart-4]);
	  p0mag =sqrt(tracks_px[0]*tracks_px[0]+tracks_py[0]*tracks_py[0]+tracks_pz[0]*tracks_pz[0]);
	  //	  tracks_px[1] = mcpart_px[nMCPart - 4] +  0.5 * (tracks_px[1] - mcpart_px[nMCPart-4]);
	  E1mag=sqrt(p1mag*p1mag+0.495*0.495);
	}

      k0Vb.SetPxPyPzE(tracks_px[0],tracks_py[0], tracks_pz[0],E0mag);



      float  phi_e=E0mag + E1mag;
      float  phi_px=tracks_px[0]+tracks_px[1];
      float  phi_py=tracks_py[0]+tracks_py[1];
      float  phi_pz=tracks_pz[0]+tracks_pz[1];

      TLorentzVector phiVb;
      phiVb.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      //      TLorentzVector eVb;
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);      

      //      eVb.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);
      //      eVb.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      if ((nTracks == 2 || tracks_trueID[2] != 4) && eVb02a.Eta() < -3.2)
	{
	tracks_px[2] = emc_epx;
	tracks_py[2] = emc_epy;
	tracks_pz[2] = emc_epz;
	p2mag = sqrt(emc_epx*emc_epx+emc_epy*emc_epy+emc_epz*emc_epz);
	}
      else
	{
      tracks_px[2] = 1.02*tracks_px[2];
      tracks_py[2] = 1.02*tracks_py[2];
      tracks_pz[2] = 1.02*tracks_pz[2];
	}

      float av_epx = (tracks_px[2]+ emc_epx)/2.0;
      float av_epy = (tracks_py[2]+ emc_epy)/2.0;
      float av_epz = (tracks_pz[2]+ emc_epz)/2.0;
      float av_eE = (p2mag+cluster_EEMC_E[0])/2.0;

      if (true)
      if (fabs(tracks_px[2]/emc_epx -1) > 0.05 || 
	  fabs(tracks_py[2]/emc_epy -1) > 0.05 || 
	  fabs(tracks_pz[2]/emc_epz -1) > 0.05)
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
      //      eVb.SetPxPyPzE(emc_epx,emc_epy,emc_epz,cluster_EEMC_E[0]);
      eVb2.SetPxPyPzE(av_epx,av_epy,av_epz,av_eE);

      TLorentzVector phiVb0;
      
      phiVb0.SetPxPyPzE(phi_px,phi_py, phi_pz,phi_e);

      //      TLorentzVector eVb0;
      //      eVb0.SetPxPyPzE(hepmcp_px[1],hepmcp_py[1],hepmcp_pz[1],hepmcp_E[1]);
      //      eVb0.SetPxPyPzE(tracks_px[2],tracks_py[2],tracks_pz[2],p2mag);

      //      if (eVb.Eta() < -3.2 || eVb.Eta() > -2.65) continue;
 
      if (eVb.Eta() > -2.65) continue;
      //tracking phi resolution gets bad above 3.2
      
      TLorentzRotation lrTrans;
      lrTrans.RotateY(+12.5e-3).Boost(sin(+12.5e-3),0,0);
      //lrTrans.RotateY(+12.5e-3);
      
      phiVb *= lrTrans;
      eVb *= lrTrans;

      
      if (fabs(phiVb.M()-1.02) > 0.04) continue;  


      phiVb0 *=lrTrans;
      //phiVb0 *=lrTrans;
				      
      k0Vb *= lrTrans;
      
      // if (!(fabs(k0Vbt.Eta()) > 1.0 && 
      // 	    fabs(k0Vbt.Eta()) < 1.5) 
      // 	  && !(fabs(k0Vbt.Eta()) < 0.2))
      h2px->Fill(k0Vb.Eta(), k0Vb.P()/sqrt(hepmcp_px[7]*hepmcp_px[7]+hepmcp_py[7]*hepmcp_py[7]+hepmcp_pz[7]*hepmcp_pz[7]));
      //      h2px->Fill(k0Vb.Eta(), k0Vb.P()/k0);
      //     h2px->Fill(eVb.Eta(), eVb.Px()/hepmcp_px[1]);

      h_px0->Fill(eVb.Px()); 
      //      h_pxb->Fill(phiVb.Px() + eVb.Px() ); 

      h_pxb->Fill(phiVb.M()); 
      
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
