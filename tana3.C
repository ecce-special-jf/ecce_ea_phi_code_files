{

  //gROOT->LoadMacro("/gpfs/mnt/gpfs02/sphenix/user/jfrantz/ecce/sartre/jlab/style/ecce-root-style/ECCEStyle.C");
  //SetECCEStyle();

event_tree->SetAlias("p0mag", "sqrt(tracks_px[0]*tracks_px[0]+tracks_py[0]*tracks_py[0]+tracks_pz[0]*tracks_pz[0])");
event_tree->SetAlias("p0magt", "sqrt(mcpart_px[nMCPart-4]*mcpart_px[nMCPart-4]+mcpart_py[nMCPart-4]*mcpart_py[nMCPart-4]+mcpart_pz[nMCPart-4]*mcpart_pz[nMCPart-4])");

 event_tree->SetAlias("p2mag","sqrt(tracks_px[2]^2+tracks_py[2]^2+tracks_pz[2]^2)");
 // event_tree->SetAlias("E0mag","sqrt(p0mag*p0mag+0.1056*0.1056)");
event_tree->SetAlias("p1mag", "sqrt(tracks_px[1]*tracks_px[1]+tracks_py[1]*tracks_py[1]+tracks_pz[1]*tracks_pz[1])");
//event_tree->SetAlias("E1mag","sqrt(p1mag*p1mag+0.1056*0.1056)");
event_tree->SetAlias("E0mag","sqrt(p0mag*p0mag+0.495*0.495)");
event_tree->SetAlias("E1mag","sqrt(p1mag*p1mag+0.495*0.495)");

 event_tree->SetAlias("emc_epx","cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_px[2]/sqrt(tracks_px[2]^2+tracks_py[2]^2))");

 event_tree->SetAlias("emc_epy","cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*(tracks_py[2]/sqrt(tracks_px[2]^2+tracks_py[2]^2))");


 // event_tree->SetAlias("emc_epy","cluster_EEMC_E[0]/cosh(cluster_EEMC_Eta[0])*sin(cluster_EEMC_Phi[0])");


 event_tree->SetAlias("phi_e","E0mag + E1mag");
 event_tree->SetAlias("phi_px","tracks_px[0]+tracks_px[1]");
 event_tree->SetAlias("phi_py","tracks_py[0]+tracks_py[1]");
 event_tree->SetAlias("phi_pz","tracks_pz[0]+tracks_pz[1]");

 event_tree->SetAlias("all_px","tracks_px[0]+tracks_px[1]+tracks_px[2]");
 event_tree->SetAlias("all_py","tracks_py[0]+tracks_py[1]+tracks_py[2]");

event_tree->SetAlias("mass_ks","sqrt((E0mag+E1mag)^2-(tracks_px[0]+tracks_px[1])^2-(tracks_py[0]+tracks_py[1])^2-(tracks_pz[0]+tracks_pz[1])^2)");
 event_tree -> SetAlias("eta_kaonst","0.5*log((p0mag+tracks_pz[0])/(p0mag-tracks_pz[0]))");
 event_tree -> SetAlias("eta_e_t","0.5*log((hepmcp_E[1]+hepmcp_pz[1])/(hepmcp_E[1]-hepmcp_pz[1]))");
 event_tree -> SetAlias("eta_kaonst2","0.5*log((p1mag+tracks_pz[1])/(p1mag-tracks_pz[1]))");
 event_tree -> SetAlias("eta_kaons_tr","0.5*log((p0magt+mcpart_pz[nMCPart-4])/(p0magt-mcpart_pz[nMCPart-4]))");
 event_tree -> SetAlias("eta_e_reco","0.5*log((p2mag+tracks_pz[2])/(p2mag-tracks_pz[2]))");
 event_tree -> SetAlias("eta_e_mct","0.5*log((mcpart_E[nMCPart-1]+mcpart_pz[nMCPart-1])/(mcpart_E[nMCPart -1]-mcpart_pz[nMCPart-1]))");


event_tree->SetAlias("t_true","hepmcp_E[2]^2-(hepmcp_px[2])^2-hepmcp_py[2]^2-hepmcp_pz[2]^2");
event_tree->SetAlias("t_true_t","-(hepmcp_px[2])^2-hepmcp_py[2]^2");
event_tree->SetAlias("q2","hepmcp_E[0]^2-(hepmcp_px[0])^2-hepmcp_py[0]^2-hepmcp_pz[0]^2");
 event_tree->SetAlias("t_r","(phi_e-hepmcp_E[0])^2-(phi_px-hepmcp_px[0])^2-(phi_py-hepmcp_py[0])^2-(phi_pz-hepmcp_pz[0])^2");
 event_tree->SetAlias("t_rr2","(phi_px+tracks_px[2])^2+(phi_py+tracks_py[2])^2");
 event_tree->SetAlias("emc_t_rr2","(phi_px+0.88*emc_epx)^2+(phi_py+0.88*emc_epy)^2");
  event_tree->SetAlias("t_rre2","(phi_px+hepmcp_px[1])^2+(phi_py+hepmcp_py[1])^2");
  event_tree->SetAlias("t_rre2_true","(hepmcp_px[5]+hepmcp_px[1])^2+(hepmcp_py[5]+hepmcp_py[1])^2");
  // event_tree->SetAlias("t_r","(hepmcp_-hepmcp_E[0])^2-(phi_px-hepmcp_px[0])^2-(phi_py-hepmcp_py[0])^2-(phi_pz-hepmcp_pz[0])^2");



}
