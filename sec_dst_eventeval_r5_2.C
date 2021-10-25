#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <eiceval/EventEvaluatorEIC.h>

#include <eccefastpidreco/ECCEFastPIDReco.h>
#include <eccefastpidreco/ECCEmRICHFastPIDMap.h>
#include <eccefastpidreco/ECCEdRICHFastPIDMap.h>
#include <eccefastpidreco/ECCEhpDIRCFastPIDMap.h>


R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libeiceval.so)
R__LOAD_LIBRARY(libECCEFastPIDReco.so)

namespace EVENT_EVALUATOR
{
  int Verbosity = 0;
  float EnergyThreshold = 0.05;
}  // namespace EVENT_EVALUATOR

void sec_dst_eventeval(const std::string fileList="pb-list.txt", int nEvents=100, std::string outputFile="event_eval.root")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);
  int verbosity = 0;

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->AddListFile(fileList);
  se->registerInputManager(hitsin);


  ECCEmRICHFastPIDMap * pidmap = new ECCEmRICHFastPIDMap();
  pidmap->Verbosity(verbosity);
  ECCEFastPIDReco * reco = new ECCEFastPIDReco(pidmap, EICPIDDefs::mRICH, "ECCEFastPIDReco-mRICH");
  se->registerSubsystem(reco);

  ECCEdRICHFastPIDMap * pidmapd = new ECCEdRICHFastPIDMap();
  pidmapd->Verbosity(verbosity);
  pidmapd->dualRICH_aerogel();
  ECCEFastPIDReco * recod = new ECCEFastPIDReco(pidmapd, EICPIDDefs::dRICH_AeroGel, "ECCEFastPIDReco-dRICH_AeroGel");
  se->registerSubsystem(recod);

  ECCEdRICHFastPIDMap * pidmapdg = new ECCEdRICHFastPIDMap();
  pidmapdg->Verbosity(verbosity);
  pidmapdg->dualRICH_C2F6();

  ECCEFastPIDReco * recodg = new ECCEFastPIDReco(pidmapd, EICPIDDefs::dRICH_Gas, "ECCEFastPIDReco-dRICH_Gas");
  se->registerSubsystem(recodg);
  
  ECCEhpDIRCFastPIDMap * pidmapD = new ECCEhpDIRCFastPIDMap();
  pidmapD->ReadMap( string(getenv("CALIBRATIONROOT")) + string("/hpDIRC/FastPID/ctr_map_p1_0.95.root"));
  ECCEFastPIDReco * recoD = new ECCEFastPIDReco(pidmap, EICPIDDefs::DIRC, "ECCEFastPIDReco-DIRC");
  se->registerSubsystem(recoD);
  

  EventEvaluatorEIC *eval = new EventEvaluatorEIC("EVENTEVALUATOR", outputFile);


  eval->set_reco_tracing_energy_thresholdMC(EVENT_EVALUATOR::EnergyThreshold);
  eval->Verbosity(EVENT_EVALUATOR::Verbosity);
  
  eval->set_do_TRACKS(true);
  eval->set_do_PROJECTIONS(true);
  eval->set_do_VERTEX(true);

  eval->set_do_HITS(true);

  eval->set_do_CEMC(true);
  eval->set_do_EEMC(true);
  eval->set_do_FEMC(true);
  eval->set_do_HCALIN(true);
  eval->set_do_HCALOUT(true);
  eval->set_do_FHCAL(true);
  eval->set_do_CLUSTERS(true);
  eval->set_do_PID_LogLikelihood(true); 

  eval->set_do_MCPARTICLES(true);
  eval->set_do_HEPMC(true);
  eval->set_do_store_event_level_info(true);

  se->registerSubsystem(eval);

  se->run(nEvents);
  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
}
