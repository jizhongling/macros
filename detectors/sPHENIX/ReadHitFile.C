#ifndef MACRO_READCLUSTERSRUNTRACKING_C
#define MACRO_READCLUSTERSRUNTRACKING_C

#include <GlobalVariables.C>

#include <G4_Magnet.C>
#include <G4_Tracking.C>

#include <g4eval/SvtxEvaluator.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libcentrality_io.so)

void ReadHitFile(const int nJob = 0, const int nEvents = 1)
{
  const int nOff = 1000+nJob;

  const string &embed_input_str1 = "/phenix/u/bogui/data/MacrosAug21/macros/detectors/sPHENIX/Aug21Hits/SvtxAug21_HitsPu50_020_1_";
  char num_field[500];
  sprintf(num_field,"%04d.root", nOff);
  string numin = num_field;
  string embed_infile1 = embed_input_str1+numin;

  const string &inputFile = Form("/phenix/spin/phnxsp01/zji/data/sphenix/output/G4sPHENIX%d.root",nJob);
  std::cout << "input file: " << inputFile << endl;

  //const string &outputFile = Form("Aug21Hits/EvalHitsPu50_020_%d.root",nOff);
  const string &outputFile = Form("/phenix/spin/phnxsp01/zji/data/sphenix/output/G4sPHENIX_g4svtx_eval%d.root",nJob);

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RANDOMSEED", nJob);
  rc->set_IntFlag("RUNNUMBER", 0);

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  MagnetInit(); // initialize fieldmap name
  TrackingInit();  // initialize tracking
  Mvtx_Clustering();
  Intt_Clustering();
  TPC_Clustering();
  Tracking_Reco();

  // explicitely configure the evaluator
  SvtxEvaluator* eval;
  const int n_maps_layer = 3;
  const int n_intt_layer = 4;
  int n_gas_layer  = 48;

  eval = new SvtxEvaluator("SVTXEVALUATOR", outputFile, "SvtxTrackMap",
                           G4MVTX::n_maps_layer,
                           G4INTT::n_intt_layer,
                           G4TPC::n_gas_layer,
                           G4MICROMEGAS::n_micromegas_layer);
  int do_default = 0;
  if(do_default)
  {
    eval->do_cluster_eval(true);
    eval->do_g4hit_eval(true);
    eval->do_hit_eval(true);  // enable to see the hits that includes the chamber physics...
    eval->do_gpoint_eval(false);
    eval->do_gtrack_eval(true);
    eval->do_eval_light(true);
  }
  else
  {
    eval->do_cluster_eval(true);
    eval->do_g4cluster_eval(false);
    eval->do_hit_eval(true);  // enable to see the hits that includes the chamber physics...
    eval->do_vertex_eval(false);
    eval->do_g4hit_eval(true);
    eval->do_gtrack_eval(false);
    eval->do_track_eval(false);
    eval->do_track_match(false);
    eval->do_gpoint_eval(false);
  }
  eval->do_eval_light(true);
  eval->set_use_initial_vertex(true);
  eval->set_use_genfit_vertex(false);
  eval->scan_for_embedded(false);  // take all tracks if false - take only embedded tracks if true
  eval->scan_for_primaries(false);
  eval->Verbosity(0);
  se->registerSubsystem(eval);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");

  //TString tstr_input(embed_infile1);
  //in->AddFile(embed_infile1);
  TString tstr_input(inputFile);
  if (tstr_input.EndsWith(".root"))
    in->AddFile(inputFile);
  else
    in->AddListFile(inputFile);

  se->registerInputManager(in);

  se->run(nEvents);
  std::cout << " Done Run, ending... " << std::endl;
  se->End();

  //se->PrintTimer();

  std::cout << " Success!! " << std::endl;
  // deleting the server shows if the memory is corrupted, if the job dies here - it is
  delete se; 
  gSystem->Exit(0);  
}

#endif
