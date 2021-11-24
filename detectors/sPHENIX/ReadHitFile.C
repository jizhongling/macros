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


void ReadHitFile(const int nJob = 1)
{

  const int nEvents = 1;
  const float pTval = 1.0;
  const int nOff = 1000+nJob;

  //  printf("Scale = %d\n",scale);

  //  char inputFile[256];
  // sprintf(inputFile,"/phenix/u/bogui/RunJune21/June21OnlyHit/SvtxJune21_HitsOnlyPu50_020_dist_1_%d.root",nOff);

  //  const string &embed_input_str1 = "DST_TRKR_RECOCLUSTER_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000001-0";
  const string &embed_input_str1 = "/phenix/u/bogui/data/MacrosAug21/macros/detectors/sPHENIX/Aug21Hits/SvtxAug21_HitsPu50_020_1_";
  char num_field[500];
  sprintf(num_field,"%04d.root", nOff);
  string numin = num_field;
  string embed_infile1 = embed_input_str1+numin;

  std::cout << "input file1: " << embed_infile1 << endl;


  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
 
  const string &inputFile = embed_input_str1;  // dummy


  char outputFile[256];
  sprintf(outputFile,"Aug21Hits/EvalHitsPu50_020_%d.root",nOff);

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RANDOMSEED", nJob);
  rc->set_IntFlag("RUNNUMER", nJob);

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(10);

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

  eval = new SvtxEvaluator("SVTXEVALUATOR", outputFile, "SvtxTrackMap", n_maps_layer, n_intt_layer, n_gas_layer);
  int do_default = 0;
  if(do_default){
    eval->do_cluster_eval(true);
    eval->do_g4hit_eval(true);
    eval->do_hit_eval(true);  // enable to see the hits that includes the chamber physics...
    eval->do_gpoint_eval(false);
    eval->do_gtrack_eval(true);
    eval->do_eval_light(true);
  }else{
    eval->do_cluster_eval(true);
    eval->do_g4cluster_eval(false);
    eval->do_hit_eval(false);  // enable to see the hits that includes the chamber physics...
    
    eval->do_vertex_eval(false);
    eval->do_g4hit_eval(false);
    eval->do_gtrack_eval(true);
    eval->do_track_eval(true);
    eval->do_track_match(true);
    
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
  
  //  TString tstr_input(inputFile);
  TString tstr_input(embed_infile1);
  in->AddFile(embed_infile1);
  /*
  if (tstr_input.EndsWith(".root"))
    in->AddFile(inputFile);
  else
    in->AddListFile(inputFile);
  */

  se->registerInputManager(in);
  
  se->run(nEvents);
  std::cout << " Done Run, ending... " << std::endl;
  se->End();

  se->PrintTimer();

  std::cout << " Success!! " << std::endl;
// deleting the server shows if the memory is corrupted, if the job dies here - it is
  delete se; 
  gSystem->Exit(0);  
  
}

#endif
