// g++ -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -L$MYINSTALL/lib -L$OFFLINE_MAIN/lib -lg4eval_io -ltrack_io `root-config --cflags --glibs` -o AnaPull AnaPull.cpp
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>

#include <g4eval/TrackEvaluationContainerv1.h>

using namespace std;

int main(int argc, const char *argv[])
{
  if(argc != 5)
  {
    cerr << "Usage: " << argv[0] << " <prefix> <start-index> <end-index> <out-file>" << endl;
    return 1;
  }
  const Int_t istart = stoi(string(argv[2]));
  const Int_t iend = stoi(string(argv[3]));

  const Int_t max_ntouch = 10;
  auto f_out = new TFile(argv[4], "RECREATE");
  TH1 *h_phi[3][max_ntouch+1];
  TH1 *h_z[3][max_ntouch+1];
  for(Int_t ilayer = 0; ilayer <= 2; ilayer++)
    for(Int_t itouch = 0; itouch <= max_ntouch; itouch++)
    {
      h_phi[ilayer][itouch] = new TH1F(Form("h_phi_li%d_ovl%d",ilayer,itouch), Form("TPC phi pull layer %d ntouch %d",ilayer,itouch), 50,-5.,5.);
      h_z[ilayer][itouch] = new TH1F(Form("h_z_li%d_ovl%d",ilayer,itouch), Form("TPC z pull layer %d ntouch %d",ilayer,itouch), 50,-5.,5.);
    }

  for(Int_t ifile = istart; ifile < iend; ifile++)
  {
    char filename[200];
    sprintf(filename, "%s%d.root", argv[1], ifile);

    TFile *f = TFile::Open(filename);
    if(!f || f->IsZombie())
    {
      cout << "Error: cannot open file " << filename << endl;
      if(f) delete f;
      continue;
    }
    cout << "Open file " << filename << endl;

    TrackEvaluationContainerv1::TrackStruct::List *v_tracks = nullptr;
    auto t_trackeval = static_cast<TTree*>(f->Get("t_trackeval"));
    t_trackeval->SetBranchAddress("tracks", &v_tracks);
    Long64_t nev = t_trackeval->GetEntries();

    TrkrDefs::cluskey training_cluskey;
    Short_t training_event, training_ntouch;
    TTree *t_training = static_cast<TTree*>(f->Get("t_training"));
    t_training->SetBranchAddress("event", &training_event);
    t_training->SetBranchAddress("cluskey", &training_cluskey);
    t_training->SetBranchAddress("ntouch", &training_ntouch);
    Long64_t training_nen = t_training->GetEntries();
    Long64_t training_ien = 0;

    for(Long64_t iev = 0; iev < nev; iev++)
    {
      t_trackeval->GetEntry(iev);
      if(!v_tracks)
      {
        cerr << "Event " << iev << " has no v_tracks!" << endl;
        continue;
      }

      map<TrkrDefs::cluskey, Short_t> m_cluskey;
      while(training_ien < training_nen)
      {
        t_training->GetEntry(training_ien);
        if(training_event < iev)
        {
          training_ien++;
          continue;
        }
        else if(training_event > iev)
        {
          break;
        }
        training_ien++;
        m_cluskey.insert(make_pair(training_cluskey, training_ntouch));
      }

      for(const auto &track : *v_tracks)
        if( track.pt > 0.3 && track.pt < 1. &&
            track.nclusters_mvtx >= 3 && track.nclusters_intt >= 1 && track.nclusters_tpc > 30 )
          for(const auto &cluster : track.clusters)
            if(cluster.phi_size >= 2)
            {
              Int_t li = -1;
              if(cluster.layer < 7)
                li = -1;
              else if(cluster.layer < 23)
                li = 0;
              else if(cluster.layer < 39)
                li = 1;
              else if(cluster.layer < 55)
                li = 2;
              Int_t ntouch = m_cluskey[cluster.key];

              Float_t dphi = (cluster.phi - cluster.truth_phi) / cluster.phi_error;
              Float_t dz = (cluster.z - cluster.truth_z) / cluster.z_error;
              if(li >= 0 && ntouch <= max_ntouch)
              {
                h_phi[li][ntouch]->Fill(dphi);
                h_z[li][ntouch]->Fill(dz);
              }
            }
    } // iev

    delete f;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
