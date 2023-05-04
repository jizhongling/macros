// g++ -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -L$MYINSTALL/lib -L$OFFLINE_MAIN/lib -lg4eval_io -ltrack_io `root-config --cflags --glibs` -o AnaHFElectron AnaHFElectron.cpp
#include <iostream>
#include <memory>
#include <array>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <g4eval/TrackEvaluationContainerv1.h>

using namespace std;

template<class T> inline constexpr T square(const T &x) { return x*x; }
template<class T> inline constexpr T get_r(const T &x, const T &y) { return sqrt(square(x) + square(y)); }

int main(int argc, const char *argv[])
{
  if(argc != 6)
  {
    cerr << "Usage: " << argv[0] << " <indir> <istart> <iend> <output>(-<ith>.root) <ith>" << endl;
    return 1;
  }

  const char *indir = argv[1];
  const Int_t istart = stoi(string(argv[2]));
  const Int_t iend = stoi(string(argv[3]));
  const char *output = argv[4];
  const Int_t ith = stoi(string(argv[5]));

  auto f_out = new TFile(Form("%s-%d.root", output, ith), "RECREATE");
  auto ntp_out = new TNtuple("ntp", "ntp of cluster energy", "trackID:pid:parent_pid:dca_r:dca_z:clus_e_cemc:clus_e_hcalin:clus_e_hcalout:clus_e_outer_cemc:clus_e_outer_hcalin:clus_e_outer_hcalout");

  for(Int_t ifile = istart; ifile < iend; ifile++)
  {
    char filename[1024];
    sprintf(filename, "%s/G4sPHENIX_g4svtx_eval-%d.root", indir, ifile);
    auto f_in = new TFile(filename);
    if(f_in->IsZombie())
    {
      cerr << "Error: cannot open file " << filename << endl;
      continue;
    }

    TrackEvaluationContainerv1::TrackStruct::List *v_tracks = nullptr;
    TTree *t_trackeval = static_cast<TTree*>(f_in->Get("t_trackeval"));
    t_trackeval->SetBranchAddress("tracks", &v_tracks);

    for(Long64_t ien = 0; ien < t_trackeval->GetEntries(); ien++)
    {
      t_trackeval->GetEntry(ien);
      if(!v_tracks)
      {
        cerr << "Entry " << ien << " has no v_tracks!" << endl;
        continue;
      }

      for(const auto &track : *v_tracks)
      {
        const Int_t nlayers = 6;
        Float_t fill_ntp[nlayers+5] = {(Float_t)track.trackID, (Float_t)track.pid, (Float_t)track.parent_pid, sqrt(square(track.x)+square(track.y)), track.z};
        for(Int_t li = 0; li < nlayers; li++)
          fill_ntp[li+5] = track.cal_cluster_e[li];
        ntp_out->Fill(fill_ntp);
      } // track
    } // ien

    delete f_in;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
