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
  auto ntp_out = new TNtuple("ntp", "ntp of cluster energy", "signal:gflavor:px:py:pz:vx:vy:vz:pcax:pcay:pcaz:dca3dxy:dca3dz:clus_e_cemc:clus_e_hcalin:clus_e_hcalout:clus_e_outer_cemc:clus_e_outer_hcalin:clus_e_outer_hcalout");

  for(Int_t ifile = istart; ifile < iend; ifile++)
  {
    char filename[1024];
    sprintf(filename, "%s/G4sPHENIX_g4svtx_eval-%d.root", indir, ifile);
    auto f_in = new TFile(filename);
    if(!f_in || f_in->IsZombie())
    {
      cerr << "Error: cannot open file " << filename << endl;
      continue;
    }

    const int nvar = 19;
    const char *str_gtrack[nvar] = {"parentID", "gflavor", "px", "py", "pz", "vx", "vy", "vz", "pcax", "pcay", "pcaz", "dca3dxy", "dca3dz", "clus_e_cemc", "clus_e_hcalin", "clus_e_hcalout", "clus_e_outer_cemc", "clus_e_outer_hcalin", "clus_e_outer_hcalout"};
    Float_t var_gtrack[nvar];
    Float_t gtrackID;
    TNtuple *ntp_gtrack = static_cast<TNtuple*>(f_in->Get("ntp_gtrack"));
    if(!ntp_gtrack) continue;
    for(int ivar = 0; ivar < nvar; ivar++)
      ntp_gtrack->SetBranchAddress(str_gtrack[ivar], &var_gtrack[ivar]);
    ntp_gtrack->SetBranchAddress("gtrackID", &gtrackID);

    for(Long64_t ien = 0; ien < ntp_gtrack->GetEntries(); ien++)
    {
      ntp_gtrack->GetEntry(ien);
      if(gtrackID < -1 || gtrackID == 1) continue;

      bool fill_this = true;
      for(int ivar = 0; ivar < nvar; ivar++)
        if(!TMath::Finite(var_gtrack[ivar]))
        {
          fill_this = false;
          break;
        }
      if(!fill_this) continue;

      if(var_gtrack[0] == 1 && var_gtrack[1] == -11)
        var_gtrack[0] = 1;
      else
        var_gtrack[0] = 0;

      ntp_out->Fill(var_gtrack);
    } // ien

    delete f_in;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
