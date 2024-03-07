// g++ -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -L$MYINSTALL/lib -L$OFFLINE_MAIN/lib -lg4eval_io -ltrack_io `root-config --cflags --glibs` -o AnaElectron AnaElectron.cpp
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
  const char *ith = argv[5];

  auto f_out = new TFile(Form("%s-%s.root", output, ith), "RECREATE");
  auto ntp_out = new TNtuple("ntp", "ntp of track and proj", "event:px:py:pz:pt:eta:phi:deltapt:deltaeta:deltaphi:siqr:siphi:sithe:six0:siy0:tpqr:tpphi:tpthe:tpx0:tpy0:charge:quality:chisq:ndf:nhits:nmaps:nintt:ntpc:nmms:ntpc1:ntpc11:ntpc2:ntpc3:nlmaps:nlintt:nltpc:nlmms:layers:vx:vy:vz:dca3dxy:dca3dxysigma:dca3dz:dca3dzsigma:pcax:pcay:pcaz:npedge:nredge:nbig:novlp:merr:msize:nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclusintt:nclusmaps:nclusmms:clus_dphi_cemc:clus_dphi_hcalin:clus_dphi_hcalout:clus_dphi_outer_cemc:clus_dphi_outer_hcalin:clus_dphi_outer_hcalout:clus_deta_cemc:clus_deta_hcalin:clus_deta_hcalout:clus_deta_outer_cemc:clus_deta_outer_hcalin:clus_deta_outer_hcalout:clus_e_cemc:clus_e_hcalin:clus_e_hcalout:clus_e_outer_cemc:clus_e_outer_hcalin:clus_e_outer_hcalout");

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

    const int nvar = 103;
    const char *str_track[nvar] = {"event", "px", "py", "pz", "pt", "eta", "phi", "deltapt", "deltaeta", "deltaphi", "siqr", "siphi", "sithe", "six0", "siy0", "tpqr", "tpphi", "tpthe", "tpx0", "tpy0", "charge", "quality", "chisq", "ndf", "nhits", "nmaps", "nintt", "ntpc", "nmms", "ntpc1", "ntpc11", "ntpc2", "ntpc3", "nlmaps", "nlintt", "nltpc", "nlmms", "layers", "vx", "vy", "vz", "dca3dxy", "dca3dxysigma", "dca3dz", "dca3dzsigma", "pcax", "pcay", "pcaz", "npedge", "nredge", "nbig", "novlp", "merr", "msize", "nhittpcall", "nhittpcin", "nhittpcmid", "nhittpcout", "nclusall", "nclustpc", "nclusintt", "nclusmaps", "nclusmms", "clus_dphi_cemc", "clus_dphi_hcalin", "clus_dphi_hcalout", "clus_dphi_outer_cemc", "clus_dphi_outer_hcalin", "clus_dphi_outer_hcalout", "clus_deta_cemc", "clus_deta_hcalin", "clus_deta_hcalout", "clus_deta_outer_cemc", "clus_deta_outer_hcalin", "clus_deta_outer_hcalout", "clus_e_cemc", "clus_e_hcalin", "clus_e_hcalout", "clus_e_outer_cemc", "clus_e_outer_hcalin", "clus_e_outer_hcalout"};
    Float_t var_track[nvar];
    TNtuple *ntp_track = static_cast<TNtuple*>(f_in->Get("ntp_track"));
    if(!ntp_track) continue;
    for(int ivar = 0; ivar < nvar; ivar++)
      ntp_track->SetBranchAddress(str_track[ivar], &var_track[ivar]);

    for(Long64_t ien = 0; ien < ntp_track->GetEntries(); ien++)
    {
      ntp_track->GetEntry(ien);

      if(!TMath::Finite(var_track[nvar-6]))
        continue;

      ntp_out->Fill(var_track);
    } // ien

    delete f_in;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
