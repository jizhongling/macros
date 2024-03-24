// g++ -Wall `root-config --cflags --glibs` -o AnaConv AnaConv.cpp
#include <iostream>
#include <string>

#include <TFile.h>
#include <TNtuple.h>

using namespace std;

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

  const Int_t nvar = 76;
  Float_t var_ntp[nvar];
  string var_name[nvar] = {"x1", "y1", "z1", "px1", "py1", "pz1", "dca3dxy1", "dca3dz1", "vposx1", "vposy1", "vposz1", "vmomx1", "vmomy1", "vmomz1", "pca_relx_1", "pca_rely_1", "pca_relz_1", "eta1", "phi1", "charge1", "tpcClusters_1", "quality1", "x2", "y2", "z2", "px2", "py2", "pz2", "dca3dxy2", "dca3dz2", "vposx2", "vposy2", "vposz2", "vmomx2", "vmomy2", "vmomz2", "pca_relx_2", "pca_rely_2", "pca_relz_2", "eta2", "phi2", "charge2", "tpcClusters_2", "quality2", "vertex_x", "vertex_y", "vertex_z", "pair_dca", "invariant_mass", "invariant_pt", "path", "runNumber", "eventNumber", "nEventTracks", "trackID1", "has_silicon1", "clus1_dphi_cemc", "clus1_dphi_hcalin", "clus1_dphi_hcalout", "clus1_deta_cemc", "clus1_deta_hcalin", "clus1_deta_hcalout", "clus1_e_cemc", "clus1_e_hcalin", "clus1_e_hcalout", "trackID2", "has_silicon2", "clus2_dphi_cemc", "clus2_dphi_hcalin", "clus2_dphi_hcalout", "clus2_deta_cemc", "clus2_deta_hcalin", "clus2_deta_hcalout", "clus2_e_cemc", "clus2_e_hcalin", "clus2_e_hcalout"};

  auto t_out = new TTree("ntp", "decay_pairs");
  for(Int_t ivar=0; ivar<nvar; ivar++)
    t_out->Branch(var_name[ivar].c_str(), &var_ntp[ivar], (var_name[ivar]+"/F").c_str());

  for(Int_t ifile = istart; ifile < iend; ifile++)
  {
    char filename[1024];
    sprintf(filename, "%s/G4sPHENIX_secvert_ntuple-%d.root", indir, ifile);
    auto f_in = new TFile(filename);
    if(!f_in || f_in->IsZombie())
    {
      cerr << "Error: cannot open file " << filename << endl;
      continue;
    }

    auto ntp_in = (TNtuple*)f_in->Get("ntp");
    for(Int_t ivar=0; ivar<nvar; ivar++)
      ntp_in->SetBranchAddress(var_name[ivar].c_str(), &var_ntp[ivar]);

    for(Long64_t ien = 0; ien < ntp_in->GetEntries(); ien++)
    {
      ntp_in->GetEntry(ien);
      if(var_ntp[0] != 0)
        t_out->Fill();
    }

    delete f_in;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
