// g++ -std=c++17 -Wall `root-config --cflags --glibs` -o AnaOmega AnaOmega.cpp
#include <iostream>

#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

using namespace std;

int main(int argc, const char *argv[])
{
  if(argc != 5)
  {
    cerr << "Usage: " << argv[0] << " <prefix> <start-index> <end-index> <out-file>" << endl;
    return 1;
  }
  const int istart = stoi(string(argv[2]));
  const int iend = stoi(string(argv[3]));

  auto f_out = new TFile(argv[4], "RECREATE");
  auto h2_m = new TH2F("h2_m", "Invariant mass without DCA;m_{#Omega};m_{#Lambda}", 2000,0.,20., 2000,0.,20);
  auto h2_m_dca = new TH2F("h2_m_dca", "Invariant mass with DCA;m_{#Omega};m_{#Lambda}", 2000,0.,20., 2000,0.,20);

  for(int ifile = istart; ifile < iend; ifile++)
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

    auto ntp_v0 = (TNtuple*)f->Get("ntp_v0");
    if(!ntp_v0 || ntp_v0->IsZombie())
    {
      cout << "Error: file " << filename << " has no ntp_v0" << endl;
      continue;
    }

    Float_t mother_pt, bch_pt, v0_pt, dau1_pt, dau2_pt, mother_eta, bch_eta, v0_eta, dau1_eta, dau2_eta, mother_m, v0_m, mother_path, v0_path, mother_decay_radius, v0_decay_radius, mother_costheta, v0_costheta, bch_dca, v0_dca, bchv0_dca, daus_dca, bch_charge, dau1_charge, dau2_charge, bch_quality, dau1_quality, dau2_quality, bch_nclus, dau1_nclus, dau2_nclus, bch_silicon, dau1_silicon, dau2_silicon, vtx_r, vtx_z;

    ntp_v0->SetBranchAddress("mother_pt", &mother_pt);
    ntp_v0->SetBranchAddress("bch_pt", &bch_pt);
    ntp_v0->SetBranchAddress("v0_pt", &v0_pt);
    ntp_v0->SetBranchAddress("dau1_pt", &dau1_pt);
    ntp_v0->SetBranchAddress("dau2_pt", &dau2_pt);
    ntp_v0->SetBranchAddress("mother_eta", &mother_eta);
    ntp_v0->SetBranchAddress("bch_eta", &bch_eta);
    ntp_v0->SetBranchAddress("v0_eta", &v0_eta);
    ntp_v0->SetBranchAddress("dau1_eta", &dau1_eta);
    ntp_v0->SetBranchAddress("dau2_eta", &dau2_eta);
    ntp_v0->SetBranchAddress("mother_m", &mother_m);
    ntp_v0->SetBranchAddress("v0_m", &v0_m);
    ntp_v0->SetBranchAddress("mother_path", &mother_path);
    ntp_v0->SetBranchAddress("v0_path", &v0_path);
    ntp_v0->SetBranchAddress("mother_decay_radius", &mother_decay_radius);
    ntp_v0->SetBranchAddress("v0_decay_radius", &v0_decay_radius);
    ntp_v0->SetBranchAddress("mother_costheta", &mother_costheta);
    ntp_v0->SetBranchAddress("v0_costheta", &v0_costheta);
    ntp_v0->SetBranchAddress("bch_dca", &bch_dca);
    ntp_v0->SetBranchAddress("v0_dca", &v0_dca);
    ntp_v0->SetBranchAddress("bchv0_dca", &bchv0_dca);
    ntp_v0->SetBranchAddress("daus_dca", &daus_dca);
    ntp_v0->SetBranchAddress("bch_charge", &bch_charge);
    ntp_v0->SetBranchAddress("dau1_charge", &dau1_charge);
    ntp_v0->SetBranchAddress("dau2_charge", &dau2_charge);
    ntp_v0->SetBranchAddress("bch_quality", &bch_quality);
    ntp_v0->SetBranchAddress("dau1_quality", &dau1_quality);
    ntp_v0->SetBranchAddress("dau2_quality", &dau2_quality);
    ntp_v0->SetBranchAddress("bch_nclus", &bch_nclus);
    ntp_v0->SetBranchAddress("dau1_nclus", &dau1_nclus);
    ntp_v0->SetBranchAddress("dau2_nclus", &dau2_nclus);
    ntp_v0->SetBranchAddress("bch_silicon", &bch_silicon);
    ntp_v0->SetBranchAddress("dau1_silicon", &dau1_silicon);
    ntp_v0->SetBranchAddress("dau2_silicon", &dau2_silicon);
    ntp_v0->SetBranchAddress("vtx_r", &vtx_r);
    ntp_v0->SetBranchAddress("vtx_z", &vtx_z);

    for(Long64_t ien; ien<ntp_v0->GetEntries(); ien++)
    {
      ntp_v0->GetEntry(ien);

      if( bch_pt>0.2 && dau1_pt>0.2 && dau2_pt>0.2 )
        h2_m->Fill(mother_m, v0_m);

      if( bch_pt>0.2 && dau1_pt>0.2 && dau2_pt>0.2 &&
          mother_path>1.25 && mother_path<30 && v0_path>1.75 && v0_path<30 && mother_path<v0_path &&
          mother_decay_radius>0.6 && v0_decay_radius>1.4 && mother_costheta>0.97 && v0_costheta>0.97 )
        h2_m_dca->Fill(mother_m, v0_m);
    }

    delete f;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
