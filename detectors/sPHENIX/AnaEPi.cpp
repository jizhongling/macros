// g++ -Wall `root-config --cflags --glibs` -o AnaEPi AnaEPi.cpp
#include <iostream>
#include <vector>
#include <algorithm>

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
  const Int_t ith = stoi(string(argv[5]));

  auto f_out = new TFile(Form("%s-%d.root", output, ith), "RECREATE");
  auto ntp_out = new TNtuple("ntp", "ntp of shower", "gflavor:etruth:ecemc:ehcalin:ehcalout");

  for(Int_t ifile = istart; ifile < iend; ifile++)
  {
    const Int_t ntypes = 3;
    const char *type[ntypes] = {"cemc", "hcalin", "hcalout"};

    TFile *f_in[ntypes];
    TNtuple *ntp_gshower[ntypes];
    Float_t event[ntypes], gparticleID[ntypes], gflavor[ntypes], etruth[ntypes], ereco[ntypes];

    bool do_fill = true;
    for(Int_t it = 0; it < ntypes; it++)
    {
      char filename[1024];
      sprintf(filename, "%s/G4sPHENIX_g4%s_eval-%d.root", indir, type[it], ifile);
      f_in[it] = new TFile(filename);
      if(f_in[it]->IsZombie())
      {
        cerr << "Error: cannot open file " << filename << endl;
        do_fill = false;
        break;
      }

      ntp_gshower[it] = static_cast<TNtuple*>(f_in[it]->Get("ntp_gshower"));
      if(!ntp_gshower[it])
      {
        do_fill = false;
        break;
      }

      ntp_gshower[it]->SetBranchAddress("event", &event[it]);
      ntp_gshower[it]->SetBranchAddress("gparticleID", &gparticleID[it]);
      ntp_gshower[it]->SetBranchAddress("gflavor", &gflavor[it]);
      ntp_gshower[it]->SetBranchAddress("ge", &etruth[it]);
      ntp_gshower[it]->SetBranchAddress("e", &ereco[it]);
    } // it
    if(!do_fill) continue;

    for(Long64_t ien = 0; ien < ntp_gshower[0]->GetEntries(); ien++)
    {
      vector<Float_t> v_labels;
      for(Int_t it = 0; it < ntypes; it++)
      {
        ntp_gshower[it]->GetEntry(ien);
        v_labels.emplace_back(event[it]);
        v_labels.emplace_back(gparticleID[it]);
        v_labels.emplace_back(gflavor[it]);
      }

      if( equal(v_labels.begin(), v_labels.begin()+3, v_labels.begin()+3, v_labels.begin()+6) &&
          equal(v_labels.begin(), v_labels.begin()+3, v_labels.begin()+6, v_labels.end()) )
      {
        Float_t fill_ntp[] = {gflavor[0], etruth[0], ereco[0], ereco[1], ereco[2]};
        ntp_out->Fill(fill_ntp);
      }
    } // ien

    for(Int_t it = 0; it < ntypes; it++)
      delete f_in[it];
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
