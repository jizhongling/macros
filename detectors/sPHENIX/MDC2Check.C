#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <fun4all/Fun4AllServer.h>
#include <frog/FROG.h>

R__LOAD_LIBRARY(libfun4all.so)

using namespace std;

void MDC2Check(const Int_t runno = 62, const Int_t index0 = 0, const Int_t nindex = 10)
{
  const Int_t ntype = 2;
  const char *prefix[ntype] = {"DST_TRUTH", "DST_TRKR_G4HIT"};
  TGraph *g_neve[ntype];

  FROG *fr = new FROG();

  for(Int_t itype = 0; itype < ntype; itype++)
  {
    g_neve[itype] = new TGraph(nindex);
    for(Int_t index = index0; index < index0 + nindex; index++)
    {
      const char *fullname = fr->location(Form("%s_sHijing_0_20fm_50kHz_bkg_0_20fm-%010d-%05d.root", prefix[itype], runno, index));
      auto f = new TFile(fullname);
      auto t = static_cast<TTree*>(f->Get("T"));
      Long64_t nentries = t->GetEntries();
      g_neve[itype]->SetPoint(index-index0, index, nentries);
      cout << fullname << "\tNEvents = " << nentries  << endl;
    } // index
  } // itype

  if(true)
  {
    auto f_out = new TFile("mdc2check.root", "RECREATE");
    for(Int_t itype = 0; itype < ntype; itype++)
      g_neve[itype]->Write();
    f_out->Close();
  }
}
