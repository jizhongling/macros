// g++ -std=c++17 -Wall -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib `root-config --cflags --glibs` -o AnaQA AnaQA.cpp
#include <iostream>
#include <map>
#include <algorithm>

#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>

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
  TH1 *h_truth = new TH1F("h_truth","Truth tracks", 200,0.,20.);
  TH1 *h_reco = new TH1F("h_reco","Reco tracks", 200,0.,20.);
  TH1 *h_all = new TH1F("h_all","All reco tracks", 40,0.,20.);
  TH1 *h_good = new TH1F("h_good","Reco tracks with correct truth association", 40,0.,20.);
  TH2 *h2_resol = new TH2F("h2_resol","Momentum resolution", 40,0.,20., 400,0.5,1.5);

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

    Float_t event;
    Float_t gtrackID, gpt, gembed;
    Float_t trackID, pt, nmaps, ntpc;

    TNtuple *ntp_gtrack = (TNtuple*)f->Get("ntp_gtrack");
    if(!ntp_gtrack || ntp_gtrack->IsZombie())
    {
      cout << "Error: file " << filename << " has no ntp_gtrack" << endl;
      continue;
    }
    Long64_t nen_truth = ntp_gtrack->GetEntries();
    ntp_gtrack->SetBranchAddress("event", &event);
    ntp_gtrack->SetBranchAddress("gtrackID", &gtrackID);
    ntp_gtrack->SetBranchAddress("gpt", &gpt);
    ntp_gtrack->SetBranchAddress("gembed", &gembed);
    ntp_gtrack->SetBranchAddress("trackID", &trackID);
    ntp_gtrack->SetBranchAddress("pt", &pt);
    ntp_gtrack->SetBranchAddress("nmaps", &nmaps);
    ntp_gtrack->SetBranchAddress("ntpc", &ntpc);

    TNtuple *ntp_track = (TNtuple*)f->Get("ntp_track");
    if(!ntp_track || ntp_track->IsZombie())
    {
      cout << "Error: file " << filename << " has no ntp_track" << endl;
      continue;
    }
    Long64_t nen_reco = ntp_track->GetEntries();
    ntp_track->SetBranchAddress("event", &event);
    ntp_track->SetBranchAddress("trackID", &trackID);
    ntp_track->SetBranchAddress("pt", &pt);
    ntp_track->SetBranchAddress("nmaps", &nmaps);
    ntp_track->SetBranchAddress("ntpc", &ntpc);
    ntp_track->SetBranchAddress("gtrackID", &gtrackID);
    ntp_track->SetBranchAddress("gpt", &gpt);
    ntp_track->SetBranchAddress("gembed", &gembed);
    ntp_track->GetEntry(ntp_track->GetEntries()-1);
    int max_event = (int)event;

    Long64_t itruth = 0, ireco = 0;
    for(int iev = 0; iev < max_event+1; iev++)
    {
      multimap<int, int> m_gid2id;
      for(; itruth < nen_truth; itruth++)
      {
        ntp_gtrack->GetEntry(itruth);
        if((int)event < iev) continue;
        if((int)event > iev) break;
        if(gtrackID > 0 && gembed > 0)
        {
          h_truth->Fill(gpt);
          if(nmaps > 2 && ntpc > 20)
          {
            h_reco->Fill(gpt);
            h2_resol->Fill(gpt, pt/gpt);
          }
        }
        if(nmaps > 2 && ntpc > 20)
          m_gid2id.insert(make_pair((int)gtrackID, (int)trackID));
      }

      multimap<int, pair<int,Float_t>> m_id2gid;
      for(; ireco < nen_reco; ireco++)
      {
        ntp_track->GetEntry(ireco);
        if((int)event < iev) continue;
        if((int)event > iev) break;
        if(nmaps > 2 && ntpc > 20)
        {
          h_all->Fill(pt);
          m_id2gid.insert(make_pair((int)trackID, make_pair((int)gtrackID, pt)));
        }
      }

      for(const auto &[id, gid_pt] : m_id2gid)
      {
        auto id_range = m_gid2id.equal_range(gid_pt.first);
        for(auto it = id_range.first; it != id_range.second; it++)
          if(it->second == id)
          {
            h_good->Fill(gid_pt.second);
            break;
          }
      }
    } // iev

    delete f;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
