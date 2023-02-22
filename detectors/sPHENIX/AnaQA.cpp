// g++ -std=c++17 -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib -lg4eval_io `root-config --cflags --glibs` -o AnaQA AnaQA.cpp
#include <iostream>
#include <tuple>
#include <map>
#include <algorithm>

#include <TDirectory.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <g4eval/TrackEvaluationContainerv1.h>

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
  auto h2_truth = new TH2F("h2_truth","Truth tracks;p_{T} (GeV);Ntrack", 200,0.,20., 100,0.,1e5);
  auto h2_reco = new TH2F("h2_reco","Reco tracks;p_{T} (GeV);Ntrack", 200,0.,20., 100,0.,1e5);
  auto h2_all = new TH2F("h2_all","All reco tracks;p_{T} (GeV);Ntrack", 40,0.,20., 100,0.,1e5);
  auto h2_good = new TH2F("h2_good","Reco tracks with correct truth association;p_{T} (GeV);Ntrack", 40,0.,20., 100,0.,1e5);
  auto h3_resol = new TH3F("h3_resol","Momentum resolution;p_{T}^{truth} (GeV);p_{T}^{reco}/p_{T}^{truth};Ntrack", 40,0.,20., 400,0.5,1.5, 100,0.,1e5);
  auto h3_eta = new TH3F("h3_eta","#eta resolution;p_{T}^{truth} (GeV);#Delta#eta;Ntrack", 40,0.,20., 100,-0.01,0.01, 100,0.,1e5);
  auto h3_phi = new TH3F("h3_phi","#phi resolution;p_{T}^{truth} (GeV);#Delta#phi;Ntrack", 40,0.,20., 100,-0.01,0.01, 100,0.,1e5);

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
    Float_t gtrackID, gpt, geta, gphi, gembed;
    Float_t trackID, pt, eta, phi, nmaps, ntpc;

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
    ntp_gtrack->SetBranchAddress("geta", &geta);
    ntp_gtrack->SetBranchAddress("gphi", &gphi);
    ntp_gtrack->SetBranchAddress("gembed", &gembed);
    ntp_gtrack->SetBranchAddress("trackID", &trackID);
    ntp_gtrack->SetBranchAddress("pt", &pt);
    ntp_gtrack->SetBranchAddress("eta", &eta);
    ntp_gtrack->SetBranchAddress("phi", &phi);
    ntp_gtrack->SetBranchAddress("nmaps", &nmaps);
    ntp_gtrack->SetBranchAddress("ntpc", &ntpc);

    ntp_gtrack->Draw("event >> h_ntrack(20,-0.5,19.5)", "gpt>0.1 && fabs(geta)<1");
    TH1 *h_ntrack = (TH1*)gDirectory->Get("h_ntrack");

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
    Int_t last_event = static_cast<Int_t>(event);

    TrackEvaluationContainerv1::TrackStruct::List *v_tracks = nullptr;
    TTree *t_trackeval = static_cast<TTree*>(f->Get("t_trackeval"));
    t_trackeval->SetBranchAddress("tracks", &v_tracks);

    Int_t max_event = min(static_cast<Int_t>(last_event), static_cast<Int_t>(t_trackeval->GetEntries()-1));

    Long64_t itruth = 0, ireco = 0;
    for(int iev = 0; iev < max_event+1; iev++)
    {
      Double_t ntrack = h_ntrack->GetBinContent(iev+1);
      multimap<int, tuple<float,float,float>> m_id;

      t_trackeval->GetEntry(iev);
      for(const auto &track : *v_tracks)
        if(track.gtrackID > 0 && track.embed > 0)
          m_id.insert(make_pair(track.trackID, make_tuple(track.truth_pt, track.truth_eta, track.truth_phi)));

      for(; itruth < nen_truth; itruth++)
      {
        ntp_gtrack->GetEntry(itruth);
        if((int)event < iev) continue;
        if((int)event > iev) break;
        if(gtrackID > 0 && gembed > 0)
          h2_truth->Fill(gpt, ntrack);
      } // itruth

      for(; ireco < nen_reco; ireco++)
      {
        ntp_track->GetEntry(ireco);
        if((int)event < iev) continue;
        if((int)event > iev) break;
        if(nmaps > 2 && ntpc > 20)
        {
          h2_all->Fill(pt, ntrack);
          const auto &it = m_id.find((int)trackID);
          if(it != m_id.end())
          {
            Float_t gpt = get<0>(it->second);
            Float_t geta = get<1>(it->second);
            Float_t gphi = get<2>(it->second);

            h2_good->Fill(pt, ntrack);
            h2_reco->Fill(gpt, ntrack);
            h3_resol->Fill(gpt, pt/gpt, ntrack);
            h3_eta->Fill(gpt, eta-geta, ntrack);
            h3_phi->Fill(gpt, phi-gphi, ntrack);
          } // matched reco
        } // all reco
      } // ireco
    } // iev

    delete h_ntrack;
    delete f;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
