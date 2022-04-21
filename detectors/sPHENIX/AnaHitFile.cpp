// g++ -Wall -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib `root-config --cflags --glibs` -o AnaHitFile AnaHitFile.cpp
#include <iostream>
#include <array>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

using namespace std;

typedef vector<vector<Float_t>> vvF;
template<class T> inline constexpr T square(const T &x) { return x*x; }
template<class T> inline constexpr T get_r(const T &x, const T &y) { return sqrt(square(x) + square(y)); }

void query(TNtuple *ntp, TString var, const char *cond, vvF &result)
{
  const Int_t nitem = var.CountChar(58) + 1;
  TSQLResult *res = ntp->Query(var, cond);
  TSQLRow *row;

  while( (row = res->Next()) )
  {
    vector<Float_t> values;
    bool store = true;

    for(Int_t i=0; i<nitem; i++)
    {
      TString field = row->GetField(i);
      Float_t value = field.Atof();
      if( !TMath::Finite(value) )
      {
        store = false;
        break;
      }
      values.emplace_back(value);
    }

    if(store)
      result.emplace_back(values);

    delete row;
  }

  delete res;
  return;
}

int main(int argc, const char *argv[])
{
  if(argc != 4)
  {
    cerr << "Usage: " << argv[0] << " <eval.root> <output>(-<ith>.root) <ith>" << endl;
    return 1;
  }

  auto f = new TFile(argv[1]);
  if(f->IsZombie())
  {
    cerr << "Error: cannot open file " << argv[1] << endl;
    return 1;
  }
  cout << "Open file " << argv[1] << endl;

  TNtuple *ntp_cluster = static_cast<TNtuple*>(f->Get("ntp_cluster"));
  TNtuple *ntp_g4cluster = static_cast<TNtuple*>(f->Get("ntp_g4cluster"));
  TNtuple *ntp_track = static_cast<TNtuple*>(f->Get("ntp_track"));

  Float_t last_event;
  ntp_cluster->SetBranchAddress("event", &last_event);
  ntp_cluster->GetEntry(ntp_cluster->GetEntries()-1);
  Int_t max_event = static_cast<Int_t>(last_event);

  const Int_t nev = 5;
  const Int_t ith = stoi(string(argv[3]));
  if(ith*nev > max_event)
  {
    delete f;
    return 0;
  }

  const size_t nd = 5;
  const size_t nc = 10;
  const size_t nb = 20;
  Short_t training_event, training_layer, training_ntouch, training_nedge, li;
  Float_t radius, center_phi, center_z, zr;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, nc> v_reco_rphi;
  array<Float_t, nc> v_reco_z;
  array<Short_t, nc> v_reco_adc;
  array<Short_t, nb> v_nreco;
  array<Float_t, nc> v_truth_rphi;
  array<Float_t, nc> v_truth_z;
  array<Short_t, nc> v_truth_adc;
  array<Short_t, nb> v_ntruth;
  Float_t track_rphi;
  Float_t track_z;

  TTree *t_training = static_cast<TTree*>(f->Get("t_training"));
  t_training->SetBranchAddress("event", &training_event);
  t_training->SetBranchAddress("layer", &training_layer);
  t_training->SetBranchAddress("ntouch", &training_ntouch);
  t_training->SetBranchAddress("nedge", &training_nedge);
  t_training->SetBranchAddress("radius", &radius);
  t_training->SetBranchAddress("phi", &center_phi);
  t_training->SetBranchAddress("z", &center_z);
  t_training->SetBranchAddress("adc", &v_adc);

  auto f_out = new TFile(Form("%s-%d.root", argv[2], ith), "RECREATE");
  auto t_out = new TTree("T", "Training data");
  t_out->Branch("layer", &li);
  t_out->Branch("ztan", &zr);
  t_out->Branch("adc", &v_adc);
  t_out->Branch("reco_rphi", &v_reco_rphi);
  t_out->Branch("reco_z", &v_reco_z);
  t_out->Branch("reco_adc", &v_reco_adc);
  t_out->Branch("nreco", &v_nreco);
  t_out->Branch("truth_rphi", &v_truth_rphi);
  t_out->Branch("truth_z", &v_truth_z);
  t_out->Branch("truth_adc", &v_truth_adc);
  t_out->Branch("ntruth", &v_ntruth);
  t_out->Branch("track_rphi", &track_rphi);
  t_out->Branch("track_z", &track_z);

  const Int_t nlayers_map = 3;
  const Int_t nlayers_intt = 4;
  const Int_t nlayers_tpc = 48;

  Long64_t ien = 0;
  Long64_t nentries = t_training->GetEntries();

  for(Int_t event = ith*nev; event < min((ith+1)*nev, max_event+1); event++)
  {
    vvF v_track;
    query(ntp_track, "px:py:pz:pcax:pcay:pcaz", Form("event==%d", event), v_track);

    for(Int_t layer = nlayers_map + nlayers_intt; layer < nlayers_map + nlayers_intt + nlayers_tpc; layer++)
    {
      if(layer < nlayers_map + nlayers_intt + nlayers_tpc/3)
        li = 0;
      else if(layer < nlayers_map + nlayers_intt + nlayers_tpc*2/3)
        li = 1;
      else
        li = 2;

      const Float_t PI = TMath::Pi();
      const Float_t width_phi[3] = {2*PI/1152, 2*PI/1536, 2*PI/2304};
      const Float_t width_z = 53. * 8. / 1000.;

      vvF v_cluster;
      vvF v_g4cluster;
      vvF v_searched_cluster;
      vvF v_searched_g4cluster;

      query(ntp_cluster, "phi:z:adc", Form("event==%d && layer==%d", event, layer), v_cluster);
      query(ntp_g4cluster, "gphi:gz:gadc:gvr", Form("event==%d && layer==%d", event, layer), v_g4cluster);

      while(ien < nentries)
      {
        t_training->GetEntry(ien);
        if(training_event < event || training_layer < layer)
        {
          ien++;
          continue;
        }
        else if(training_event > event || training_layer > layer)
        {
          break;
        }
        ien++;

        zr = center_z / radius;
        v_reco_rphi.fill(0.);
        v_reco_z.fill(0.);
        v_reco_adc.fill(0);
        v_nreco.fill(0);
        v_truth_rphi.fill(0.);
        v_truth_z.fill(0.);
        v_truth_adc.fill(0);
        v_ntruth.fill(0);
        track_rphi = 9999.;
        track_z = 9999.;
        size_t counter;

        counter = 0;
        for(const auto &cluster : v_cluster)
        {
          int iphi_diff = round((cluster[0] - center_phi) / width_phi[li]);
          int iz_diff = round((cluster[1] - center_z) / width_z);
          if( abs(iphi_diff) <= nd && abs(iz_diff) <= nd &&
              v_adc[(iphi_diff+nd)*(2*nd+1)+(iz_diff+nd)] > 0 &&
              find(v_searched_cluster.begin(), v_searched_cluster.end(), cluster) == v_searched_cluster.end() )
          {
            size_t ic = min(counter++, nc - 1);
            v_reco_rphi[ic] = radius * (cluster[0] - center_phi);
            v_reco_z[ic] = cluster[1] - center_z;
            v_reco_adc[ic] = static_cast<Short_t>(cluster[2]);
            size_t ib = min(static_cast<size_t>(floor(cluster[2]/(400./nb))), nb - 1);
            v_nreco[ib]++;
            v_searched_cluster.emplace_back(cluster);
          }
        }

        counter = 0;
        for(const auto &g4cluster : v_g4cluster)
        {
          int iphi_diff = round((g4cluster[0] - center_phi) / width_phi[li]);
          int iz_diff = round((g4cluster[1] - center_z) / width_z);
          if( g4cluster[3] < 25. &&
              abs(iphi_diff) <= nd && abs(iz_diff) <= nd &&
              v_adc[(iphi_diff+nd)*(2*nd+1)+(iz_diff+nd)] > 0 &&
              find(v_searched_g4cluster.begin(), v_searched_g4cluster.end(), g4cluster) == v_searched_g4cluster.end() )
          {
            size_t ic = min(counter++, nc - 1);
            v_truth_rphi[ic] = radius * (g4cluster[0] - center_phi);
            v_truth_z[ic] = g4cluster[1] - center_z;
            v_truth_adc[ic] = static_cast<Short_t>(g4cluster[2]);
            size_t ib = min(static_cast<size_t>(floor(g4cluster[2]/(400./nb))), nb - 1);
            v_ntruth[ib]++;
            v_searched_g4cluster.emplace_back(g4cluster);
          }
        }

        if(training_ntouch > 0)
        {
          float min_dist2 = 9999.;
          for(const auto &track : v_track)
          {
            const auto px = track[0];
            const auto py = track[1];
            const auto pz = track[2];
            const auto x = track[3];
            const auto y = track[4];
            const auto z = track[5];

            const auto r = get_r(x, y);
            const auto dr = radius - r;
            const auto drdt = (x * px + y * py) / r;
            const auto dxdr = px / drdt;
            const auto dydr = py / drdt;
            const auto dzdr = pz / drdt;

            const auto trk_x = x + dr * dxdr;
            const auto trk_y = y + dr * dydr;
            const auto trk_z = z + dr * dzdr;
            const auto trk_r = get_r(trk_x, trk_y);
            const auto trk_phi = atan2(trk_y, trk_x);
            const auto trk_rphi = trk_r * trk_phi;
            const auto center_rphi = radius * center_phi;

            float dist2 = square(trk_rphi - center_rphi) + square(trk_z - center_z);
            if(dist2 < min_dist2)
            {
              min_dist2 = dist2;
              track_rphi = trk_rphi - center_rphi;
              track_z = trk_z - center_z;
            }
          }
        }

        t_out->Fill();
      } // t_training
    } // layer
  } // event

  f_out->Write();
  f_out->Close();
  delete f;
  return 0;
}
