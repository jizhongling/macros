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

  Float_t last_event;
  ntp_cluster->SetBranchAddress("event", &last_event);
  ntp_cluster->GetEntry(ntp_cluster->GetEntries()-1);
  Int_t max_event = static_cast<Int_t>(last_event);

  const Int_t nev = 10;
  const Int_t ith = stoi(string(argv[3]));
  if(ith*nev > max_event)
  {
    delete f;
    return 0;
  }

  const size_t nd = 5;
  const size_t nc = 10;
  const size_t nb = 20;
  Short_t training_event, training_layer, li;
  Float_t center_phi, center_z, zr;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, nc> v_reco_phi;
  array<Float_t, nc> v_reco_z;
  array<Short_t, nc> v_reco_adc;
  array<Short_t, nb> v_nreco;
  array<Float_t, nc> v_truth_phi;
  array<Float_t, nc> v_truth_z;
  array<Float_t, nc> v_truth_phicov;
  array<Float_t, nc> v_truth_zcov;
  array<Short_t, nc> v_truth_adc;
  array<Short_t, nb> v_ntruth;

  TTree *t_training = static_cast<TTree*>(f->Get("t_training"));
  t_training->SetBranchAddress("event", &training_event);
  t_training->SetBranchAddress("layer", &training_layer);
  t_training->SetBranchAddress("phi", &center_phi);
  t_training->SetBranchAddress("z", &center_z);
  t_training->SetBranchAddress("adc", &v_adc);

  auto f_out = new TFile(Form("%s-%d.root", argv[2], ith), "RECREATE");
  auto t_out = new TTree("T", "Training data");
  t_out->Branch("layer", &li);
  t_out->Branch("ztan", &zr);
  t_out->Branch("adc", &v_adc);
  t_out->Branch("reco_phi", &v_reco_phi);
  t_out->Branch("reco_z", &v_reco_z);
  t_out->Branch("reco_adc", &v_reco_adc);
  t_out->Branch("nreco", &v_nreco);
  t_out->Branch("truth_phi", &v_truth_phi);
  t_out->Branch("truth_z", &v_truth_z);
  t_out->Branch("truth_phicov", &v_truth_phicov);
  t_out->Branch("truth_zcov", &v_truth_zcov);
  t_out->Branch("truth_adc", &v_truth_adc);
  t_out->Branch("ntruth", &v_ntruth);

  const Int_t nlayers_map = 3;
  const Int_t nlayers_intt = 4;
  const Int_t nlayers_tpc = 48;

  Long64_t ien = 0;
  Long64_t nentries = t_training->GetEntries();

  for(Int_t event = ith*nev; event < min((ith+1)*nev, max_event+1); event++)
    for(Int_t layer = nlayers_map + nlayers_intt; layer < nlayers_map + nlayers_intt + nlayers_tpc; layer++)
    {
      if(layer < nlayers_map + nlayers_intt + nlayers_tpc/3)
        li = 0;
      else if(layer < nlayers_map + nlayers_intt + nlayers_tpc*2/3)
        li = 1;
      else
        li = 2;

      const Float_t PI = TMath::Pi();
      const Float_t radius[3] = {(30.+40.)/2., (40.+60.)/2., (60.+77.)/2.};
      const Float_t width_phi[3] = {2*PI/1152, 2*PI/1536, 2*PI/2304};
      const Float_t width_z = 53. * 8. / 1000.;
      const Float_t region_phi = nd * width_phi[li];
      const Float_t region_z = nd * width_z;

      vvF v_cluster;
      vvF v_g4cluster;
      vvF v_searched_cluster;
      vvF v_searched_g4cluster;

      query(ntp_cluster, "phi:z:adc", Form("event==%d && layer==%d", event, layer), v_cluster);
      query(ntp_g4cluster, "gphi:gz:gadc:gvr:ephi:ez", Form("event==%d && layer==%d", event, layer), v_g4cluster);

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

        zr = center_z / radius[li];
        v_reco_phi.fill(0.);
        v_reco_z.fill(0.);
        v_reco_adc.fill(0);
        v_nreco.fill(0);
        v_truth_phi.fill(0.);
        v_truth_z.fill(0.);
        v_truth_phicov.fill(0.);
        v_truth_zcov.fill(0.);
        v_truth_adc.fill(0);
        v_ntruth.fill(0);
        size_t counter;

        counter = 0;
        for(const auto &cluster : v_cluster)
          if( fabs(cluster[0] - center_phi) < region_phi &&
              fabs(cluster[1] - center_z) < region_z &&
              find(v_searched_cluster.begin(), v_searched_cluster.end(), cluster) == v_searched_cluster.end() )
          {
            size_t ic = min(counter++, nc - 1);
            v_reco_phi[ic] = (cluster[0] - center_phi) / width_phi[li];
            v_reco_z[ic] = (cluster[1] - center_z) / width_z;
            v_reco_adc[ic] = static_cast<Short_t>(cluster[2]);
            size_t ib = min(static_cast<size_t>(floor(cluster[2]/(400./nb))), nb - 1);
            v_nreco[ib]++;
            v_searched_cluster.emplace_back(cluster);
          }

        counter = 0;
        for(const auto &g4cluster : v_g4cluster)
          if( g4cluster[2] < 1000. && g4cluster[3] < 25. &&
              g4cluster[4] < 1. && g4cluster[5] < 1. &&
              fabs(g4cluster[0] - center_phi) < region_phi &&
              fabs(g4cluster[1] - center_z) < region_z &&
              find(v_searched_g4cluster.begin(), v_searched_g4cluster.end(), g4cluster) == v_searched_g4cluster.end() )
          {
            size_t ic = min(counter++, nc - 1);
            v_truth_phi[ic] = (g4cluster[0] - center_phi) / width_phi[li];
            v_truth_z[ic] = (g4cluster[1] - center_z) / width_z;
            v_truth_phicov[ic] = g4cluster[4] * g4cluster[4];
            v_truth_zcov[ic] = g4cluster[5] * g4cluster[5];
            v_truth_adc[ic] = static_cast<Short_t>(g4cluster[2]);
            size_t ib = min(static_cast<size_t>(floor(g4cluster[2]/(400./nb))), nb - 1);
            v_ntruth[ib]++;
            v_searched_g4cluster.emplace_back(g4cluster);
          }

        t_out->Fill();
      } // t_training
    } // event, layer

  f_out->Write();
  f_out->Close();
  delete f;
  return 0;
}
