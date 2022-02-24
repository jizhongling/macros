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

  TNtuple *ntp_hit = static_cast<TNtuple*>(f->Get("ntp_hit"));
  TNtuple *ntp_cluster = static_cast<TNtuple*>(f->Get("ntp_cluster"));
  TNtuple *ntp_g4cluster = static_cast<TNtuple*>(f->Get("ntp_g4cluster"));

  Float_t last_event;
  ntp_cluster->SetBranchAddress("event", &last_event);
  ntp_cluster->GetEntry(ntp_cluster->GetEntries()-1);
  Int_t max_event = static_cast<Int_t>(last_event);

  const Int_t nev = 5;
  const Int_t ith = stoi(string(argv[3]));
  if(ith*nev >= max_event)
  {
    delete f;
    return 0;
  }

  const size_t nd = 5;
  const size_t nc = 10;
  const size_t nb = 20;
  Short_t li;
  Float_t zr;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, nc> v_reco_phi;
  array<Float_t, nc> v_reco_z;
  array<Short_t, nc> v_reco_adc;
  array<Short_t, nb> v_nreco;
  array<Float_t, nc> v_truth_phi;
  array<Float_t, nc> v_truth_z;
  array<Float_t, nc> v_truth_edep;
  array<Short_t, nb> v_ntruth;

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
  t_out->Branch("truth_edep", &v_truth_edep);
  t_out->Branch("ntruth", &v_ntruth);

  const Int_t nlayers_map = 3;
  const Int_t nlayers_intt = 4;
  const Int_t nlayers_tpc = 48;

  for(Int_t event = ith*nev; event < min((ith+1)*nev, max_event); event++)
    for(Int_t layer = nlayers_map + nlayers_intt; layer < nlayers_map + nlayers_intt + nlayers_tpc; layer++)
    {
      if(layer < nlayers_map + nlayers_intt + nlayers_tpc/3)
        li = 0;
      else if(layer < nlayers_map + nlayers_intt + nlayers_tpc*2/3)
        li = 1;
      else
        li = 2;

      vvF v_hit;
      vvF v_cluster;
      vvF v_g4cluster;
      vvF v_searched;

      query(ntp_hit, "phi:z:adc", Form("event==%d && layer==%d", event, layer), v_hit);
      query(ntp_cluster, "phi:z:adc", Form("event==%d && layer==%d", event, layer), v_cluster);
      query(ntp_g4cluster, "gphi:gz:gedep:gprimary", Form("event==%d && layer==%d", event, layer), v_g4cluster);

      for(const auto &cluster : v_cluster)
        if( find(v_searched.begin(), v_searched.end(), cluster) == v_searched.end() )
        {
          const Float_t PI = TMath::Pi();
          const Float_t radius[3] = {(30.+40.)/2., (40.+60.)/2., (60.+77.)/2.};
          const Float_t width_phi[3] = {2*PI/1152, 2*PI/1536, 2*PI/2304};
          const Float_t width_z = 53. * 8. / 1000.;
          const Float_t region_phi = nd * width_phi[li];
          const Float_t region_z = nd * width_z;

          zr = cluster[1] / radius[li];
          v_adc.fill(0);
          v_reco_phi.fill(0.);
          v_reco_z.fill(0.);
          v_reco_adc.fill(0);
          v_nreco.fill(0);
          v_truth_phi.fill(0.);
          v_truth_z.fill(0.);
          v_truth_edep.fill(0.);
          v_ntruth.fill(0);
          size_t counter;

          for(const auto &hit : v_hit)
            if( fabs(hit[0] - cluster[0]) < region_phi &&
                fabs(hit[1] - cluster[1]) < region_z )
            {
              Float_t bin_phi = round((hit[0] - cluster[0]) / width_phi[li]);
              Float_t bin_z = round((hit[1] - cluster[1]) / width_z);
              size_t bin_i = static_cast<size_t>((bin_phi+nd)*(2*nd+1) + (bin_z+nd));
              v_adc[bin_i] = static_cast<Short_t>(hit[2]);
            }

          counter = 0;
          for(const auto &searched : v_cluster)
            if( fabs(searched[0] - cluster[0]) < region_phi &&
                fabs(searched[1] - cluster[1]) < region_z )
            {
              size_t ic = min(counter++, nc - 1);
              v_reco_phi[ic] = (searched[0] - cluster[0]) / width_phi[li];
              v_reco_z[ic] = (searched[1] - cluster[1]) / width_z;
              v_reco_adc[ic] = static_cast<Short_t>(searched[2]);
              size_t ib = min(static_cast<size_t>(floor(searched[2]/(400./nb))), nb - 1);
              v_nreco[ib]++;
              v_searched.emplace_back(searched);
            }

          counter = 0;
          for(const auto &g4cluster : v_g4cluster)
            if( g4cluster[3] > 0.5 &&
                fabs(g4cluster[0] - cluster[0]) < region_phi &&
                fabs(g4cluster[1] - cluster[1]) < region_z )
            {
              size_t ic = min(counter++, nc - 1);
              v_truth_phi[ic] = (g4cluster[0] - cluster[0]) / width_phi[li];
              v_truth_z[ic] = (g4cluster[1] - cluster[1]) / width_z;
              v_truth_edep[ic] = g4cluster[2];
              size_t ib = min(static_cast<size_t>(floor(g4cluster[2]/(4e-6/nb))), nb - 1);
              v_ntruth[ib]++;
            }

          t_out->Fill();
        } // cluster
    } // event, layer

  f_out->Write();
  f_out->Close();
  delete f;
  return 0;
}
