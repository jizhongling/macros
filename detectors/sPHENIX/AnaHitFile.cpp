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
  const Short_t nitem = var.CountChar(58) + 1;
  TSQLResult *res = ntp->Query(var, cond);
  TSQLRow *row;

  while( (row = res->Next()) )
  {
    vector<Float_t> values;
    bool store = true;

    for(Short_t i=0; i<nitem; i++)
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
  if(argc != 3)
  {
    cerr << "Usage: " << argv[0] << " <eval.root> <output.root>" << endl;
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
  TNtuple *ntp_hit = static_cast<TNtuple*>(f->Get("ntp_hit"));
  TNtuple *ntp_g4cluster = static_cast<TNtuple*>(f->Get("ntp_g4cluster"));

  Float_t last_event;
  ntp_cluster->SetBranchAddress("event", &last_event);
  ntp_cluster->GetEntry(ntp_cluster->GetEntries()-1);

  const size_t nd = 5;
  const size_t nin = (2*nd+1)*(2*nd+1) + 2;  // adc + layer, z/r*1000
  array<Short_t, nin> v_in;
  vector<Short_t> v_truth_bin;
  Short_t ntruth;

  auto f_out = new TFile(argv[2], "RECREATE");
  auto t_out = new TTree("T", "Training data");
  t_out->Branch("input", &v_in);
  t_out->Branch("truth_bin", &v_truth_bin);
  t_out->Branch("ntruth", &ntruth);

  const Short_t nlayers_map = 3;
  const Short_t nlayers_intt = 4;
  const Short_t nlayers_tpc = 48;

  for(Short_t event = 0; event < static_cast<Short_t>(last_event); event++)
    for(Short_t layer = nlayers_map + nlayers_intt; layer < nlayers_map + nlayers_intt + nlayers_tpc; layer++)
    {
      Short_t li;
      if(layer < nlayers_map + nlayers_intt + nlayers_tpc/3)
        li = 0;
      else if(layer < nlayers_map + nlayers_intt + nlayers_tpc*2/3)
        li = 1;
      else
        li = 2;

      vvF v_cluster;
      vvF v_searched;
      vvF v_hit;
      vvF v_g4cluster;

      query(ntp_cluster, "phi:z", Form("event==%hd&&layer==%hd", event, layer), v_cluster);
      query(ntp_hit, "phi:z:adc", Form("event==%hd&&layer==%hd", event, layer), v_hit);
      query(ntp_g4cluster, "gphi:gz", Form("event==%hd&&layer==%hd", event, layer), v_g4cluster);

      for(const auto &cluster : v_cluster)
        if( find(v_searched.begin(), v_searched.end(), cluster) == v_searched.end() )
        {
          const Float_t PI = TMath::Pi();
          const Float_t radius[3] = {(30.+40.)/2., (40.+60.)/2., (60.+77.)/2.};
          const Float_t width_phi[3] = {2*PI/1152, 2*PI/1536, 2*PI/2304};
          const Float_t width_z = 53. * 8. / 1000.;
          const Float_t region_phi = nd * width_phi[li];
          const Float_t region_z = nd * width_z;

          v_in.fill(0);
          v_truth_bin.clear();
          ntruth = 0;

          for(const auto &searched : v_cluster)
            if( fabs(searched[0] - cluster[0]) < region_phi && fabs(searched[1] - cluster[1]) < region_z )
              v_searched.emplace_back(searched);

          for(const auto &hit : v_hit)
            if( fabs(hit[0] - cluster[0]) < region_phi && fabs(hit[1] - cluster[1]) < region_z )
            {
              Float_t bin_phi = round((hit[0] - cluster[0]) / width_phi[li]);
              Float_t bin_z = round((hit[1] - cluster[1]) / width_z);
              size_t bin_i = static_cast<size_t>((bin_phi+nd)*(2*nd+1) + (bin_z+nd));
              v_in[bin_i] = static_cast<Short_t>(hit[2]);
            }
          v_in[nin-2] = layer;
          v_in[nin-1] = static_cast<Short_t>(round(cluster[1] / radius[li] * 1000.));

          for(const auto &g4cluster : v_g4cluster)
            if( fabs(g4cluster[0] - cluster[0]) < region_phi && fabs(g4cluster[1] - cluster[1]) < region_z )
            {
              Float_t bin_phi = round((g4cluster[0] - cluster[0]) / width_phi[li]);
              Float_t bin_z = round((g4cluster[1] - cluster[1]) / width_z);
              Short_t bin_i = static_cast<Short_t>((bin_phi+nd)*(2*nd+1) + (bin_z+nd));
              v_truth_bin.emplace_back(bin_i);
            }
          ntruth = v_truth_bin.size();

          t_out->Fill();
        } // cluster
    } // event, layer

  f_out->Write();
  f_out->Close();
  delete f;
  return 0;
}
