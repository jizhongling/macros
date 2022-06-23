// g++ -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -L$MYINSTALL/lib -L$OFFLINE_MAIN/lib -lg4eval_io -ltrack_io `root-config --cflags --glibs` -o AnaHitFile AnaHitFile.cpp
#include <iostream>
#include <memory>
#include <array>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <g4eval/TrackEvaluationContainerv1.h>

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

  Float_t last_event;
  TNtuple *ntp_cluster = static_cast<TNtuple*>(f->Get("ntp_cluster"));
  ntp_cluster->SetBranchAddress("event", &last_event);
  ntp_cluster->GetEntry(ntp_cluster->GetEntries()-1);

  Float_t g4cluster_event, g4cluster_ID, g4cluster_embed, g4cluster_x, g4cluster_y, g4cluster_z, g4cluster_r;
  TNtuple *ntp_g4cluster = static_cast<TNtuple*>(f->Get("ntp_g4cluster"));
  ntp_g4cluster->SetBranchAddress("event", &g4cluster_event);
  ntp_g4cluster->SetBranchAddress("gtrackID", &g4cluster_ID);
  ntp_g4cluster->SetBranchAddress("gembed", &g4cluster_embed);
  ntp_g4cluster->SetBranchAddress("gx", &g4cluster_x);
  ntp_g4cluster->SetBranchAddress("gy", &g4cluster_y);
  ntp_g4cluster->SetBranchAddress("gz", &g4cluster_z);
  ntp_g4cluster->SetBranchAddress("gr", &g4cluster_r);

  Float_t gtrack_event, gtrack_ID, gtrack_embed, gtrack_px, gtrack_py, gtrack_pz;
  TNtuple *ntp_gtrack = static_cast<TNtuple*>(f->Get("ntp_gtrack"));
  ntp_gtrack->SetBranchAddress("event", &gtrack_event);
  ntp_gtrack->SetBranchAddress("gtrackID", &gtrack_ID);
  ntp_gtrack->SetBranchAddress("gembed", &gtrack_embed);
  ntp_gtrack->SetBranchAddress("gpx", &gtrack_px);
  ntp_gtrack->SetBranchAddress("gpy", &gtrack_py);
  ntp_gtrack->SetBranchAddress("gpz", &gtrack_pz);

  TrackEvaluationContainerv1::TrackStruct::List *v_tracks = nullptr;
  TTree *t_trackeval = static_cast<TTree*>(f->Get("t_trackeval"));
  t_trackeval->SetBranchAddress("tracks", &v_tracks);

  auto v_gtracks = make_unique<TrackEvaluationContainerv1::TrackStruct::List>();

  Int_t max_event = min(static_cast<Int_t>(last_event), static_cast<Int_t>(t_trackeval->GetEntries()-1));

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
  TrkrDefs::cluskey training_cluskey;
  Short_t training_event, training_layer, training_ntouch, training_nedge, li;
  Float_t radius, center_phi, center_z, zr;
  Float_t track_phi, track_z;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, nc> v_reco_phi;
  array<Float_t, nc> v_reco_z;
  array<Short_t, nc> v_reco_adc;
  array<Short_t, nb> v_nreco;
  array<Float_t, nc> v_truth_phi;
  array<Float_t, nc> v_truth_z;
  array<Short_t, nc> v_truth_adc;
  array<Short_t, nb> v_ntruth;

  TTree *t_training = static_cast<TTree*>(f->Get("t_training"));
  t_training->SetBranchAddress("event", &training_event);
  t_training->SetBranchAddress("cluskey", &training_cluskey);
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
  t_out->Branch("reco_phi", &v_reco_phi);
  t_out->Branch("reco_z", &v_reco_z);
  t_out->Branch("reco_adc", &v_reco_adc);
  t_out->Branch("nreco", &v_nreco);
  t_out->Branch("truth_phi", &v_truth_phi);
  t_out->Branch("truth_z", &v_truth_z);
  t_out->Branch("truth_adc", &v_truth_adc);
  t_out->Branch("ntruth", &v_ntruth);
  t_out->Branch("track_phi", &track_phi);
  t_out->Branch("track_z", &track_z);
  t_out->Branch("ntouch", &training_ntouch);

  const Int_t nlayers_map = 3;
  const Int_t nlayers_intt = 4;
  const Int_t nlayers_tpc = 48;

  Long64_t ien_training = 0;
  Long64_t nen_training = t_training->GetEntries();
  Long64_t ien_g4cluster = 0;
  Long64_t nen_g4cluster = ntp_g4cluster->GetEntries();
  Long64_t ien_gtrack = 0;
  Long64_t nen_gtrack = ntp_gtrack->GetEntries();

  for(Int_t event = ith*nev; event < min((ith+1)*nev, max_event+1); event++)
  {
    t_trackeval->GetEntry(event);
    if(!v_tracks)
    {
      cerr << "Event " << event << " has no v_tracks!" << endl;
      continue;
    }

    while(ien_gtrack < nen_gtrack)
    {
      ntp_gtrack->GetEntry(ien_gtrack);
      if(gtrack_event < event)
      {
        ien_gtrack++;
        continue;
      }
      else if(gtrack_event > event)
      {
        break;
      }
      ien_gtrack++;
      if(gtrack_ID <= 0 || gtrack_embed <= 0)
        continue;

      TrackFitUtils::position_vector_t xy_pts;
      TrackFitUtils::position_vector_t rz_pts;
      while(ien_g4cluster < nen_g4cluster)
      {
        ntp_g4cluster->GetEntry(ien_g4cluster);
        if(g4cluster_event < event)
        {
          ien_g4cluster++;
          continue;
        }
        else if(g4cluster_event > event)
        {
          break;
        }
        ien_g4cluster++;
        if(g4cluster_ID == gtrack_ID)
        {
          xy_pts.emplace_back(make_pair(g4cluster_x, g4cluster_y));
          rz_pts.emplace_back(make_pair(g4cluster_r, g4cluster_z));
        }
      } // g4cluster
      if(xy_pts.size() < 3)
        continue;

      auto [cir_R, cir_X0, cir_Y0] = TrackFitUtils::circle_fit_by_taubin(xy_pts);
      auto [lin_k, lin_b] = TrackFitUtils::line_fit(rz_pts);
      Double_t cir_D = sqrt(cir_X0*cir_X0 + cir_Y0*cir_Y0);

      TrackEvaluationContainerv1::TrackStruct gtrack;
      gtrack.px = gtrack_px;
      gtrack.py = gtrack_py;
      gtrack.pz = gtrack_pz;
      gtrack.x = (1 - cir_R / cir_D) * cir_X0;
      gtrack.y = (1 - cir_R / cir_D) * cir_Y0;
      gtrack.z = lin_b;
      gtrack.R = cir_R;
      gtrack.X0 = cir_X0;
      gtrack.Y0 = cir_Y0;
      v_gtracks->emplace_back(gtrack);
    } // gtrack

    vector<TrkrDefs::cluskey> v_cluskey;
    for(const auto &track : *v_tracks)
      for(const auto &cluster : track.clusters)
        v_cluskey.emplace_back(cluster.key);

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
      query(ntp_g4cluster, "gphi:gz:gadc:gprimary:gembed", Form("event==%d && layer==%d", event, layer), v_g4cluster);

      while(ien_training < nen_training)
      {
        t_training->GetEntry(ien_training);
        if(training_event < event || training_layer < layer)
        {
          ien_training++;
          continue;
        }
        else if(training_event > event || training_layer > layer)
        {
          break;
        }
        ien_training++;
        if( find(v_cluskey.begin(), v_cluskey.end(), training_cluskey) == v_cluskey.end() )
          continue;

        zr = center_z / radius;
        v_reco_phi.fill(0.);
        v_reco_z.fill(0.);
        v_reco_adc.fill(0);
        v_nreco.fill(0);
        v_truth_phi.fill(0.);
        v_truth_z.fill(0.);
        v_truth_adc.fill(0);
        v_ntruth.fill(0);
        track_phi = 9999.;
        track_z = 9999.;
        size_t counter;

        counter = 0;
        for(const auto &cluster : v_cluster)
        {
          Float_t phi_diff = (cluster[0] - center_phi) / width_phi[li];
          Float_t z_diff = (cluster[1] - center_z) / width_z;
          Int_t iphi_diff = round(phi_diff);
          Int_t iz_diff = round(z_diff);
          if( abs(iphi_diff) <= nd && abs(iz_diff) <= nd &&
              v_adc[(iphi_diff+nd)*(2*nd+1)+(iz_diff+nd)] > 0 &&
              find(v_searched_cluster.begin(), v_searched_cluster.end(), cluster) == v_searched_cluster.end() )
          {
            size_t ic = min(counter++, nc - 1);
            v_reco_phi[ic] = phi_diff;
            v_reco_z[ic] = z_diff;
            v_reco_adc[ic] = static_cast<Short_t>(cluster[2]);
            size_t ib = min(static_cast<size_t>(floor(cluster[2]/(400./nb))), nb - 1);
            v_nreco[ib]++;
            v_searched_cluster.emplace_back(cluster);
          }
        }

        counter = 0;
        for(const auto &g4cluster : v_g4cluster)
        {
          Float_t phi_diff = (g4cluster[0] - center_phi) / width_phi[li];
          Float_t z_diff = (g4cluster[1] - center_z) / width_z;
          Int_t iphi_diff = round(phi_diff);
          Int_t iz_diff = round(z_diff);
          if( g4cluster[3] > 0 && g4cluster[4] > 0 &&
              abs(iphi_diff) <= nd && abs(iz_diff) <= nd &&
              v_adc[(iphi_diff+nd)*(2*nd+1)+(iz_diff+nd)] > 0 &&
              find(v_searched_g4cluster.begin(), v_searched_g4cluster.end(), g4cluster) == v_searched_g4cluster.end() )
          {
            size_t ic = min(counter++, nc - 1);
            v_truth_phi[ic] = phi_diff;
            v_truth_z[ic] = z_diff;
            v_truth_adc[ic] = static_cast<Short_t>(g4cluster[2]);
            size_t ib = min(static_cast<size_t>(floor(g4cluster[2]/(400./nb))), nb - 1);
            v_ntruth[ib]++;
            v_searched_g4cluster.emplace_back(g4cluster);
          }
        }

        Float_t min_dist2 = 9999.;
        for(const auto &track : *v_gtracks)
        {
          Float_t trk_px = track.px;
          Float_t trk_py = track.py;
          Float_t trk_pz = track.pz;
          Float_t dca_x = track.x;
          Float_t dca_y = track.y;
          Float_t dca_z = track.z;
          Float_t cir_R = track.R;
          Float_t cir_X0 = track.X0;
          Float_t cir_Y0 = track.Y0;

          // Circle intersection
          // math.stackexchange.com/questions/256100/how-can-i-find-the-points-at-which-two-circles-intersect?newreg=b9dc98e45b514173ae85b6dbaf4d2508
          Float_t cir_D = sqrt(cir_X0*cir_X0 + cir_Y0*cir_Y0);
          if(cir_D < fabs(cir_R - radius) || cir_D > fabs(cir_R + radius)) continue;
          Float_t cir_a = (cir_R*cir_R - radius*radius) / (cir_D*cir_D);
          Float_t cir_b = sqrt(2 * (cir_R*cir_R + radius*radius) / (cir_D*cir_D) - square((cir_R*cir_R - radius*radius) / (cir_D*cir_D)) - 1);
          Float_t cir_sign = trk_px*cir_Y0 - trk_py*cir_X0 > 0 ? 1 : -1;
          Float_t sec_x = (1 - cir_a) * cir_X0/2 + cir_sign * cir_b * cir_Y0/2;
          Float_t sec_y = (1 - cir_a) * cir_Y0/2 - cir_sign * cir_b * cir_X0/2;

          Float_t sec_phi = atan2(sec_x-cir_X0, sec_y-cir_Y0);
          Float_t dca_phi = atan2(dca_x-cir_X0, dca_y-cir_Y0);
          Float_t dphi = sec_phi - dca_phi;
          while(dphi >= PI) dphi -= 2*PI;
          while(dphi < -PI) dphi += 2*PI;
          Float_t dt = fabs(dphi * cir_R) / sqrt(trk_px*trk_px + trk_py*trk_py);

          Float_t trk_phi = atan2(sec_y, sec_x);
          Float_t trk_z = dca_z + trk_pz * dt;

          Float_t phi_diff = (trk_phi - center_phi) / width_phi[li];
          Float_t z_diff = (trk_z - center_z) / width_z;
          Float_t dist2 = square(phi_diff) + square(z_diff);
          if(dist2 < min_dist2)
          {
            min_dist2 = dist2;
            track_phi = phi_diff;
            track_z = z_diff;
          }
        } // v_gtracks

        if(min_dist2 < 1000.)
          t_out->Fill();
      } // t_training
    } // layer
  } // event

  f_out->Write();
  f_out->Close();
  delete f;
  return 0;
}
