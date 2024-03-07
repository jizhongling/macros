// g++ -Wall -I$MYINSTALL/include -I$OFFLINE_MAIN/include -L$MYINSTALL/lib -L$OFFLINE_MAIN/lib -lg4eval_io -ltrack_io `root-config --cflags --glibs` -o AnaEEMass AnaEEMass.cpp
#include <iostream>
#include <set>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH3.h>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <g4eval/TrackEvaluationContainerv1.h>

using namespace std;

template<class T> inline constexpr T square(const T &x) { return x*x; }
template<class T> inline constexpr T get_r(const T &x, const T &y) { return sqrt(square(x) + square(y)); }
template <class T> inline ROOT::Math::PxPyPzEVector get_pe(const T &a, const T &b, const T &c)
{
  T d = sqrt(a*a + b*b + c*c);
  return ROOT::Math::PxPyPzEVector(a, b, c, d);
}

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
  const char *ith = argv[5];

  auto f_out = new TFile(Form("%s-%s.root", output, ith), "RECREATE");
  auto h3_minv = new TH3F("h3_minv", "Invariant mass;z_{1} (cm);z_{2} (cm);m_{inv} (GeV)", 100,-50.,50., 100,-50.,50., 200,0.,0.4);
  auto h3_e_ratio = new TH3F("h3_e_ratio", "Energy ratio;etabin;hcalin/cemc;cemc/mom", 80,-0.5,79.5, 20,0.,0.2, 300,0.,3.);

  for(Int_t ifile = istart; ifile < iend; ifile++)
  {
    char filename[1024];
    sprintf(filename, "%s/G4sPHENIX_g4svtx_eval-%d.root", indir, ifile);
    auto f_in = new TFile(filename);
    if(!f_in || f_in->IsZombie())
    {
      cerr << "Error: cannot open file " << filename << endl;
      continue;
    }

    TrackEvaluationContainerv1::TrackStruct::List *v_tracks = nullptr;
    TTree *t_trackeval = static_cast<TTree*>(f_in->Get("t_trackeval"));
    t_trackeval->SetBranchAddress("tracks", &v_tracks);

    for(Long64_t event = 0; event < t_trackeval->GetEntries(); event++)
    {
      t_trackeval->GetEntry(event);
      if(!v_tracks)
      {
        cerr << "Event " << event << " has no v_tracks!" << endl;
        continue;
      }

      set<size_t> filled;
      for(size_t itrk = 0; itrk < v_tracks->size(); itrk++)
      {
        Int_t charge1 = v_tracks->at(itrk).charge;
        Float_t cir1_R = v_tracks->at(itrk).R;
        Float_t cir1_X0 = v_tracks->at(itrk).X0;
        Float_t cir1_Y0 = v_tracks->at(itrk).Y0;
        Float_t trk1_x = v_tracks->at(itrk).x;
        Float_t trk1_y = v_tracks->at(itrk).y;
        Float_t trk1_z = v_tracks->at(itrk).z;
        Float_t trk1_px = v_tracks->at(itrk).px;
        Float_t trk1_py = v_tracks->at(itrk).py;
        Float_t trk1_pz = v_tracks->at(itrk).pz;
        Float_t trk1_pt = v_tracks->at(itrk).pt;
        Float_t trk1_p = v_tracks->at(itrk).p;
        auto trk1_pe = get_pe<Float_t>(trk1_px, trk1_py, trk1_pz);

        for(size_t jtrk = itrk+1; jtrk < v_tracks->size(); jtrk++)
        {
          Int_t charge2 = v_tracks->at(jtrk).charge;
          if(charge1*charge2 != -1) continue;
          Float_t cir2_R = v_tracks->at(jtrk).R;
          Float_t cir2_X0 = v_tracks->at(jtrk).X0;
          Float_t cir2_Y0 = v_tracks->at(jtrk).Y0;
          Float_t trk2_x = v_tracks->at(jtrk).x;
          Float_t trk2_y = v_tracks->at(jtrk).y;
          Float_t trk2_z = v_tracks->at(jtrk).z;
          Float_t trk2_px = v_tracks->at(jtrk).px;
          Float_t trk2_py = v_tracks->at(jtrk).py;
          Float_t trk2_pz = v_tracks->at(jtrk).pz;
          Float_t trk2_pt = v_tracks->at(jtrk).pt;
          Float_t trk2_p = v_tracks->at(jtrk).p;
          auto trk2_pe = get_pe<Float_t>(trk2_px, trk2_py, trk2_pz);

          // Circle intersection
          // math.stackexchange.com/questions/256100/how-can-i-find-the-points-at-which-two-circles-intersect?newreg=b9dc98e45b514173ae85b6dbaf4d2508
          Float_t cir_D = get_r(cir1_X0-cir2_X0, cir1_Y0-cir2_Y0);
          if(cir_D < fabs(cir1_R - cir2_R) || cir_D > fabs(cir1_R + cir2_R)) continue;
          Float_t cir_a = (cir1_R*cir1_R - cir2_R*cir2_R) / (cir_D*cir_D);
          Float_t cir_b = sqrt(2 * (cir1_R*cir1_R + cir2_R*cir2_R) / (cir_D*cir_D) - square((cir1_R*cir1_R - cir2_R*cir2_R) / (cir_D*cir_D)) - 1);
          Float_t cir_sign = trk1_px*(cir2_Y0-cir1_Y0) - trk1_py*(cir2_X0-cir1_X0) > 0 ? 1 : -1;
          Float_t sec_x = (cir1_X0 + cir2_X0)/2 + cir_a * (cir2_X0 - cir1_X0)/2 + cir_sign * cir_b * (cir2_Y0 - cir1_Y0)/2;
          Float_t sec_y = (cir1_Y0 + cir2_Y0)/2 + cir_a * (cir2_Y0 - cir1_Y0)/2 - cir_sign * cir_b * (cir2_X0 - cir1_X0)/2;
          Float_t sec_r = get_r(sec_x, sec_y);

          const Float_t PI = TMath::Pi();
          Float_t sec1_phi = atan2(sec_x-cir1_X0, sec_y-cir1_Y0);
          Float_t dca1_phi = atan2(trk1_x-cir1_X0, trk1_y-cir1_Y0);
          Float_t dphi1 = sec1_phi - dca1_phi;
          while(dphi1 >= PI) dphi1 -= 2*PI;
          while(dphi1 < -PI) dphi1 += 2*PI;
          Float_t dt1 = fabs(dphi1 * cir1_R) / trk1_pt;
          Float_t sec1_z = trk1_z + trk1_pz * dt1;

          Float_t sec2_phi = atan2(sec_x-cir2_X0, sec_y-cir2_Y0);
          Float_t dca2_phi = atan2(trk2_x-cir2_X0, trk2_y-cir2_Y0);
          Float_t dphi2 = sec2_phi - dca2_phi;
          while(dphi2 >= PI) dphi2 -= 2*PI;
          while(dphi2 < -PI) dphi2 += 2*PI;
          Float_t dt2 = fabs(dphi2 * cir2_R) / trk2_pt;
          Float_t sec2_z = trk2_z + trk2_pz * dt2;

          Float_t minv = (trk1_pe+trk2_pe).M();

          Int_t trk1_etabin = v_tracks->at(itrk).cal_etabin[0];
          Float_t trk1_cemc_e = v_tracks->at(itrk).cal_cluster_e[0];
          Int_t trk1_hcalin_e = v_tracks->at(itrk).cal_cluster_e[1];
          Int_t trk2_etabin = v_tracks->at(jtrk).cal_etabin[0];
          Float_t trk2_cemc_e = v_tracks->at(jtrk).cal_cluster_e[0];
          Int_t trk2_hcalin_e = v_tracks->at(jtrk).cal_cluster_e[1];

          if( filled.count(itrk) == 0 && filled.count(jtrk) == 0 && 
              trk1_hcalin_e/trk1_cemc_e < 0.2 && trk2_hcalin_e/trk2_cemc_e < 0.2 &&
              (fabs(sec_r - 12) < 4 || fabs(sec_r - 35) < 1) )
          {
            h3_minv->Fill(sec1_z, sec2_z, minv);
            h3_e_ratio->Fill(trk1_etabin, trk1_hcalin_e/trk1_cemc_e, trk1_cemc_e/trk1_p);
            h3_e_ratio->Fill(trk2_etabin, trk2_hcalin_e/trk2_cemc_e, trk2_cemc_e/trk2_p);
            filled.insert(itrk);
            filled.insert(jtrk);
          }
        } // jtrk
      } // itrk
    } // event

    delete f_in;
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
