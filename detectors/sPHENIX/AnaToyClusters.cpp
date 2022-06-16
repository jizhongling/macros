// g++ -Wall -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib `root-config --cflags --glibs` -o AnaToyClusters AnaToyClusters.cpp
#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <utility>
#include <algorithm>

#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TNtuple.h>

using namespace std;

typedef vector<vector<Short_t>> vvS;
template<class T> inline constexpr T square(const T &x) { return x*x; }

Short_t maxHalfSizePhi = 5;
Short_t maxHalfSizeZ = 5;
Short_t phibins = 11;
Short_t zbins = 11;

struct ihit
{
  Short_t iphi = 0;
  Short_t iz = 0;
  Short_t adc = 0;
  Short_t edge = 0;
};

void remove_hit(Int_t adc, Int_t phibin, Int_t zbin, multimap<Short_t, ihit> &all_hit_map, vvS &adcval)
{
  typedef multimap<Short_t, ihit>::iterator hit_iterator;
  pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
  hit_iterator it = iterpair.first;
  for (; it != iterpair.second; ++it) {
    if (it->second.iphi == phibin && it->second.iz == zbin) { 
      all_hit_map.erase(it);
      break;
    }
  }
  adcval[phibin][zbin] = 0;
}

void remove_hits(vector<ihit> &ihit_list, multimap<Short_t, ihit> &all_hit_map, vvS &adcval)
{
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    Short_t adc    = iter->adc; 
    Short_t phibin = iter->iphi;
    Short_t zbin   = iter->iz;
    remove_hit(adc,phibin,zbin,all_hit_map,adcval);
  }
}

void find_z_range(Int_t phibin, Int_t zbin, const vvS &adcval, Int_t& zdown, Int_t& zup, Int_t &touch, Int_t &edge)
{
  const Int_t FitRangeZ = (Int_t) maxHalfSizeZ;
  const Int_t NZBinsMax = (Int_t) zbins;
  zup = 0;
  zdown = 0;
  for(Int_t iz=0; iz< FitRangeZ; iz++){
    Int_t cz = zbin + iz;

    if(cz <= 0 || cz >= NZBinsMax){
      // zup = iz;
      edge++;
      break; // truncate edge
    }

    if(adcval[phibin][cz] <= 0) {
      break;
    }
    //check local minima and break at minimum.
    if(cz<NZBinsMax-4){//make sure we stay clear from the edge
      if(adcval[phibin][cz]+adcval[phibin][cz+1] < 
          adcval[phibin][cz+2]+adcval[phibin][cz+3]){//rising again
        zup = iz+1;
        touch++;
        break;
      }
    }
    zup = iz;
  }
  for(Int_t iz=0; iz< FitRangeZ; iz++){
    Int_t cz = zbin - iz;
    if(cz <= 0 || cz >= NZBinsMax){
      //      zdown = iz;
      edge++;
      break; // truncate edge
    }
    if(adcval[phibin][cz] <= 0) {
      break;
    }
    if(cz>4){//make sure we stay clear from the edge
      if(adcval[phibin][cz]+adcval[phibin][cz-1] < 
          adcval[phibin][cz-2]+adcval[phibin][cz-3]){//rising again
        zdown = iz+1;
        touch++;
        break;
      }
    }
    zdown = iz;
  }
  return;
}

void find_phi_range(Int_t phibin, Int_t zbin, const vvS &adcval, Int_t& phidown, Int_t& phiup, Int_t &touch, Int_t &edge)
{
  Int_t FitRangePHI = (Int_t) maxHalfSizePhi;
  Int_t NPhiBinsMax = (Int_t) phibins;
  phidown = 0;
  phiup = 0;
  for(Int_t iphi=0; iphi< FitRangePHI; iphi++){
    Int_t cphi = phibin + iphi;
    if(cphi < 0 || cphi >= NPhiBinsMax){
      // phiup = iphi;
      edge++;
      break; // truncate edge
    }

    //break when below minimum
    if(adcval[cphi][zbin] <= 0) {
      // phiup = iphi;
      break;
    }
    //check local minima and break at minimum.
    if(cphi<NPhiBinsMax-4){//make sure we stay clear from the edge
      if(adcval[cphi][zbin]+adcval[cphi+1][zbin] < 
          adcval[cphi+2][zbin]+adcval[cphi+3][zbin]){//rising again
        phiup = iphi+1;
        touch++;
        break;
      }
    }
    phiup = iphi;
  }

  for(Int_t iphi=0; iphi< FitRangePHI; iphi++){
    Int_t cphi = phibin - iphi;
    if(cphi < 0 || cphi >= NPhiBinsMax){
      // phidown = iphi;
      edge++;
      break; // truncate edge
    }

    if(adcval[cphi][zbin] <= 0) {
      //phidown = iphi;
      break;
    }
    if(cphi>4){//make sure we stay clear from the edge
      if(adcval[cphi][zbin]+adcval[cphi-1][zbin] < 
          adcval[cphi-2][zbin]+adcval[cphi-3][zbin]){//rising again
        phidown = iphi+1;
        touch++;
        break;
      }
    }
    phidown = iphi;
  }
  return;
}

void get_cluster(Int_t phibin, Int_t zbin, const vvS &adcval, vector<ihit> &ihit_list, Int_t &touch, Int_t &edge)
{
  // search along phi at the peak in z

  Int_t zup =0;
  Int_t zdown =0;
  find_z_range(phibin, zbin, adcval, zdown, zup, touch, edge);
  //now we have the z extent of the cluster, go find the phi edges

  for(Int_t iz=zbin - zdown ; iz<= zbin + zup; iz++){
    Int_t phiup = 0;
    Int_t phidown = 0;
    find_phi_range(phibin, iz, adcval, phidown, phiup, touch, edge);
    for (Int_t iphi = phibin - phidown; iphi <= (phibin + phiup); iphi++){
      if(adcval[iphi][iz]>0){
        ihit hit;
        hit.iphi = iphi;
        hit.iz = iz;
        hit.adc = adcval[iphi][iz];
        if(touch>0){
          if((iphi == (phibin - phidown))||
              (iphi == (phibin + phiup))){
            hit.edge = 1;
          }
        }
        ihit_list.push_back(hit);
      }
    }
  }
  return;
}

tuple<Float_t, Float_t, Int_t> calc_cluster_parameter(const vector<ihit> &ihit_list)
{
  Int_t z_sum = 0.;
  Int_t phi_sum = 0.;
  Int_t adc_sum = 0.;

  Int_t clus_size = ihit_list.size();
  if(clus_size == 1) return make_tuple(0, 0, 0);

  // loop over the hits in this cluster
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    Int_t iphi = iter->iphi;
    Int_t iz   = iter->iz;
    Int_t adc = iter->adc; 
    if (adc <= 0) continue;

    phi_sum += iphi * adc;
    z_sum += iz * adc;
    adc_sum += adc;
  }
  if (adc_sum < 10){
    return make_tuple(0, 0, 0);  // skip obvious noise "clusters"
  }  

  // This is the global position
  Float_t clusphi = (Float_t)phi_sum / adc_sum;
  Float_t clusz = (Float_t)z_sum / adc_sum;
  return make_tuple(clusphi, clusz, adc_sum);
}

int main(int argc, const char *argv[])
{
  if(argc != 4)
  {
    cerr << "Usage: " << argv[0] << " <nev> <output>(-<ith>.root) <ith>" << endl;
    return 1;
  }

  const Int_t nev = stoi(string(argv[1]));
  const Int_t ith = stoi(string(argv[3]));

  const Int_t nd = 5;
  const Int_t nc = 10;
  const Int_t nb = 20;
  Short_t training_ntouch, li;
  Float_t zr;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, nc> v_reco_rphi;
  array<Float_t, nc> v_reco_z;
  array<Short_t, nc> v_reco_adc;
  array<Short_t, nb> v_nreco;
  array<Float_t, nc> v_truth_rphi;
  array<Float_t, nc> v_truth_z;
  array<Short_t, nc> v_truth_adc;
  array<Short_t, nb> v_ntruth;

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
  t_out->Branch("ntouch", &training_ntouch);

  auto rnd = new TRandom3(0);

  for(Int_t iev=0; iev<nev; iev++)
  {
    li = 0;
    zr = 0.;
    v_adc.fill(0);
    v_reco_rphi.fill(0.);
    v_reco_z.fill(0.);
    v_reco_adc.fill(0);
    v_nreco.fill(0);
    v_truth_rphi.fill(0.);
    v_truth_z.fill(0.);
    v_truth_adc.fill(0);
    v_ntruth.fill(0);
    training_ntouch = 0;

    vvS adcval(phibins, vector<Short_t>(zbins, 0));
    multimap<Short_t, ihit> all_hit_map;

    Float_t truth_peak[2];
    for(Int_t ic=0; ic<2; ic++)
    {
      v_truth_rphi[ic] = rnd->Gaus(0, 0.8);
      v_truth_z[ic] = rnd->Gaus(0, 0.8);
      truth_peak[ic] = rnd->Gaus(60, 10);
    }
    v_ntruth[0] = 2;

    for(Int_t i=-nd; i<=nd; i++)
      for(Int_t j=-nd; j<=nd; j++)
      {
        Short_t adc = 0;
        for(Int_t ic=0; ic<2; ic++)
        {
          Float_t dist = sqrt( square(i - v_truth_rphi[ic]) + square(j - v_truth_z[ic]) );
          Short_t gadc = TMath::Abs( truth_peak[ic] * TMath::Gaus(dist, 0, 1, kTRUE) * (1 + rnd->Gaus(0, 0.08)) );
          if(gadc > 0)
            v_truth_adc[ic] += static_cast<Short_t>(round(gadc));
          adc += static_cast<Short_t>(round(TMath::Abs( gadc * (1 + rnd->Gaus(0, 0.08)) )));
        } // ic
        if(adc > 0)
        {
          ihit thisHit;
          thisHit.iphi = i + nd;
          thisHit.iz = j + nd;
          thisHit.adc = adc;
          thisHit.edge = 0;
          all_hit_map.insert(make_pair(adc, thisHit));
          adcval[i+nd][j+nd] = adc;
          v_adc[(i+nd)*(2*nd+1)+(j+nd)] = adc;
        } // adc > 0
      } // i, j

    Int_t iclus = 0;
    while(all_hit_map.size() > 0 && iclus < nc)
    {
      auto iter = all_hit_map.rbegin();
      ihit hiHit = iter->second;
      Int_t iphi = hiHit.iphi;
      Int_t iz = hiHit.iz;

      // put all hits in the all_hit_map (sorted by adc)
      // start with highest adc hit
      // -> cluster around it and get vector of hits
      vector<ihit> ihit_list;
      Int_t ntouch = 0;
      Int_t nedge = 0;
      get_cluster(iphi, iz, adcval, ihit_list, ntouch, nedge);
      training_ntouch += ntouch;

      // -> calculate cluster parameters
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      auto [clusphi, clusz, clusadc] = calc_cluster_parameter(ihit_list);
      v_reco_rphi[iclus] = clusphi - nd;
      v_reco_z[iclus] = clusz - nd;
      v_reco_adc[iclus] = clusadc;

      remove_hits(ihit_list, all_hit_map, adcval);
      ihit_list.clear();
      iclus++;
    }
    v_nreco[0] = iclus;
    training_ntouch /= iclus;

    t_out->Fill();
  } // iev

  f_out->Write();
  f_out->Close();
  return 0;
}
