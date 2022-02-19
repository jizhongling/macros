// g++ -Wall -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib `root-config --cflags --glibs` -lc10 -ltorch_cpu -ltorch -o LoadTorchScript LoadTorchScript.cpp
#include <iostream>
#include <memory>
#include <array>
#include <vector>
#include <algorithm>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>

// One-stop header
#include <torch/script.h>

using namespace std;

int main(int argc, const char *argv[])
{
  if(argc != 5)
  {
    cerr << "Usage: " << argv[0] << " <script-module> <type> <root-file> <ientry>" << endl;
    return 1;
  }

  torch::jit::script::Module module;
  try
  {
    // Deserialize the ScriptModule from a file using torch::jit::load()
    module = torch::jit::load(argv[1]);
  }
  catch(const c10::Error &e)
  {
    cerr << "Error: cannot load model " << argv[1] << endl;
    return 1;
  }
  const int type = stoi(string(argv[2]));
  cout << "Load model " << argv[1] << " with type " << type << endl;

  auto f = new TFile(argv[3]);
  if(f->IsZombie())
  {
    cerr << "Error: cannot open file " << argv[3] << endl;
    return 1;
  }
  cout << "Open file " << argv[3] << endl;

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

  TTree *T = static_cast<TTree*>(f->Get("T"));
  T->SetBranchAddress("layer", &li);
  T->SetBranchAddress("ztan", &zr);
  T->SetBranchAddress("adc", &v_adc);
  T->SetBranchAddress("reco_phi", &v_reco_phi);
  T->SetBranchAddress("reco_z", &v_reco_z);
  T->SetBranchAddress("reco_adc", &v_reco_adc);
  T->SetBranchAddress("nreco", &v_nreco);
  T->SetBranchAddress("truth_phi", &v_truth_phi);
  T->SetBranchAddress("truth_z", &v_truth_z);
  T->SetBranchAddress("truth_edep", &v_truth_edep);
  T->SetBranchAddress("ntruth", &v_ntruth);

  T->GetEntry(stoi(string(argv[4])));

  // Create a vector of inputs
  vector<torch::jit::IValue> inputs;
  inputs.emplace_back(torch::stack({
        torch::from_blob(vector<Float_t>(v_adc.begin(), v_adc.end()).data(), {1, 2*nd+1, 2*nd+1}, torch::kFloat32).sub(75).clamp_min(0),
        torch::full({1, 2*nd+1, 2*nd+1}, li, torch::kFloat32),
        torch::full({1, 2*nd+1, 2*nd+1}, zr, torch::kFloat32)
        }, 1));

  // Execute the model and turn its output into a tensor
  at::Tensor output = module.forward(inputs).toTensor();
  cout << "NN predictions: " << output.slice(/*dim=*/1, /*start=*/0, /*end=*/type).slice(2, 0, type) << endl
    << "NN first (phi, z): " << output.data_ptr<float>()[0] << ", " << output.data_ptr<float>()[type] << endl;

  at::Tensor truth_phi = torch::from_blob(vector<Short_t>(v_ntruth.begin(), v_ntruth.end()).data(), {1, nb}, torch::kFloat32);
  at::Tensor truth_z = torch::from_blob(vector<Float_t>(v_truth_z.begin(), v_truth_z.end()).data(), {1, nc}, torch::kFloat32);
  cout << "Truth (phi,z): ";
  for(size_t i=0; i<nc; i++)
    cout << "(" << v_truth_phi[i] << ", " << v_truth_z[i] << "), ";
  cout << endl;

  return 0;
}
