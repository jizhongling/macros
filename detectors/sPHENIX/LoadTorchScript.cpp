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
  if(argc != 4)
  {
    cerr << "Usage: " << argv[0] << " <script-module> <root-file> <ientry>" << endl;
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
  cout << "Load model " << argv[1] << endl;

  auto f = new TFile(argv[2]);
  if(f->IsZombie())
  {
    cerr << "Error: cannot open file " << argv[2] << endl;
    return 1;
  }
  cout << "Open file " << argv[2] << endl;

  const size_t nd = 10;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, (2*nd+1)*(2*nd+1)> v_gedep;
  Short_t li;
  Float_t zr;
  vector<Float_t> *v_reco_phi = 0;
  vector<Float_t> *v_reco_z = 0;
  Short_t nreco;
  vector<Float_t> *v_truth_phi = 0;
  vector<Float_t> *v_truth_z = 0;
  Short_t ntruth;

  TTree *T = static_cast<TTree*>(f->Get("T"));
  T->SetBranchAddress("adc", &v_adc);
  T->SetBranchAddress("gedep", &v_gedep);
  T->SetBranchAddress("layer", &li);
  T->SetBranchAddress("ztan", &zr);
  T->SetBranchAddress("reco_phi", &v_reco_phi);
  T->SetBranchAddress("reco_z", &v_reco_z);
  T->SetBranchAddress("nreco", &nreco);
  T->SetBranchAddress("truth_phi", &v_truth_phi);
  T->SetBranchAddress("truth_z", &v_truth_z);
  T->SetBranchAddress("ntruth", &ntruth);

  T->GetEntry(stoi(string(argv[3])));

  // Create a vector of inputs
  vector<torch::jit::IValue> inputs;
  inputs.push_back(torch::stack({
        torch::from_blob(vector<Float_t>(v_adc.begin(), v_adc.end()).data(), {1, 11, 11}, torch::kFloat32).sub(75).clamp_min(0),
        torch::full({1, 11, 11}, li, torch::kFloat32),
        torch::full({1, 11, 11}, zr, torch::kFloat32)
        }, 1));

  // Execute the model and turn its output into a tensor
  at::Tensor output = module.forward(inputs).toTensor();
  cout << output.slice(/*dim=*/1, /*start=*/0, /*end=*/2).slice(2, 0, 2) << endl;
  cout << "NN phi position (3 clusters): " << output.data_ptr<float>()[0] << ", " << output.data_ptr<float>()[3] << endl;

  at::Tensor target = torch::from_blob(vector<Float_t>(v_gedep.begin(), v_gedep.end()).data(), {1, 21, 21}, torch::kFloat32).mul(1e6).clamp_max(10);
  cout << target.slice(/*dim=*/1, /*start=*/0, /*end=*/2).slice(2, 0, 2) << endl;
  cout << "Target phi position (3 clusters): " << target.data_ptr<float>()[0] << ", " << target.data_ptr<float>()[3] << endl;

  return 0;
}
