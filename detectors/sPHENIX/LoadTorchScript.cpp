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
  if(argc != 6)
  {
    cerr << "Usage: " << argv[0] << " <script-module> <type> <nout> <root-file> <ientry>" << endl;
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
    cerr << "Error: Cannot load model " << argv[1] << endl;
    return 1;
  }
  const int type = stoi(string(argv[2]));
  const int nout = stoi(string(argv[3]));
  cout << "Load model " << argv[1] << " with type " << type << " and nout " << nout << endl;

  auto f = new TFile(argv[4]);
  if(f->IsZombie())
  {
    cerr << "Error: Cannot open file " << argv[4] << endl;
    return 1;
  }
  cout << "Open file " << argv[4] << endl;

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
  array<Float_t, nc> v_truth_phisize;
  array<Float_t, nc> v_truth_zsize;
  array<Short_t, nc> v_truth_adc;
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
  T->SetBranchAddress("truth_phisize", &v_truth_phisize);
  T->SetBranchAddress("truth_zsize", &v_truth_zsize);
  T->SetBranchAddress("truth_adc", &v_truth_adc);
  T->SetBranchAddress("ntruth", &v_ntruth);

  T->GetEntry(stoi(string(argv[5])));

  try
  {
    // Create a vector of inputs
    vector<torch::jit::IValue> inputs;
    inputs.emplace_back(torch::stack({
          torch::from_blob(vector<Float_t>(v_adc.begin(), v_adc.end()).data(), {1, 2*nd+1, 2*nd+1}, torch::kFloat32).sub(75).clamp_min(0),
          torch::full({1, 2*nd+1, 2*nd+1}, li, torch::kFloat32),
          torch::full({1, 2*nd+1, 2*nd+1}, zr, torch::kFloat32)
          }, 1));

    // Execute the model and turn its output into a tensor
    at::Tensor output = module.forward(inputs).toTensor();
    if(type == 0)
    {
      int ntruth = 0;
      for(size_t i=5; i<nb; i++)
        ntruth += v_ntruth[i];
      cout << "NN predictions: " << output.argmax(1)[0].item<int>() << endl
        << "Truth clusters: " << ntruth << endl;
    }
    else if(type == 1)
    {
      output = output.round().clamp(-(float)nd, (float)nd);
      cout << "NN (phi,z): ";
      for(int i=0; i<nout; i++)
        cout << "(" << output[0][0][i].item<float>() << ", " << output[0][1][i].item<float>() << "), ";
      cout << endl;

      cout << "Truth (phi,z): ";
      for(size_t i=0; i<nc; i++)
        cout << "(" << v_truth_phi[i] << ", " << v_truth_z[i] << "), ";
      cout << endl;
    }
    else if(type <= 4)
    {
      cout << "NN predictions: ";
      for(int i=0; i<nout; i++)
        cout << output[0][i].item<float>() << ", ";
      cout << endl;
    }
  }
  catch(const c10::Error &e)
  {
    cerr << "Error: Failed to execute NN modules" << endl;
    return 1;
  }

  return 0;
}
