// g++ -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib -lc10 -ltorch_cpu -ltorch -o LoadTorchScript LoadTorchScript.cpp
#include <iostream>
#include <memory>

// One-stop header
#include <torch/script.h>

using namespace std;

int main(int argc, const char *argv[])
{
  if(argc != 2)
  {
    cerr << "Usage: " << argv[0] << " <path-to-exported-script-module>" << endl;
    return 1;
  }

  torch::jit::script::Module module;
  try
  {
    // Deserialize the ScriptModule from a file using torch::jit::load()
    module = torch::jit::load(argv[1]);
  }
  catch(const c10::Error& e)
  {
    cerr << "Error: cannot load model " << argv[1] << endl;
    return 1;
  }
  cout << "Load model " << argv[1] << endl;

  // Create a vector of inputs
  vector<torch::jit::IValue> inputs;
  inputs.push_back(torch::ones({1, 25}));

  // Execute the model and turn its output into a tensor
  at::Tensor output = module.forward(inputs).toTensor();
  cout << output.slice(/*dim=*/1, /*start=*/0, /*end=*/2) << endl;

  return 0;
}
