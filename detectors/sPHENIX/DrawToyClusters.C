void DrawToyClusters(const Int_t iev = 0)
{
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

  auto f = new TFile("/phenix/spin/phnxsp01/zji/data/sphenix/histos/training-0-0.root");
  TTree *t = static_cast<TTree*>(f->Get("T"));
  t->SetBranchAddress("layer", &li);
  t->SetBranchAddress("ztan", &zr);
  t->SetBranchAddress("adc", &v_adc);
  t->SetBranchAddress("reco_rphi", &v_reco_rphi);
  t->SetBranchAddress("reco_z", &v_reco_z);
  t->SetBranchAddress("reco_adc", &v_reco_adc);
  t->SetBranchAddress("nreco", &v_nreco);
  t->SetBranchAddress("truth_rphi", &v_truth_rphi);
  t->SetBranchAddress("truth_z", &v_truth_z);
  t->SetBranchAddress("truth_adc", &v_truth_adc);
  t->SetBranchAddress("ntruth", &v_ntruth);
  t->SetBranchAddress("ntouch", &training_ntouch);

  auto h2_adc = new TH2F("h2_adc","ADC", 2*nd+1,-nd-0.5,nd+0.5, 2*nd+1,-nd-0.5,nd+0.5);

  //for(Int_t iev=0; iev<1000; iev++)
  {
    t->GetEntry(iev);
    for(Int_t i=-nd; i<=nd; i++)
      for(Int_t j=-nd; j<=nd; j++)
        h2_adc->Fill(i, j, v_adc[(i+nd)*(2*nd+1)+(j+nd)]);
  }

  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gStyle->SetOptStat(0);
  h2_adc->GetXaxis()->SetTitle("#phi");
  h2_adc->GetYaxis()->SetTitle("z");
  h2_adc->Draw("COLZ");
  c->Print("results/ToyADC.pdf");
}
