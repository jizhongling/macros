void DrawToyClusters(const Int_t iev = 0)
{
  const Int_t nd = 5;
  const Int_t nc = 10;
  const Int_t nb = 20;
  Short_t training_ntouch, li;
  Float_t zr;
  array<Short_t, (2*nd+1)*(2*nd+1)> v_adc;
  array<Float_t, nc> v_reco_phi;
  array<Float_t, nc> v_reco_z;
  array<Short_t, nc> v_reco_adc;
  array<Short_t, nb> v_nreco;
  array<Float_t, nc> v_truth_phi;
  array<Float_t, nc> v_truth_z;
  array<Short_t, nc> v_truth_adc;
  array<Short_t, nb> v_ntruth;

  auto f = new TFile("/phenix/spin/phnxsp01/zji/data/sphenix/histos/training.root");
  TTree *t = static_cast<TTree*>(f->Get("T"));
  t->SetBranchAddress("layer", &li);
  t->SetBranchAddress("ztan", &zr);
  t->SetBranchAddress("adc", &v_adc);
  t->SetBranchAddress("reco_phi", &v_reco_phi);
  t->SetBranchAddress("reco_z", &v_reco_z);
  t->SetBranchAddress("reco_adc", &v_reco_adc);
  t->SetBranchAddress("nreco", &v_nreco);
  t->SetBranchAddress("truth_phi", &v_truth_phi);
  t->SetBranchAddress("truth_z", &v_truth_z);
  t->SetBranchAddress("truth_adc", &v_truth_adc);
  t->SetBranchAddress("ntruth", &v_ntruth);
  t->SetBranchAddress("ntouch", &training_ntouch);

  auto h2_phiz = new TH2F("h2_phiz","#phi vs z", 200,-1,1, 200,-1,1);
  auto h2_adc = new TH2F("h2_adc","ADC", 2*nd+1,-nd-0.5,nd+0.5, 2*nd+1,-nd-0.5,nd+0.5);

  for(Long64_t iev=0; iev<t->GetEntries(); iev++)
  {
    t->GetEntry(iev);

    for(Int_t ic=0; ic<nc; ic++)
      if(v_truth_adc[ic] > 0)
        h2_phiz->Fill(v_truth_phi[ic], v_truth_z[ic]);

    for(Int_t i=-nd; i<=nd; i++)
      for(Int_t j=-nd; j<=nd; j++)
        if(v_adc[(i+nd)*(2*nd+1)+(j+nd)] > 0)
          h2_adc->Fill(i, j, v_adc[(i+nd)*(2*nd+1)+(j+nd)]);
  }

  auto h_phi = h2_phiz->ProjectionX("h_phi");
  auto h_z = h2_phiz->ProjectionY("h_z");

  auto h_adc_phi = h2_adc->ProjectionX("h_adc_phi");
  auto h_adc_z = h2_adc->ProjectionY("h_adc_z");

  TCanvas *c = new TCanvas("c", "c", 1800, 1200);
  c->Divide(3, 2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  c->cd(1);
  gPad->SetLogz();
  h2_phiz->GetXaxis()->SetTitle("#phi");
  h2_phiz->GetYaxis()->SetTitle("z");
  h2_phiz->Draw("COLZ");

  c->cd(2);
  h_phi->GetXaxis()->SetTitle("#phi");
  h_phi->Fit("gaus", "Q");

  c->cd(3);
  h_z->GetXaxis()->SetTitle("z");
  h_z->Fit("gaus", "Q");

  c->cd(4);
  gPad->SetLogz();
  h2_adc->GetXaxis()->SetTitle("#phi");
  h2_adc->GetYaxis()->SetTitle("z");
  h2_adc->Draw("COLZ");

  c->cd(5);
  h_adc_phi->GetXaxis()->SetTitle("#phi");
  h_adc_phi->Fit("gaus", "Q");

  c->cd(6);
  h_adc_z->GetXaxis()->SetTitle("z");
  h_adc_z->Fit("gaus", "Q");

  c->Print("results/ToyADC.pdf");
}
