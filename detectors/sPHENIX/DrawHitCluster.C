void DrawHitCluster()
{
  const Int_t nh = 2;
  const char *name[nh] = {"hit", "cluster"};
  const char *opt[nh] = {"COL", "SAME"};

  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *f = new TFile("/phenix/spin/phnxsp01/zji/data/sphenix/output/G4sPHENIX_g4svtx_eval0.root");
  for(Int_t ih=0; ih<nh; ih++)
  {
    TNtuple *ntp = (TNtuple*)f->Get(Form("ntp_%s",name[ih]));
    Float_t event, layer, adc, phi, z;
    ntp->SetBranchAddress("event", &event);
    ntp->SetBranchAddress("layer", &layer);
    ntp->SetBranchAddress("adc", &adc);
    ntp->SetBranchAddress("phi", &phi);
    ntp->SetBranchAddress("z", &z);

    TH2 *h2 = new TH2F(Form("h2_%s",name[ih]), "", 100,10.,40., 100,-0.8,0.);
    for(Int_t ien=0; ien<ntp->GetEntries(); ien++)
    {
      ntp->GetEntry(ien);
      if(event==0 && layer==10)
        h2->Fill(z, phi, adc);
    }
    h2->GetXaxis()->SetTitle("z (cm)");
    h2->GetYaxis()->SetTitle("#phi (rad)");
    h2->GetYaxis()->SetTitleOffset(0.95);
    h2->Draw(opt[ih]);
  }

  c->Print("results/HitCluster.pdf");
}
