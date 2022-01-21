void DrawHitCluster()
{
  const Int_t nh = 2;
  const char *name[nh] = {"hit", "cluster"};
  TH2 *h2[nh*2];

  TCanvas *c = new TCanvas("c", "c", 600, 600);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TFile *f = new TFile("/phenix/spin/phnxsp01/zji/data/sphenix/output/G4sPHENIX_g4svtx_eval.root");
  for(Int_t ih=0; ih<nh; ih++)
  {
    TNtuple *ntp = (TNtuple*)f->Get(Form("ntp_%s",name[ih]));
    Float_t event, layer, adc, phi, z, gphi, gz;
    ntp->SetBranchAddress("event", &event);
    ntp->SetBranchAddress("layer", &layer);
    ntp->SetBranchAddress("adc", &adc);
    ntp->SetBranchAddress("phi", &phi);
    ntp->SetBranchAddress("z", &z);
    ntp->SetBranchAddress("gphi", &gphi);
    ntp->SetBranchAddress("gz", &gz);

    h2[ih] = new TH2F(Form("h2_%s",name[ih]), "", 500,-50.,50., 500,-2.5,2.5);
    h2[ih+nh] = new TH2F(Form("h2_g4%s",name[ih]), "", 500,-50.,50., 500,-2.5,2.5);
    for(Int_t ien=0; ien<ntp->GetEntries(); ien++)
    {
      ntp->GetEntry(ien);
      if(event==0 && layer==10)
      {
        h2[ih]->Fill(z, phi, adc);
        h2[ih+nh]->Fill(gz, gphi, adc);
      }
    }
  }

  h2[0]->GetXaxis()->SetTitle("z (cm)");
  h2[0]->GetYaxis()->SetTitle("#phi (rad)");
  h2[0]->GetYaxis()->SetTitleOffset(0.95);
  h2[0]->Draw("COL");
  h2[1]->Draw("CONT SAME");
  h2[3]->Draw("SAME");
  c->Print("results/HitCluster.pdf");
}
