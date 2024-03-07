void DrawEEMass()
{
  const char *dataset[3] = {"ph", "auau", "auau-cuts"};

  TCanvas *c0 = new TCanvas("c0", "c0", 3*600, 2*600);
  c0->Divide(3, 2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  Int_t ipad0 = 1;

  for(Int_t id=0; id<3; id++)
  {
    auto f = new TFile(Form("results/ee-minv-%s.root", dataset[id]));
    auto h2_minv_r = (TH2*)f->Get("h_minv_r");

    mcd(0, id+1);
    TH1 *h_minv = h2_minv_r->ProjectionX("h_minv");
    h_minv->SetTitle(dataset[id]);
    aset(h_minv);
    style(h_minv, 20, 2);
    h_minv->GetXaxis()->SetNdivisions(5);
    h_minv->DrawCopy();
    delete h_minv;

    mcd(0, id+4);
    TH1 *h_pcar = h2_minv_r->ProjectionY("h_pcar");
    h_pcar->SetTitle(dataset[id]);
    aset(h_pcar, "PCA_{r} (cm)","", 0.,6.);
    style(h_pcar, 20, 2);
    h_pcar->GetXaxis()->SetNdivisions(5);
    h_pcar->DrawCopy();
    delete h_pcar;
  }

  c0->Print("results/EEMass.pdf");
}
