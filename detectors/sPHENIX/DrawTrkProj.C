void DrawTrkProj()
{
  const char *dataset[2] = {"trk-proj-auau", "hf-electron-auau"};
  const char *type[3] = {"cemc", "hcalin", "hcalout"};

  TCanvas *c_tmp = new TCanvas("c_tmp", "c_tmp", 600, 600);
  TCanvas *c0 = new TCanvas("c0", "c0", 4*600, 3*600);
  c0->Divide(4, 3);
  TCanvas *c1 = new TCanvas("c1", "c1", 2*600, 1*600);
  c1->Divide(2, 1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  Int_t ipad0 = 1;
  Int_t ipad1 = 1;

  for(Int_t id=0; id<2; id++)
  {
    c_tmp->cd();
    auto f = new TFile(Form("results/%s.root", dataset[id]));
    auto ntp = (TNtuple*)f->Get("ntp");
    if(id == 0)
      for(Int_t it=0; it<3; it++)
      {
        ntp->Draw(Form("clus_dphi_%s >> h_dphi_%s", type[it], type[it]));
        ntp->Draw(Form("clus_deta_%s >> h_deta_%s", type[it], type[it]));
        ntp->Draw(Form("clus_dphi_%s:clus_dphi_outer_%s >> h_dphi_diff_%s", type[it], type[it], type[it]));
        ntp->Draw(Form("clus_deta_%s:clus_deta_outer_%s >> h_deta_diff_%s", type[it], type[it], type[it]));
      }
    ntp->Draw("clus_e_cemc/sqrt(pt*pt+pz*pz) >> h_e_p", "fabs(clus_e_cemc/sqrt(pt*pt+pz*pz)-1)<0.5");

    if(id == 0)
      for(Int_t it=0; it<3; it++)
    {
      mcd(0, ipad0++);
      ((TH1*)gDirectory->Get(Form("h_dphi_%s", type[it])))->Draw();

      mcd(0, ipad0++);
      ((TH1*)gDirectory->Get(Form("h_deta_%s", type[it])))->Draw();

      mcd(0, ipad0++);
      ((TH2*)gDirectory->Get(Form("h_dphi_diff_%s", type[it])))->Draw("COLZ");

      mcd(0, ipad0++);
      ((TH2*)gDirectory->Get(Form("h_deta_diff_%s", type[it])))->Draw("COLZ");
    }

    mcd(1, ipad1++);
    TH1 *h_e_p = (TH1*)gDirectory->Get("h_e_p");
    h_e_p->SetTitle(id?"Electron":"All");
    h_e_p->GetXaxis()->SetTitle("clus_e/mom");
    h_e_p->Draw();
  }

  c0->Print("results/Trk-Proj.pdf(");
  c1->Print("results/Trk-Proj.pdf)");
}
