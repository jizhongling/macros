void DrawPull()
{
  const Int_t max_ntouch = 4;
  TH1 *h_phi[3][max_ntouch+1];
  TH1 *h_z[3][max_ntouch+1];
  TGraphErrors *g_sigma[3];

  auto f = new TFile("results/pull.root");
  for(Int_t ilayer = 0; ilayer <= 2; ilayer++)
  {
    for(Int_t itouch = 0; itouch <= max_ntouch; itouch++)
    {
      h_phi[ilayer][itouch] = (TH1*)f->Get(Form("h_phi_li%d_ovl%d",ilayer,itouch));
      h_z[ilayer][itouch] = (TH1*)f->Get(Form("h_z_li%d_ovl%d",ilayer,itouch));
    }
    g_sigma[ilayer] = new TGraphErrors(max_ntouch+1);
  }


  auto fn_fit = new TF1("fn_fit", "gaus(0)+pol2(3)", -5.,5.);

  auto c0 = new TCanvas("c0", "c0", 600*3, 600*2);
  auto c1 = new TCanvas("c1", "c1", 600*3, 600);
  c0->Divide(3, 2);
  c1->Divide(3, 1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  for(Int_t ilayer = 0; ilayer <= 2; ilayer++)
  {
    for(Int_t itouch = 0; itouch <= max_ntouch; itouch++)
    {
      Double_t max = h_phi[ilayer][itouch]->GetMaximum();
      Double_t par[6] = {max,0.,1., 0.,0.,0.};
      fn_fit->SetParameters(par);
      h_phi[ilayer][itouch]->SetLineColor(itouch+1);
      fn_fit->SetLineColor(itouch+1);
      h_phi[ilayer][itouch]->Fit(fn_fit, "RQ0");
      c0->cd(ilayer+1);
      h_phi[ilayer][itouch]->DrawCopy(itouch==0?"":"SAME");
      fn_fit->DrawCopy("SAME");
      if(itouch > 0)
      {
        c0->cd(ilayer+4);
        h_phi[ilayer][itouch]->DrawCopy(itouch==1?"":"SAME");
        fn_fit->DrawCopy("SAME");
      }

      Double_t sigma = fn_fit->GetParameter(2);
      Double_t esigma = fn_fit->GetParError(2);
      g_sigma[ilayer]->SetPoint(itouch, itouch, sigma);
      g_sigma[ilayer]->SetPointError(itouch, 0., esigma);
    }

    c1->cd(ilayer+1);
    g_sigma[ilayer]->SetTitle(Form("Layer %d",ilayer));
    aset(g_sigma[ilayer], "ntouch","#sigma", -1.,5., 0.,0., 0.5,0.5);
    style(g_sigma[ilayer], 20, 1);
    g_sigma[ilayer]->Draw("AP");
  }

  c0->Print("results/pull-fit.pdf");
  c1->Print("results/pull-sigma.pdf");
}
