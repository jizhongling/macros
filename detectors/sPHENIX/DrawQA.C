void DrawQA()
{
  auto f = new TFile("results/qa.root");
  auto h2_truth = (TH2*)f->Get("h2_truth");
  auto h2_reco = (TH2*)f->Get("h2_reco");
  auto h2_all = (TH2*)f->Get("h2_all");
  auto h2_good = (TH2*)f->Get("h2_good");
  auto h3_resol = (TH3*)f->Get("h3_resol");

  auto h_truth = h2_truth->ProjectionX("h_truth", 30, -1);
  auto h_reco = h2_reco->ProjectionX("h_reco", 30, -1);
  auto h_all = h2_all->ProjectionX("h_all", 30, -1);
  auto h_good = h2_good->ProjectionX("h_good", 30, -1);

  TGraphAsymmErrors *g_eff = new TGraphAsymmErrors(h_reco, h_truth);
  TGraphAsymmErrors *g_purity = new TGraphAsymmErrors(h_good, h_all);

  TCanvas *c0 = new TCanvas("c0","c0", 1200,600);
  c0->Divide(2, 1);
  gStyle->SetOptStat(0);

  c0->cd(1);
  gPad->SetLogx();
  g_eff->SetTitle("Efficiency");
  g_eff->GetXaxis()->SetTitle("p_{T} (GeV)");
  g_eff->GetYaxis()->SetTitle("Eff");
  g_eff->GetXaxis()->SetRangeUser(0.3, 20.);
  g_eff->GetYaxis()->SetRangeUser(0., 1.1);
  g_eff->Draw("AP");

  c0->cd(2);
  g_purity->SetTitle("Purity");
  g_purity->GetXaxis()->SetTitle("p_{T} (GeV)");
  g_purity->GetYaxis()->SetTitle("Purity");
  g_purity->GetYaxis()->SetRangeUser(0., 1.1);
  g_purity->Draw("AP");

  c0->Print("results/eff-purity.pdf");

  TGraphErrors *g_resol = new TGraphErrors(40);
  TF1 *fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.9, 1.1);

  TCanvas *c1 = new TCanvas("c1","c1", 4800,3000);
  c1->Divide(8, 5);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  for(int ipt = 0; ipt < 40; ipt++)
  {
    c1->cd(ipt+1);
    TH1 *h_resol = h3_resol->ProjectionY("h_resol", ipt+1,ipt+1, 30,-1);
    h_resol->SetTitle(Form("p_{T}: %.1f-%.1f GeV", 0.5*ipt,0.5*(ipt+1)));
    h_resol->GetXaxis()->SetTitle("pt/gpt");
    const Double_t max = h_resol->GetMaximum();
    Double_t par[] = {max,1.,0.1, 0.,0.,0.};
    fn_fit->SetParameters(par);
    h_resol->Fit(fn_fit, "RQ0");
    const Double_t sigma = TMath::Abs(fn_fit->GetParameter(2));
    const Double_t esigma = fn_fit->GetParError(2);
    h_resol->DrawCopy();
    fn_fit->DrawCopy("SAME");
    g_resol->SetPoint(ipt, 0.5*(ipt+0.5), sigma);
    g_resol->SetPointError(ipt, 0., esigma);
    delete h_resol;
  }
  c1->Print("results/resolution-histo.pdf");

  TCanvas *c2 = new TCanvas("c2","c2", 600,600);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  g_resol->GetYaxis()->SetRangeUser(0., 0.04);
  g_resol->SetTitle("Momentum resolution");
  g_resol->GetXaxis()->SetTitle("p_{T} (GeV)");
  g_resol->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T}");
  g_resol->Draw("AP");
  c2->Print("results/resolution.pdf");
}
