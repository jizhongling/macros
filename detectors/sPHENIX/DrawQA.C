void PlotMulRes(TFile *f, const char *hname, const char *res_name, Double_t mean0, Double_t sigma0)
{
  auto h3_res = (TH3*)f->Get(Form("h3_%s", hname));
  auto g_res = new TGraphErrors(4);
  auto fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", mean0-sigma0*5, mean0+sigma0*5);

  auto c = new TCanvas("c","c", 3*600,2*600);
  c->Divide(3, 2);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  for(Int_t imul = 0; imul < 4; imul++)
  {
    c->cd(imul+1);
    auto h_res = h3_res->ProjectionY("h_res", 2,6, imul*20,(imul+1)*20);
    h_res->SetTitle(Form("ntrack: %dk--%dk", imul*20,(imul+1)*20));
    h_res->GetXaxis()->SetTitle(res_name);
    const Double_t max = h_res->GetMaximum();
    Double_t par[] = {max,mean0,sigma0, 0.,0.,0.};
    fn_fit->SetParameters(par);
    h_res->Fit(fn_fit, "RQ0");
    const Double_t sigma = TMath::Abs(fn_fit->GetParameter(2));
    const Double_t esigma = fn_fit->GetParError(2);
    h_res->DrawCopy();
    fn_fit->DrawCopy("SAME");
    g_res->SetPoint(imul, (imul+0.5)*20, sigma);
    g_res->SetPointError(imul, 0., esigma);
    delete h_res;
  }

  c->cd(5);
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  g_res->SetTitle(res_name);
  aset(g_res, "Multiplicity #times 10^{3}",res_name, 0.,1e5, 0.,sigma0*1.5);
  style(g_res, 20, 1);
  g_res->Draw("AP");

  c->Print(Form("results/resolution-%s.pdf", hname));
  delete c;
  delete fn_fit;
  delete g_res;
}

void DrawQA()
{
  const Int_t ntype = 1;
  const char *qa_type[] = {"", "-nn"};
  const char *event_type = "";
  const Int_t ntrack = 0;

  auto c0 = new TCanvas("c0","c0", 2*600,1*600);
  c0->Divide(2, 1);

  auto c1 = new TCanvas("c1","c1", 8*600,5*600);
  c1->Divide(8, 5);

  auto c2 = new TCanvas("c2","c2", 1*600,1*600);

  for(Int_t id=0; id<ntype; id++)
  {
    auto f = new TFile(Form("results/qa%s.root", qa_type[id]));
    auto h2_truth = (TH2*)f->Get("h2_truth");
    auto h2_reco = (TH2*)f->Get("h2_reco");
    auto h2_all = (TH2*)f->Get("h2_all");
    auto h2_good = (TH2*)f->Get("h2_good");
    auto h3_resol = (TH3*)f->Get("h3_resol");

    auto h_truth = h2_truth->ProjectionX("h_truth", ntrack, -1);
    auto h_reco = h2_reco->ProjectionX("h_reco", ntrack, -1);
    auto h_all = h2_all->ProjectionX("h_all", ntrack, -1);
    auto h_good = h2_good->ProjectionX("h_good", ntrack, -1);

    auto g_eff = new TGraphAsymmErrors(h_reco, h_truth);
    auto g_purity = new TGraphAsymmErrors(h_good, h_all);

    gStyle->SetOptStat(0);

    c0->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetLogx();
    g_eff->SetTitle("Efficiency");
    aset(g_eff, "p_{T} (GeV)","Eff", 0.3,20., 0.,1.1);
    style(g_eff, 24+id, 1+id);
    g_eff->Draw(id==0?"AP":"P");

    c0->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    g_purity->SetTitle("Purity");
    aset(g_purity, "p_{T} (GeV)","Purity", 0.,20., 0.,1.1);
    style(g_purity, 24+id, 1+id);
    g_purity->Draw(id==0?"AP":"P");

    auto g_resol = new TGraphErrors(40);
    auto fn_fit = new TF1("fn_fit", "gaus(0) + pol2(3)", 0.9, 1.1);

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    for(Int_t ipt = 0; ipt < 40; ipt++)
    {
      c1->cd(ipt+1);
      auto h_resol = h3_resol->ProjectionY("h_resol", ipt+1,ipt+1, 30,-1);
      h_resol->SetTitle(Form("p_{T}: %.1f--%.1f GeV", 0.5*ipt,0.5*(ipt+1)));
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

    c2->cd();
    gStyle->SetOptStat(0);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    g_resol->SetTitle("Momentum resolution");
    aset(g_resol, "p_{T} (GeV)","#sigma(p_{T})/p_{T}", 0.,20., 0.,0.04);
    style(g_resol, 24+id, 1+id);
    g_resol->Draw(id==0?"AP":"P");

    PlotMulRes(f, "resol", "#sigma(p_{T})/p_{T}", 1., 0.01);
    PlotMulRes(f, "eta", "#sigma(#eta)", 0., 0.002);
    PlotMulRes(f, "phi", "#sigma(#phi)", 0., 0.002);
  }

  c0->Print(Form("results/eff-purity%s.pdf", event_type));
  c1->Print(Form("results/resolution-histo%s.pdf", event_type));
  c2->Print(Form("results/resolution%s.pdf", event_type));
}
