void DrawRes()
{
  const char *dataset[3] = {"simple", "pp", "auau"};
  const char *type[3] = {"p_gp", "e_p", "e_p_conv"};

  TCanvas *c_tmp = new TCanvas("c_tmp", "c_tmp", 600, 600);
  TCanvas *c0 = new TCanvas("c0", "c0", 3*600, 3*600);
  c0->Divide(3, 3);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  Int_t ipad = 1;

  for(Int_t id=0; id<3; id++)
  {
    c_tmp->cd();
    auto f = new TFile(Form("results/hf-electron-%s.root", dataset[id]));
    auto ntp = (TNtuple*)f->Get("ntp");
    ntp->Draw("sqrt(pt*pt+pz*pz)/sqrt(gpt*gpt+gpz*gpz) >> h_p_gp", "fabs(sqrt(pt*pt+pz*pz)/sqrt(gpt*gpt+gpz*gpz)-1)<0.5");
    ntp->Draw("clus_e_cemc/sqrt(pt*pt+pz*pz) >> h_e_p", "fabs(clus_e_cemc/sqrt(pt*pt+pz*pz)-1)<0.5");
    if(id == 2)
      ntp->Draw("clus_e_cemc/sqrt(pt*pt+pz*pz) >> h_e_p_conv", "fabs(clus_e_cemc/sqrt(pt*pt+pz*pz)-1)<0.5 && gparentPID==22");

    for(Int_t it=0; it<3; it++)
    {
      mcd(0, ipad++);
      TH1* h_res = (TH1*)gDirectory->Get(Form("h_%s", type[it]));
      if(!h_res) continue;
      h_res->SetTitle(dataset[id]);
      TString tmp = type[it];
      aset(h_res, tmp.Replace(1,1,"/").Data(),"");
      h_res->Draw();
    }
  }
  c0->Print("results/HFRes.pdf");
}
