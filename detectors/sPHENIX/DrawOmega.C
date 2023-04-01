#include "DivideFunctions.h"

void DrawOmega()
{
  auto c0 = new TCanvas("c0", "c0", 3*600, 2*600);
  auto c1 = new TCanvas("c1", "c1", 600, 600);
  c0->Divide(3, 2);

  auto f = new TFile("/phenix/spin/phnxsp01/zji/data/sphenix/output/KFP.root");

  auto ntp_truth = (TNtuple*)f->Get("ntp_truth");
  auto ntp_kfp = (TNtuple*)f->Get("ntp_kfp");

  TString str_truth = "gpt_omega>0.5 && gpt_omega<10 && fabs(geta_omega)<0.8 && gpt_lambda > 0.5 && gpt_lambda<10 && fabs(geta_lambda)<0.8";
  TString str_dca = " && Lambda0_track_1_pT>0.2 && Lambda0_track_2_pT>0.2 && track_3_pT>0.2 && track_1_track_2_DCA<1.5 && Omegaminus_decayLength<Lambda0_decayLength && Omegaminus_DIRA>0.97 && Lambda0_DIRA>0.97";
  TString str_omega_mass = " && Omegaminus_mass>1.5 && Omegaminus_mass<1.9";
  TString str_omega_sig = " && Omegaminus_mass>1.64 && Omegaminus_mass<1.71";
  TString str_omega_bg = " && Omegaminus_mass>1.71 && Omegaminus_mass<1.78";
  TString str_lambda_mass = " && Lambda0_mass>1 && Lambda0_mass<1.3";
  TString str_lambda_sig = " && Lambda0_mass>1.09 && Lambda0_mass<1.14";

  for(Int_t ith=0; ith<2; ith++)
  {
    TString str_type = ith ? " with truth" : " without truth";
    TString str_dist2 = ith ? " && dist2_omega<0.1 && dist2_lambda<0.1" : "";

    c1->cd();
    ntp_truth->Draw(Form("gpt_omega >> h_gpt_%d(10,0.5,10.5)", ith), str_truth);
    ntp_kfp->Draw(Form("Omegaminus_pT >> h_pt_sig_%d(10,0.5,10.5)", ith), str_truth + str_dist2 + str_dca + str_lambda_sig + str_omega_sig);
    ntp_kfp->Draw(Form("Omegaminus_pT >> h_pt_bg_%d(10,0.5,10.5)", ith), str_truth + str_dist2 + str_dca + str_lambda_sig + str_omega_bg);
    ntp_kfp->Draw(Form("Omegaminus_mass >> h_minv_omega_%d(200,1.5,1.9)", ith), str_truth + str_dist2 + str_dca + str_lambda_sig + str_omega_mass);
    ntp_kfp->Draw(Form("Lambda0_mass >> h_minv_lambda_%d(200,1,1.3)", ith), str_truth + str_dist2 + str_dca + str_omega_sig + str_lambda_mass);

    auto h_gpt = (TH1*)gDirectory->Get(Form("h_gpt_%d", ith));
    auto h_pt_sig = (TH1*)gDirectory->Get(Form("h_pt_sig_%d", ith));
    auto h_pt_bg = (TH1*)gDirectory->Get(Form("h_pt_bg_%d", ith));
    auto h_minv_omega = (TH1*)gDirectory->Get(Form("h_minv_omega_%d", ith));
    auto h_minv_lambda = (TH1*)gDirectory->Get(Form("h_minv_lambda_%d", ith));
    h_pt_sig->Add(h_pt_bg, -1.);

    c0->cd(1+ith*3);
    h_minv_omega->SetTitle("#Omega mass" + str_type);
    h_minv_omega->GetXaxis()->SetTitle("m (GeV)");
    h_minv_omega->Draw();

    c0->cd(2+ith*3);
    h_minv_lambda->SetTitle("#Lambda mass" + str_type);
    h_minv_lambda->GetXaxis()->SetTitle("m (GeV)");
    h_minv_lambda->Draw();

    c0->cd(3+ith*3);
    auto g_eff = DivideHisto(h_pt_sig, h_gpt);
    g_eff->SetTitle("#Omega efficiency" + str_type);
    g_eff->GetXaxis()->SetTitle("p_{T} (GeV)");
    g_eff->Draw("AP");
  }

  c0->Print("results/omega-eff.pdf");
}
