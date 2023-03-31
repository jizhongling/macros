#include "DivideFunctions.h"

void DrawOmega()
{
  auto f = new TFile("/phenix/spin/phnxsp01/zji/data/sphenix/output/KFP.root");

  auto ntp_truth = (TNtuple*)f->Get("ntp_truth");
  auto ntp_kfp = (TNtuple*)f->Get("ntp_kfp");

  ntp_truth->Draw("gpt_omega >> h_gpt(10,0.5,10.5)", "gpt_omega>0.5 && gpt_omega<10 && fabs(geta_omega)<0.8 && gpt_lambda > 0.5 && gpt_lambda<10 && fabs(geta_lambda)<0.8");
  ntp_kfp->Draw("Omegaminus_pT >> h_pt(10,0.5,10.5)", "gpt_omega>0.5 && gpt_omega<10 && fabs(geta_omega)<0.8 && gpt_lambda > 0.5 && gpt_lambda<10 && fabs(geta_lambda)<0.8 && dist2_omega<0.1 && dist2_lambda<0.1 && Omegaminus_DIRA>0.97 && Lambda0_DIRA>0.97 && Omegaminus_decayLength<Lambda0_decayLength && Lambda0_track_1_pT>0.2 && Lambda0_track_2_pT>0.2 && track_3_pT>0.2 && track_1_track_2_DCA<1.5 && Lambda0_mass>1.09 && Lambda0_mass<1.14 && Omegaminus_mass>1.64 && Omegaminus_mass<1.71");
  ntp_kfp->Draw("Omegaminus_pT >> h_pt_bg(10,0.5,10.5)", "gpt_omega>0.5 && gpt_omega<10 && fabs(geta_omega)<0.8 && gpt_lambda > 0.5 && gpt_lambda<10 && fabs(geta_lambda)<0.8 && dist2_omega<0.1 && dist2_lambda<0.1 && Omegaminus_DIRA>0.97 && Lambda0_DIRA>0.97 && Omegaminus_decayLength<Lambda0_decayLength && Lambda0_track_1_pT>0.2 && Lambda0_track_2_pT>0.2 && track_3_pT>0.2 && track_1_track_2_DCA<1.5 && Lambda0_mass>1.09 && Lambda0_mass<1.14 && Omegaminus_mass>1.71 && Omegaminus_mass<1.78");
  ntp_kfp->Draw("Omegaminus_mass >> h_minv_omega(200,1.5,1.9)", "gpt_omega>0.5 && gpt_omega<10 && fabs(geta_omega)<0.8 && gpt_lambda > 0.5 && gpt_lambda<10 && fabs(geta_lambda)<0.8 && dist2_omega<0.1 && dist2_lambda<0.1 && Omegaminus_DIRA>0.97 && Lambda0_DIRA>0.97 && Omegaminus_decayLength<Lambda0_decayLength && Lambda0_track_1_pT>0.2 && Lambda0_track_2_pT>0.2 && track_3_pT>0.2 && track_1_track_2_DCA<1.5 && Lambda0_mass>1.09 && Lambda0_mass<1.14 && Omegaminus_mass>1.5 && Omegaminus_mass<1.9");
  ntp_kfp->Draw("Lambda0_mass >> h_minv_lambda(200,1,1.3)", "gpt_omega>0.5 && gpt_omega<10 && fabs(geta_omega)<0.8 && gpt_lambda > 0.5 && gpt_lambda<10 && fabs(geta_lambda)<0.8 && dist2_omega<0.1 && dist2_lambda<0.1 && Omegaminus_DIRA>0.97 && Lambda0_DIRA>0.97 && Omegaminus_decayLength<Lambda0_decayLength && Lambda0_track_1_pT>0.2 && Lambda0_track_2_pT>0.2 && track_3_pT>0.2 && track_1_track_2_DCA<1.5 && Omegaminus_mass>1.64 && Omegaminus_mass<1.71");

  auto h_gpt = (TH1*)gDirectory->Get("h_gpt");
  auto h_pt = (TH1*)gDirectory->Get("h_pt");
  auto h_pt_bg = (TH1*)gDirectory->Get("h_pt_bg");
  auto h_minv_omega = (TH1*)gDirectory->Get("h_minv_omega");
  auto h_minv_lambda = (TH1*)gDirectory->Get("h_minv_lambda");
  h_pt->Add(h_pt_bg, -1.);

  TCanvas *c0 = new TCanvas("c0", "c0", 3*600, 600);
  c0->Divide(3, 1);

  c0->cd(1);
  h_minv_omega->SetTitle("#Omega mass");
  h_minv_omega->GetXaxis()->SetTitle("m (GeV)");
  h_minv_omega->Draw();

  c0->cd(2);
  h_minv_lambda->SetTitle("#Lambda mass");
  h_minv_lambda->GetXaxis()->SetTitle("m (GeV)");
  h_minv_lambda->Draw();

  c0->cd(3);
  auto g_eff = DivideHisto(h_pt, h_gpt);
  g_eff->SetTitle("#Omega efficiency");
  g_eff->GetXaxis()->SetTitle("p_{T} (GeV)");
  g_eff->Draw("AP");

  c0->Print("results/omega-eff.pdf");
}
