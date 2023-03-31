// g++ -std=c++17 -Wall `root-config --cflags --glibs` -o AnaOmega AnaOmega.cpp
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>

using namespace std;

template<class T> inline T square(T x)
{
  return x*x;
}

template<class T> inline T dist2(T eta0, T phi0, T eta1, T phi1)
{
  return square(eta0-eta1) + square(phi0-phi1);
}

int main(int argc, const char *argv[])
{
  if(argc != 5)
  {
    cerr << "Usage: " << argv[0] << " <dir> <start-index> <end-index> <out-file>" << endl;
    return 1;
  }
  const int istart = stoi(string(argv[2]));
  const int iend = stoi(string(argv[3]));

  auto f_out = new TFile(argv[4], "RECREATE");
  auto ntp_truth = new TNtuple("ntp_truth", "Truth tree", "gpt_omega:geta_omega:gphi_omega:gpt_lambda:geta_lambda:gphi_lambda");
  auto ntp_kfp = new TNtuple("ntp_kfp", "KFParticle decay tree", "gpt_omega:geta_omega:gphi_omega:gpt_lambda:geta_lambda:gphi_lambda:dist2_omega:dist2_lambda:Omegaminus_mass:Omegaminus_massErr:Omegaminus_decayTime:Omegaminus_decayTimeErr:Omegaminus_decayLength:Omegaminus_decayLengthErr:Omegaminus_DIRA:Omegaminus_FDchi2:Omegaminus_IP:Omegaminus_IPchi2:Omegaminus_IPErr:Omegaminus_IP_xy:Omegaminus_x:Omegaminus_y:Omegaminus_z:Omegaminus_px:Omegaminus_py:Omegaminus_pz:Omegaminus_pE:Omegaminus_p:Omegaminus_pErr:Omegaminus_pT:Omegaminus_pTErr:Omegaminus_pseudorapidity:Omegaminus_rapidity:Omegaminus_theta:Omegaminus_phi:Omegaminus_vertex_volume:Omegaminus_chi2:Lambda0_mass:Lambda0_massErr:Lambda0_decayTime:Lambda0_decayTimeErr:Lambda0_decayLength:Lambda0_decayLengthErr:Lambda0_DIRA:Lambda0_FDchi2:Lambda0_IP:Lambda0_IPchi2:Lambda0_IPErr:Lambda0_IP_xy:Lambda0_x:Lambda0_y:Lambda0_z:Lambda0_px:Lambda0_py:Lambda0_pz:Lambda0_pE:Lambda0_p:Lambda0_pErr:Lambda0_pT:Lambda0_pTErr:Lambda0_pseudorapidity:Lambda0_rapidity:Lambda0_theta:Lambda0_phi:Lambda0_vertex_volume:Lambda0_chi2:Lambda0_track_1_mass:Lambda0_track_1_IP:Lambda0_track_1_IPchi2:Lambda0_track_1_IPErr:Lambda0_track_1_IP_xy:Lambda0_track_1_x:Lambda0_track_1_y:Lambda0_track_1_z:Lambda0_track_1_px:Lambda0_track_1_py:Lambda0_track_1_pz:Lambda0_track_1_pE:Lambda0_track_1_p:Lambda0_track_1_pErr:Lambda0_track_1_pT:Lambda0_track_1_pTErr:Lambda0_track_1_jT:Lambda0_track_1_pseudorapidity:Lambda0_track_1_rapidity:Lambda0_track_1_theta:Lambda0_track_1_phi:Lambda0_track_1_chi2:Lambda0_track_2_mass:Lambda0_track_2_IP:Lambda0_track_2_IPchi2:Lambda0_track_2_IPErr:Lambda0_track_2_IP_xy:Lambda0_track_2_x:Lambda0_track_2_y:Lambda0_track_2_z:Lambda0_track_2_px:Lambda0_track_2_py:Lambda0_track_2_pz:Lambda0_track_2_pE:Lambda0_track_2_p:Lambda0_track_2_pErr:Lambda0_track_2_pT:Lambda0_track_2_pTErr:Lambda0_track_2_jT:Lambda0_track_2_pseudorapidity:Lambda0_track_2_rapidity:Lambda0_track_2_theta:Lambda0_track_2_phi:Lambda0_track_2_chi2:track_3_mass:track_3_IP:track_3_IPchi2:track_3_IPErr:track_3_IP_xy:track_3_x:track_3_y:track_3_z:track_3_px:track_3_py:track_3_pz:track_3_pE:track_3_p:track_3_pErr:track_3_pT:track_3_pTErr:track_3_jT:track_3_pseudorapidity:track_3_rapidity:track_3_theta:track_3_phi:track_3_chi2:track_1_track_2_DCA:track_1_track_3_DCA:track_2_track_3_DCA:primary_vertex_x:primary_vertex_y:primary_vertex_z:primary_vertex_volume:primary_vertex_chi2:secondary_vertex_mass_pionPID");

  for(int ifile = istart; ifile < iend; ifile++)
  {
    const char *fname[2] = {"G4sPHENIX_g4svtx_eval", "KFParticleOutput_Omega-_reconstruction"};
    TFile *f[2];
    bool badfile = false;
    for(int i=0; i<2; i++)
    {
      char filename[200];
      sprintf(filename, "%s/%s-%d.root", argv[1], fname[i], ifile);
      f[i] = TFile::Open(filename);
      if(!f[i] || f[i]->IsZombie())
      {
        cout << "Error: cannot open file " << filename << endl;
        if(f[i]) delete f[i];
        badfile = true;
        break;
      }
    }
    if(badfile) continue;

    auto ntp_gtrack = (TNtuple*)f[0]->Get("ntp_gtrack");
    if(!ntp_gtrack || ntp_gtrack->IsZombie())
    {
      cout << "Error: file has no ntp_gtrack" << endl;
      continue;
    }
    Float_t event, gtrackID, gflavor, gpt, geta, gphi;
    ntp_gtrack->SetBranchAddress("event", &event);
    ntp_gtrack->SetBranchAddress("gtrackID", &gtrackID);
    ntp_gtrack->SetBranchAddress("gflavor", &gflavor);
    ntp_gtrack->SetBranchAddress("gpt", &gpt);
    ntp_gtrack->SetBranchAddress("geta", &geta);
    ntp_gtrack->SetBranchAddress("gphi", &gphi);

    auto DecayTree = (TTree*)f[1]->Get("DecayTree");
    if(!DecayTree || DecayTree->IsZombie())
    {
      cout << "Error: file has no DecayTree" << endl;
      continue;
    }
    Int_t Omegaminus_charge, Omegaminus_nDoF, Omegaminus_PDG_ID, Lambda0_charge, Lambda0_nDoF, Lambda0_PDG_ID, Lambda0_track_1_charge, Lambda0_track_1_nDoF, Lambda0_track_1_track_ID, Lambda0_track_1_PDG_ID, Lambda0_track_2_charge, Lambda0_track_2_nDoF, Lambda0_track_2_track_ID, Lambda0_track_2_PDG_ID, track_3_charge, track_3_nDoF, track_3_track_ID, track_3_PDG_ID, primary_vertex_nTracks, primary_vertex_nDoF, nPrimaryVertices, nEventTracks, runNumber, eventNumber;
    Float_t Omegaminus_mass, Omegaminus_massErr, Omegaminus_decayTime, Omegaminus_decayTimeErr, Omegaminus_decayLength, Omegaminus_decayLengthErr, Omegaminus_DIRA, Omegaminus_FDchi2, Omegaminus_IP, Omegaminus_IPchi2, Omegaminus_IPErr, Omegaminus_IP_xy, Omegaminus_x, Omegaminus_y, Omegaminus_z, Omegaminus_px, Omegaminus_py, Omegaminus_pz, Omegaminus_pE, Omegaminus_p, Omegaminus_pErr, Omegaminus_pT, Omegaminus_pTErr, Omegaminus_pseudorapidity, Omegaminus_rapidity, Omegaminus_theta, Omegaminus_phi, Omegaminus_vertex_volume, Omegaminus_chi2, Lambda0_mass, Lambda0_massErr, Lambda0_decayTime, Lambda0_decayTimeErr, Lambda0_decayLength, Lambda0_decayLengthErr, Lambda0_DIRA, Lambda0_FDchi2, Lambda0_IP, Lambda0_IPchi2, Lambda0_IPErr, Lambda0_IP_xy, Lambda0_x, Lambda0_y, Lambda0_z, Lambda0_px, Lambda0_py, Lambda0_pz, Lambda0_pE, Lambda0_p, Lambda0_pErr, Lambda0_pT, Lambda0_pTErr, Lambda0_pseudorapidity, Lambda0_rapidity, Lambda0_theta, Lambda0_phi, Lambda0_vertex_volume, Lambda0_chi2, Lambda0_track_1_mass, Lambda0_track_1_IP, Lambda0_track_1_IPchi2, Lambda0_track_1_IPErr, Lambda0_track_1_IP_xy, Lambda0_track_1_x, Lambda0_track_1_y, Lambda0_track_1_z, Lambda0_track_1_px, Lambda0_track_1_py, Lambda0_track_1_pz, Lambda0_track_1_pE, Lambda0_track_1_p, Lambda0_track_1_pErr, Lambda0_track_1_pT, Lambda0_track_1_pTErr, Lambda0_track_1_jT, Lambda0_track_1_pseudorapidity, Lambda0_track_1_rapidity, Lambda0_track_1_theta, Lambda0_track_1_phi, Lambda0_track_1_chi2, Lambda0_track_2_mass, Lambda0_track_2_IP, Lambda0_track_2_IPchi2, Lambda0_track_2_IPErr, Lambda0_track_2_IP_xy, Lambda0_track_2_x, Lambda0_track_2_y, Lambda0_track_2_z, Lambda0_track_2_px, Lambda0_track_2_py, Lambda0_track_2_pz, Lambda0_track_2_pE, Lambda0_track_2_p, Lambda0_track_2_pErr, Lambda0_track_2_pT, Lambda0_track_2_pTErr, Lambda0_track_2_jT, Lambda0_track_2_pseudorapidity, Lambda0_track_2_rapidity, Lambda0_track_2_theta, Lambda0_track_2_phi, Lambda0_track_2_chi2, track_3_mass, track_3_IP, track_3_IPchi2, track_3_IPErr, track_3_IP_xy, track_3_x, track_3_y, track_3_z, track_3_px, track_3_py, track_3_pz, track_3_pE, track_3_p, track_3_pErr, track_3_pT, track_3_pTErr, track_3_jT, track_3_pseudorapidity, track_3_rapidity, track_3_theta, track_3_phi, track_3_chi2, track_1_track_2_DCA, track_1_track_3_DCA, track_2_track_3_DCA, primary_vertex_x, primary_vertex_y, primary_vertex_z, primary_vertex_volume, primary_vertex_chi2, secondary_vertex_mass_pionPID;
    DecayTree->SetBranchAddress("Omegaminus_mass", &Omegaminus_mass);
    DecayTree->SetBranchAddress("Omegaminus_massErr", &Omegaminus_massErr);
    DecayTree->SetBranchAddress("Omegaminus_decayTime", &Omegaminus_decayTime);
    DecayTree->SetBranchAddress("Omegaminus_decayTimeErr", &Omegaminus_decayTimeErr);
    DecayTree->SetBranchAddress("Omegaminus_decayLength", &Omegaminus_decayLength);
    DecayTree->SetBranchAddress("Omegaminus_decayLengthErr", &Omegaminus_decayLengthErr);
    DecayTree->SetBranchAddress("Omegaminus_DIRA", &Omegaminus_DIRA);
    DecayTree->SetBranchAddress("Omegaminus_FDchi2", &Omegaminus_FDchi2);
    DecayTree->SetBranchAddress("Omegaminus_IP", &Omegaminus_IP);
    DecayTree->SetBranchAddress("Omegaminus_IPchi2", &Omegaminus_IPchi2);
    DecayTree->SetBranchAddress("Omegaminus_IPErr", &Omegaminus_IPErr);
    DecayTree->SetBranchAddress("Omegaminus_IP_xy", &Omegaminus_IP_xy);
    DecayTree->SetBranchAddress("Omegaminus_x", &Omegaminus_x);
    DecayTree->SetBranchAddress("Omegaminus_y", &Omegaminus_y);
    DecayTree->SetBranchAddress("Omegaminus_z", &Omegaminus_z);
    DecayTree->SetBranchAddress("Omegaminus_px", &Omegaminus_px);
    DecayTree->SetBranchAddress("Omegaminus_py", &Omegaminus_py);
    DecayTree->SetBranchAddress("Omegaminus_pz", &Omegaminus_pz);
    DecayTree->SetBranchAddress("Omegaminus_pE", &Omegaminus_pE);
    DecayTree->SetBranchAddress("Omegaminus_p", &Omegaminus_p);
    DecayTree->SetBranchAddress("Omegaminus_pErr", &Omegaminus_pErr);
    DecayTree->SetBranchAddress("Omegaminus_pT", &Omegaminus_pT);
    DecayTree->SetBranchAddress("Omegaminus_pTErr", &Omegaminus_pTErr);
    DecayTree->SetBranchAddress("Omegaminus_charge", &Omegaminus_charge);
    DecayTree->SetBranchAddress("Omegaminus_pseudorapidity", &Omegaminus_pseudorapidity);
    DecayTree->SetBranchAddress("Omegaminus_rapidity", &Omegaminus_rapidity);
    DecayTree->SetBranchAddress("Omegaminus_theta", &Omegaminus_theta);
    DecayTree->SetBranchAddress("Omegaminus_phi", &Omegaminus_phi);
    DecayTree->SetBranchAddress("Omegaminus_vertex_volume", &Omegaminus_vertex_volume);
    DecayTree->SetBranchAddress("Omegaminus_chi2", &Omegaminus_chi2);
    DecayTree->SetBranchAddress("Omegaminus_nDoF", &Omegaminus_nDoF);
    DecayTree->SetBranchAddress("Omegaminus_PDG_ID", &Omegaminus_PDG_ID);
    DecayTree->SetBranchAddress("Lambda0_mass", &Lambda0_mass);
    DecayTree->SetBranchAddress("Lambda0_massErr", &Lambda0_massErr);
    DecayTree->SetBranchAddress("Lambda0_decayTime", &Lambda0_decayTime);
    DecayTree->SetBranchAddress("Lambda0_decayTimeErr", &Lambda0_decayTimeErr);
    DecayTree->SetBranchAddress("Lambda0_decayLength", &Lambda0_decayLength);
    DecayTree->SetBranchAddress("Lambda0_decayLengthErr", &Lambda0_decayLengthErr);
    DecayTree->SetBranchAddress("Lambda0_DIRA", &Lambda0_DIRA);
    DecayTree->SetBranchAddress("Lambda0_FDchi2", &Lambda0_FDchi2);
    DecayTree->SetBranchAddress("Lambda0_IP", &Lambda0_IP);
    DecayTree->SetBranchAddress("Lambda0_IPchi2", &Lambda0_IPchi2);
    DecayTree->SetBranchAddress("Lambda0_IPErr", &Lambda0_IPErr);
    DecayTree->SetBranchAddress("Lambda0_IP_xy", &Lambda0_IP_xy);
    DecayTree->SetBranchAddress("Lambda0_x", &Lambda0_x);
    DecayTree->SetBranchAddress("Lambda0_y", &Lambda0_y);
    DecayTree->SetBranchAddress("Lambda0_z", &Lambda0_z);
    DecayTree->SetBranchAddress("Lambda0_px", &Lambda0_px);
    DecayTree->SetBranchAddress("Lambda0_py", &Lambda0_py);
    DecayTree->SetBranchAddress("Lambda0_pz", &Lambda0_pz);
    DecayTree->SetBranchAddress("Lambda0_pE", &Lambda0_pE);
    DecayTree->SetBranchAddress("Lambda0_p", &Lambda0_p);
    DecayTree->SetBranchAddress("Lambda0_pErr", &Lambda0_pErr);
    DecayTree->SetBranchAddress("Lambda0_pT", &Lambda0_pT);
    DecayTree->SetBranchAddress("Lambda0_pTErr", &Lambda0_pTErr);
    DecayTree->SetBranchAddress("Lambda0_charge", &Lambda0_charge);
    DecayTree->SetBranchAddress("Lambda0_pseudorapidity", &Lambda0_pseudorapidity);
    DecayTree->SetBranchAddress("Lambda0_rapidity", &Lambda0_rapidity);
    DecayTree->SetBranchAddress("Lambda0_theta", &Lambda0_theta);
    DecayTree->SetBranchAddress("Lambda0_phi", &Lambda0_phi);
    DecayTree->SetBranchAddress("Lambda0_vertex_volume", &Lambda0_vertex_volume);
    DecayTree->SetBranchAddress("Lambda0_chi2", &Lambda0_chi2);
    DecayTree->SetBranchAddress("Lambda0_nDoF", &Lambda0_nDoF);
    DecayTree->SetBranchAddress("Lambda0_PDG_ID", &Lambda0_PDG_ID);
    DecayTree->SetBranchAddress("Lambda0_track_1_mass", &Lambda0_track_1_mass);
    DecayTree->SetBranchAddress("Lambda0_track_1_IP", &Lambda0_track_1_IP);
    DecayTree->SetBranchAddress("Lambda0_track_1_IPchi2", &Lambda0_track_1_IPchi2);
    DecayTree->SetBranchAddress("Lambda0_track_1_IPErr", &Lambda0_track_1_IPErr);
    DecayTree->SetBranchAddress("Lambda0_track_1_IP_xy", &Lambda0_track_1_IP_xy);
    DecayTree->SetBranchAddress("Lambda0_track_1_x", &Lambda0_track_1_x);
    DecayTree->SetBranchAddress("Lambda0_track_1_y", &Lambda0_track_1_y);
    DecayTree->SetBranchAddress("Lambda0_track_1_z", &Lambda0_track_1_z);
    DecayTree->SetBranchAddress("Lambda0_track_1_px", &Lambda0_track_1_px);
    DecayTree->SetBranchAddress("Lambda0_track_1_py", &Lambda0_track_1_py);
    DecayTree->SetBranchAddress("Lambda0_track_1_pz", &Lambda0_track_1_pz);
    DecayTree->SetBranchAddress("Lambda0_track_1_pE", &Lambda0_track_1_pE);
    DecayTree->SetBranchAddress("Lambda0_track_1_p", &Lambda0_track_1_p);
    DecayTree->SetBranchAddress("Lambda0_track_1_pErr", &Lambda0_track_1_pErr);
    DecayTree->SetBranchAddress("Lambda0_track_1_pT", &Lambda0_track_1_pT);
    DecayTree->SetBranchAddress("Lambda0_track_1_pTErr", &Lambda0_track_1_pTErr);
    DecayTree->SetBranchAddress("Lambda0_track_1_jT", &Lambda0_track_1_jT);
    DecayTree->SetBranchAddress("Lambda0_track_1_charge", &Lambda0_track_1_charge);
    DecayTree->SetBranchAddress("Lambda0_track_1_pseudorapidity", &Lambda0_track_1_pseudorapidity);
    DecayTree->SetBranchAddress("Lambda0_track_1_rapidity", &Lambda0_track_1_rapidity);
    DecayTree->SetBranchAddress("Lambda0_track_1_theta", &Lambda0_track_1_theta);
    DecayTree->SetBranchAddress("Lambda0_track_1_phi", &Lambda0_track_1_phi);
    DecayTree->SetBranchAddress("Lambda0_track_1_chi2", &Lambda0_track_1_chi2);
    DecayTree->SetBranchAddress("Lambda0_track_1_nDoF", &Lambda0_track_1_nDoF);
    DecayTree->SetBranchAddress("Lambda0_track_1_track_ID", &Lambda0_track_1_track_ID);
    DecayTree->SetBranchAddress("Lambda0_track_1_PDG_ID", &Lambda0_track_1_PDG_ID);
    DecayTree->SetBranchAddress("Lambda0_track_2_mass", &Lambda0_track_2_mass);
    DecayTree->SetBranchAddress("Lambda0_track_2_IP", &Lambda0_track_2_IP);
    DecayTree->SetBranchAddress("Lambda0_track_2_IPchi2", &Lambda0_track_2_IPchi2);
    DecayTree->SetBranchAddress("Lambda0_track_2_IPErr", &Lambda0_track_2_IPErr);
    DecayTree->SetBranchAddress("Lambda0_track_2_IP_xy", &Lambda0_track_2_IP_xy);
    DecayTree->SetBranchAddress("Lambda0_track_2_x", &Lambda0_track_2_x);
    DecayTree->SetBranchAddress("Lambda0_track_2_y", &Lambda0_track_2_y);
    DecayTree->SetBranchAddress("Lambda0_track_2_z", &Lambda0_track_2_z);
    DecayTree->SetBranchAddress("Lambda0_track_2_px", &Lambda0_track_2_px);
    DecayTree->SetBranchAddress("Lambda0_track_2_py", &Lambda0_track_2_py);
    DecayTree->SetBranchAddress("Lambda0_track_2_pz", &Lambda0_track_2_pz);
    DecayTree->SetBranchAddress("Lambda0_track_2_pE", &Lambda0_track_2_pE);
    DecayTree->SetBranchAddress("Lambda0_track_2_p", &Lambda0_track_2_p);
    DecayTree->SetBranchAddress("Lambda0_track_2_pErr", &Lambda0_track_2_pErr);
    DecayTree->SetBranchAddress("Lambda0_track_2_pT", &Lambda0_track_2_pT);
    DecayTree->SetBranchAddress("Lambda0_track_2_pTErr", &Lambda0_track_2_pTErr);
    DecayTree->SetBranchAddress("Lambda0_track_2_jT", &Lambda0_track_2_jT);
    DecayTree->SetBranchAddress("Lambda0_track_2_charge", &Lambda0_track_2_charge);
    DecayTree->SetBranchAddress("Lambda0_track_2_pseudorapidity", &Lambda0_track_2_pseudorapidity);
    DecayTree->SetBranchAddress("Lambda0_track_2_rapidity", &Lambda0_track_2_rapidity);
    DecayTree->SetBranchAddress("Lambda0_track_2_theta", &Lambda0_track_2_theta);
    DecayTree->SetBranchAddress("Lambda0_track_2_phi", &Lambda0_track_2_phi);
    DecayTree->SetBranchAddress("Lambda0_track_2_chi2", &Lambda0_track_2_chi2);
    DecayTree->SetBranchAddress("Lambda0_track_2_nDoF", &Lambda0_track_2_nDoF);
    DecayTree->SetBranchAddress("Lambda0_track_2_track_ID", &Lambda0_track_2_track_ID);
    DecayTree->SetBranchAddress("Lambda0_track_2_PDG_ID", &Lambda0_track_2_PDG_ID);
    DecayTree->SetBranchAddress("track_3_mass", &track_3_mass);
    DecayTree->SetBranchAddress("track_3_IP", &track_3_IP);
    DecayTree->SetBranchAddress("track_3_IPchi2", &track_3_IPchi2);
    DecayTree->SetBranchAddress("track_3_IPErr", &track_3_IPErr);
    DecayTree->SetBranchAddress("track_3_IP_xy", &track_3_IP_xy);
    DecayTree->SetBranchAddress("track_3_x", &track_3_x);
    DecayTree->SetBranchAddress("track_3_y", &track_3_y);
    DecayTree->SetBranchAddress("track_3_z", &track_3_z);
    DecayTree->SetBranchAddress("track_3_px", &track_3_px);
    DecayTree->SetBranchAddress("track_3_py", &track_3_py);
    DecayTree->SetBranchAddress("track_3_pz", &track_3_pz);
    DecayTree->SetBranchAddress("track_3_pE", &track_3_pE);
    DecayTree->SetBranchAddress("track_3_p", &track_3_p);
    DecayTree->SetBranchAddress("track_3_pErr", &track_3_pErr);
    DecayTree->SetBranchAddress("track_3_pT", &track_3_pT);
    DecayTree->SetBranchAddress("track_3_pTErr", &track_3_pTErr);
    DecayTree->SetBranchAddress("track_3_jT", &track_3_jT);
    DecayTree->SetBranchAddress("track_3_charge", &track_3_charge);
    DecayTree->SetBranchAddress("track_3_pseudorapidity", &track_3_pseudorapidity);
    DecayTree->SetBranchAddress("track_3_rapidity", &track_3_rapidity);
    DecayTree->SetBranchAddress("track_3_theta", &track_3_theta);
    DecayTree->SetBranchAddress("track_3_phi", &track_3_phi);
    DecayTree->SetBranchAddress("track_3_chi2", &track_3_chi2);
    DecayTree->SetBranchAddress("track_3_nDoF", &track_3_nDoF);
    DecayTree->SetBranchAddress("track_3_track_ID", &track_3_track_ID);
    DecayTree->SetBranchAddress("track_3_PDG_ID", &track_3_PDG_ID);
    DecayTree->SetBranchAddress("track_1_track_2_DCA", &track_1_track_2_DCA);
    DecayTree->SetBranchAddress("track_1_track_3_DCA", &track_1_track_3_DCA);
    DecayTree->SetBranchAddress("track_2_track_3_DCA", &track_2_track_3_DCA);
    DecayTree->SetBranchAddress("primary_vertex_x", &primary_vertex_x);
    DecayTree->SetBranchAddress("primary_vertex_y", &primary_vertex_y);
    DecayTree->SetBranchAddress("primary_vertex_z", &primary_vertex_z);
    DecayTree->SetBranchAddress("primary_vertex_nTracks", &primary_vertex_nTracks);
    DecayTree->SetBranchAddress("primary_vertex_volume", &primary_vertex_volume);
    DecayTree->SetBranchAddress("primary_vertex_chi2", &primary_vertex_chi2);
    DecayTree->SetBranchAddress("primary_vertex_nDoF", &primary_vertex_nDoF);
    DecayTree->SetBranchAddress("secondary_vertex_mass_pionPID", &secondary_vertex_mass_pionPID);
    DecayTree->SetBranchAddress("nPrimaryVertices", &nPrimaryVertices);
    DecayTree->SetBranchAddress("nEventTracks", &nEventTracks);
    DecayTree->SetBranchAddress("runNumber", &runNumber);
    DecayTree->SetBranchAddress("eventNumber", &eventNumber);

    Long64_t ien_gtr = 0, ien_kfp = 0;
    Int_t this_event = -1;
    Float_t gpt_omega, geta_omega, gphi_omega, gpt_lambda, geta_lambda, gphi_lambda;
    bool has_omega, has_lambda, has_proton, has_kaon, filled;
    for(; ien_gtr<ntp_gtrack->GetEntries(); ien_gtr++)
    {
      ntp_gtrack->GetEntry(ien_gtr);
      if((Int_t)event != this_event)
      {
        this_event = event;
        has_omega = false;
        has_lambda = false;
        has_proton = false;
        has_kaon = false;
        filled = false;
      }
      else if(filled) continue;

      if(gflavor==3334 && gtrackID==1)
      {
        gpt_omega = gpt;
        geta_omega = geta;
        gphi_omega = gphi;
        has_omega = true;
      }
      else if(gflavor==3122 && gtrackID>-3)
      {
        gpt_lambda = gpt;
        geta_lambda = geta;
        gphi_lambda = gphi;
        has_lambda = true;
      }
      else if(gflavor==2212)
        has_proton = true;
      else if(gflavor==-321)
        has_kaon = true;

      if( !has_omega || !has_lambda || !has_proton || !has_kaon )
        continue;

      Float_t fill_truth[] = {gpt_omega, geta_omega, gphi_omega, gpt_lambda, geta_lambda, gphi_lambda};
      ntp_truth->Fill(fill_truth);
      filled = true;

      for(; ien_kfp<DecayTree->GetEntries(); ien_kfp++)
      {
        DecayTree->GetEntry(ien_kfp);
        if(eventNumber-1 < this_event)
          continue;
        else if(eventNumber-1 > this_event)
          break;

        Float_t dist2_omega = dist2(Omegaminus_pseudorapidity, Omegaminus_phi, geta_omega, gphi_omega);
        Float_t dist2_lambda = dist2(Lambda0_pseudorapidity, Lambda0_phi, geta_lambda, gphi_lambda);
        Float_t fill_kfp[] = {gpt_omega, geta_omega, gphi_omega, gpt_lambda, geta_lambda, gphi_lambda, dist2_omega, dist2_lambda, Omegaminus_mass, Omegaminus_massErr, Omegaminus_decayTime, Omegaminus_decayTimeErr, Omegaminus_decayLength, Omegaminus_decayLengthErr, Omegaminus_DIRA, Omegaminus_FDchi2, Omegaminus_IP, Omegaminus_IPchi2, Omegaminus_IPErr, Omegaminus_IP_xy, Omegaminus_x, Omegaminus_y, Omegaminus_z, Omegaminus_px, Omegaminus_py, Omegaminus_pz, Omegaminus_pE, Omegaminus_p, Omegaminus_pErr, Omegaminus_pT, Omegaminus_pTErr, Omegaminus_pseudorapidity, Omegaminus_rapidity, Omegaminus_theta, Omegaminus_phi, Omegaminus_vertex_volume, Omegaminus_chi2, Lambda0_mass, Lambda0_massErr, Lambda0_decayTime, Lambda0_decayTimeErr, Lambda0_decayLength, Lambda0_decayLengthErr, Lambda0_DIRA, Lambda0_FDchi2, Lambda0_IP, Lambda0_IPchi2, Lambda0_IPErr, Lambda0_IP_xy, Lambda0_x, Lambda0_y, Lambda0_z, Lambda0_px, Lambda0_py, Lambda0_pz, Lambda0_pE, Lambda0_p, Lambda0_pErr, Lambda0_pT, Lambda0_pTErr, Lambda0_pseudorapidity, Lambda0_rapidity, Lambda0_theta, Lambda0_phi, Lambda0_vertex_volume, Lambda0_chi2, Lambda0_track_1_mass, Lambda0_track_1_IP, Lambda0_track_1_IPchi2, Lambda0_track_1_IPErr, Lambda0_track_1_IP_xy, Lambda0_track_1_x, Lambda0_track_1_y, Lambda0_track_1_z, Lambda0_track_1_px, Lambda0_track_1_py, Lambda0_track_1_pz, Lambda0_track_1_pE, Lambda0_track_1_p, Lambda0_track_1_pErr, Lambda0_track_1_pT, Lambda0_track_1_pTErr, Lambda0_track_1_jT, Lambda0_track_1_pseudorapidity, Lambda0_track_1_rapidity, Lambda0_track_1_theta, Lambda0_track_1_phi, Lambda0_track_1_chi2, Lambda0_track_2_mass, Lambda0_track_2_IP, Lambda0_track_2_IPchi2, Lambda0_track_2_IPErr, Lambda0_track_2_IP_xy, Lambda0_track_2_x, Lambda0_track_2_y, Lambda0_track_2_z, Lambda0_track_2_px, Lambda0_track_2_py, Lambda0_track_2_pz, Lambda0_track_2_pE, Lambda0_track_2_p, Lambda0_track_2_pErr, Lambda0_track_2_pT, Lambda0_track_2_pTErr, Lambda0_track_2_jT, Lambda0_track_2_pseudorapidity, Lambda0_track_2_rapidity, Lambda0_track_2_theta, Lambda0_track_2_phi, Lambda0_track_2_chi2, track_3_mass, track_3_IP, track_3_IPchi2, track_3_IPErr, track_3_IP_xy, track_3_x, track_3_y, track_3_z, track_3_px, track_3_py, track_3_pz, track_3_pE, track_3_p, track_3_pErr, track_3_pT, track_3_pTErr, track_3_jT, track_3_pseudorapidity, track_3_rapidity, track_3_theta, track_3_phi, track_3_chi2, track_1_track_2_DCA, track_1_track_3_DCA, track_2_track_3_DCA, primary_vertex_x, primary_vertex_y, primary_vertex_z, primary_vertex_volume, primary_vertex_chi2, secondary_vertex_mass_pionPID};
        ntp_kfp->Fill(fill_kfp);
      } // ien_kfp
    } // ien_gtr

    for(int i=0; i<2; i++)
      delete f[i];
  } // ifile

  f_out->Write();
  f_out->Close();
  return 0;
}
