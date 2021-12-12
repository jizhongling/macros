// g++ -Wall -I$OFFLINE_MAIN/include -L$OFFLINE_MAIN/lib `root-config --cflags --glibs` -o AnaHitFile AnaHitFile.cpp
#include <iostream>

#include <TMath.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

using namespace std;

void query(TFile *f, const char *tname, TString var, const char *cond, vector<vector<Float_t>> &result)
{
  const Int_t nitem = var.CountChar(58) + 1;
  vector<Float_t> values(nitem);
  result.clear();

  TNtuple *ntp = static_cast<TNtuple*>(f->Get(tname));
  TSQLResult *res = ntp->Query(var, cond);
  TSQLRow *row;

  while( (row = res->Next()) )
  {
    bool store = true;
    values.clear();

    for(Int_t i=0; i<nitem; i++)
    {
      TString field = row->GetField(i);
      Float_t value = field.Atof();
      if( !TMath::Finite(value) )
      {
        store = false;
        break;
      }
      values.push_back(value);
    }

    if(store)
      result.push_back(values);

    delete row;
  }

  delete res;
}

int main(int argc, const char *argv[])
{
  if(argc != 2)
  {
    cerr << "Usage: " << argv[0] << " <path-to-eval-file>" << endl;
    return 1;
  }

  TFile *f = new TFile(argv[1]);
  if(f->IsZombie())
  {
    cerr << "Error: cannot open file " << argv[1] << endl;
    return 1;
  }
  cout << "Open file " << argv[1] << endl;

  vector<vector<Float_t>> v_g4cluster;
  vector<vector<Float_t>> v_hit;
  vector<vector<Float_t>>::iterator it_g4cluster;
  vector<vector<Float_t>>::iterator it_hit;

  Int_t layer = 35;
  query(f, "ntp_g4cluster", "event:gphi:gz", Form("event==0&&layer==%d", layer), v_g4cluster);
  for(it_g4cluster = v_g4cluster.begin(); it_g4cluster != v_g4cluster.end(); it_g4cluster++)
  {
    query(f, "ntp_hit", "phi:z:phibin:zbin", Form("event==%d&&layer==%d", Int_t(it_g4cluster->at(0)), layer), v_hit);
    for(it_hit = v_hit.begin(); it_hit != v_hit.end(); it_hit++)
      if( fabs(it_g4cluster->at(1) - it_hit->at(0)) < 0.05 && fabs(it_g4cluster->at(2) - it_hit->at(1)) < 0.1 )
        cout << it_hit->at(2) << " | " << it_hit->at(3) << endl;
  }

  delete f;
  return 0;
}
