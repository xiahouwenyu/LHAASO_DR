#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>

#include"inc_astro/astro.h"
#include "inc_astro/papi.h"

using namespace std;
int main(int argc, char * argv[])
{ 

  double zen_src, azi_src, dec_src, ra_src;
  double mjd;
  ra_src  = 0.0;
  double dec_width=1.0;

  TFile *fout = new TFile("data/zenc.root","recreate");
  TH1D *h_zen[110];

  for(int i_dec=0;i_dec<110;i_dec++){
        cout<<i_dec<<endl;
        h_zen[i_dec]=new TH1D(Form("dec_%02d",i_dec),Form("dec_%02d",i_dec),90,0,90);
        dec_src=-25.+i_dec*dec_width;
        double bin=86400.-236.;
        for (int i=0; i<bin; i++){
          mjd=58949.000000+1.0/86400.*i+0.5/86400.;
          papi::eqm2hcs(mjd,0.,ra_src*papi::degrad,dec_src*papi::degrad,zen_src,azi_src);
          zen_src=zen_src*papi::raddeg;
          h_zen[i_dec]->Fill(zen_src);
        }
        for(int i=60;i<90;i++){
            h_zen[i_dec]->SetBinContent(i+1,0);
        }
    }
    fout->cd();
    for(int i=0;i<110;i++){
        h_zen[i]->Write();
    }
    fout->Close();

  return 0;     
}
 

