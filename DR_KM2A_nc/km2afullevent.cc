#define km2afullevent_cxx
// #include "km2afullevent.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include "dataCut/EventCut.h"
#include "inc_astro/papi.h"
#include <TF1.h>
#include <TMath.h>

using namespace std;

namespace gv3{
    const int dec_bin=110;
    double dec_begin=-25.;
    double dec_width=1.;
    float E_min =  0;
    float E_max =  4.;
    const int  zen_bin =  50;
    const int    E_bin = 40;
    float E_binwid=(E_max-E_min)*1.0/E_bin;

    float psf_p0[EventCut::Nhit_bin]={1.,1.,1.,1.,1.,1., 1., 1., 1., 1., 1., 1., 1., 1.};
    float psf_p1[EventCut::Nhit_bin]={0.5808,  0.1433,   0.08849,    0.07374,     0.02708,    0.01834, 0.01415, 0.01795, 0.7076, 0.8052, 0.005663, 1, 0.07808, 0.07808};
    float psf_p2[EventCut::Nhit_bin]={0.2298,     1,         1,            1,           1,           1,     1,      0.2349, 0.01488, 0.009653, 0.6416, 0.0277, 0.003531, 0.003531};
    float psf_p3[EventCut::Nhit_bin]={0.8686,  0.2425,    0.1742,   0.1117,       0.07669,   0.05229, 0.04402, 0.02936, 0.0399, 0.02972, 0.02163, 0.0003678, 5.25e-6, 5.25e-6};
    
    float area = TMath::Pi() * 100000 * 100000;

    Double_t timearray3 =  710.746617;
    Double_t timearray2 = 215.895577;
    Double_t timearray1 = 289.596898;
    Double_t timeall = timearray3 + timearray2 + timearray1;
}

void km2afullevent::Loop(TH1D ***PSF_dec_nh, TF1  ***PSF_dec_nh_fit, TH1D ***EnSig_dec_nh)
{
    if (km2afullevent::fChain == 0) return;

    std::cout<<" *** Read weight from normalizing zenith angle distribution ..."<<std::endl;
    double Num_E_mczen[gv3::E_bin][gv3::zen_bin];
    TFile *ftime = TFile::Open("./data/zenc.root");

    TH1D  *Time[gv3::dec_bin];
    for (int i=0;i<gv3::dec_bin;i++)
        Time[i] = (TH1D *) ftime->Get(Form("dec_%02d", i));

    TFile *simuall=TFile::Open("/data/home/cwy/Science/DR_KM2A/km2atot.root");
    TH2D *hmap=(TH2D*)simuall->Get("hmap");
    double Nall[40][90];
    for(int i=0;i<40;i++){
        for(int j=0;j<90;j++){
            Nall[i][j]=hmap->GetBinContent(i+1,j+1);
        }
    }
    TF1 *spectrum=new TF1("spectrum","1.0*pow(10.,-11)*pow(x,-3.09)",0.0001,10000.);
    TF1 *spectrum_SCI=new TF1("spectrum_SCI","8.2*pow(10.,-14)*pow(x/10.,-2.90-0.19*log10(x/10.))",0.0001,10000.);
    TF1* ff1=new TF1("ff1","cos(x)",0,papi::twopi);

    float ratio[gv3::dec_bin][gv3::E_bin][gv3::zen_bin];
    for (int idx_dec=0;idx_dec<gv3::dec_bin;idx_dec++){
        std::cout<<idx_dec<<std::endl;
        for(int idx_E=0;idx_E<gv3::E_bin;idx_E++){
            double E_lo = pow(10.,gv3::E_binwid*idx_E+gv3::E_min);
            double E_hi = pow(10.,gv3::E_binwid*(idx_E+1)+gv3::E_min);
            double E_center = pow(10.,gv3::E_binwid*(idx_E+0.5)+gv3::E_min);
            for(int idx_zen=0;idx_zen<gv3::zen_bin;idx_zen++){
                ratio[idx_dec][idx_E][idx_zen]=(spectrum->Integral(E_lo,E_hi))*(ff1->Eval((idx_zen+0.5)*papi::degrad))*gv3::area*(Time[idx_dec]->GetBinContent(idx_zen+1));
            }
        }
    }

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    // std::ofstream outputFile("example3.txt"); 
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if ( jentry%100000==0 ) {std::cout<<"  simu entry = "<<jentry<<" "<<ientry<<std::endl;
        std::cout<<Redge<<std::endl;}
        // outputFile<<rec_Eage<<endl;
        
        int inhbin = EventCut::EventSelection(this);
        if(inhbin<0) continue;
        int idx_E     = int( log10(E/1000)/gv3::E_binwid );
        int idx_mczen = int( theta*papi::raddeg);
        double distance=papi::dir2ang(rec_theta,rec_phi,theta,phi)*papi::raddeg;
            for(int idecbin=0;idecbin<gv3::dec_bin;idecbin++){
                if(Nall[idx_E][idx_mczen]<=0.){
                    // std::cout<<"Warning: Unexpected Nall!"<<std::endl;
                    continue;
                }
                PSF_dec_nh[idecbin][inhbin]->Fill(distance, gv3::timearray3/gv3::timeall*ratio[idecbin][idx_E][idx_mczen]/Nall[idx_E][idx_mczen]);
                EnSig_dec_nh[idecbin][inhbin]->Fill(log10(E/1000), gv3::timearray3/gv3::timeall*ratio[idecbin][idx_E][idx_mczen]/Nall[idx_E][idx_mczen]);
                // EnSig_dec_nh[idecbin][inhbin]->Fill(log10(E/1000), ratio[idecbin][idx_E][idx_mczen]/Nall[idx_E][idx_mczen]);
            }
        }

        for (int i=0;i<gv3::dec_bin;i++){                            
            for (int j=0;j<EventCut::Nhit_bin;j++){
                // if (PSF_dec_nh[i][j]->GetSumOfWeights()>0){
                PSF_dec_nh[i][j]->Scale(1./PSF_dec_nh[i][j]->GetSumOfWeights());
                double max = PSF_dec_nh[i][j]->GetBinContent(PSF_dec_nh[i][j]->GetMaximumBin());
                double sigma = PSF_dec_nh[i][j]->GetRMS();
                PSF_dec_nh_fit[i][j]->SetParameter(0, max);
                PSF_dec_nh_fit[i][j]->SetParameter(1, 0.95);
                PSF_dec_nh_fit[i][j]->SetParameter(2, sigma);
                PSF_dec_nh_fit[i][j]->SetParameter(3, 2.*sigma);

                PSF_dec_nh_fit[i][j]->SetParLimits(0, 0.,max);
                PSF_dec_nh_fit[i][j]->SetParLimits(1, 0.,1.);
                PSF_dec_nh_fit[i][j]->SetParLimits(2, 0.,3);
                PSF_dec_nh_fit[i][j]->SetParLimits(3, 0.,3);

                PSF_dec_nh_fit[i][j]->FixParameter(0,gv3::psf_p0[j]);
                PSF_dec_nh_fit[i][j]->FixParameter(1,gv3::psf_p1[j]);
                PSF_dec_nh_fit[i][j]->FixParameter(2,gv3::psf_p2[j]);
                PSF_dec_nh_fit[i][j]->FixParameter(3,gv3::psf_p3[j]);
                // cout<<gv3::psf_p0[j]<<", "<<gv3::psf_p1[j]<<", "<<gv3::psf_p2[j]<<", "<<gv3::psf_p3[j]<<endl;
                // cout<<PSF_dec_nh_fit[i][j]->GetParameter(0)<< ", "<<PSF_dec_nh_fit[i][j]->GetParameter(1)<< ", "<<PSF_dec_nh_fit[i][j]->GetParameter(2)<<", "<< PSF_dec_nh_fit[i][j]->GetParameter(3)<<endl;
                PSF_dec_nh[i][j]->Fit(PSF_dec_nh_fit[i][j], "Q",0,10.);
                // }
            }
        }

        ftime->Close();
    }