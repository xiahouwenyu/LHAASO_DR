#define events_cxx
#include "events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

namespace gv{
    const int dec_bin=110;
    double dec_begin=-25.;
    double dec_width=1.;
    float E_min =  0;
    float E_max =  4.;
    const int  zen_bin =  70;
    float zen_binwidth =  70./zen_bin;
    const int    E_bin = 40;
    float E_binwid=(E_max-E_min)*1.0/E_bin;
    float area=100000.*100000.*3.1415926;
    float livetime[3]={296.086412,216.648228,340.581586};
    const float N_MC[4]={20.*pow(10.,7),20.*pow(10.,6),20.*pow(10.,5),20.*pow(10.,4)};
    const int Erec_bin=14;
    float cut_E[Erec_bin+1]={0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4};
    float psf_p0[Erec_bin]={1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    float psf_p1[Erec_bin]={0.,0.455768 , 0.0500514, 0.0286098, 0.0435536, 0.0133566, 0.0120976,0.016009 , 0.076813 , 0.925462 , 0.0767159, 0.60586  , 0.130434 ,0.947004};
    float psf_p2[Erec_bin]={0.,1.69255e-08, 2.85972e+00, 2.71855e+00, 1.75280e+00, 5.00000e+00,3.23699e+00, 8.20615e-01, 9.05739e-02, 2.10865e-02, 1.80398e-02,3.04091e-04, 1.60081e-02, 1.34311e-05};
    float psf_p3[Erec_bin]={0.,4.75946e-01, 2.45181e-01, 1.91231e-01, 1.25895e-01, 7.70719e-02,5.52710e-02, 4.37323e-02, 2.43867e-02, 2.10837e-02, 3.46485e-08,1.66127e-02, 2.19763e-05, 8.38775e-03};

}

int Cut_E(float Erec){
    for (int i=0;i<gv::Erec_bin;i++){
        if ( Erec < gv::cut_E[i+1] && Erec > gv::cut_E[i] ){
            //std::cout<<Erec <<" "<< i <<std::endl;
            return i;
        }
    }
    return -99;
}



void events::Loop(TH1D ***PSF_dec_nh, TH1D ***EnSig_dec_nh,int phase)
{
    if (fChain == 0) return;

    std::cout<<" *** Read weight from normalizing zenith angle distribution ..."<<std::endl;
    Long64_t nentries = fChain->GetEntriesFast();
    TFile *ftime = TFile::Open("zenc.root");
    TH1D  *Time[gv::dec_bin];
    for (int i=0;i<gv::dec_bin;i++){
        Time[i] = (TH1D *) ftime->Get(Form("dec_%02d", i));
    }

    TF1* F_theta=new TF1("F_theta","sin(x)*cos(x)",0,70*papi::degrad);
    double theta_all=F_theta->Integral(0,70*papi::degrad);

    double Num_mczen[gv::zen_bin];
    for(int idx_zen=0;idx_zen<gv::zen_bin;idx_zen++){
        Num_mczen[idx_zen]= F_theta->Integral(idx_zen*gv::zen_binwidth*papi::degrad,(idx_zen+1)*gv::zen_binwidth*papi::degrad)/theta_all;
        std::cout<< idx_zen<<" "<< Num_mczen[idx_zen]<<std::endl;
    }
    


    TF1 *spectrum=new TF1("spectrum","1.0*pow(10.,-12)*pow(x,-2.)",0.1,100000.);
    TF1 *sp_init=new TF1("spectrum_init","8.84*pow(10.,-7)*pow(x/20.,-2.)",0.1,100000.);
    TF1* ff1=new TF1("ff1","cos(x)",0,papi::twopi);

    TH2D *for_test =new TH2D("theta_E","theta_E",gv::zen_bin,0,70,gv::E_bin,0,4);

    double weight[gv::dec_bin][gv::E_bin][gv::zen_bin];
    for(int idx_dec=0;idx_dec<gv::dec_bin;idx_dec++){
        for(int idx_E=0;idx_E<gv::E_bin;idx_E++){
            double E_lo = pow(10.,gv::E_binwid*idx_E+gv::E_min);
            double E_hi = pow(10.,gv::E_binwid*(idx_E+1.)+gv::E_min);
            double ratio= spectrum->Integral(E_lo,E_hi)/spectrum->Integral(pow(10.,idx_E/10) ,pow(10.,1+idx_E/10));
            //std::cout<< idx_E<<" "<<spectrum->Integral(E_lo,E_hi)/spectrum->Integral(pow(10.,idx_E/10) ,pow(10.,1+idx_E/10))<<std::endl;
            for(int idx_zen=0;idx_zen<gv::zen_bin;idx_zen++){
                //weight[idx_dec][idx_E][idx_zen]=(spectrum->Integral(E_lo,E_hi))*(ff1->Eval((idx_zen+0.5)*gv::zen_binwidth*papi::degrad))*gv::area*(Time[idx_dec]->GetBinContent(idx_zen+1))*gv::livetime[phase] /  ( gv::N_MC[idx_E/10] * ratio * Num_mczen[idx_zen])/(853.31623) ;
                weight[idx_dec][idx_E][idx_zen]=(spectrum->Integral(E_lo,E_hi))*(ff1->Eval((idx_zen+0.5)*gv::zen_binwidth*papi::degrad))*gv::area*(Time[idx_dec]->GetBinContent(idx_zen+1))*gv::livetime[phase] /  ( gv::N_MC[idx_E/10] * ratio * Num_mczen[idx_zen])/(gv::livetime[0]+gv::livetime[1]+gv::livetime[2]) ;
                //weight[idx_dec][idx_E][idx_zen]=(spectrum->Integral(E_lo,E_hi))*(F_theta->Integral(0*gv::zen_binwidth*papi::degrad,70*gv::zen_binwidth*papi::degrad))*(Time[idx_dec]->GetBinContent(idx_zen+1)) /  ( sp_init->Integral(E_lo,E_hi)*sin((idx_zen+0.5)*papi::degrad));

                if(idx_dec==57){
                    for_test->SetBinContent(idx_zen,idx_E,(spectrum->Integral(E_lo,E_hi))*(ff1->Eval((idx_zen+0.5)*gv::zen_binwidth*papi::degrad))*gv::area*(Time[idx_dec]->GetBinContent(idx_zen+1))*gv::livetime[phase]);
                }
            }
        }
    }
    TFile *file_test = TFile::Open(Form("tst%d.root",phase),"RECREATE");
    file_test->cd();
    for_test->Write();
    file_test->Close();
    



    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if(jentry%100000==0) std::cout<<jentry<<std::endl;
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if(rec_theta*papi::raddeg > 50) continue;
        int idx_E     = int( (prim_E-gv::E_min)/gv::E_binwid );
        double E_Low = pow(10.,idx_E*gv::E_binwid);
        double E_Hig = pow(10.,(idx_E+1)*gv::E_binwid);
        int idx_mczen = int(prim_theta/(gv::zen_binwidth)*papi::raddeg);
        if(idx_mczen >=gv::zen_bin || idx_mczen<0 ) continue;
        int inhbin = Cut_E(rec_E);
        if(inhbin<0) continue;

        for (int idecbin=0;idecbin<gv::dec_bin;idecbin++){
            EnSig_dec_nh[idecbin][inhbin] ->Fill( prim_E , weight[idecbin][idx_E][idx_mczen]);
            PSF_dec_nh[idecbin][inhbin]->Fill(psf_r, weight[idecbin][idx_E][idx_mczen]);
            //double weight_new= spectrum->Integral( E_Low,E_Hig)*Time[idecbin]->GetBinContent(idx_mczen+1)*57.3*0.5*sin(70*papi::degrad)*sin(70*papi::degrad) /(sp_init->Integral(E_Low,E_Hig ))/sin(70*papi::degrad);
            //    weight[idx_dec][idx_E][idx_zen]=(spectrum->Integral(E_lo,E_hi))*(F_theta->Integral(0*gv::zen_binwidth*papi::degrad,70*gv::zen_binwidth*papi::degrad))*(Time[idx_dec]->GetBinContent(idx_zen+1)) /  ( sp_init->Integral(E_lo,E_hi)*sin((idx_zen+0.5)*papi::degrad));
            //std::cout<<prim_E <<" "<< E_Low<<" "<< E_Hig <<" "<< inhbin <<" " <<" "<<weight_new<<" "<< weight[idecbin][idx_E][idx_mczen]<<""<<std::endl;
            //EnSig_dec_nh[idecbin][inhbin] ->Fill( prim_E , weight_new);
            //PSF_dec_nh[idecbin][inhbin]->Fill(psf_r, weight_new);
        }
    }
}


void events::Fit(TH1D ***PSF_dec_nh,TF1  ***PSF_dec_nh_fit)
{

    for (int i=0;i<gv::dec_bin;i++)
    {
        for (int j=0;j<gv::Erec_bin;j++)
        {
            if (PSF_dec_nh[i][j]->GetSumOfWeights()>0)
            {
                PSF_dec_nh[i][j]->Scale(1./PSF_dec_nh[i][j]->GetSumOfWeights());
                double max = PSF_dec_nh[i][j]->GetBinContent(PSF_dec_nh[i][j]->GetMaximumBin());
                double sigma = PSF_dec_nh[i][j]->GetRMS();
                //PSF_dec_nh_fit[i][j]->SetParameter(0, max);
                //PSF_dec_nh_fit[i][j]->SetParameter(1, 0.30);
                //PSF_dec_nh_fit[i][j]->SetParameter(2, 0.30);
                //PSF_dec_nh_fit[i][j]->SetParameter(3, 1*sigma);
                //PSF_dec_nh_fit[i][j]->SetParameter(4, 3.*sigma);
                //PSF_dec_nh_fit[i][j]->SetParameter(5, 5.*sigma);
                //PSF_dec_nh_fit[i][j]->SetParameter(6, 0.30);

                PSF_dec_nh_fit[i][j]->SetParameter(0, max);
                PSF_dec_nh_fit[i][j]->SetParameter(1, 0.95);
                PSF_dec_nh_fit[i][j]->SetParameter(2, sigma);
                PSF_dec_nh_fit[i][j]->SetParameter(3, 2.*sigma);

                //PSF_dec_nh_fit[i][j]->FixParameter(0,gv::psf_p0[j]);
                //PSF_dec_nh_fit[i][j]->FixParameter(1,gv::psf_p1[j]);
                //PSF_dec_nh_fit[i][j]->FixParameter(2,gv::psf_p2[j]);
                //PSF_dec_nh_fit[i][j]->FixParameter(3,gv::psf_p3[j]);

                PSF_dec_nh[i][j]->Fit(PSF_dec_nh_fit[i][j], "Q",0,5.);
            }
        }
    }

}
