/***********************************
 *         GM,Xiang  &SC,Hu         *
 *            @2021.2.1             *
 *               IHEP               *
 *             Response             *
 ***********************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TBranch.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TString.h>
#include <TParameter.h>
#include <TCut.h>
#include <TF1.h>
#include <TString.h>
#include "inc/papi.h"
#include "events.h"
//#include "Rec.h"


using namespace std;
struct AnalysisBin{
        TCut cuts_;
            Int_t id;
};

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cout<<"USAGE: "<<argv[0]<<" output.root !!!\n"<<endl;
        exit(0);
    }

    string prefix = "root://eos01.ihep.ac.cn";

    //define the efficiency array
    //static float ***Nsrc;
    int NUM_E=14;
    float data_count,bkg_count;
    const long nside=pow(2,10);
    const long npix=12*nside*nside;
    //Nsrc = new float**[npix];
    double ha,dec,deltat,time_begin;
    Double_t mjd;
    long livetime=0;
    double mjd_begin,mjd_end;
    double simdec,lowerEdge,upperEdge;
    
    int dec_bin=110;
    double dec_begin=-25.;
    double dec_width=1.;
    float cut_nhit[NUM_E+1]={0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4};

    TH1D ***PSF_dec_nh     = new TH1D**[dec_bin];
    TH1D ***EnSig_dec_nh   = new TH1D**[dec_bin];
    TF1  ***PSF_dec_nh_fit = new TF1 **[dec_bin];
    for (int i=0;i<dec_bin;i++){
        PSF_dec_nh[i]     = new TH1D*[NUM_E];
        EnSig_dec_nh[i]   = new TH1D*[NUM_E];
        PSF_dec_nh_fit[i] = new TF1 *[NUM_E];
        for (int j=0;j<NUM_E;j++){
            PSF_dec_nh[i][j] = new TH1D(Form("PSF_dec%d_nh%d", i, j), Form("PSF_dec%d_nh%d", i, j), 100, 0, 5);
            EnSig_dec_nh[i][j] = new TH1D(Form("EnSig_dec%d_nh%d", i, j), Form("EnSig_dec%d_nh%d", i, j), 40, 0, 4);
            //PSF_dec_nh_fit[i][j] = new TF1(Form("PSF_dec%d_nh%d_fit", i, j), "[0]*x*(([1]*exp(-(x*((x/2)/[3]))))+([2]*exp(-(x*((x/2)/[4]))))+([6]*exp(-(x*((x/2)/[5])))))", 0, 10);
            PSF_dec_nh_fit[i][j] = new TF1(Form("PSF_dec%d_nh%d_fit", i, j), "[0]*(x*(([1]*exp(-(x*((x/2)/[2]))))+((1-[1])*exp(-(x*((x/2)/[3]))))))", 0, 5);
            //PSF_dec_nh_fit[i][j] = new TF1(Form("PSF_dec%d_nh%d_fit", i, j), "[0]*x*(exp(-(x*((x/2)/[2]))))", 0, 5);
        }
    }

    // input, Gamma-ray simulation data
    //TFile *fin  = TFile::Open("GammaMC-12array.root");
    //TTree *tin = (TTree*)fin->Get("Rec");
    TFile *fin  = TFile::Open("MC12.root");
    TTree *tin = (TTree*)fin->Get("events");
    //Rec genresponse1(tin);
    events genresponse1(tin);
    genresponse1.Loop(PSF_dec_nh, EnSig_dec_nh,0);

    fin  = TFile::Open("MC34.root");
    tin = (TTree*)fin->Get("events");
    events genresponse2(tin);
    genresponse2.Loop(PSF_dec_nh, EnSig_dec_nh,1);

    fin  = TFile::Open("MCfull.root");
    tin = (TTree*)fin->Get("events");
    events genresponse3(tin);
    genresponse3.Loop(PSF_dec_nh, EnSig_dec_nh,2);

    genresponse3.Fit(PSF_dec_nh,PSF_dec_nh_fit);

    // output : detector_response.root
    TF1 *fspectrum = new TF1("LogLogSpectrum", "log10([0])-([1]*x)", 0, 4);
    fspectrum->SetParameter(0, 1.e-12);
    fspectrum->SetParameter(1, 2.);

    TFile *fout = TFile::Open(Form("%s/%s", prefix.data(), argv[1]),"RECREATE");
    TTree *DecBins=new TTree("DecBins","DecBins");
    DecBins->Branch("simdec",&simdec,"simdec/D");
    DecBins->Branch("lowerEdge",&lowerEdge,"lowerEdge/D");
    DecBins->Branch("upperEdge",&upperEdge,"upperEdge/D");
    for(int i=0;i<dec_bin;i++){
        simdec=dec_begin+i*dec_width;
        lowerEdge=dec_begin+(i-0.5)*dec_width;
        upperEdge=dec_begin+(i+0.5)*dec_width;
        DecBins->Fill();
    }

    AnalysisBin nbin;
    int id;
    TTree *AnalysisBinsTree=new TTree("AnalysisBins","AnalysisBins");
    AnalysisBinsTree->Branch("cuts", "TCut", &nbin.cuts_);
    AnalysisBinsTree->Branch("id", &id,"id/I");

    for(int i=0;i<NUM_E;i++){
        nbin.cuts_=Form("(log10(Erec)>=%.1f) && (log10(Erec)<%.1f)",cut_nhit[i],cut_nhit[i+1]);
        id=Int_t(i);
        AnalysisBinsTree->Fill();
    }


    fout->cd();
    fspectrum->Write();
    DecBins->Write();
    AnalysisBinsTree->Write();
    for(int i=0;i<dec_bin;i++){
        for(int i_E=0;i_E<NUM_E;i_E++){
            fout->mkdir(Form("dec_%02d/nh_%02d",i,i_E));
            fout->cd(Form("dec_%02d/nh_%02d",i,i_E));
            PSF_dec_nh[i][i_E]->Write();
            EnSig_dec_nh[i][i_E]->Write();
            PSF_dec_nh_fit[i][i_E]->Write();
        }
        //data->Write();
    }
    fout->Close();
    fin->Close();

    /*for (int i=0;i<dec_bin;i++)
        for (int j=0;j<NUM_E;j++){
            delete []PSF_dec_nh[i][j];
            delete []EnSig_dec_nh[i][j];
            delete []PSF_dec_nh_fit[i][j];
        }
    for (int i=0;i<dec_bin;i++){
        delete []PSF_dec_nh[i];
        delete []EnSig_dec_nh[i];
        delete []PSF_dec_nh_fit[i];
    }
    delete []PSF_dec_nh;
    delete []EnSig_dec_nh;
    delete []PSF_dec_nh_fit;*/

    cout<<"response over"<<endl;
    return 0;
}
