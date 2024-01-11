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
#include <TString.h>
#include <TCut.h>
#include <TParameter.h>
#include "dataCut/EventCut.h"
#include "inc_astro/papi.h"
#include <TH1D.h>
#include <TF1.h>


using namespace std;
struct AnalysisBin {
    TCut cuts_;
    Int_t id;
};

int main(int argc, char *argv[])
{   
    static const int Nhit_bin = 14;
    const int Ebin[Nhit_bin+1]   =  {0.6, 0.8, 1,  1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 10000};
    const float R_Cut[Nhit_bin] =  {-5.11,-5.11,-5.11, -5.24, -5.95, -6.08, -2.34, -2.35, -2.36, -2.36, -2.36, -2.36,-2.36,-2.36};
    
    if(argc < 3)
    {
        cout<<"USAGE: "<<argv[0]<<" inputroot output.root !!!\n"<<endl;
        exit(0);
    }

    string prefix = "root://eos01.ihep.ac.cn";

    double simdec,lowerEdge,upperEdge;
    
    int dec_bin=110;
    double dec_begin=-25.;
    double dec_width=1.;

    TH1D ***PSF_dec_nh     = new TH1D**[dec_bin];
    TH1D ***EnSig_dec_nh   = new TH1D**[dec_bin];
    TF1  ***PSF_dec_nh_fit = new TF1 **[dec_bin];
    for (int i=0;i<dec_bin;i++){
        PSF_dec_nh[i]     = new TH1D*[Nhit_bin];
        EnSig_dec_nh[i]   = new TH1D*[Nhit_bin];
        PSF_dec_nh_fit[i] = new TF1 *[Nhit_bin];
        for (int j=0;j<Nhit_bin;j++){
            PSF_dec_nh[i][j] = new TH1D(Form("PSF_dec%d_nh%d", i, j), Form("PSF_dec%d_nh%d", i, j), 200, 0, 10);
            EnSig_dec_nh[i][j] = new TH1D(Form("EnSig_dec%d_nh%d", i, j), Form("EnSig_dec%d_nh%d", i, j), 40, 0, 4);
            PSF_dec_nh_fit[i][j] = new TF1(Form("PSF_dec%d_nh%d_fit", i, j), "[0]*(x*(([1]*exp(-(x*((x/2)/[2]))))+((1-[1])*exp(-(x*((x/2)/[3]))))))", 0, 10);
        }
    }

    TFile *fin3  = TFile::Open(argv[3]);
    std::cout<<argv[3]<<std::endl;
    if (!fin3){
        cerr<<"Error : Can not open "<<argv[3]<<". Returned!"<<endl;
        return -1;
    }   
    if (fin3->IsZombie()){
        cerr<<"Error : "<<argv[3]<<"is Zombie. Returned!"<<endl;
        return -1;
    }   
    if (fin3->GetEND()<10000){
        cerr<<"Error : "<<argv[3]<<"is small. Returned!"<<endl;
        return -1;
    }   
    TTree *tin3 = (TTree*)fin3->Get("events");
    if (!tin3) {
        cerr<<"Error : "<<argv[3]<<"has no tree named events"<<endl;
        return -1;
    }

    KM2Afliter genresponse3(tin3);
    cout<<"detector3 loop ~~~";
    genresponse3.Loop(PSF_dec_nh, PSF_dec_nh_fit, EnSig_dec_nh, 3);
    // fin3->Close();
    
    // input, Gamma-ray simulation data
    //TFile *fin  = TFile::Open(Form("%s/%s", prefix.data(), argv[1]));
    TFile *fin  = TFile::Open(argv[1]);
    std::cout<<argv[1]<<std::endl;
    if (!fin){
        cerr<<"Error : Can not open "<<argv[1]<<". Returned!"<<endl;
        return -1;
    }   
    if (fin->IsZombie()){
        cerr<<"Error : "<<argv[1]<<"is Zombie. Returned!"<<endl;
        return -1;
    }   
    if (fin->GetEND()<10000){
        cerr<<"Error : "<<argv[1]<<"is small. Returned!"<<endl;
        return -1;
    }   

    TTree *tin = (TTree*)fin->Get("events");
    if (!tin) {
        cerr<<"Error : "<<argv[1]<<"has no tree named events"<<endl;
        return -1;
    }

    KM2Afliter genresponse1(tin);
    cout<<"detector1 loop ~~~";
    genresponse1.Loop(PSF_dec_nh, PSF_dec_nh_fit, EnSig_dec_nh, 2);
    // fin->Close();
    
    TFile *fin2  = TFile::Open(argv[2]);
    std::cout<<argv[2]<<std::endl;
    if (!fin2){
        cerr<<"Error : Can not open "<<argv[2]<<". Returned!"<<endl;
        return -1;
    }   
    if (fin2->IsZombie()){
        cerr<<"Error : "<<argv[2]<<"is Zombie. Returned!"<<endl;
        return -1;
    }   
    if (fin2->GetEND()<10000){
        cerr<<"Error : "<<argv[2]<<"is small. Returned!"<<endl;
        return -1;
    }   
    TTree *tin2 = (TTree*)fin2->Get("events");
    if (!tin2) {
        cerr<<"Error : "<<argv[2]<<"has no tree named events"<<endl;
        return -1;
    }

    KM2Afliter genresponse2(tin2);
    cout<<"detector2 loop ~~~";
    genresponse2.Loop(PSF_dec_nh, PSF_dec_nh_fit, EnSig_dec_nh, 1);
    // fin2->Close();






    TFile *fout = TFile::Open(Form("%s", argv[4]),"RECREATE");
    // output : detector_response.root
    TF1 *fspectrum = new TF1("LogLogSpectrum", "log10([0])-([1]*x)", -3, 6);
    fspectrum->SetParameter(0, 1.e-12);
    fspectrum->SetParameter(1, 2);
    //TF1 *fspectrum = new TF1("LogLogSpectrum", "log10([0])+(-[1]-[2]*(x-1))*(x-1)", -3, 6);
    //fspectrum->SetParameter(0, 8.2e-14);
    //fspectrum->SetParameter(1, 2.90);
    //fspectrum->SetParameter(2, 0.19);

    cout<<"debug outfile!";
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

    cout<<"debug outfile2!";
    AnalysisBin nbin;
    int id;
    TTree *AnalysisBinsTree=new TTree("AnalysisBins","AnalysisBins");
    AnalysisBinsTree->Branch("cuts", "TCut", &nbin.cuts_);
    AnalysisBinsTree->Branch("id", &id,"id/I");

    for(int i=0;i<Nhit_bin;i++){
        nbin.cuts_=Form("(NuW2>=%d) && (NuW2<%d) &&(R<%.2f)",Ebin[i],Ebin[i+1],R_Cut[i]);
        id=Int_t(i);
        AnalysisBinsTree->Fill();
    }

    cout<<"debug outfile3!";
    fout->cd();
    fspectrum->Write();
    DecBins->Write();
    AnalysisBinsTree->Write();

    for(int i=0;i<dec_bin;i++){
        for(int i_E=0;i_E<Nhit_bin;i_E++){
            fout->mkdir(Form("dec_%02d/nh_%02d",i,i_E));
            fout->cd(Form("dec_%02d/nh_%02d",i,i_E));
            // PSF_dec_nh[i][i_E]->Write();
            EnSig_dec_nh[i][i_E]->Write();
            PSF_dec_nh_fit[i][i_E]->Write();
        }
        //data->Write();
    }
    fout->Close();
    fin3->Close();
    fin2->Close();
    fin->Close();
    
    cout<<"response over"<<endl;
    return 0;
}
