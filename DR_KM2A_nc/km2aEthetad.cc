#include <TTree.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>

void km2aEthetad(){
    int totol = 11110000 * 20;
    int area = 3.1415926*100000*100000;

    TFile *fileout = new TFile("./km2atot.root","RECREATE");
    TF2 *toudian = new TF2("toudian", "20*8.84e-7*pow(x/20,-2)*sin(pi*y/180)*cos(pi*y/180)*3.1415926*100000*100000/(57.3*0.5*pow(sin(pi*70/180),2))");
    TH2D *hmap = new TH2D("hmap", "hmap", 40, 0, 4, 90, 0, 90);
    for (int i=1; i<=40; i++){
        for (int j=1; j<=90; j++){
            if (j<=70){
                hmap->SetBinContent(i,j,toudian->Integral(pow(10,0.1*(i-1)),pow(10, 0.1*i), j-1,j));
            }else{
                hmap->SetBinContent(i,j,0);
            }
        }
    }
    // hmap->Draw("colz");
    // cout<<hmap->GetSum()<<endl;
    fileout->cd();
    hmap->Write();
    fileout->Close();


    // TFile *fileout = new TFile("./km2atot_xgm.root","RECREATE");
    // float N_MC[4]={20.*1e7, 20.*1e6, 20.*1e5, 20.*1e4};
    // TF1* F_theta=new TF1("F_theta","sin(x)*cos(x)",0,70*3.1415926/180);
    // TF1 *spectrum=new TF1("spectrum","1.0*pow(10.,-12)*pow(x,-2.)",0.1,100000.);
    
    // double theta_all=F_theta->Integral(0,70*3.1415926/180);
    // double Num_mczen[70];
    // for(int idx_zen=0;idx_zen<70;idx_zen++){
    //     Num_mczen[idx_zen]= F_theta->Integral(idx_zen*1*3.1415926/180,(idx_zen+1)*1*3.1415926/180)/theta_all;
    //     std::cout<< idx_zen<<" "<< Num_mczen[idx_zen]<<std::endl;
    // }
    

    // TH2D *hmap = new TH2D("hmap", "hmap", 40, 0, 4, 90, 0, 90);
    // for (int i=0; i<40; i++){
    //     for (int j=0; j<90; j++){
    //         if (j<70){
    //             double ratio= spectrum->Integral(pow(10,0.1*(i-1)),pow(10, 0.1*i))/spectrum->Integral(pow(10.,i/10) ,pow(10.,1+i/10));
    //             hmap->SetBinContent(i,j,N_MC[i/10] * ratio * Num_mczen[j]);
    //         }else{
    //             hmap->SetBinContent(i,j,0);
    //         }
    //     }
    // }

    // std::cout<<hmap->GetSum()<<std::endl;
    // fileout->cd();
    // hmap->Write();
    // fileout->Close();
}