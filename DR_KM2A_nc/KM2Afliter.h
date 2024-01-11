//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 10 18:44:24 2024 by ROOT version 6.24/06
// from TTree events/mc event after cut
// found on file: MCfull.root
//////////////////////////////////////////////////////////

#ifndef KM2Afliter_h
#define KM2Afliter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

#include <stdio.h>
#include <stdlib.h>

// Header file for the classes stored in the TTree if any.

class KM2Afliter {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         prim_theta;
   Float_t         rec_theta;
   Float_t         prim_E;
   Float_t         rec_E;
   Float_t         psf_r;
   Double_t        w_sed;

   // List of branches
   TBranch        *b_prim_theta;   //!
   TBranch        *b_rec_theta;   //!
   TBranch        *b_prim_E;   //!
   TBranch        *b_rec_E;   //!
   TBranch        *b_psf_r;   //!
   TBranch        *b_w_sed;   //!

   KM2Afliter(TTree *tree=0);
   virtual ~KM2Afliter();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TH1D ***PSF_dec_nh, TF1  ***PSF_dec_nh_fit, TH1D ***EnSig_dec_nh, int detflag);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef KM2Afliter_cxx
KM2Afliter::KM2Afliter(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MCfull.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MCfull.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
}

KM2Afliter::~KM2Afliter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t KM2Afliter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t KM2Afliter::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void KM2Afliter::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("prim_theta", &prim_theta, &b_prim_theta);
   fChain->SetBranchAddress("rec_theta", &rec_theta, &b_rec_theta);
   fChain->SetBranchAddress("prim_E", &prim_E, &b_prim_E);
   fChain->SetBranchAddress("rec_E", &rec_E, &b_rec_E);
   fChain->SetBranchAddress("psf_r", &psf_r, &b_psf_r);
   fChain->SetBranchAddress("w_sed", &w_sed, &b_w_sed);
   Notify();
}

Bool_t KM2Afliter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void KM2Afliter::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t KM2Afliter::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef KM2Afliter_cxx
