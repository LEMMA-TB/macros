// *****************************************************
//
// author: Alessandra Cappati
//         03/05/2019
// 
// usage: specify the input files (Data and MC), 
//        the output directory, and other options 
//        at the end of the script
//
// run with:
//        root -l -b -q emittanceUnc_bootstrap.C++
//
// *****************************************************

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TAxis.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TRandom3.h"
#include "TPaveText.h"

#include "FillBeamInfo.h"


using namespace std;



Double_t getemittance(vector<Double_t> xv, vector<Double_t> xpv){

  UInt_t centering=1; // default: 1, only other reasonable value: 0
	
  UInt_t n_points = xv.size();

  Double_t emittance=0.;

  Double_t x=0.;
  Double_t x2=0.;
  Double_t xp=0.;
  Double_t xp2=0.;
  Double_t xxp=0.;
  for(UInt_t i=0;i<xv.size();i++){
    x+=(1000000*xv[i]);                // [nm]
    x2+=(1000000*xv[i]*1000000*xv[i]); // [nm x nm]
    xp+=(xpv[i]);                      // [rad]
    xp2+=(xpv[i]*xpv[i]);              // [rad x rad]
    xxp+=(1000000*xv[i]*xpv[i]);       // [nm x rad]
  }
  x*=1./float(n_points);
  x2*=1./float(n_points);
  xp*=1./float(n_points);
  xp2*=1./float(n_points);
  xxp*=1./float(n_points);

  if( centering==0 ){
    emittance=TMath::Sqrt(x2*xp2-xxp*xxp); // [nm x rad]
  }else{
    emittance=((x2-x*x)*(xp2-xp*xp)-(xxp-x*xp)*(xxp-x*xp));
    emittance=TMath::Sqrt(emittance);	  
  }

  return emittance;

}




void emittanceUnc_bootstrap(){


  // ----------------------------------
  TString inputFileName = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/reco-mupmum-Be6cm-GausGaus.root";
  TString label = "MC";

  //double zEndTarget = 10.*(457.9+3.-84.6);   // [mm] - dataset: AUGUST 2018    Be target 6 cm
  double zEndTarget = 10.*(460.93+3.-82.78); // [mm] - dataset: SEPTEMBER 2018 Be target 6 cm and C target 6cm
  //double zEndTarget = 10.*(460.93+1.-82.78); // [mm] - dataset: SEPTEMBER 2018 C  target 2 cm 

  // ----------------------------------




  bool isMC = false;                // for Data
  if(label == "MC"){ isMC = true; } // for MC


  Double_t chi2m;
  Double_t x_pos_mum[12];
  Double_t x_pos_mum_err[12];
  Double_t z_x_pos_mum[12];
  Double_t x_pos_DT_mum[8];
  Double_t z_pos_DT_mum[8];
  Double_t p_mum;
  Double_t p_mup;
  Double_t chi2p;
  Double_t x_pos_mup[12];
  Double_t x_pos_mup_err[12];
  Double_t z_x_pos_mup[12];
  Double_t x_pos_DT_mup[8];
  Double_t z_pos_DT_mup[8];
  Int_t    subdet[100];
  Int_t    itrack[100];
  Double_t xh[100];
  Double_t yh[100];
  Double_t zh[100];
  Int_t    nhits;
  Double_t Calo_EnDep[25];
  Int_t    event_type;
  Double_t gen_pos_mum[12]; // used for MC only
  Double_t gen_pos_mup[12]; // used for MC only
  Double_t gen_vtx_mum[7];  // used for MC only
  Double_t gen_vtx_mup[7];  // used for MC only

  TFile* inputFile = new TFile(inputFileName);
  TTree* inputTree = (TTree*)inputFile->Get("lemma");

  inputTree->SetBranchAddress("chi2m",	        &chi2m);	     
  inputTree->SetBranchAddress("x_pos_mum",      &x_pos_mum[0]); 
  inputTree->SetBranchAddress("x_pos_mum_err",  &x_pos_mum_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mum",    &z_x_pos_mum[0]);
  inputTree->SetBranchAddress("x_pos_DT_mum",   &x_pos_DT_mum[0]);
  inputTree->SetBranchAddress("z_pos_DT_mum",   &z_pos_DT_mum[0]);
  inputTree->SetBranchAddress("p_mum",          &p_mum);	     
  inputTree->SetBranchAddress("p_mup",          &p_mup);	     
  inputTree->SetBranchAddress("chi2p",          &chi2p);	       
  inputTree->SetBranchAddress("x_pos_mup",      &x_pos_mup[0]); 
  inputTree->SetBranchAddress("x_pos_mup_err",  &x_pos_mup_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mup",    &z_x_pos_mup[0]);
  inputTree->SetBranchAddress("x_pos_DT_mup",   &x_pos_DT_mup[0]);
  inputTree->SetBranchAddress("z_pos_DT_mup",   &z_pos_DT_mup[0]);
  inputTree->SetBranchAddress("subdet",         &subdet[0]);   
  inputTree->SetBranchAddress("itrack",         &itrack[0]);   
  inputTree->SetBranchAddress("xh",             &xh[0]);	     
  inputTree->SetBranchAddress("yh",             &yh[0]);	     
  inputTree->SetBranchAddress("zh",             &zh[0]);	     
  inputTree->SetBranchAddress("nhits",          &nhits);	     
  inputTree->SetBranchAddress("Calo_EnDep",     &Calo_EnDep[0]);
  inputTree->SetBranchAddress("event_type",     &event_type);   
  if(isMC){
    inputTree->SetBranchAddress("gen_pos_mum", &gen_pos_mum[0]);
    inputTree->SetBranchAddress("gen_pos_mup", &gen_pos_mup[0]); 
    inputTree->SetBranchAddress("gen_vtx_mum", &gen_vtx_mum[0]);
    inputTree->SetBranchAddress("gen_vtx_mup", &gen_vtx_mup[0]); 
  }    


  // define Histos 
  TH1F* hist_emittanceValue_mup = new TH1F("hist_emittanceValue_mup","hist_emittanceValue_mup",100,5.,20.);
  TH1F* hist_emittanceValue_mum = new TH1F("hist_emittanceValue_mum","hist_emittanceValue_mum",100,5.,20.);


  // ---------------------------
  // --- get tree entries
  Long64_t entries = inputTree->GetEntries();
  cout<<" Reading input file ..."<<endl;



  // repeat 1000 times the procedure
  for(Int_t j=0; j<1000; j++){

    if(j%10 == 0) { cout<<" Running iteration n "<<j<<endl; }


    // def vectors for emittance estimate
    vector<double> vec_emittance_x_mup_MC; 
    vec_emittance_x_mup_MC.clear();
    vector<double> vec_emittance_xprime_mup_MC;
    vec_emittance_xprime_mup_MC.clear();

    vector<double> vec_emittance_x_mum_MC ;
    vec_emittance_x_mum_MC.clear();
    vector<double> vec_emittance_xprime_mum_MC;
    vec_emittance_xprime_mum_MC.clear();


  
    // --- define Trandom variable
    TRandom3* qwerty = new TRandom3(0); 

    // --- for cycle from 0 to 1000 forming the subsample
    for(Int_t z=0; z<1000; z++){

      inputTree->GetEntry(entries*qwerty->Rndm());


      if( p_mup > 0. && p_mum > 0. ) {


        // --- ONLY FOR MC 
        // -----------------------------
        // --- emittance 
        // -----------------------------
        // x  = x  @det30 - x  incoming e+
        // x' = x' @det30 - x' incoming e+
        if(isMC){
          // estimate is done on a reference plane 
          Double_t z_ref  = zEndTarget; // [mm] z ref (end of the target)
        

          // --- e+ incoming
          // px of e+ = Cx mu- * En of mu- + Cx mu+ * En of mu+ = px of mu- + px of mu+
          Double_t px_eplus = gen_vtx_mum[3]*gen_vtx_mum[6] + gen_vtx_mup[3]*gen_vtx_mup[6];
          Double_t py_eplus = gen_vtx_mum[4]*gen_vtx_mum[6] + gen_vtx_mup[4]*gen_vtx_mup[6];
          Double_t pz_eplus = gen_vtx_mum[5]*gen_vtx_mum[6] + gen_vtx_mup[5]*gen_vtx_mup[6];
          // x  of e+ (extrapolation on reference plane)  = x_vtx - (z_vtx - z_ref)*(px_e+/pz_e+)
          // N.B. gen_vtx_mup[0] = gen_vtx_mum[0] and gen_vtx_mup[2] = gen_vtx_mum[2] 
          Double_t x_atZref_eplus = gen_vtx_mup[0] - (gen_vtx_mup[2] - z_ref)*(px_eplus/pz_eplus);
          // x' of e+ (extrapolation on reference plane)  = px / p  (the direction remain the same as at vtx)
          Double_t x_prime_atZref_eplus = px_eplus / sqrt(px_eplus*px_eplus + py_eplus*py_eplus + pz_eplus*pz_eplus);

       

          // --- mu+
          // x  (extrapolation on reference plane): x_det30_atZref_mup = x_onDet30 - (z_onDet30 - z_ref)* px_onDet30 / pz_onDet30
          Double_t x_det30_atZref_mup = gen_pos_mup[0] - (gen_pos_mup[2] - z_ref)*(gen_pos_mup[3]/gen_pos_mup[5]); 
          // x' (extrapolation on reference plane): the direction remain the same as on det30 or at vtx
          Double_t pTot_genLev_mup = sqrt(gen_pos_mup[3]*gen_pos_mup[3] + gen_pos_mup[4]*gen_pos_mup[4] + gen_pos_mup[5]*gen_pos_mup[5]);
          Double_t x_prime_ondet30_mup = gen_pos_mup[3] / pTot_genLev_mup; 
         
          // emittance of mu+
          Double_t x_emittance_mup       = x_det30_atZref_mup  - x_atZref_eplus;
          Double_t x_prime_emittance_mup = x_prime_ondet30_mup - x_prime_atZref_eplus;
          // --- fill vectors
          vec_emittance_x_mup_MC     .push_back(x_emittance_mup);
          vec_emittance_xprime_mup_MC.push_back(x_prime_emittance_mup);



          // --- mu-
          // x  (extrapolation on reference plane): x_det30_atZref_mum = x_onDet30 - (z_onDet30 - z_ref)* px_onDet30 / pz_onDet30
          Double_t x_det30_atZref_mum = gen_pos_mum[0] - (gen_pos_mum[2] - z_ref)*(gen_pos_mum[3]/gen_pos_mum[5]); 
          // x' (extrapolation on reference plane): the direction remain the same as on det30 or at vtx
          Double_t pTot_genLev_mum = sqrt(gen_pos_mum[3]*gen_pos_mum[3] + gen_pos_mum[4]*gen_pos_mum[4] + gen_pos_mum[5]*gen_pos_mum[5]);
          Double_t x_prime_ondet30_mum = gen_pos_mum[3] / pTot_genLev_mum; 

          // emittance of mu-
          Double_t x_emittance_mum       = x_det30_atZref_mum  - x_atZref_eplus;
          Double_t x_prime_emittance_mum = x_prime_ondet30_mum - x_prime_atZref_eplus;
          // --- fill vectors
          vec_emittance_x_mum_MC     .push_back(x_emittance_mum);
          vec_emittance_xprime_mum_MC.push_back(x_prime_emittance_mum);


        }// end if isMC        
                      

      } // end if (p_mup > 0. && p_mum > 0.)

    }//end for cycle from 0 to 1000 forming the subsample

  
  
    // ---------------------------------
    //  compute emittance from vectors
    // ---------------------------------

    hist_emittanceValue_mup->Fill(getemittance(vec_emittance_x_mup_MC, vec_emittance_xprime_mup_MC)); 
    hist_emittanceValue_mum->Fill(getemittance(vec_emittance_x_mum_MC, vec_emittance_xprime_mum_MC)); 


  }//end for cycle from 0 to 1000 repeating the procedure 1000 times


  
  TCanvas* c_emittanceValue_mup = new TCanvas("c_emittanceValue_mup","c_emittanceValue_mup");
  c_emittanceValue_mup->cd();
  hist_emittanceValue_mup->Draw();
  c_emittanceValue_mup->SaveAs("emittanceUnc_mup.png");


  TCanvas* c_emittanceValue_mum = new TCanvas("c_emittanceValue_mum","c_emittanceValue_mum");
  c_emittanceValue_mum->cd();
  hist_emittanceValue_mum->Draw();
  c_emittanceValue_mum->SaveAs("emittanceUnc_mum.png");



}
