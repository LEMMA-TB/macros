// *****************************************************
//
// author: Alessandra Cappati
//         04/04/2019
// 
// usage: specify the input files (Data and MC), 
//        the output directory, and other options 
//        at the end of the script
//
// run with:
//        root -l -b -q plotEmittance_graph.C++
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



void fillVecGetEmittance(TString inputFileName, TString label, double zEndTarget, double &emittanceValue_mup, double &emittanceValue_mum){

  bool isMC = false;                // for Data
  if(label == "MC"){ isMC = true; } // for MC


  Double_t chi2Si5MuM;
  Double_t x_pos_mum[12];
  Double_t x_pos_mum_err[12];
  Double_t z_x_pos_mum[12];
  Double_t x_pos_DT_mum[8];
  Double_t z_pos_DT_mum[8];
  Double_t p_mum;
  Double_t p_mup;
  Double_t chi2Si5MuP;
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

  inputTree->SetBranchAddress("chi2Si5MuM",     &chi2Si5MuM);	     
  inputTree->SetBranchAddress("x_pos_mum",      &x_pos_mum[0]); 
  inputTree->SetBranchAddress("x_pos_mum_err",  &x_pos_mum_err[0]);
  inputTree->SetBranchAddress("z_x_pos_mum",    &z_x_pos_mum[0]);
  inputTree->SetBranchAddress("x_pos_DT_mum",   &x_pos_DT_mum[0]);
  inputTree->SetBranchAddress("z_pos_DT_mum",   &z_pos_DT_mum[0]);
  inputTree->SetBranchAddress("p_mum",          &p_mum);	     
  inputTree->SetBranchAddress("p_mup",          &p_mup);	     
  inputTree->SetBranchAddress("chi2Si5MuP",     &chi2Si5MuP);	       
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


  // def vectors for emittance estimate

  
  // MC vectors
  vector<double> vec_emittance_x_mup_MC;
  vector<double> vec_emittance_xprime_mup_MC;

  vector<double> vec_emittance_x_mum_MC;
  vector<double> vec_emittance_xprime_mum_MC;


  
  // ---------------------------
  // --- loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);


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
        Double_t x_emittance_mup = x_det30_atZref_mup - x_atZref_eplus;
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
        Double_t x_emittance_mum = x_det30_atZref_mum - x_atZref_eplus;
        Double_t x_prime_emittance_mum = x_prime_ondet30_mum - x_prime_atZref_eplus;
        // --- fill vectors
        vec_emittance_x_mum_MC     .push_back(x_emittance_mum);
        vec_emittance_xprime_mum_MC.push_back(x_prime_emittance_mum);


      }// end if isMC        
                      

    } // end if (p_mup > 0. && p_mum > 0.)

  }//end over tree entries 

  
  
  // ---------------------------------
  //  compute emittance from vectors
  // ---------------------------------

  emittanceValue_mup = getemittance(vec_emittance_x_mup_MC, vec_emittance_xprime_mup_MC); 
  emittanceValue_mum = getemittance(vec_emittance_x_mum_MC, vec_emittance_xprime_mum_MC); 


}



// main function 
void plotEmittance_graph(){


  TString inputPath = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/mc_test/";
  double zEndTarget = 10.*(457.9+3.-84.6);   // [mm] - dataset: AUGUST 2018    Be target 6 cm


  // --- file 1
  double emittValue_mup_1 = -99.;
  double emittValue_mum_1 = -99.;
  fillVecGetEmittance(inputPath + "reco-mupmum-2k-1.root" , "MC", zEndTarget, emittValue_mup_1, emittValue_mum_1);
  cout<<"--- File 1: "<<endl;
  cout<<"emittance mup: "<<emittValue_mup_1<<" nm x rad"<<endl;
  cout<<"emittance mum: "<<emittValue_mum_1<<" nm x rad"<<endl;

  // --- file 2
  double emittValue_mup_2 = -99.;
  double emittValue_mum_2 = -99.;
  fillVecGetEmittance(inputPath + "reco-mupmum-2k-2.root" , "MC", zEndTarget, emittValue_mup_2, emittValue_mum_2);
  cout<<"--- File 2: "<<endl;
  cout<<"emittance mup: "<<emittValue_mup_2<<" nm x rad"<<endl;
  cout<<"emittance mum: "<<emittValue_mum_2<<" nm x rad"<<endl;

  // --- file 3
  double emittValue_mup_3 = -99.;
  double emittValue_mum_3 = -99.;
  fillVecGetEmittance(inputPath + "reco-mupmum-2k-3.root" , "MC", zEndTarget, emittValue_mup_3, emittValue_mum_3);
  cout<<"--- File 3: "<<endl;
  cout<<"emittance mup: "<<emittValue_mup_3<<" nm x rad"<<endl;
  cout<<"emittance mum: "<<emittValue_mum_3<<" nm x rad"<<endl;

  // --- file 4
  double emittValue_mup_4 = -99.;
  double emittValue_mum_4 = -99.;
  fillVecGetEmittance(inputPath + "reco-mupmum-2k-4.root" , "MC", zEndTarget, emittValue_mup_4, emittValue_mum_4);
  cout<<"--- File 4: "<<endl;
  cout<<"emittance mup: "<<emittValue_mup_4<<" nm x rad"<<endl;
  cout<<"emittance mum: "<<emittValue_mum_4<<" nm x rad"<<endl;

  // --- file 5
  double emittValue_mup_5 = -99.;
  double emittValue_mum_5 = -99.;
  fillVecGetEmittance(inputPath + "reco-mupmum-2k-5.root" , "MC", zEndTarget, emittValue_mup_5, emittValue_mum_5);
  cout<<"--- File 5: "<<endl;
  cout<<"emittance mup: "<<emittValue_mup_5<<" nm x rad"<<endl;
  cout<<"emittance mum: "<<emittValue_mum_5<<" nm x rad"<<endl;

  // --- file 12345 (it contains all the previous ones)
  double emittValue_mup_12345 = -99.;
  double emittValue_mum_12345 = -99.;
  fillVecGetEmittance(inputPath + "reco-mupmum-2k-12345.root" , "MC", zEndTarget, emittValue_mup_12345, emittValue_mum_12345);
  cout<<"--- File 12345: "<<endl;
  cout<<"emittance mup: "<<emittValue_mup_12345<<" nm x rad"<<endl;
  cout<<"emittance mum: "<<emittValue_mum_12345<<" nm x rad"<<endl;


  // -----------
  // --- Graph
  // -----------
  double vec_file[5] = {1.,2.,3.,4.,5.};
  double vec_mup[5]  = {emittValue_mup_1, emittValue_mup_2, emittValue_mup_3, emittValue_mup_4, emittValue_mup_5}; 
  double vec_mum[5]  = {emittValue_mum_1, emittValue_mum_2, emittValue_mum_3, emittValue_mum_4, emittValue_mum_5}; 


  TGraph* g_mup = new TGraph(5, vec_file, vec_mup);
  TGraph* g_mum = new TGraph(5, vec_file, vec_mum);


  TCanvas* c_emittance = new TCanvas("c_emittance","c_emittance",800,500);
  c_emittance->Divide(2,1);
  c_emittance->cd(1);
  g_mup->SetTitle("emittance #mu^{+}");
  g_mup->GetXaxis()->SetTitle("file");
  g_mup->GetYaxis()->SetTitle("#epsilon [nm x rad]");
  g_mup->GetYaxis()->SetTitleOffset(1.5);
  g_mup->SetMarkerStyle(20); 
  g_mup->SetMarkerStyle(20); 
  g_mup->SetMarkerColor(kRed);
  g_mup->GetHistogram()->SetMaximum(15.);
  g_mup->GetHistogram()->SetMinimum(10.);
  g_mup->Draw("AP");
  TLine* l_mup = new TLine(0.7,emittValue_mup_12345,5.3,emittValue_mup_12345); //emittance value of merged file
  l_mup->SetLineColor(kRed);
  l_mup->SetLineStyle(2);
  l_mup->Draw();
  c_emittance->cd(2);
  g_mum->SetTitle("emittance #mu^{-}");
  g_mum->GetXaxis()->SetTitle("file");
  g_mum->GetYaxis()->SetTitle("#epsilon [nm x rad]");
  g_mum->GetYaxis()->SetTitleOffset(1.5);
  g_mum->SetMarkerStyle(20); 
  g_mum->SetMarkerColor(kBlue);
  g_mum->GetHistogram()->SetMaximum(15.);
  g_mum->GetHistogram()->SetMinimum(10.);
  g_mum->Draw("AP");
  TLine* l_mum = new TLine(0.7,emittValue_mum_12345,5.3,emittValue_mum_12345); //emittance value of merged file
  l_mum->SetLineColor(kBlue);
  l_mum->SetLineStyle(2);
  l_mum->Draw();
  c_emittance->SaveAs("emittance_files.png");

}
