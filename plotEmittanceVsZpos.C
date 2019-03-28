// *****************************************************
//
// author: Alessandra Cappati
//         28/03/2019
// 
// usage: specify the input files (Data and MC), 
//        the output directory, and other options 
//        at the end of the script
//
// run with:
//        root -l -b -q plotEmittanceVsZpos.C++
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
    emittance=TMath::Sqrt(emittance);	   // [nm x rad]
  }

  return emittance;

}





// doTheHistos function: read root file and do histos 
void doTheHistos(TString inputFileName, TString label, double zEndTarget, double lTarget, TString plotOutputPath){

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


  // --- def vectors for emittance estimate

  // DATA vectors  

  // by now we don' t have complete estimate of the emittance on data.. 


  // MC vectors
  vector<double> vec_emittance_x_mup_MC_a;
  vector<double> vec_emittance_xprime_mup_MC_a;
  vector<double> vec_emittance_x_mup_MC_b;
  vector<double> vec_emittance_xprime_mup_MC_b;
  vector<double> vec_emittance_x_mup_MC_c;
  vector<double> vec_emittance_xprime_mup_MC_c;
  


  vector<double> vec_emittance_x_mum_MC_a;
  vector<double> vec_emittance_xprime_mum_MC_a;
  vector<double> vec_emittance_x_mum_MC_b;
  vector<double> vec_emittance_xprime_mum_MC_b;
  vector<double> vec_emittance_x_mum_MC_c;
  vector<double> vec_emittance_xprime_mum_MC_c;



  // ------------------------------------------------
  // estimate is done on a reference plane 
  Double_t z_ref  = zEndTarget; // [mm] z ref (end of the target)

  // --- divide by 3 the target length for computing emittance vs vertex position
  Double_t ltarget_by3 = lTarget/3.; 

  // reference points:  z_ref_3  <  z_ref_2  <  z_ref_1  <  z_ref
  Double_t z_ref_1 = z_ref -   ltarget_by3;
  Double_t z_ref_2 = z_ref - 2*ltarget_by3;
  Double_t z_ref_3 = z_ref - 3*ltarget_by3;
  // -------------------------------------------------


  // counters
  int count1_mup =0;
  int count2_mup =0;
  int count3_mup =0;
  int count1_mum =0;
  int count2_mum =0;
  int count3_mum =0;



  // ---------------------------
  // --- loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; z++){

    inputTree->GetEntry(z);



    // Call function to fill positron tracks info --------------------------
    // ONLY FOR DATA
    if(!isMC){
      // Max 10 positrons tracks allowed
      Double_t xpatz[10], thatz[10];
      Double_t zpos(zEndTarget); //Z position exit face of the Be target
      Int_t nxbe = ReturnBeamInfo(subdet, xh, zh, nhits, zpos, xpatz, thatz);  // number of positrons
      //hist_npos_Data->Fill(nxbe);

      for(Int_t j=0; j<nxbe; j++){
        //hist_xbe_positrons_Data->Fill(xpatz[j]); //position [mm]
        //hist_the_positrons_Data->Fill(thatz[j]); //angle [rad]
  
        if(j==1){
          
        }else if(j>1){

        }
      }
    }//end if !isMC
    // ----------------------------------------------------------------------


    
    // --- condition for candidate events
    if( p_mup > 0. && p_mum > 0. ) {


      // --- ONLY FOR MC 
      // -----------------------------
      // --- emittance in x 1D and 2D histos
      // -----------------------------
      // x  = x  @det30 - x  incoming e+
      // x' = x' @det30 - x' incoming e+
      if(isMC){       


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

     
        if(gen_vtx_mup[2] > z_ref_3  &&  gen_vtx_mup[2] <= z_ref_2){
          vec_emittance_x_mup_MC_a     .push_back(x_emittance_mup);
          vec_emittance_xprime_mup_MC_a.push_back(x_prime_emittance_mup);
          count1_mup++;
        }
        else if(gen_vtx_mup[2] > z_ref_2  &&  gen_vtx_mup[2] <= z_ref_1){
          vec_emittance_x_mup_MC_b     .push_back(x_emittance_mup);
          vec_emittance_xprime_mup_MC_b.push_back(x_prime_emittance_mup);
          count2_mup++;
        }
        else if(gen_vtx_mup[2] > z_ref_1  &&  gen_vtx_mup[2] < z_ref){
          vec_emittance_x_mup_MC_c     .push_back(x_emittance_mup);
          vec_emittance_xprime_mup_MC_c.push_back(x_prime_emittance_mup);
          count3_mup++;
        }



        // --- mu-
        // x  (extrapolation on reference plane): x_det30_atZref_mum = x_onDet30 - (z_onDet30 - z_ref)* px_onDet30 / pz_onDet30
        Double_t x_det30_atZref_mum = gen_pos_mum[0] - (gen_pos_mum[2] - z_ref)*(gen_pos_mum[3]/gen_pos_mum[5]); 
        // x' (extrapolation on reference plane): the direction remain the same as on det30 or at vtx
        Double_t pTot_genLev_mum = sqrt(gen_pos_mum[3]*gen_pos_mum[3] + gen_pos_mum[4]*gen_pos_mum[4] + gen_pos_mum[5]*gen_pos_mum[5]);
        Double_t x_prime_ondet30_mum = gen_pos_mum[3] / pTot_genLev_mum; 


        // emittance of mu-
        Double_t x_emittance_mum       = x_det30_atZref_mum  - x_atZref_eplus;
        Double_t x_prime_emittance_mum = x_prime_ondet30_mum - x_prime_atZref_eplus;
        
        
        if(gen_vtx_mum[2] > z_ref_3  &&  gen_vtx_mum[2] <= z_ref_2){
          vec_emittance_x_mum_MC_a     .push_back(x_emittance_mum);
          vec_emittance_xprime_mum_MC_a.push_back(x_prime_emittance_mum);
          count1_mum++;
        }
        else if(gen_vtx_mum[2] > z_ref_2  &&  gen_vtx_mum[2] <= z_ref_1){
          vec_emittance_x_mum_MC_b     .push_back(x_emittance_mum);
          vec_emittance_xprime_mum_MC_b.push_back(x_prime_emittance_mum);
          count2_mum++;
        }
        else if(gen_vtx_mum[2] > z_ref_1  &&  gen_vtx_mum[2] < z_ref){
          vec_emittance_x_mum_MC_c     .push_back(x_emittance_mum);
          vec_emittance_xprime_mum_MC_c.push_back(x_prime_emittance_mum);
          count3_mum++;
        }


        
      }        
                


      

    } // end if (p_mup > 0. && p_mum > 0.)

  }//end over tree entries 



  // ------------------
  //  plot histos
  // ------------------

  gStyle->SetOptStat(1);
 
  

  // beam info 
  // DATA ONLY
  if(!isMC){
  
    
  }
  // end beam info (DATA ONLY)


  // MC HISTOS
  if(isMC){

    
    // --- z pos vector
    Double_t vec_z_pos_mean[3];
    vec_z_pos_mean[0] = z_ref_3 + (z_ref_2 - z_ref_3)/2.;
    vec_z_pos_mean[1] = z_ref_2 + (z_ref_1 - z_ref_2)/2.;
    vec_z_pos_mean[2] = z_ref_1 + (z_ref   - z_ref_1)/2.;
    

    // --- mu + 
    Double_t  vec_emittanceValue_mup[3];
    vec_emittanceValue_mup[0] = getemittance(vec_emittance_x_mup_MC_a, vec_emittance_xprime_mup_MC_a);
    vec_emittanceValue_mup[1] = getemittance(vec_emittance_x_mup_MC_b, vec_emittance_xprime_mup_MC_b);
    vec_emittanceValue_mup[2] = getemittance(vec_emittance_x_mup_MC_c, vec_emittance_xprime_mup_MC_c);

        
    TCanvas* c_emittance_mup = new TCanvas("c_emittance_mup","c_emittance_mup");
    c_emittance_mup->cd();
    TGraph* g_mup = new TGraph(3,vec_z_pos_mean,vec_emittanceValue_mup);
    g_mup->SetTitle("emittance #mu^{+}");
    g_mup->SetMarkerStyle(20);
    g_mup->SetMarkerColor(kRed);
    g_mup->GetXaxis()->SetTitle("z [mm]");
    g_mup->GetYaxis()->SetTitle("Emittance [nm x rad]");
    g_mup->Draw("AP");    
    // TLine *line_0_mup = new TLine(z_ref_3,7.5,z_ref_3,10.);
    // line_0_mup->SetLineColor(kRed);
    // line_0_mup->SetLineStyle(1);
    // line_0_mup->Draw();
    // TLine *line_1_mup = new TLine(z_ref_2,7.5,z_ref_2,10);
    // line_1_mup->SetLineColor(kRed);
    // line_1_mup->SetLineStyle(2);
    // line_1_mup->Draw();
    // TLine *line_2_mup = new TLine(z_ref_1,7.5,z_ref_1,10.);
    // line_2_mup->SetLineColor(kRed);
    // line_2_mup->SetLineStyle(2);
    // line_2_mup->Draw();
    // TLine *line_3_mup = new TLine(z_ref,7.5,z_ref,10.);
    // line_3_mup->SetLineColor(kRed);
    // line_3_mup->SetLineStyle(1);
    // line_3_mup->Draw();
    c_emittance_mup->Update();
    c_emittance_mup->SaveAs((plotOutputPath + "/" + c_emittance_mup->GetName() + ".png"));


    cout<<vec_z_pos_mean[0]<<" "<<vec_z_pos_mean[1]<<" "<<vec_z_pos_mean[2]<<endl;
    cout<<vec_emittanceValue_mup[0]<<" "<<vec_emittanceValue_mup[1]<<" "<<vec_emittanceValue_mup[2]<<endl;
    cout<<count1_mup<<" "<<count2_mup<<" "<<count3_mup<<" "<<endl;






    // --- mu -
    Double_t  vec_emittanceValue_mum[3];
    vec_emittanceValue_mum[0] = getemittance(vec_emittance_x_mum_MC_a, vec_emittance_xprime_mum_MC_a);
    vec_emittanceValue_mum[1] = getemittance(vec_emittance_x_mum_MC_b, vec_emittance_xprime_mum_MC_b);
    vec_emittanceValue_mum[2] = getemittance(vec_emittance_x_mum_MC_c, vec_emittance_xprime_mum_MC_c);

        
    TCanvas* c_emittance_mum = new TCanvas("c_emittance_mum","c_emittance_mum");
    c_emittance_mum->cd();
    TGraph* g_mum = new TGraph(3,vec_z_pos_mean,vec_emittanceValue_mum);
    g_mum->SetTitle("emittance #mu^{-}");
    g_mum->SetMarkerStyle(20);
    g_mum->SetMarkerColor(kBlue);
    g_mum->GetXaxis()->SetTitle("z [mm]");
    g_mum->GetYaxis()->SetTitle("Emittance [nm x rad]");
    g_mum->Draw("AP");    
    // TLine *line_0_mum = new TLine(z_ref_3,7.5,z_ref_3,10.);
    // line_0_mum->SetLineColor(kRed);
    // line_0_mum->SetLineStyle(1);
    // line_0_mum->Draw();
    // TLine *line_1_mum = new TLine(z_ref_2,7.5,z_ref_2,10);
    // line_1_mum->SetLineColor(kRed);
    // line_1_mum->SetLineStyle(2);
    // line_1_mum->Draw();
    // TLine *line_2_mum = new TLine(z_ref_1,7.5,z_ref_1,10.);
    // line_2_mum->SetLineColor(kRed);
    // line_2_mum->SetLineStyle(2);
    // line_2_mum->Draw();
    // TLine *line_3_mum = new TLine(z_ref,7.5,z_ref,10.);
    // line_3_mum->SetLineColor(kRed);
    // line_3_mum->SetLineStyle(1);
    // line_3_mum->Draw();
    c_emittance_mum->Update();
    c_emittance_mum->SaveAs((plotOutputPath + "/" + c_emittance_mum->GetName() + ".png"));


    cout<<vec_z_pos_mean[0]<<" "<<vec_z_pos_mean[1]<<" "<<vec_z_pos_mean[2]<<endl;
    cout<<vec_emittanceValue_mum[0]<<" "<<vec_emittanceValue_mum[1]<<" "<<vec_emittanceValue_mum[2]<<endl;
    cout<<count1_mum<<" "<<count2_mum<<" "<<count3_mum<<" "<<endl;




  
  } //end plot mc plots

  cout<<" Plots done! =) "<<endl; 


}





// main function 
void plotEmittanceVsZpos(){

  // define input files 
  TString inputFile_Data_Aug2018_Be6cm = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/reco-333to352.root"; 
  TString inputFile_MC_Aug2018_Be6cm   = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/reco-mupmum.root";   
  TString inputFile_MC_Sep2018_Be6cm   = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/reco-mupmum-Be6cm.root";
  TString inputFile_MC_Sep2018_C6cm    = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/reco-mupmum-C6cm.root";
  TString inputFile_MC_Sep2018_C2cm    = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/reco-mupmum-C2cm.root";
  
  
  // define output path and make output directory 

  //TString plotOutputPath = "190328_EmittanceVsZpos_August2018_targetBe6cm_DATA";
  //TString plotOutputPath = "190328_EmittanceVsZpos_August2018_targetBe6cm_MC";
  //TString plotOutputPath = "190328_EmittanceVsZpos_September2018_targetBe6cm_MC";
  //TString plotOutputPath = "190328_EmittanceVsZpos_September2018_targetC6cm_MC";
  TString plotOutputPath = "190328_EmittanceVsZpos_September2018_targetC2cm_MC";
  gSystem->Exec(("mkdir -p "+plotOutputPath));



  // choose type of target
  //double zEndTarget = 10.*(457.9+3.-84.6);   // [mm] - dataset: AUGUST 2018    Be target 6 cm
  //double zEndTarget = 10.*(460.93+3.-82.78); // [mm] - dataset: SEPTEMBER 2018 Be target 6 cm and C target 6cm
  double zEndTarget = 10.*(460.93+1.-82.78); // [mm] - dataset: SEPTEMBER 2018 C  target 2 cm 


  // choose type of target
  //double lTarget = 60.; // [mm] - target 6 cm
  double lTarget = 20.; // [mm] - target 2 cm



  // --- call do the histos function
  // arguments: input file, label for data or MC

  //doTheHistos(inputFile_Data_Aug2018_Be6cm, "DATA", zEndTarget, lTarget, plotOutputPath);
  //doTheHistos(inputFile_MC_Aug2018_Be6cm,   "MC",   zEndTarget, lTarget, plotOutputPath);
  //doTheHistos(inputFile_MC_Sep2018_Be6cm,   "MC",   zEndTarget, lTarget, plotOutputPath);
  //doTheHistos(inputFile_MC_Sep2018_C6cm,    "MC",   zEndTarget, lTarget, plotOutputPath); 
  doTheHistos(inputFile_MC_Sep2018_C2cm,    "MC",   zEndTarget, lTarget, plotOutputPath);

}
