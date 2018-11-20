// *****************************************************
//
// author: Alessandra Cappati
//         13/11/2018
// 
// usage: specify the input files (Data and MC) and the 
//        output directory at the end of the script
//
// run with:
//        root -l -b -q plotVariables_BeamInfo.C++
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

#include "FillBeamInfo.h"


using namespace std;



//Extrapolate track function: track points are refitted
Double_t extrapolate_track_x(Double_t z0 ,Double_t x_pos_mum[12], Double_t x_pos_mum_err[12], Double_t z_x_pos_mum[12], Double_t& x_ext, Double_t& x_ext_err, Double_t& dx_on_dz_ext, Double_t& dx_on_dz_ext_err){

   //Magnetic field (box)
   Double_t zM= 17193-846;
   Double_t z1= zM-1000.;
   Double_t z2= zM+1000.;
   Double_t B = 1.7476;

   Double_t chi2 = 99999;

   //3 cases: z > z_magnet_out=z2, z < z_magnet_entry=z1, z1<z<z2 

   if (z0>z2) {

      Int_t npoints=0;

      for (Int_t k=0; k<12; k++) {
 
         if (z_x_pos_mum[k]>z2) npoints++;

         }

      if (npoints < 2) return chi2;
   
      Double_t z_x_pos_mum_err[12];
      for (Int_t k=0; k<12; k++) {z_x_pos_mum_err[k]=0;}

      //Change TGraph with TGraphErrors to use errors
      //TGraphErrors* graph = new TGraphErrors(npoints, z_x_pos_mum, x_pos_mum, z_x_pos_mum_err, x_pos_mum_err);

      TGraph* graph = new TGraph(npoints, z_x_pos_mum, x_pos_mum);

      //Linear fit using points after magnet

      TF1* line = new TF1("line", "[0]+[1]*(x-[2])");

      Double_t theta_min = 0.03;
      Double_t theta_max = 0.07;

      line->FixParameter(2,z0);


      graph->Fit(line,"Q");


      chi2 =  line->GetChisquare()/line->GetNDF();

      x_ext = line->GetParameter(0);
      x_ext_err = line->GetParError(0);

      dx_on_dz_ext = line->GetParameter(1);
      dx_on_dz_ext_err = line->GetParError(1);

      delete graph;
      delete line;

   } else if (z0 < z1)
   {

      Int_t npoints=0;

      for (Int_t k=0; k<12; k++) {
 
         if ((z_x_pos_mum[k]<-0.001 || z_x_pos_mum[k]>0.001)  && z_x_pos_mum[k] < 22700) npoints++;

         }    

      if (npoints < 2) return chi2;

      Double_t z_x_pos_mum_err[12];
      for (Int_t k=0; k<12; k++) {z_x_pos_mum_err[k]=0;}

      //Change TGraph with TGraphErrors to use errors
      //TGraphErrors* graph = new TGraphErrors(npoints, z_x_pos_mum, x_pos_mum, z_x_pos_mum_err, x_pos_mum_err);

      TGraph* graph = new TGraph(npoints, z_x_pos_mum, x_pos_mum);

      //Fit to the complete track: line + parabula + line, 3 free parameters

      TF1* trajectory = new TF1("trajectory", "(x<[3])*([0]+[1]*(x-[3]))+(x>[4])*([0]-[2]*([4]-[3])*([4]-[3])+([1]+2*[2]*([4]-[3]))*(x-[3]))+(x>[3])*(x<[4])*([0]+[1]*(x-[3])+[2]*(x-[3])*(x-[3]))");

      trajectory->FixParameter(3,z1);
      trajectory->FixParameter(4,z2);
  
      graph->Fit(trajectory,"Q");

      Double_t R = 1./(2.*trajectory->GetParameter(2));
      Double_t p = -B/(1e+9/TMath::C()) * R;
  
      chi2 = trajectory->GetChisquare()/trajectory->GetNDF();

      x_ext = trajectory->Eval(z0);
      x_ext_err = trajectory->GetParError(0);

      dx_on_dz_ext = trajectory->GetParameter(1);
      dx_on_dz_ext_err = trajectory->GetParError(1);


      delete graph;
      delete trajectory;

   } else 
   {
   cout<< "Extrapolation not defined" << endl;
   }

   return chi2;

}


// doTheHistos function: read root file and do histos 
void doTheHistos(TString inputFileName, TString label){

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
  Double_t gen_pos_mum[7]; // used for MC only
  Double_t gen_pos_mup[7]; // used for MC only

  TFile* inputFile = TFile::Open(inputFileName);
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
  }    

  // def histos 
  TH1F* hist_pMuPlus     = new TH1F("hist_pMuPlus",    "hist_pMuPlus",    30,16000.,30000.);
  TH1F* hist_pMuMinus    = new TH1F("hist_pMuMinus",   "hist_pMuMinus",   30,16000.,30000.);
  TH1F* hist_pTot        = new TH1F("hist_pTot",       "hist_pTot",       80,24800.,56800.);  
  TH1F* hist_pTot_smear03_bias000 = new TH1F("hist_pTot_smear03_bias000", "hist_pTot_smear03_bias000", 80,24800.,56800.); // used for MC only
  TH1F* hist_pTot_smear03_bias099 = new TH1F("hist_pTot_smear03_bias099", "hist_pTot_smear03_bias099", 80,24800.,56800.); // used for MC only
  TH1F* hist_pTot_smear03_bias101 = new TH1F("hist_pTot_smear03_bias101", "hist_pTot_smear03_bias101", 80,24800.,56800.); // used for MC only
  TH1F* hist_chi2MuPlus  = new TH1F("hist_chi2MuPlus", "hist_chi2MuPlus", 50,0.,10.);
  TH1F* hist_chi2MuMinus = new TH1F("hist_chi2MuMinus","hist_chi2MuMinus",50,0.,10.);
  // TH1F* hist_ThetaMuPlus  = new TH1F("hist_ThetaMuPlus","hist_ThetaMuPlus",10,0.,10.);    //angle in bending plane
  // TH1F* hist_ThetaMuMinus = new TH1F("hist_ThetaMuMinus","hist_ThetaMuMinus",10,0.,10.);  //angle in bending plane

  TH1F* hist_xh_det10_MuPlus = new TH1F("hist_xh_det10_MuPlus","hist_xh_det10_MuPlus",25,   -9.5,   9.5);
  TH1F* hist_xh_det20_MuPlus = new TH1F("hist_xh_det20_MuPlus","hist_xh_det20_MuPlus",25,   -9.5,   9.5);	  
  TH1F* hist_xh_det30_MuPlus = new TH1F("hist_xh_det30_MuPlus","hist_xh_det30_MuPlus",25,  -46.5,  46.5);	  
  TH1F* hist_xh_det31_MuPlus = new TH1F("hist_xh_det31_MuPlus","hist_xh_det31_MuPlus",25,  -46.5,  46.5);	  
  TH1F* hist_xh_det32_MuPlus = new TH1F("hist_xh_det32_MuPlus","hist_xh_det32_MuPlus",25,   40.0, 120.0);	  
  TH1F* hist_xh_det33_MuPlus = new TH1F("hist_xh_det33_MuPlus","hist_xh_det33_MuPlus",25, -120.0, -40.0);	  
  TH1F* hist_xh_det34_MuPlus = new TH1F("hist_xh_det34_MuPlus","hist_xh_det34_MuPlus",25,   93.5, 186.5);	  
  TH1F* hist_xh_det35_MuPlus = new TH1F("hist_xh_det35_MuPlus","hist_xh_det35_MuPlus",25, -186.5, -93.5);
  TH1F* hist_xh_det36_MuPlus = new TH1F("hist_xh_det36_MuPlus","hist_xh_det36_MuPlus",25,  150.0, 330.0);
  TH1F* hist_xh_det37_MuPlus = new TH1F("hist_xh_det37_MuPlus","hist_xh_det37_MuPlus",25, -330.0,-150.0);

  TH1F* hist_xh_det10_MuMinus = new TH1F("hist_xh_det10_MuMinus","hist_xh_det10_MuMinus",25,   -9.5,   9.5);
  TH1F* hist_xh_det20_MuMinus = new TH1F("hist_xh_det20_MuMinus","hist_xh_det20_MuMinus",25,   -9.5,   9.5);
  TH1F* hist_xh_det30_MuMinus = new TH1F("hist_xh_det30_MuMinus","hist_xh_det30_MuMinus",25,  -46.5,  46.5);
  TH1F* hist_xh_det31_MuMinus = new TH1F("hist_xh_det31_MuMinus","hist_xh_det31_MuMinus",25,  -46.5,  46.5);
  TH1F* hist_xh_det32_MuMinus = new TH1F("hist_xh_det32_MuMinus","hist_xh_det32_MuMinus",25,   40.0, 120.0);
  TH1F* hist_xh_det33_MuMinus = new TH1F("hist_xh_det33_MuMinus","hist_xh_det33_MuMinus",25, -120.0, -40.0);
  TH1F* hist_xh_det34_MuMinus = new TH1F("hist_xh_det34_MuMinus","hist_xh_det34_MuMinus",25,   93.5, 186.5);
  TH1F* hist_xh_det35_MuMinus = new TH1F("hist_xh_det35_MuMinus","hist_xh_det35_MuMinus",25, -186.5, -93.5);
  TH1F* hist_xh_det36_MuMinus = new TH1F("hist_xh_det36_MuMinus","hist_xh_det36_MuMinus",25,  150.0, 330.0);
  TH1F* hist_xh_det37_MuMinus = new TH1F("hist_xh_det37_MuMinus","hist_xh_det37_MuMinus",25, -330.0,-150.0);

  TH1F* hist_xh_det62_MuPlus  = new TH1F("hist_xh_det62_MuPlus", "hist_xh_det62_MuPlus", 50,-600.,600.);
  TH1F* hist_xh_det61_MuMinus = new TH1F("hist_xh_det61_MuMinus","hist_xh_det61_MuMinus",50,-600.,600.);

  TH1F* hist_xext_MuMinus = new TH1F("hist_xext_MuMinus","hist_xext_MuMinus",20,-50.,50.);
  TH1F* hist_xext_MuPlus  = new TH1F("hist_xext_MuPlus", "hist_xext_MuPlus", 20,-50.,50.);

  TH1F* hist_theta_xz_mup    = new TH1F("hist_theta_xz_mup", "hist_theta_xz_mup", 100,-0.025,0.025); //[rad]
  TH1F* hist_theta_xz_mum    = new TH1F("hist_theta_xz_mum", "hist_theta_xz_mum", 100,-0.025,0.025); //[rad]
  TH1F* hist_InvMass_mupmum  = new TH1F("hist_InvMass_mupmum", "hist_InvMass_mupmum", 100,100.,300.);

  TH1F* hist_xcross = new TH1F("hist_xcross", "hist_xcross", 15,-30.,30.);    // [mm]
  TH1F* hist_zcross = new TH1F("hist_zcross", "hist_zcross", 50,2000.,7000.); // [mm]

  TH1F* hist_npos = new TH1F("hist_npos", "N positrons", 11, -0.5, 10.5);
  TH1F* hist_xbe_positrons = new TH1F("hist_xbe_positrons", "Positron: Be exit point (mm)", 100, -15.0, 15.0);
  TH1F* hist_the_positrons = new TH1F("hist_the_positrons", "Positron: theta exit (urad)", 100, -500, 500);

  // loop over tree entries 
  Long64_t entries = inputTree->GetEntries();
  for(Long64_t z=0; z<entries; ++z){

    inputTree->GetEntry(z);

    
    // --- condition for candidate events
    if( p_mup > 0. && p_mum > 0. ) {


      // Call function to fill positron tracks info --------------------------
      // Max 10 positrons tracks allowed
      Double_t xpatz[10], thatz[10];
      Double_t zpos(454.9+3); //Z position exit face of the Be target
      Int_t nxbe = ReturnBeamInfo(subdet, xh, zh, nhits, zpos, xpatz, thatz);
      hist_npos->Fill(nxbe);
      for(Int_t j=0; j<nxbe; j++){
        hist_xbe_positrons->Fill(xpatz[j]);
        hist_the_positrons->Fill(thatz[j]*1e6);
      }
      // ----------------------------------------------------------------------


      // --- fill momentum histos 
      hist_pMuPlus->Fill(p_mup);      //momentum for mu plus
      hist_chi2MuPlus->Fill(chi2p);   //chi2 for mu plus tracks

      hist_pMuMinus->Fill(p_mum);     //momentum for mu minus
      hist_chi2MuMinus->Fill(chi2m);  //chi2 for mu minus tracks

      hist_pTot->Fill(p_mum + p_mup); //total momentum 

      // total momentum with smearing
      // -- gen_pos_mum[3] = px (MC)
      // -- gen_pos_mum[4] = py (MC)
      // -- gen_pos_mum[5] = pz (MC)
      if(isMC){  
        TRandom3* r0= new TRandom3(0);
        Float_t pSum0=0;
        pSum0+= sqrt(gen_pos_mum[3]*gen_pos_mum[3] + gen_pos_mum[4]*gen_pos_mum[4] + gen_pos_mum[5]*gen_pos_mum[5])*(1+r0->Gaus(0.,0.03));
        pSum0+= sqrt(gen_pos_mup[3]*gen_pos_mup[3] + gen_pos_mup[4]*gen_pos_mup[4] + gen_pos_mup[5]*gen_pos_mup[5])*(1+r0->Gaus(0.,0.03));
        hist_pTot_smear03_bias000->Fill(pSum0);

        TRandom3* r1= new TRandom3(0);
        Float_t pSum1=0;
        pSum1+= sqrt(gen_pos_mum[3]*gen_pos_mum[3] + gen_pos_mum[4]*gen_pos_mum[4] + gen_pos_mum[5]*gen_pos_mum[5])*(r1->Gaus(0.99,0.03));
        pSum1+= sqrt(gen_pos_mup[3]*gen_pos_mup[3] + gen_pos_mup[4]*gen_pos_mup[4] + gen_pos_mup[5]*gen_pos_mup[5])*(r1->Gaus(0.99,0.03));
        hist_pTot_smear03_bias099->Fill(pSum1);
        
        TRandom3* r2= new TRandom3(0);
        Float_t pSum2=0;
        pSum2+= sqrt(gen_pos_mum[3]*gen_pos_mum[3] + gen_pos_mum[4]*gen_pos_mum[4] + gen_pos_mum[5]*gen_pos_mum[5])*(r2->Gaus(1.01,0.03));
        pSum2+= sqrt(gen_pos_mup[3]*gen_pos_mup[3] + gen_pos_mup[4]*gen_pos_mup[4] + gen_pos_mup[5]*gen_pos_mup[5])*(r2->Gaus(1.01,0.03));
        hist_pTot_smear03_bias101->Fill(pSum2);
      }

      // histos for DTs
      hist_xh_det62_MuPlus->Fill(x_pos_DT_mup[0]);
      hist_xh_det61_MuMinus->Fill(x_pos_DT_mum[0]);


      // use extrapolate_track_x function 
      Double_t z0 = 10.*(569.5-84.6); // subdet 30
      Double_t x_ext_mup     = -9999;
      Double_t x_ext_mup_err = -9999;
      Double_t x_ext_mum     = -9999;
      Double_t x_ext_mum_err = -9999;
      Double_t dx_on_dz_ext_mup     = -9999;
      Double_t dx_on_dz_ext_mup_err = -9999;
      Double_t dx_on_dz_ext_mum     = -9999;
      Double_t dx_on_dz_ext_mum_err = -9999;
      Double_t chi2_mup = 99999;
      Double_t chi2_mum = 99999;
      //inputs are: z0, x_pos_mum, x_pos_mum_err, z_pos_mum
      //outputs are saved in:
      //x_ext (extrapolated x position at z=z0), 
      //x_ext_err, 
      //dx_on_dz_ext (extrapolated dx/dz at z=z0), 
      //dx_on_dz_ext_err
      chi2_mum = extrapolate_track_x(z0, x_pos_mum, x_pos_mum_err, z_x_pos_mum, x_ext_mum, x_ext_mum_err, dx_on_dz_ext_mum, dx_on_dz_ext_mum_err);
      if( chi2_mum < 9999. ){
        // cout << x_ext << " " << x_ext_err << endl;
        hist_xext_MuMinus->Fill(x_ext_mum);
      }
      chi2_mup = extrapolate_track_x(z0, x_pos_mup, x_pos_mup_err, z_x_pos_mup, x_ext_mup, x_ext_mup_err, dx_on_dz_ext_mup, dx_on_dz_ext_mup_err);
      if( chi2_mup < 9999. ){
        // cout << x_ext << " " << x_ext_err << endl;
        hist_xext_MuPlus->Fill(x_ext_mup);
      }


      // ---- x and z on the target (primary vertex of mu+ mu- production)
      // derived assuming: 
      // -   x_mup(z) = x_ext_mup + dx_on_dz_ext_mup * (z-z0)
      // -   x_mum(z) = x_ext_mum + dx_on_dz_ext_mum * (z-z0) 
      // -   x_mup(z_cross) = x_mum(z_cross)  
      Double_t z_cross = z0 + (x_ext_mum - x_ext_mup)/(dx_on_dz_ext_mup - dx_on_dz_ext_mum);
      Double_t x_cross = x_ext_mup + dx_on_dz_ext_mup * (z_cross - z0);

      // if MC we use gen level quantities
      // - on det 30: gen_pos_mum[0]=xh[i]; gen_pos_mum[1]=yh[i]; gen_pos_mum[2]=zh[i];
      //              gen_pos_mup[0]=xh[i]; gen_pos_mup[1]=yh[i]; gen_pos_mup[2]=zh[i];
      // - on det 31: gen_pos_mum[3]=xh[i]; gen_pos_mum[4]=yh[i]; gen_pos_mum[5]=zh[i];
      //              gen_pos_mup[3]=xh[i]; gen_pos_mup[4]=yh[i]; gen_pos_mup[5]=zh[i];
      if(isMC){
        Double_t gen_dx_on_dz_ext_mup = (gen_pos_mup[3] - gen_pos_mup[0])/(gen_pos_mup[5] - gen_pos_mup[2]);
        Double_t gen_dx_on_dz_ext_mum = (gen_pos_mum[3] - gen_pos_mum[0])/(gen_pos_mum[5] - gen_pos_mum[2]);
        z_cross = z0 + (gen_pos_mum[0] - gen_pos_mup[0])/(gen_dx_on_dz_ext_mup - gen_dx_on_dz_ext_mum);
      }

      hist_xcross->Fill(x_cross);
      hist_zcross->Fill(z_cross);
      // ----


      // --- mu+ mu- invariant mass
      Double_t theta_xz_mup = atan(dx_on_dz_ext_mup);
      Double_t theta_xz_mum = atan(dx_on_dz_ext_mum);

      Double_t restMass_mu = 105.66; // [MeV]

      Double_t E_mup = sqrt((restMass_mu*restMass_mu) + (p_mup*p_mup));
      Double_t E_mum = sqrt((restMass_mu*restMass_mu) + (p_mum*p_mum));

      Double_t invMass_mupmum = sqrt(2*restMass_mu*restMass_mu + 2*(E_mup*E_mum - p_mup*p_mum*cos(theta_xz_mup-theta_xz_mum)));
      
      hist_theta_xz_mup->Fill(theta_xz_mup);
      hist_theta_xz_mum->Fill(theta_xz_mum);
      hist_InvMass_mupmum->Fill(invMass_mupmum);
      // -----


      // histos for hits in Si det
      for(int i=0;i<100;i++){

        if(xh[i] == -9990 ){continue;}   // skip empty events

        // mu plus events
        if(itrack[i] == -1){             // itrack = -1 for mu+  (like PDG ID)
	  if(subdet[i] == 10) {hist_xh_det10_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 20) {hist_xh_det20_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 30) {hist_xh_det30_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 31) {hist_xh_det31_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 32) {hist_xh_det32_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 33) {hist_xh_det33_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 34) {hist_xh_det34_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 35) {hist_xh_det35_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 36) {hist_xh_det36_MuPlus->Fill(xh[i]);}
	  if(subdet[i] == 37) {hist_xh_det37_MuPlus->Fill(xh[i]);}
        }
        // mu minus events
        else if(itrack[i] == 1){         // itrack = +1 for mu-  (like PDG ID)
	  if(subdet[i] == 10) {hist_xh_det10_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 20) {hist_xh_det20_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 30) {hist_xh_det30_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 31) {hist_xh_det31_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 32) {hist_xh_det32_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 33) {hist_xh_det33_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 34) {hist_xh_det34_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 35) {hist_xh_det35_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 36) {hist_xh_det36_MuMinus->Fill(xh[i]);}
	  if(subdet[i] == 37) {hist_xh_det37_MuMinus->Fill(xh[i]);}
        }
        else continue;

      } // end loop over i

    } // end if (p_mup > 0. && p_mum > 0.)

  }//end over tree entries 



  // save histos into a root file 
  TString outFileName = "plotVariables_" + label + ".root";
  TFile* fOutHistos = new TFile(outFileName,"recreate");
  fOutHistos->cd();

  hist_pMuPlus->Write(hist_pMuPlus->GetName());
  hist_pMuMinus->Write(hist_pMuMinus->GetName());
  hist_pTot->Write(hist_pTot->GetName());  
  if(isMC){
    hist_pTot_smear03_bias000->Write(hist_pTot_smear03_bias000->GetName());
    hist_pTot_smear03_bias099->Write(hist_pTot_smear03_bias099->GetName());
    hist_pTot_smear03_bias101->Write(hist_pTot_smear03_bias101->GetName());
  }    
  hist_chi2MuPlus->Write(hist_chi2MuPlus->GetName());
  hist_chi2MuMinus->Write(hist_chi2MuMinus->GetName());

  hist_theta_xz_mup->Write(hist_theta_xz_mup->GetName());
  hist_theta_xz_mum->Write(hist_theta_xz_mum->GetName());
  hist_InvMass_mupmum->Write(hist_InvMass_mupmum->GetName());

  hist_xcross->Write(hist_xcross->GetName());
  hist_zcross->Write(hist_zcross->GetName());
                                
  hist_xh_det10_MuPlus->Write(hist_xh_det10_MuPlus->GetName());
  hist_xh_det20_MuPlus->Write(hist_xh_det20_MuPlus->GetName());
  hist_xh_det30_MuPlus->Write(hist_xh_det30_MuPlus->GetName());
  hist_xh_det31_MuPlus->Write(hist_xh_det31_MuPlus->GetName());
  hist_xh_det32_MuPlus->Write(hist_xh_det32_MuPlus->GetName());
  hist_xh_det33_MuPlus->Write(hist_xh_det33_MuPlus->GetName());
  hist_xh_det34_MuPlus->Write(hist_xh_det34_MuPlus->GetName());
  hist_xh_det35_MuPlus->Write(hist_xh_det35_MuPlus->GetName());
  hist_xh_det36_MuPlus->Write(hist_xh_det36_MuPlus->GetName());
  hist_xh_det37_MuPlus->Write(hist_xh_det37_MuPlus->GetName());
                              
  hist_xh_det10_MuMinus->Write(hist_xh_det10_MuMinus->GetName()); 
  hist_xh_det20_MuMinus->Write(hist_xh_det20_MuMinus->GetName()); 
  hist_xh_det30_MuMinus->Write(hist_xh_det30_MuMinus->GetName()); 
  hist_xh_det31_MuMinus->Write(hist_xh_det31_MuMinus->GetName()); 
  hist_xh_det32_MuMinus->Write(hist_xh_det32_MuMinus->GetName()); 
  hist_xh_det33_MuMinus->Write(hist_xh_det33_MuMinus->GetName()); 
  hist_xh_det34_MuMinus->Write(hist_xh_det34_MuMinus->GetName()); 
  hist_xh_det35_MuMinus->Write(hist_xh_det35_MuMinus->GetName()); 
  hist_xh_det36_MuMinus->Write(hist_xh_det36_MuMinus->GetName()); 
  hist_xh_det37_MuMinus->Write(hist_xh_det37_MuMinus->GetName()); 
                              
  hist_xh_det62_MuPlus->Write(hist_xh_det62_MuPlus->GetName());  
  hist_xh_det61_MuMinus->Write(hist_xh_det61_MuMinus->GetName()); 
                              
  hist_xext_MuMinus->Write(hist_xext_MuMinus->GetName());
  hist_xext_MuPlus->Write(hist_xext_MuPlus->GetName()); 

  // save some positron beam info
  hist_npos->Write(hist_npos->GetName());
  hist_xbe_positrons->Write(hist_xbe_positrons->GetName());
  hist_the_positrons->Write(hist_the_positrons->GetName());
  //

  fOutHistos->Close();
  delete fOutHistos;

  cout<<" root file filled and created!"<<endl;

}



// data MC comparison function
void dataMCComparison(TString plotDataMCOutputPath){
  
  // read data file 
  TFile *inFile_Data = TFile::Open("plotVariables_DATA.root");

  // beam info
  TH1F* hist_npos = (TH1F*)inFile_Data->Get("hist_npos");
  TH1F* hist_xbe_positrons = (TH1F*)inFile_Data->Get("hist_xbe_positrons");
  TH1F* hist_the_positrons = (TH1F*)inFile_Data->Get("hist_the_positrons");
  // 

  TH1F* hist_pMuPlus_Data     = (TH1F*)inFile_Data->Get("hist_pMuPlus");
  TH1F* hist_pMuMinus_Data    = (TH1F*)inFile_Data->Get("hist_pMuMinus");
  TH1F* hist_pTot_Data        = (TH1F*)inFile_Data->Get("hist_pTot");      
  TH1F* hist_chi2MuPlus_Data  = (TH1F*)inFile_Data->Get("hist_chi2MuPlus");
  TH1F* hist_chi2MuMinus_Data = (TH1F*)inFile_Data->Get("hist_chi2MuMinus");

  TH1F* hist_theta_xz_mup_Data   = (TH1F*)inFile_Data->Get("hist_theta_xz_mup");
  TH1F* hist_theta_xz_mum_Data   = (TH1F*)inFile_Data->Get("hist_theta_xz_mum");
  TH1F* hist_InvMass_mupmum_Data = (TH1F*)inFile_Data->Get("hist_InvMass_mupmum");

  TH1F* hist_xcross_Data = (TH1F*)inFile_Data->Get("hist_xcross"); 
  TH1F* hist_zcross_Data = (TH1F*)inFile_Data->Get("hist_zcross");
                                
  TH1F* hist_xh_det10_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det10_MuPlus");
  TH1F* hist_xh_det20_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det20_MuPlus");
  TH1F* hist_xh_det30_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det30_MuPlus");
  TH1F* hist_xh_det31_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det31_MuPlus");
  TH1F* hist_xh_det32_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det32_MuPlus");
  TH1F* hist_xh_det33_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det33_MuPlus");
  TH1F* hist_xh_det34_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det34_MuPlus");
  TH1F* hist_xh_det35_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det35_MuPlus");
  TH1F* hist_xh_det36_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det36_MuPlus");
  TH1F* hist_xh_det37_MuPlus_Data = (TH1F*)inFile_Data->Get("hist_xh_det37_MuPlus");
                              
  TH1F* hist_xh_det10_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det10_MuMinus"); 
  TH1F* hist_xh_det20_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det20_MuMinus"); 
  TH1F* hist_xh_det30_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det30_MuMinus"); 
  TH1F* hist_xh_det31_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det31_MuMinus"); 
  TH1F* hist_xh_det32_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det32_MuMinus"); 
  TH1F* hist_xh_det33_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det33_MuMinus"); 
  TH1F* hist_xh_det34_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det34_MuMinus"); 
  TH1F* hist_xh_det35_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det35_MuMinus"); 
  TH1F* hist_xh_det36_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det36_MuMinus"); 
  TH1F* hist_xh_det37_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det37_MuMinus"); 
                              
  TH1F* hist_xh_det62_MuPlus_Data  = (TH1F*)inFile_Data->Get("hist_xh_det62_MuPlus");  
  TH1F* hist_xh_det61_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xh_det61_MuMinus"); 
                              
  TH1F* hist_xext_MuMinus_Data = (TH1F*)inFile_Data->Get("hist_xext_MuMinus");
  TH1F* hist_xext_MuPlus_Data  = (TH1F*)inFile_Data->Get("hist_xext_MuPlus"); 


  // read MC file 
  TFile *inFile_MC = TFile::Open("plotVariables_MC.root");

  TH1F* hist_pMuPlus_MC     = (TH1F*)inFile_MC->Get("hist_pMuPlus");
  TH1F* hist_pMuMinus_MC    = (TH1F*)inFile_MC->Get("hist_pMuMinus");
  TH1F* hist_pTot_MC        = (TH1F*)inFile_MC->Get("hist_pTot");      
  TH1F* hist_pTot_smear03_bias000_MC = (TH1F*)inFile_MC->Get("hist_pTot_smear03_bias000");      
  TH1F* hist_pTot_smear03_bias099_MC = (TH1F*)inFile_MC->Get("hist_pTot_smear03_bias099");    
  TH1F* hist_pTot_smear03_bias101_MC = (TH1F*)inFile_MC->Get("hist_pTot_smear03_bias101");        
  TH1F* hist_chi2MuPlus_MC  = (TH1F*)inFile_MC->Get("hist_chi2MuPlus");
  TH1F* hist_chi2MuMinus_MC = (TH1F*)inFile_MC->Get("hist_chi2MuMinus");

  TH1F* hist_theta_xz_mup_MC   = (TH1F*)inFile_MC->Get("hist_theta_xz_mup");
  TH1F* hist_theta_xz_mum_MC   = (TH1F*)inFile_MC->Get("hist_theta_xz_mum");
  TH1F* hist_InvMass_mupmum_MC = (TH1F*)inFile_MC->Get("hist_InvMass_mupmum");

  TH1F* hist_xcross_MC = (TH1F*)inFile_MC->Get("hist_xcross"); 
  TH1F* hist_zcross_MC = (TH1F*)inFile_MC->Get("hist_zcross");
                                 
  TH1F* hist_xh_det10_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det10_MuPlus");
  TH1F* hist_xh_det20_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det20_MuPlus");
  TH1F* hist_xh_det30_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det30_MuPlus");
  TH1F* hist_xh_det31_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det31_MuPlus");
  TH1F* hist_xh_det32_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det32_MuPlus");
  TH1F* hist_xh_det33_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det33_MuPlus");
  TH1F* hist_xh_det34_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det34_MuPlus");
  TH1F* hist_xh_det35_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det35_MuPlus");
  TH1F* hist_xh_det36_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det36_MuPlus");
  TH1F* hist_xh_det37_MuPlus_MC = (TH1F*)inFile_MC->Get("hist_xh_det37_MuPlus");
                              
  TH1F* hist_xh_det10_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det10_MuMinus"); 
  TH1F* hist_xh_det20_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det20_MuMinus"); 
  TH1F* hist_xh_det30_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det30_MuMinus"); 
  TH1F* hist_xh_det31_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det31_MuMinus"); 
  TH1F* hist_xh_det32_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det32_MuMinus"); 
  TH1F* hist_xh_det33_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det33_MuMinus"); 
  TH1F* hist_xh_det34_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det34_MuMinus"); 
  TH1F* hist_xh_det35_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det35_MuMinus"); 
  TH1F* hist_xh_det36_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det36_MuMinus"); 
  TH1F* hist_xh_det37_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det37_MuMinus"); 
                              
  TH1F* hist_xh_det62_MuPlus_MC  = (TH1F*)inFile_MC->Get("hist_xh_det62_MuPlus");  
  TH1F* hist_xh_det61_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xh_det61_MuMinus"); 
                              
  TH1F* hist_xext_MuMinus_MC = (TH1F*)inFile_MC->Get("hist_xext_MuMinus");
  TH1F* hist_xext_MuPlus_MC  = (TH1F*)inFile_MC->Get("hist_xext_MuPlus"); 


  gStyle->SetOptStat(0);
 
  //plot histos

  // beam info 
  TCanvas * c_pos = new TCanvas("c_pos","c_pos",800, 1200);
  c_pos->Divide(1,3);
  c_pos->cd(1); hist_npos->Draw();
  c_pos->cd(2); hist_xbe_positrons->Draw();
  c_pos->cd(3); hist_the_positrons->Draw();
  c_pos->SaveAs((plotDataMCOutputPath + "/" + c_pos->GetName() + ".png"));



  // pMuPlus plot
  TCanvas* c_pMuPlus = new TCanvas("c_pMuPlus","c_pMuPlus");
  c_pMuPlus->cd();
  hist_pMuPlus_MC->SetTitle("p #mu^{+}");
  hist_pMuPlus_MC->GetXaxis()->SetTitle("p #mu^{+} [MeV]");
  hist_pMuPlus_MC->GetYaxis()->SetTitle("events");
  hist_pMuPlus_MC->SetLineColor(kRed);
  hist_pMuPlus_MC->SetFillColor(kRed-10);
  hist_pMuPlus_MC->Scale(hist_pMuPlus_Data->Integral() / hist_pMuPlus_MC->Integral()); //normalize MC to Data
  hist_pMuPlus_MC->SetMaximum(1.5 * max(hist_pMuPlus_MC->GetMaximum(),hist_pMuPlus_Data->GetMaximum()));
  hist_pMuPlus_MC->Draw("hist");
  hist_pMuPlus_Data->SetMarkerStyle(20);
  hist_pMuPlus_Data->SetMarkerColor(kRed);
  hist_pMuPlus_Data->SetLineColor(kBlack);
  hist_pMuPlus_Data->Draw("samepe");
  TLegend* l_pMuPlus = new TLegend(0.72,0.67,0.98,0.97);
  l_pMuPlus->AddEntry(hist_pMuPlus_MC,"MC #mu^{+}","f");
  l_pMuPlus->AddEntry((TObject*)0,Form("entries: %.2f",hist_pMuPlus_MC->GetEntries()),"");
  l_pMuPlus->AddEntry((TObject*)0,Form("mean: %.2f",hist_pMuPlus_MC->GetMean()),"");
  l_pMuPlus->AddEntry(hist_pMuPlus_Data, "Data #mu^{+}", "pl");
  l_pMuPlus->AddEntry((TObject*)0,Form("entries: %.2f",hist_pMuPlus_Data->GetEntries()),"");
  l_pMuPlus->AddEntry((TObject*)0,Form("mean: %.2f",hist_pMuPlus_Data->GetMean()),"");
  l_pMuPlus->SetFillColor(kWhite);
  l_pMuPlus->SetLineColor(kBlack);
  l_pMuPlus->SetTextFont(43);
  l_pMuPlus->SetTextSize(20);
  l_pMuPlus->Draw();
  c_pMuPlus->Update();
  c_pMuPlus->SaveAs((plotDataMCOutputPath + "/" + c_pMuPlus->GetName() + ".png"));

  // pMuMinus plot
  TCanvas* c_pMuMinus = new TCanvas("c_pMuMinus","c_pMuMinus");
  c_pMuMinus->cd();
  hist_pMuMinus_MC->SetTitle("p #mu^{-}");
  hist_pMuMinus_MC->GetXaxis()->SetTitle("p #mu^{-} [MeV]");
  hist_pMuMinus_MC->GetYaxis()->SetTitle("events");
  hist_pMuMinus_MC->SetLineColor(kBlue);
  hist_pMuMinus_MC->SetFillColor(kBlue-10);
  hist_pMuMinus_MC->Scale(hist_pMuMinus_Data->Integral() / hist_pMuMinus_MC->Integral()); //normalize MC to Data
  hist_pMuMinus_MC->SetMaximum(1.5 * max(hist_pMuMinus_MC->GetMaximum(),hist_pMuMinus_Data->GetMaximum()));
  hist_pMuMinus_MC->Draw("hist");
  hist_pMuMinus_Data->SetMarkerStyle(20);
  hist_pMuMinus_Data->SetMarkerColor(kBlue);
  hist_pMuMinus_Data->SetLineColor(kBlack);
  hist_pMuMinus_Data->Draw("samepe");
  TLegend* l_pMuMinus = new TLegend(0.72,0.67,0.98,0.97);
  l_pMuMinus->AddEntry(hist_pMuMinus_MC,"MC #mu^{-}","f");
  l_pMuMinus->AddEntry((TObject*)0,Form("entries: %.2f",hist_pMuMinus_MC->GetEntries()),"");
  l_pMuMinus->AddEntry((TObject*)0,Form("mean: %.2f",hist_pMuMinus_MC->GetMean()),"");
  l_pMuMinus->AddEntry(hist_pMuMinus_Data, "Data #mu^{-}", "pl");
  l_pMuMinus->AddEntry((TObject*)0,Form("entries: %.2f",hist_pMuMinus_Data->GetEntries()),"");
  l_pMuMinus->AddEntry((TObject*)0,Form("mean: %.2f",hist_pMuMinus_Data->GetMean()),"");
  l_pMuMinus->SetFillColor(kWhite);
  l_pMuMinus->SetLineColor(kBlack);
  l_pMuMinus->SetTextFont(43);
  l_pMuMinus->SetTextSize(20);
  l_pMuMinus->Draw();
  c_pMuMinus->Update();
  c_pMuMinus->SaveAs((plotDataMCOutputPath + "/" + c_pMuMinus->GetName() + ".png"));
 
  TCanvas* c_pTot = new TCanvas("c_pTot","c_pTot");
  c_pTot->cd();
  hist_pTot_MC->SetTitle("p #mu^{+} + p #mu^{-}");
  hist_pTot_MC->GetXaxis()->SetTitle("p #mu^{+} + p #mu^{-} [MeV]");
  hist_pTot_MC->GetYaxis()->SetTitle("events");
  hist_pTot_MC->SetLineColor(kViolet+1);
  hist_pTot_MC->SetLineWidth(2);
  //hist_pTot_MC->SetLineColor(kViolet-6);
  //hist_pTot_MC->SetFillColorAlpha(kViolet-4, 0.571);
  hist_pTot_MC->Scale(hist_pTot_Data->Integral() / hist_pTot_MC->Integral()); //normalize MC to Data
  hist_pTot_smear03_bias000_MC->SetLineColor(kOrange+7);
  hist_pTot_smear03_bias000_MC->SetLineWidth(2);
  //hist_pTot_smear03_bias000_MC->SetLineColor(kOrange+7);
  //hist_pTot_smear03_bias000_MC->SetFillColorAlpha(kOrange-3, 0.571);
  hist_pTot_smear03_bias000_MC->Scale(hist_pTot_Data->Integral() / hist_pTot_smear03_bias000_MC->Integral()); //normalize MC to Data
  hist_pTot_Data->SetMarkerStyle(20);
  hist_pTot_Data->SetMarkerColor(kBlack);
  hist_pTot_Data->SetLineColor(kBlack);
  hist_pTot_MC->SetMaximum(1.2 * max(max(hist_pTot_MC->GetMaximum(),hist_pTot_smear03_bias000_MC->GetMaximum()),hist_pTot_Data->GetMaximum()));
  hist_pTot_MC->Draw("hist");
  hist_pTot_smear03_bias000_MC->Draw("histsame");
  hist_pTot_Data->Draw("samepe");
  TLegend* l_pTot = new TLegend(0.76,0.55,0.98,0.96);
  l_pTot->AddEntry(hist_pTot_MC,"MC","f");
  l_pTot->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_MC->GetEntries()),"");
  l_pTot->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_MC->GetMean()),"");  
  l_pTot->AddEntry(hist_pTot_smear03_bias000_MC,"MC smear","f");
  l_pTot->AddEntry((TObject*)0,"gauss(0.00,0.03)","");
  l_pTot->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_smear03_bias000_MC->GetEntries()),"");
  l_pTot->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_smear03_bias000_MC->GetMean()),"");
  l_pTot->AddEntry(hist_pTot_Data, "Data", "pl");
  l_pTot->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_Data->GetEntries()),"");
  l_pTot->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_Data->GetMean()),"");
  l_pTot->SetFillColor(kWhite);
  l_pTot->SetLineColor(kBlack);
  l_pTot->SetTextFont(43);
  l_pTot->SetTextSize(14);
  l_pTot->Draw();
  c_pTot->Update();
  c_pTot->SaveAs((plotDataMCOutputPath + "/" + c_pTot->GetName() + ".png"));


  TCanvas* c_pTot_smear = new TCanvas("c_pTot_smear","c_pTot_smear");
  c_pTot_smear->cd();
  hist_pTot_smear03_bias099_MC->SetTitle("p #mu^{+} + p #mu^{-}");
  hist_pTot_smear03_bias099_MC->GetXaxis()->SetTitle("p #mu^{+} + p #mu^{-} [MeV]");
  hist_pTot_smear03_bias099_MC->GetYaxis()->SetTitle("events");
  hist_pTot_smear03_bias099_MC->SetLineColor(kGreen+1);
  hist_pTot_smear03_bias099_MC->SetLineWidth(2);
  //hist_pTot_smear03_bias099_MC->SetLineColor(kGreen+3);
  //hist_pTot_smear03_bias099_MC->SetFillColorAlpha(kGreen-9, 0.571);
  hist_pTot_smear03_bias099_MC->Scale(hist_pTot_Data->Integral() / hist_pTot_smear03_bias099_MC->Integral()); //normalize MC to Data
  hist_pTot_smear03_bias101_MC->SetLineColor(kAzure-3);
  hist_pTot_smear03_bias101_MC->SetLineWidth(2);
  //hist_pTot_smear03_bias101_MC->SetLineColor(kAzure+2);
  //hist_pTot_smear03_bias101_MC->SetFillColorAlpha(kAzure+6, 0.571);
  hist_pTot_smear03_bias101_MC->Scale(hist_pTot_Data->Integral() / hist_pTot_smear03_bias101_MC->Integral()); //normalize MC to Data
  hist_pTot_smear03_bias000_MC->SetLineColor(kOrange+7);
  hist_pTot_smear03_bias000_MC->SetLineWidth(2);
  //hist_pTot_smear03_bias000_MC->SetLineColor(kOrange+7);
  //hist_pTot_smear03_bias000_MC->SetFillColorAlpha(kOrange-3, 0.571);
  hist_pTot_smear03_bias000_MC->Scale(hist_pTot_Data->Integral() / hist_pTot_smear03_bias000_MC->Integral()); //normalize MC to Data
  hist_pTot_Data->SetMarkerStyle(20);
  hist_pTot_Data->SetMarkerColor(kBlack);
  hist_pTot_Data->SetLineColor(kBlack);
  hist_pTot_smear03_bias099_MC->SetMaximum(1.2 * max(max(hist_pTot_smear03_bias099_MC->GetMaximum(),hist_pTot_smear03_bias101_MC->GetMaximum()),max(hist_pTot_Data->GetMaximum(),hist_pTot_smear03_bias000_MC->GetMaximum())));
  hist_pTot_smear03_bias099_MC->Draw("histsame");
  hist_pTot_smear03_bias101_MC->Draw("histsame");
  hist_pTot_smear03_bias000_MC->Draw("histsame");
  hist_pTot_Data->Draw("samepe");
  TLegend* l_pTot_smear = new TLegend(0.76,0.36,0.98,0.97);
  l_pTot_smear->AddEntry(hist_pTot_smear03_bias099_MC,"MC smear","f");
  l_pTot_smear->AddEntry((TObject*)0,"gauss(0.99,0.03)","");
  l_pTot_smear->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_smear03_bias099_MC->GetEntries()),"");
  l_pTot_smear->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_smear03_bias099_MC->GetMean()),"");
  l_pTot_smear->AddEntry(hist_pTot_smear03_bias101_MC,"MC smear","f");
  l_pTot_smear->AddEntry((TObject*)0,"gauss(1.01,0.03)","");
  l_pTot_smear->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_smear03_bias101_MC->GetEntries()),"");
  l_pTot_smear->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_smear03_bias101_MC->GetMean()),"");
  l_pTot_smear->AddEntry(hist_pTot_smear03_bias000_MC,"MC smear","f");
  l_pTot_smear->AddEntry((TObject*)0,"gauss(0.00,0.03)","");
  l_pTot_smear->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_smear03_bias000_MC->GetEntries()),"");
  l_pTot_smear->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_smear03_bias000_MC->GetMean()),"");
  l_pTot_smear->AddEntry(hist_pTot_Data, "Data", "pl");
  l_pTot_smear->AddEntry((TObject*)0,Form("entries: %.2f",hist_pTot_Data->GetEntries()),"");
  l_pTot_smear->AddEntry((TObject*)0,Form("mean: %.2f",hist_pTot_Data->GetMean()),"");
  l_pTot_smear->SetFillColor(kWhite);
  l_pTot_smear->SetLineColor(kBlack);
  l_pTot_smear->SetTextFont(43);
  l_pTot_smear->SetTextSize(14);
  l_pTot_smear->Draw();
  c_pTot_smear->Update();
  c_pTot_smear->SaveAs((plotDataMCOutputPath + "/" + c_pTot_smear->GetName() + ".png"));
 
 

  TCanvas* c_chi2MuPlus = new TCanvas("c_chi2MuPlus","c_chi2MuPlus");
  c_chi2MuPlus->cd();
  hist_chi2MuPlus_MC->SetTitle("#Chi^{2} #mu^{+}");
  hist_chi2MuPlus_MC->GetXaxis()->SetTitle("#Chi^{2} #mu^{+}");
  hist_chi2MuPlus_MC->GetYaxis()->SetTitle("events");
  hist_chi2MuPlus_MC->SetLineColor(kRed);   
  hist_chi2MuPlus_MC->SetFillColor(kRed-10);
  hist_chi2MuPlus_MC->Scale(hist_chi2MuPlus_Data->Integral() / hist_chi2MuPlus_MC->Integral()); //normalize MC to Data
  hist_chi2MuPlus_MC->SetMaximum(1.2 * max(hist_chi2MuPlus_MC->GetMaximum(),hist_chi2MuPlus_Data->GetMaximum()));
  hist_chi2MuPlus_MC->Draw("hist");
  hist_chi2MuPlus_Data->SetMarkerStyle(20);  
  hist_chi2MuPlus_Data->SetMarkerColor(kRed);
  hist_chi2MuPlus_Data->SetLineColor(kBlack);
  hist_chi2MuPlus_Data->Draw("samepe");
  TLegend* l_chi2MuPlus = new TLegend(0.72,0.67,0.98,0.97);
  l_chi2MuPlus->AddEntry(hist_chi2MuPlus_MC,"MC","f");
  l_chi2MuPlus->AddEntry((TObject*)0,Form("entries: %.2f",hist_chi2MuPlus_MC->GetEntries()),"");
  l_chi2MuPlus->AddEntry((TObject*)0,Form("mean: %.2f",hist_chi2MuPlus_MC->GetMean()),"");
  l_chi2MuPlus->AddEntry(hist_chi2MuPlus_Data, "Data", "pl");
  l_chi2MuPlus->AddEntry((TObject*)0,Form("entries: %.2f",hist_chi2MuPlus_Data->GetEntries()),"");
  l_chi2MuPlus->AddEntry((TObject*)0,Form("mean: %.2f",hist_chi2MuPlus_Data->GetMean()),"");
  l_chi2MuPlus->SetFillColor(kWhite);
  l_chi2MuPlus->SetLineColor(kBlack);
  l_chi2MuPlus->SetTextFont(43);
  l_chi2MuPlus->SetTextSize(20);
  l_chi2MuPlus->Draw();
  c_chi2MuPlus->Update();
  c_chi2MuPlus->SaveAs((plotDataMCOutputPath + "/" + c_chi2MuPlus->GetName() + ".png"));
 

  TCanvas* c_chi2MuMinus = new TCanvas("c_chi2MuMinus","c_chi2MuMinus");
  c_chi2MuMinus->cd();
  hist_chi2MuMinus_MC->SetTitle("#Chi^{2} #mu^{-}");
  hist_chi2MuMinus_MC->GetXaxis()->SetTitle("#Chi^{2} #mu^{-}");
  hist_chi2MuMinus_MC->GetYaxis()->SetTitle("events");
  hist_chi2MuMinus_MC->SetLineColor(kBlue);   
  hist_chi2MuMinus_MC->SetFillColor(kBlue-10);
  hist_chi2MuMinus_MC->Scale(hist_chi2MuMinus_Data->Integral() / hist_chi2MuMinus_MC->Integral()); //normalize MC to Data
  hist_chi2MuMinus_MC->SetMaximum(1.2 * max(hist_chi2MuMinus_MC->GetMaximum(),hist_chi2MuMinus_Data->GetMaximum()));
  hist_chi2MuMinus_MC->Draw("hist");
  hist_chi2MuMinus_Data->SetMarkerStyle(20);  
  hist_chi2MuMinus_Data->SetMarkerColor(kBlue);
  hist_chi2MuMinus_Data->SetLineColor(kBlack);
  hist_chi2MuMinus_Data->Draw("samepe");
  TLegend* l_chi2MuMinus = new TLegend(0.72,0.67,0.98,0.97);
  l_chi2MuMinus->AddEntry(hist_chi2MuMinus_MC,"MC","f");
  l_chi2MuMinus->AddEntry((TObject*)0,Form("entries: %.2f",hist_chi2MuMinus_MC->GetEntries()),"");
  l_chi2MuMinus->AddEntry((TObject*)0,Form("mean: %.2f",hist_chi2MuMinus_MC->GetMean()),"");
  l_chi2MuMinus->AddEntry(hist_chi2MuMinus_Data, "Data", "pl");
  l_chi2MuMinus->AddEntry((TObject*)0,Form("entries: %.2f",hist_chi2MuMinus_Data->GetEntries()),"");
  l_chi2MuMinus->AddEntry((TObject*)0,Form("mean: %.2f",hist_chi2MuMinus_Data->GetMean()),"");
  l_chi2MuMinus->SetFillColor(kWhite);
  l_chi2MuMinus->SetLineColor(kBlack);
  l_chi2MuMinus->SetTextFont(43);
  l_chi2MuMinus->SetTextSize(20);
  l_chi2MuMinus->Draw();
  c_chi2MuMinus->Update();
  c_chi2MuMinus->SaveAs((plotDataMCOutputPath + "/" + c_chi2MuMinus->GetName() + ".png"));

  
  TCanvas* c_theta_xz_mup = new TCanvas("c_theta_xz_mup","c_theta_xz_mup");
  c_theta_xz_mup->cd();
  hist_theta_xz_mup_MC->SetTitle("#theta xz #mu^{+}");
  hist_theta_xz_mup_MC->GetXaxis()->SetTitle("#theta [rad]");
  hist_theta_xz_mup_MC->GetYaxis()->SetTitle("events");
  hist_theta_xz_mup_MC->SetLineColor(kRed);   
  hist_theta_xz_mup_MC->SetFillColor(kRed-10);
  hist_theta_xz_mup_MC->Scale(hist_theta_xz_mup_Data->Integral() / hist_theta_xz_mup_MC->Integral()); //normalize MC to Data
  hist_theta_xz_mup_MC->SetMaximum(1.2 * max(hist_theta_xz_mup_MC->GetMaximum(),hist_theta_xz_mup_Data->GetMaximum()));
  hist_theta_xz_mup_MC->Draw("hist");
  hist_theta_xz_mup_Data->SetMarkerStyle(20);  
  hist_theta_xz_mup_Data->SetMarkerColor(kBlack);
  hist_theta_xz_mup_Data->SetLineColor(kBlack);
  hist_theta_xz_mup_Data->Draw("samepe");
  TLegend* l_theta_xz_mup = new TLegend(0.75,0.67,0.98,0.95);
  l_theta_xz_mup->AddEntry(hist_theta_xz_mup_MC,"MC","f");
  l_theta_xz_mup->AddEntry((TObject*)0,Form("entries: %.2f",hist_theta_xz_mup_MC->GetEntries()),"");
  l_theta_xz_mup->AddEntry((TObject*)0,Form("mean: %.2f",hist_theta_xz_mup_MC->GetMean()),"");
  l_theta_xz_mup->AddEntry(hist_theta_xz_mup_Data, "Data", "pl");
  l_theta_xz_mup->AddEntry((TObject*)0,Form("entries: %.2f",hist_theta_xz_mup_Data->GetEntries()),"");
  l_theta_xz_mup->AddEntry((TObject*)0,Form("mean: %.2f",hist_theta_xz_mup_Data->GetMean()),"");
  l_theta_xz_mup->SetFillColor(kWhite);
  l_theta_xz_mup->SetLineColor(kBlack);
  l_theta_xz_mup->SetTextFont(43);
  l_theta_xz_mup->SetTextSize(16);
  l_theta_xz_mup->Draw();
  c_theta_xz_mup->Update();
  c_theta_xz_mup->SaveAs((plotDataMCOutputPath + "/" + c_theta_xz_mup->GetName() + ".png"));


  TCanvas* c_theta_xz_mum = new TCanvas("c_theta_xz_mum","c_theta_xz_mum");
  c_theta_xz_mum->cd();
  hist_theta_xz_mum_MC->SetTitle("#theta xz #mu^{-}");
  hist_theta_xz_mum_MC->GetXaxis()->SetTitle("#theta [rad]");
  hist_theta_xz_mum_MC->GetYaxis()->SetTitle("events");
  hist_theta_xz_mum_MC->SetLineColor(kBlue);   
  hist_theta_xz_mum_MC->SetFillColor(kBlue-10);
  hist_theta_xz_mum_MC->Scale(hist_theta_xz_mum_Data->Integral() / hist_theta_xz_mum_MC->Integral()); //normalize MC to Data
  hist_theta_xz_mum_MC->SetMaximum(1.2 * max(hist_theta_xz_mum_MC->GetMaximum(),hist_theta_xz_mum_Data->GetMaximum()));
  hist_theta_xz_mum_MC->Draw("hist");
  hist_theta_xz_mum_Data->SetMarkerStyle(20);  
  hist_theta_xz_mum_Data->SetMarkerColor(kBlack);
  hist_theta_xz_mum_Data->SetLineColor(kBlack);
  hist_theta_xz_mum_Data->Draw("samepe");
  TLegend* l_theta_xz_mum = new TLegend(0.75,0.67,0.98,0.95);
  l_theta_xz_mum->AddEntry(hist_theta_xz_mum_MC,"MC","f");
  l_theta_xz_mum->AddEntry((TObject*)0,Form("entries: %.2f",hist_theta_xz_mum_MC->GetEntries()),"");
  l_theta_xz_mum->AddEntry((TObject*)0,Form("mean: %.2f",hist_theta_xz_mum_MC->GetMean()),"");
  l_theta_xz_mum->AddEntry(hist_theta_xz_mum_Data, "Data", "pl");
  l_theta_xz_mum->AddEntry((TObject*)0,Form("entries: %.2f",hist_theta_xz_mum_Data->GetEntries()),"");
  l_theta_xz_mum->AddEntry((TObject*)0,Form("mean: %.2f",hist_theta_xz_mum_Data->GetMean()),"");
  l_theta_xz_mum->SetFillColor(kWhite);
  l_theta_xz_mum->SetLineColor(kBlack);
  l_theta_xz_mum->SetTextFont(43);
  l_theta_xz_mum->SetTextSize(16);
  l_theta_xz_mum->Draw();
  c_theta_xz_mum->Update();
  c_theta_xz_mum->SaveAs((plotDataMCOutputPath + "/" + c_theta_xz_mum->GetName() + ".png"));



  TCanvas* c_InvMass_mupmum = new TCanvas("c_InvMass_mupmum","c_InvMass_mupmum");
  c_InvMass_mupmum->cd();
  hist_InvMass_mupmum_MC->SetTitle("Invariant Mass #mu^{+} #mu^{-}");
  hist_InvMass_mupmum_MC->GetXaxis()->SetTitle("m #mu^{+} #mu^{-} [MeV]");
  hist_InvMass_mupmum_MC->GetYaxis()->SetTitle("events");
  hist_InvMass_mupmum_MC->SetLineColor(kGreen+2);   
  hist_InvMass_mupmum_MC->SetFillColor(kGreen-9);
  hist_InvMass_mupmum_MC->Scale(hist_InvMass_mupmum_Data->Integral() / hist_InvMass_mupmum_MC->Integral()); //normalize MC to Data
  hist_InvMass_mupmum_MC->SetMaximum(1.2 * max(hist_InvMass_mupmum_MC->GetMaximum(),hist_InvMass_mupmum_Data->GetMaximum()));
  hist_InvMass_mupmum_MC->Draw("hist");
  hist_InvMass_mupmum_Data->SetMarkerStyle(20);  
  hist_InvMass_mupmum_Data->SetMarkerColor(kBlack);
  hist_InvMass_mupmum_Data->SetLineColor(kBlack);
  hist_InvMass_mupmum_Data->Draw("samepe");
  TLegend* l_InvMass_mupmum = new TLegend(0.75,0.67,0.98,0.95);
  l_InvMass_mupmum->AddEntry(hist_InvMass_mupmum_MC,"MC","f");
  l_InvMass_mupmum->AddEntry((TObject*)0,Form("entries: %.2f",hist_InvMass_mupmum_MC->GetEntries()),"");
  l_InvMass_mupmum->AddEntry((TObject*)0,Form("mean: %.2f",hist_InvMass_mupmum_MC->GetMean()),"");
  l_InvMass_mupmum->AddEntry(hist_InvMass_mupmum_Data, "Data", "pl");
  l_InvMass_mupmum->AddEntry((TObject*)0,Form("entries: %.2f",hist_InvMass_mupmum_Data->GetEntries()),"");
  l_InvMass_mupmum->AddEntry((TObject*)0,Form("mean: %.2f",hist_InvMass_mupmum_Data->GetMean()),"");
  l_InvMass_mupmum->SetFillColor(kWhite);
  l_InvMass_mupmum->SetLineColor(kBlack);
  l_InvMass_mupmum->SetTextFont(43);
  l_InvMass_mupmum->SetTextSize(16);
  l_InvMass_mupmum->Draw();
  c_InvMass_mupmum->Update();
  c_InvMass_mupmum->SaveAs((plotDataMCOutputPath + "/" + c_InvMass_mupmum->GetName() + ".png"));


  TCanvas* c_xcross = new TCanvas("c_xcross","c_xcross");
  c_xcross->cd();
  hist_xcross_MC->SetTitle("x Primary Vertex");
  hist_xcross_MC->GetXaxis()->SetTitle("x cross [mm]");
  hist_xcross_MC->GetYaxis()->SetTitle("events");
  hist_xcross_MC->SetLineColor(kOrange+7);   
  hist_xcross_MC->SetFillColor(kOrange-3);
  hist_xcross_MC->Scale(hist_xcross_Data->Integral() / hist_xcross_MC->Integral()); //normalize MC to Data
  hist_xcross_MC->SetMaximum(1.2 * max(hist_xcross_MC->GetMaximum(),hist_xcross_Data->GetMaximum()));
  hist_xcross_MC->Draw("hist");
  hist_xcross_Data->SetMarkerStyle(20);  
  hist_xcross_Data->SetMarkerColor(kBlack);
  hist_xcross_Data->SetLineColor(kBlack);
  hist_xcross_Data->Draw("samepe");
  TLegend* l_xcross = new TLegend(0.75,0.67,0.98,0.95);
  l_xcross->AddEntry(hist_xcross_MC,"MC","f");
  l_xcross->AddEntry((TObject*)0,Form("entries: %.2f",hist_xcross_MC->GetEntries()),"");
  l_xcross->AddEntry((TObject*)0,Form("mean: %.2f",hist_xcross_MC->GetMean()),"");
  l_xcross->AddEntry(hist_xcross_Data, "Data", "pl");
  l_xcross->AddEntry((TObject*)0,Form("entries: %.2f",hist_xcross_Data->GetEntries()),"");
  l_xcross->AddEntry((TObject*)0,Form("mean: %.2f",hist_xcross_Data->GetMean()),"");
  l_xcross->SetFillColor(kWhite);
  l_xcross->SetLineColor(kBlack);
  l_xcross->SetTextFont(43);
  l_xcross->SetTextSize(16);
  l_xcross->Draw();
  c_xcross->Update();
  c_xcross->SaveAs((plotDataMCOutputPath + "/" + c_xcross->GetName() + ".png"));


  TCanvas* c_zcross = new TCanvas("c_zcross","c_zcross");
  c_zcross->cd();
  hist_zcross_MC->SetTitle("z Primary Vertex");
  hist_zcross_MC->GetXaxis()->SetTitle("z cross [mm]");
  hist_zcross_MC->GetYaxis()->SetTitle("events");
  hist_zcross_MC->SetLineColor(kOrange+7);   
  hist_zcross_MC->SetFillColor(kOrange-3);
  hist_zcross_MC->Scale(hist_zcross_Data->Integral() / hist_zcross_MC->Integral()); //normalize MC to Data
  hist_zcross_MC->SetMaximum(1.2 * max(hist_zcross_MC->GetMaximum(),hist_zcross_Data->GetMaximum()));
  hist_zcross_MC->Draw("hist");
  hist_zcross_Data->SetMarkerStyle(20);  
  hist_zcross_Data->SetMarkerColor(kBlack);
  hist_zcross_Data->SetLineColor(kBlack);
  hist_zcross_Data->Draw("samepe");
  TLegend* l_zcross = new TLegend(0.75,0.67,0.98,0.95);
  l_zcross->AddEntry(hist_zcross_MC,"MC","f");
  l_zcross->AddEntry((TObject*)0,Form("entries: %.2f",hist_zcross_MC->GetEntries()),"");
  l_zcross->AddEntry((TObject*)0,Form("mean: %.2f",hist_zcross_MC->GetMean()),"");
  l_zcross->AddEntry(hist_zcross_Data, "Data", "pl");
  l_zcross->AddEntry((TObject*)0,Form("entries: %.2f",hist_zcross_Data->GetEntries()),"");
  l_zcross->AddEntry((TObject*)0,Form("mean: %.2f",hist_zcross_Data->GetMean()),"");
  l_zcross->SetFillColor(kWhite);
  l_zcross->SetLineColor(kBlack);
  l_zcross->SetTextFont(43);
  l_zcross->SetTextSize(16);
  l_zcross->Draw();
  c_zcross->Update();
  TLine *line_zcross_centre = new TLine(3763.3,gPad->GetUymin(),3763.3,gPad->GetUymax());
  line_zcross_centre->SetLineColor(kRed);
  line_zcross_centre->SetLineStyle(2);
  line_zcross_centre->Draw();
  c_zcross->SaveAs((plotDataMCOutputPath + "/" + c_zcross->GetName() + ".png"));
 

  // xh in det30
  TCanvas* c_det30 = new TCanvas("c_det30","c_det30");
  c_det30->cd();
  hist_xh_det30_MuMinus_MC->SetTitle("xh in det30");
  hist_xh_det30_MuMinus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det30_MuMinus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det30_MuMinus_MC->SetLineColor(kBlue);
  //hist_xh_det30_MuMinus_MC->SetFillColorAlpha(kBlue-10, 0.571);
  hist_xh_det30_MuMinus_MC->SetFillColor(kBlue-10);
  hist_xh_det30_MuMinus_MC->Scale(hist_xh_det30_MuMinus_Data->Integral() / hist_xh_det30_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det30_MuPlus_MC->SetTitle("");
  hist_xh_det30_MuPlus_MC->SetLineColor(kRed);
  //hist_xh_det30_MuPlus_MC->SetFillColorAlpha(kRed-10, 0.571);
  hist_xh_det30_MuPlus_MC->SetFillColor(kRed-10);
  hist_xh_det30_MuPlus_MC->Scale(hist_xh_det30_MuPlus_Data->Integral() / hist_xh_det30_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det30_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det30_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det30_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det30_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det30_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det30_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det30_MuMinus_MC->SetMaximum(1.5 * max(max(hist_xh_det30_MuMinus_MC->GetMaximum(),hist_xh_det30_MuMinus_Data->GetMaximum()),max(hist_xh_det30_MuPlus_MC->GetMaximum(),hist_xh_det30_MuPlus_Data->GetMaximum())));
  hist_xh_det30_MuMinus_MC->Draw("hist");  
  hist_xh_det30_MuPlus_MC->Draw("histsame");
  hist_xh_det30_MuMinus_Data->Draw("samepe"); 
  hist_xh_det30_MuPlus_Data->Draw("samepe");
  TLegend* l_det30 = new TLegend(0.79,0.49,0.98,0.97);
  l_det30->AddEntry(hist_xh_det30_MuMinus_MC,"MC #mu^{-}","f");
  l_det30->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det30_MuMinus_MC->GetEntries()),"");
  l_det30->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det30_MuMinus_MC->GetMean()),"");
  l_det30->AddEntry(hist_xh_det30_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det30->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det30_MuMinus_Data->GetEntries()),"");
  l_det30->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det30_MuMinus_Data->GetMean()),"");
  l_det30->AddEntry(hist_xh_det30_MuPlus_MC,"MC #mu^{+}","f");
  l_det30->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det30_MuPlus_MC->GetEntries()),"");
  l_det30->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det30_MuPlus_MC->GetMean()),"");
  l_det30->AddEntry(hist_xh_det30_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det30->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det30_MuPlus_Data->GetEntries()),"");
  l_det30->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det30_MuPlus_Data->GetMean()),"");
  l_det30->SetFillColor(kWhite);
  l_det30->SetLineColor(kBlack);
  l_det30->SetTextFont(43);
  l_det30->SetTextSize(14);
  l_det30->Draw();
  c_det30->Update();
  c_det30->SaveAs((plotDataMCOutputPath + "/" + c_det30->GetName() + ".png"));  


  // xh in det31
  TCanvas* c_det31 = new TCanvas("c_det31","c_det31");
  c_det31->cd();
  hist_xh_det31_MuMinus_MC->SetTitle("xh in det31");
  hist_xh_det31_MuMinus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det31_MuMinus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det31_MuMinus_MC->SetLineColor(kBlue);
  //hist_xh_det31_MuMinus_MC->SetFillColorAlpha(kBlue-10, 0.571);
  hist_xh_det31_MuMinus_MC->SetFillColor(kBlue-10);
  hist_xh_det31_MuMinus_MC->Scale(hist_xh_det31_MuMinus_Data->Integral() / hist_xh_det31_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det31_MuPlus_MC->SetTitle("");
  hist_xh_det31_MuPlus_MC->SetLineColor(kRed);
  //hist_xh_det31_MuPlus_MC->SetFillColorAlpha(kRed-10, 0.571);
  hist_xh_det31_MuPlus_MC->SetFillColor(kRed-10);
  hist_xh_det31_MuPlus_MC->Scale(hist_xh_det31_MuPlus_Data->Integral() / hist_xh_det31_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det31_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det31_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det31_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det31_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det31_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det31_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det31_MuMinus_MC->SetMaximum(1.5 * max(max(hist_xh_det31_MuMinus_MC->GetMaximum(),hist_xh_det31_MuMinus_Data->GetMaximum()),max(hist_xh_det31_MuPlus_MC->GetMaximum(),hist_xh_det31_MuPlus_Data->GetMaximum())));
  hist_xh_det31_MuMinus_MC->Draw("hist");  
  hist_xh_det31_MuPlus_MC->Draw("histsame");
  hist_xh_det31_MuMinus_Data->Draw("samepe"); 
  hist_xh_det31_MuPlus_Data->Draw("samepe");
  TLegend* l_det31 = new TLegend(0.79,0.49,0.98,0.97);
  l_det31->AddEntry(hist_xh_det31_MuMinus_MC,"MC #mu^{-}","f");
  l_det31->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det31_MuMinus_MC->GetEntries()),"");
  l_det31->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det31_MuMinus_MC->GetMean()),"");
  l_det31->AddEntry(hist_xh_det31_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det31->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det31_MuMinus_Data->GetEntries()),"");
  l_det31->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det31_MuMinus_Data->GetMean()),"");
  l_det31->AddEntry(hist_xh_det31_MuPlus_MC,"MC #mu^{+}","f");
  l_det31->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det31_MuPlus_MC->GetEntries()),"");
  l_det31->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det31_MuPlus_MC->GetMean()),"");
  l_det31->AddEntry(hist_xh_det31_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det31->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det31_MuPlus_Data->GetEntries()),"");
  l_det31->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det31_MuPlus_Data->GetMean()),"");
  l_det31->SetFillColor(kWhite);
  l_det31->SetLineColor(kBlack);
  l_det31->SetTextFont(43);
  l_det31->SetTextSize(14);
  l_det31->Draw();
  c_det31->Update();
  c_det31->SaveAs((plotDataMCOutputPath + "/" + c_det31->GetName() + ".png"));  

  
  // xh in det32
  TCanvas* c_det32 = new TCanvas("c_det32","c_det32");
  c_det32->cd();
  hist_xh_det32_MuMinus_MC->SetTitle("xh in det32");
  hist_xh_det32_MuMinus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det32_MuMinus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det32_MuMinus_MC->SetLineColor(kBlue);
  hist_xh_det32_MuMinus_MC->SetFillColor(kBlue-10);
  hist_xh_det32_MuMinus_MC->Scale(hist_xh_det32_MuMinus_Data->Integral() / hist_xh_det32_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det32_MuPlus_MC->SetTitle("");
  hist_xh_det32_MuPlus_MC->SetLineColor(kRed);
  hist_xh_det32_MuPlus_MC->SetFillColor(kRed-10);
  //hist_xh_det32_MuPlus_MC->Scale(hist_xh_det32_MuPlus_Data->Integral() / hist_xh_det32_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det32_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det32_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det32_MuMinus_Data->SetLineColor(kBlack);  
  hist_xh_det32_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det32_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det32_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det32_MuMinus_MC->SetMaximum(1.5 * max(max(hist_xh_det32_MuMinus_MC->GetMaximum(),hist_xh_det32_MuMinus_Data->GetMaximum()),max(hist_xh_det32_MuPlus_MC->GetMaximum(),hist_xh_det32_MuPlus_Data->GetMaximum())));
  hist_xh_det32_MuMinus_MC->Draw("hist");
  hist_xh_det32_MuPlus_MC->Draw("histsame");
  hist_xh_det32_MuMinus_Data->Draw("samepe");
  hist_xh_det32_MuPlus_Data->Draw("samepe");
  TLegend* l_det32 = new TLegend(0.79,0.49,0.98,0.97); 
  l_det32->AddEntry(hist_xh_det32_MuMinus_MC,"MC #mu^{-}","f");
  l_det32->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det32_MuMinus_MC->GetEntries()),"");
  l_det32->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det32_MuMinus_MC->GetMean()),"");
  l_det32->AddEntry(hist_xh_det32_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det32->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det32_MuMinus_Data->GetEntries()),"");
  l_det32->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det32_MuMinus_Data->GetMean()),"");
  l_det32->AddEntry(hist_xh_det32_MuPlus_MC,"MC #mu^{+}","f");
  l_det32->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det32_MuPlus_MC->GetEntries()),"");
  l_det32->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det32_MuPlus_MC->GetMean()),"");
  l_det32->AddEntry(hist_xh_det32_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det32->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det32_MuPlus_Data->GetEntries()),"");
  l_det32->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det32_MuPlus_Data->GetMean()),"");
  l_det32->SetFillColor(kWhite);
  l_det32->SetLineColor(kBlack);
  l_det32->SetTextFont(43);
  l_det32->SetTextSize(14);
  l_det32->Draw();
  c_det32->Update();
  c_det32->SaveAs((plotDataMCOutputPath + "/" + c_det32->GetName() + ".png"));   

  // xh in det33
  TCanvas* c_det33 = new TCanvas("c_det33","c_det33");
  c_det33->cd();
  hist_xh_det33_MuPlus_MC->SetTitle("xh in det33");
  hist_xh_det33_MuPlus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det33_MuPlus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det33_MuPlus_MC->SetLineColor(kRed);
  hist_xh_det33_MuPlus_MC->SetFillColor(kRed-10);
  hist_xh_det33_MuPlus_MC->Scale(hist_xh_det33_MuPlus_Data->Integral() / hist_xh_det33_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det33_MuMinus_MC->SetTitle("");
  hist_xh_det33_MuMinus_MC->SetLineColor(kBlue);
  hist_xh_det33_MuMinus_MC->SetFillColor(kBlue-10);
  //hist_xh_det33_MuMinus_MC->Scale(hist_xh_det33_MuMinus_Data->Integral() / hist_xh_det33_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det33_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det33_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det33_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det33_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det33_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det33_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det33_MuPlus_MC->SetMaximum(1.5 * max(max(hist_xh_det33_MuMinus_MC->GetMaximum(),hist_xh_det33_MuMinus_Data->GetMaximum()),max(hist_xh_det33_MuPlus_MC->GetMaximum(),hist_xh_det33_MuPlus_Data->GetMaximum())));
  hist_xh_det33_MuPlus_MC->Draw("hist");
  hist_xh_det33_MuMinus_MC->Draw("histsame");
  hist_xh_det33_MuPlus_Data->Draw("samepe");
  hist_xh_det33_MuMinus_Data->Draw("samepe");
  TLegend* l_det33 = new TLegend(0.12,0.49,0.31,0.97);
  l_det33->AddEntry(hist_xh_det33_MuPlus_MC,"MC #mu^{+}","f");
  l_det33->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det33_MuPlus_MC->GetEntries()),"");
  l_det33->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det33_MuPlus_MC->GetMean()),"");
  l_det33->AddEntry(hist_xh_det33_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det33->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det33_MuPlus_Data->GetEntries()),"");
  l_det33->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det33_MuPlus_Data->GetMean()),"");
  l_det33->AddEntry(hist_xh_det33_MuMinus_MC,"MC #mu^{-}","f");
  l_det33->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det33_MuMinus_MC->GetEntries()),"");
  l_det33->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det33_MuMinus_MC->GetMean()),"");
  l_det33->AddEntry(hist_xh_det33_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det33->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det33_MuMinus_Data->GetEntries()),"");
  l_det33->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det33_MuMinus_Data->GetMean()),"");
  l_det33->SetFillColor(kWhite);
  l_det33->SetLineColor(kBlack);
  l_det33->SetTextFont(43);
  l_det33->SetTextSize(14);
  l_det33->Draw();
  c_det33->Update();
  c_det33->SaveAs((plotDataMCOutputPath + "/" + c_det33->GetName() + ".png"));   

  // xh in det34
  TCanvas* c_det34 = new TCanvas("c_det34","c_det34");
  c_det34->cd();
  hist_xh_det34_MuMinus_MC->SetTitle("xh in det34");
  hist_xh_det34_MuMinus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det34_MuMinus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det34_MuMinus_MC->SetLineColor(kBlue);
  hist_xh_det34_MuMinus_MC->SetFillColor(kBlue-10);
  hist_xh_det34_MuMinus_MC->Scale(hist_xh_det34_MuMinus_Data->Integral() / hist_xh_det34_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det34_MuPlus_MC->SetTitle("");
  hist_xh_det34_MuPlus_MC->SetLineColor(kRed);
  hist_xh_det34_MuPlus_MC->SetFillColor(kRed-10);
  //hist_xh_det34_MuPlus_MC->Scale(hist_xh_det34_MuPlus_Data->Integral() / hist_xh_det34_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det34_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det34_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det34_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det34_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det34_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det34_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det34_MuMinus_MC->SetMaximum(1.5 * max(max(hist_xh_det34_MuMinus_MC->GetMaximum(),hist_xh_det34_MuMinus_Data->GetMaximum()),max(hist_xh_det34_MuPlus_MC->GetMaximum(),hist_xh_det34_MuPlus_Data->GetMaximum())));
  hist_xh_det34_MuMinus_MC->Draw("hist");
  hist_xh_det34_MuPlus_MC->Draw("histsame");
  hist_xh_det34_MuMinus_Data->Draw("samepe");
  hist_xh_det34_MuPlus_Data->Draw("samepe");
  TLegend* l_det34 = new TLegend(0.79,0.49,0.98,0.97);
  l_det34->AddEntry(hist_xh_det34_MuMinus_MC,"MC #mu^{-}","f");
  l_det34->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det34_MuMinus_MC->GetEntries()),"");
  l_det34->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det34_MuMinus_MC->GetMean()),"");
  l_det34->AddEntry(hist_xh_det34_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det34->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det34_MuMinus_Data->GetEntries()),"");
  l_det34->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det34_MuMinus_Data->GetMean()),"");
  l_det34->AddEntry(hist_xh_det34_MuPlus_MC,"MC #mu^{+}","f");
  l_det34->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det34_MuPlus_MC->GetEntries()),"");
  l_det34->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det34_MuPlus_MC->GetMean()),"");
  l_det34->AddEntry(hist_xh_det34_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det34->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det34_MuPlus_Data->GetEntries()),"");
  l_det34->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det34_MuPlus_Data->GetMean()),"");
  l_det34->SetFillColor(kWhite);
  l_det34->SetLineColor(kBlack);
  l_det34->SetTextFont(43);
  l_det34->SetTextSize(14);
  l_det34->Draw();
  c_det34->Update();
  c_det34->SaveAs((plotDataMCOutputPath + "/" + c_det34->GetName() + ".png")); 

  // xh in det35
  TCanvas* c_det35 = new TCanvas("c_det35","c_det35");
  c_det35->cd();
  hist_xh_det35_MuPlus_MC->SetTitle("xh in det35");
  hist_xh_det35_MuPlus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det35_MuPlus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det35_MuPlus_MC->SetLineColor(kRed);
  hist_xh_det35_MuPlus_MC->SetFillColor(kRed-10);
  hist_xh_det35_MuPlus_MC->Scale(hist_xh_det35_MuPlus_Data->Integral() / hist_xh_det35_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det35_MuMinus_MC->SetTitle("");
  hist_xh_det35_MuMinus_MC->SetLineColor(kBlue);
  hist_xh_det35_MuMinus_MC->SetFillColor(kBlue-10);
  //hist_xh_det35_MuMinus_MC->Scale(hist_xh_det35_MuMinus_Data->Integral() / hist_xh_det35_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det35_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det35_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det35_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det35_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det35_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det35_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det35_MuPlus_MC->SetMaximum(1.5 * max(max(hist_xh_det35_MuMinus_MC->GetMaximum(),hist_xh_det35_MuMinus_Data->GetMaximum()),max(hist_xh_det35_MuPlus_MC->GetMaximum(),hist_xh_det35_MuPlus_Data->GetMaximum())));
  hist_xh_det35_MuPlus_MC->Draw("hist");
  hist_xh_det35_MuMinus_MC->Draw("histsame");
  hist_xh_det35_MuPlus_Data->Draw("samepe");
  hist_xh_det35_MuMinus_Data->Draw("samepe");
  TLegend* l_det35 = new TLegend(0.12,0.49,0.31,0.97);
  l_det35->AddEntry(hist_xh_det35_MuPlus_MC,"MC #mu^{+}","f");
  l_det35->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det35_MuPlus_MC->GetEntries()),"");
  l_det35->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det35_MuPlus_MC->GetMean()),"");
  l_det35->AddEntry(hist_xh_det35_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det35->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det35_MuPlus_Data->GetEntries()),"");
  l_det35->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det35_MuPlus_Data->GetMean()),"");
  l_det35->AddEntry(hist_xh_det35_MuMinus_MC,"MC #mu^{-}","f");
  l_det35->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det35_MuMinus_MC->GetEntries()),"");
  l_det35->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det35_MuMinus_MC->GetMean()),"");
  l_det35->AddEntry(hist_xh_det35_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det35->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det35_MuMinus_Data->GetEntries()),"");
  l_det35->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det35_MuMinus_Data->GetMean()),"");
  l_det35->SetFillColor(kWhite);
  l_det35->SetLineColor(kBlack);
  l_det35->SetTextFont(43);
  l_det35->SetTextSize(14);
  l_det35->Draw();
  c_det35->Update();
  c_det35->SaveAs((plotDataMCOutputPath + "/" + c_det35->GetName() + ".png"));  

  // xh in det36
  TCanvas* c_det36 = new TCanvas("c_det36","c_det36");
  c_det36->cd();
  hist_xh_det36_MuMinus_MC->SetTitle("xh in det36");
  hist_xh_det36_MuMinus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det36_MuMinus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det36_MuMinus_MC->SetLineColor(kBlue);
  hist_xh_det36_MuMinus_MC->SetFillColor(kBlue-10);
  hist_xh_det36_MuMinus_MC->Scale(hist_xh_det36_MuMinus_Data->Integral() / hist_xh_det36_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det36_MuPlus_MC->SetTitle("");
  hist_xh_det36_MuPlus_MC->SetLineColor(kRed);
  hist_xh_det36_MuPlus_MC->SetFillColor(kRed-10);
  //hist_xh_det36_MuPlus_MC->Scale(hist_xh_det36_MuPlus_Data->Integral() / hist_xh_det36_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det36_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det36_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det36_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det36_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det36_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det36_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det36_MuMinus_MC->SetMaximum(1.5 * max(max(hist_xh_det36_MuMinus_MC->GetMaximum(),hist_xh_det36_MuMinus_Data->GetMaximum()),max(hist_xh_det36_MuPlus_MC->GetMaximum(),hist_xh_det36_MuPlus_Data->GetMaximum())));
  hist_xh_det36_MuMinus_MC->Draw("hist");
  hist_xh_det36_MuPlus_MC->Draw("histsame");
  hist_xh_det36_MuMinus_Data->Draw("samepe");
  hist_xh_det36_MuPlus_Data->Draw("samepe");
  TLegend* l_det36 = new TLegend(0.79,0.49,0.98,0.97);
  l_det36->AddEntry(hist_xh_det36_MuMinus_MC,"MC #mu^{-}","f");
  l_det36->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det36_MuMinus_MC->GetEntries()),"");
  l_det36->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det36_MuMinus_MC->GetMean()),"");
  l_det36->AddEntry(hist_xh_det36_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det36->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det36_MuMinus_Data->GetEntries()),"");
  l_det36->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det36_MuMinus_Data->GetMean()),"");
  l_det36->AddEntry(hist_xh_det36_MuPlus_MC,"MC #mu^{+}","f");
  l_det36->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det36_MuPlus_MC->GetEntries()),"");
  l_det36->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det36_MuPlus_MC->GetMean()),"");
  l_det36->AddEntry(hist_xh_det36_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det36->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det36_MuPlus_Data->GetEntries()),"");
  l_det36->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det36_MuPlus_Data->GetMean()),"");
  l_det36->SetFillColor(kWhite);
  l_det36->SetLineColor(kBlack);
  l_det36->SetTextFont(43);
  l_det36->SetTextSize(14);
  l_det36->Draw();
  c_det36->Update();
  c_det36->SaveAs((plotDataMCOutputPath + "/" + c_det36->GetName() + ".png")); 

  // xh in det37
  TCanvas* c_det37 = new TCanvas("c_det37","c_det37");
  c_det37->cd();
  hist_xh_det37_MuPlus_MC->SetTitle("xh in det37");
  hist_xh_det37_MuPlus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det37_MuPlus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det37_MuPlus_MC->SetLineColor(kRed);
  hist_xh_det37_MuPlus_MC->SetFillColor(kRed-10);
  hist_xh_det37_MuPlus_MC->Scale(hist_xh_det37_MuPlus_Data->Integral() / hist_xh_det37_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det37_MuMinus_MC->SetTitle("");
  hist_xh_det37_MuMinus_MC->SetLineColor(kBlue);
  hist_xh_det37_MuMinus_MC->SetFillColor(kBlue-10);
  //hist_xh_det37_MuMinus_MC->Scale(hist_xh_det37_MuMinus_Data->Integral() / hist_xh_det37_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det37_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det37_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det37_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det37_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det37_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det37_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det37_MuPlus_MC->SetMaximum(1.5 * max(max(hist_xh_det37_MuMinus_MC->GetMaximum(),hist_xh_det37_MuMinus_Data->GetMaximum()),max(hist_xh_det37_MuPlus_MC->GetMaximum(),hist_xh_det37_MuPlus_Data->GetMaximum())));
  hist_xh_det37_MuPlus_MC->Draw("hist");
  hist_xh_det37_MuMinus_MC->Draw("histsame");
  hist_xh_det37_MuPlus_Data->Draw("samepe");
  hist_xh_det37_MuMinus_Data->Draw("samepe");
  TLegend* l_det37 = new TLegend(0.12,0.49,0.31,0.97);
  l_det37->AddEntry(hist_xh_det37_MuPlus_MC,"MC #mu^{+}","f");
  l_det37->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det37_MuPlus_MC->GetEntries()),"");
  l_det37->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det37_MuPlus_MC->GetMean()),"");
  l_det37->AddEntry(hist_xh_det37_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det37->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det37_MuPlus_Data->GetEntries()),"");
  l_det37->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det37_MuPlus_Data->GetMean()),"");
  l_det37->AddEntry(hist_xh_det37_MuMinus_MC,"MC #mu^{-}","f");
  l_det37->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det37_MuMinus_MC->GetEntries()),"");
  l_det37->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det37_MuMinus_MC->GetMean()),"");
  l_det37->AddEntry(hist_xh_det37_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det37->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det37_MuMinus_Data->GetEntries()),"");
  l_det37->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det37_MuMinus_Data->GetMean()),"");
  l_det37->SetFillColor(kWhite);
  l_det37->SetLineColor(kBlack);
  l_det37->SetTextFont(43);
  l_det37->SetTextSize(14);
  l_det37->Draw();
  c_det37->Update();
  c_det37->SaveAs((plotDataMCOutputPath + "/" + c_det37->GetName() + ".png"));     


  TCanvas* c_det6x = new TCanvas("c_det6x","c_det6x");
  c_det6x->cd();
  hist_xh_det62_MuPlus_MC->SetTitle("xh in DTs");   
  hist_xh_det62_MuPlus_MC->GetXaxis()->SetTitle("mm");
  hist_xh_det62_MuPlus_MC->GetYaxis()->SetTitle("events");
  hist_xh_det62_MuPlus_MC->SetLineColor(kRed); 
  hist_xh_det62_MuPlus_MC->SetFillColor(kRed-10);
  hist_xh_det62_MuPlus_MC->Scale(hist_xh_det62_MuPlus_Data->Integral() / hist_xh_det62_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xh_det61_MuMinus_MC->SetTitle("");	   
  hist_xh_det61_MuMinus_MC->SetLineColor(kBlue);   
  hist_xh_det61_MuMinus_MC->SetFillColor(kBlue-10);
  hist_xh_det61_MuMinus_MC->Scale(hist_xh_det61_MuMinus_Data->Integral() / hist_xh_det61_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xh_det62_MuPlus_Data->SetMarkerStyle(20);  
  hist_xh_det62_MuPlus_Data->SetMarkerColor(kRed);
  hist_xh_det62_MuPlus_Data->SetLineColor(kBlack);
  hist_xh_det61_MuMinus_Data->SetMarkerStyle(20);  
  hist_xh_det61_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xh_det61_MuMinus_Data->SetLineColor(kBlack);
  hist_xh_det62_MuPlus_MC->SetMaximum(1.5 * max(max(hist_xh_det61_MuMinus_MC->GetMaximum(),hist_xh_det61_MuMinus_Data->GetMaximum()),max(hist_xh_det62_MuPlus_MC->GetMaximum(),hist_xh_det62_MuPlus_Data->GetMaximum())));
  hist_xh_det62_MuPlus_MC->Draw("hist");
  hist_xh_det61_MuMinus_MC->Draw("samehist");
  hist_xh_det62_MuPlus_Data->Draw("samepe");
  hist_xh_det61_MuMinus_Data->Draw("samepe");
  TLegend* l_det6x = new TLegend(0.80,0.47,0.98,0.95);
  l_det6x->AddEntry(hist_xh_det62_MuPlus_MC,"MC #mu^{+}","f");
  l_det6x->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det62_MuPlus_MC->GetEntries()),"");
  l_det6x->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det62_MuPlus_MC->GetMean()),"");
  l_det6x->AddEntry(hist_xh_det62_MuPlus_Data, "Data #mu^{+}", "pl");
  l_det6x->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det62_MuPlus_Data->GetEntries()),"");
  l_det6x->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det62_MuPlus_Data->GetMean()),"");
  l_det6x->AddEntry(hist_xh_det61_MuMinus_MC,"MC #mu^{-}","f");
  l_det6x->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det61_MuMinus_MC->GetEntries()),"");
  l_det6x->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det61_MuMinus_MC->GetMean()),"");
  l_det6x->AddEntry(hist_xh_det61_MuMinus_Data, "Data #mu^{-}", "pl");
  l_det6x->AddEntry((TObject*)0,Form("entries: %.2f",hist_xh_det61_MuMinus_Data->GetEntries()),"");
  l_det6x->AddEntry((TObject*)0,Form("mean: %.2f",hist_xh_det61_MuMinus_Data->GetMean()),"");
  l_det6x->SetFillColor(kWhite);
  l_det6x->SetLineColor(kBlack);
  l_det6x->SetTextFont(43);
  l_det6x->SetTextSize(14);
  l_det6x->Draw();
  c_det6x->Update();
  TLine *line_det6xA = new TLine(261.5,gPad->GetUymin(),261.5,gPad->GetUymax());
  line_det6xA->SetLineColor(kBlack);
  line_det6xA->SetLineStyle(1);
  line_det6xA->Draw();
  // TLine *line_det6xB = new TLine(261.5 + 735.,gPad->GetUymin(),261.5 + 735.,gPad->GetUymax());
  // line_det6xB->SetLineColor(kBlack);
  // line_det6xB->SetLineStyle(2);
  // line_det6xB->Draw();
  TLine *line_det6xC = new TLine(-947.4 + 735.,gPad->GetUymin(),-947.4 + 735.,gPad->GetUymax());
  line_det6xC->SetLineColor(kBlack);
  line_det6xC->SetLineStyle(2);
  line_det6xC->Draw();
  // TLine *line_det6xD = new TLine(-947.4,gPad->GetUymin(),-947.4,gPad->GetUymax());
  // line_det6xD->SetLineColor(kBlack);
  // line_det6xD->SetLineStyle(1);
  // line_det6xD->Draw();
  c_det6x->Update();
  c_det6x->SaveAs((plotDataMCOutputPath + "/" + c_det6x->GetName() + ".png"));   


  TCanvas* c_xext = new TCanvas("c_xext","c_xext");
  c_xext->cd();
  hist_xext_MuPlus_MC->SetTitle("x ext");
  hist_xext_MuPlus_MC->GetXaxis()->SetTitle("mm");
  hist_xext_MuPlus_MC->GetYaxis()->SetTitle("events");
  hist_xext_MuPlus_MC->SetLineColor(kRed);
  hist_xext_MuPlus_MC->SetFillColor(kRed-10);
  hist_xext_MuPlus_MC->Scale(hist_xext_MuPlus_Data->Integral() / hist_xext_MuPlus_MC->Integral()); //normalize MC to Data
  hist_xext_MuMinus_MC->SetTitle("");
  hist_xext_MuMinus_MC->SetLineColor(kBlue);
  hist_xext_MuMinus_MC->SetFillColorAlpha(kBlue-10, 0.571); // color with transparency - https://root.cern.ch/doc/master/classTColor.html
  hist_xext_MuMinus_MC->Scale(hist_xext_MuMinus_Data->Integral() / hist_xext_MuMinus_MC->Integral()); //normalize MC to Data
  hist_xext_MuPlus_Data->SetMarkerStyle(20);  
  hist_xext_MuPlus_Data->SetMarkerColor(kRed);
  hist_xext_MuPlus_Data->SetLineColor(kBlack);
  hist_xext_MuMinus_Data->SetMarkerStyle(20);  
  hist_xext_MuMinus_Data->SetMarkerColor(kBlue);
  hist_xext_MuMinus_Data->SetLineColor(kBlack);
  hist_xext_MuPlus_MC->SetMaximum(1.5 * max(max(hist_xext_MuMinus_MC->GetMaximum(),hist_xext_MuMinus_Data->GetMaximum()),max(hist_xext_MuPlus_MC->GetMaximum(),hist_xext_MuPlus_Data->GetMaximum())));
  hist_xext_MuPlus_MC->Draw("hist");
  hist_xext_MuMinus_MC->Draw("histsame");
  hist_xext_MuPlus_Data->Draw("samepe");
  hist_xext_MuMinus_Data->Draw("samepe");
  TLegend* l_xext = new TLegend(0.72,0.47,0.98,0.97);
  l_xext->AddEntry(hist_xext_MuPlus_MC,"MC #mu^{+}","f");
  l_xext->AddEntry((TObject*)0,Form("entries: %.2f",hist_xext_MuPlus_MC->GetEntries()),"");
  l_xext->AddEntry((TObject*)0,Form("mean: %.2f",hist_xext_MuPlus_MC->GetMean()),"");
  l_xext->AddEntry(hist_xext_MuPlus_Data, "Data #mu^{+}", "pl");
  l_xext->AddEntry((TObject*)0,Form("entries: %.2f",hist_xext_MuPlus_Data->GetEntries()),"");
  l_xext->AddEntry((TObject*)0,Form("mean: %.2f",hist_xext_MuPlus_Data->GetMean()),"");
  l_xext->AddEntry(hist_xext_MuMinus_MC,"MC #mu^{-}","f");
  l_xext->AddEntry((TObject*)0,Form("entries: %.2f",hist_xext_MuMinus_MC->GetEntries()),"");
  l_xext->AddEntry((TObject*)0,Form("mean: %.2f",hist_xext_MuMinus_MC->GetMean()),"");
  l_xext->AddEntry(hist_xext_MuMinus_Data, "Data #mu^{-}", "pl");
  l_xext->AddEntry((TObject*)0,Form("entries: %.2f",hist_xext_MuMinus_Data->GetEntries()),"");
  l_xext->AddEntry((TObject*)0,Form("mean: %.2f",hist_xext_MuMinus_Data->GetMean()),"");
  l_xext->SetFillColor(kWhite);
  l_xext->SetLineColor(kBlack);
  l_xext->SetTextFont(43);
  l_xext->SetTextSize(20);
  l_xext->Draw();
  c_xext->Update();
  c_xext->SaveAs((plotDataMCOutputPath + "/" + c_xext->GetName() + ".png"));


  cout<<" Plots done! =) "<<endl; 


}


// main function 
void plotVariables_BeamInfo(){

  // define input files 
  TString inputFile_Data = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/reco-333to352.root";
  TString inputFile_MC   = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/reco-mupmum.root"; 

  // define output path and make output directory for data/MC comparison
  TString plotDataMCOutputPath = "181119_LemmaVariables_DataMCComparison_reco-333to352";
  gSystem->Exec(("mkdir -p "+plotDataMCOutputPath));

  // call do the histos function
  // arguments: input file, label for data or MC
  doTheHistos(inputFile_Data, "DATA");
  doTheHistos(inputFile_MC, "MC");

  // call data/MC comparison function
  dataMCComparison(plotDataMCOutputPath);


}
