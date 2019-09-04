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

Int_t nRecoEventsTot = 61;




//Extrapolate track function: track points are refitted
Double_t extrapolate_track_x(TString inputFileName, Double_t z0, Double_t x_pos_mum[12], Double_t x_pos_mum_err[12], Double_t z_x_pos_mum[12], Double_t& x_ext, Double_t& x_ext_err, Double_t& dx_on_dz_ext, Double_t& dx_on_dz_ext_err){

  //Magnetic field (box)
  Double_t zM=0.,B=0.;
  if ( inputFileName.Contains("aug18") ){
    zM=17193-846;
    B=1.7476;
  }
  if( inputFileName.Contains("sep18") ){
    zM=10.*(1573.79+0.5*(1773.79-1573.79)-82.78);
    B=2.01;
  }
  Double_t z1=zM-1000.;
  Double_t z2=zM+1000.;
  if( B==0. ) cout << "B undefined in extrapolate_track_x !" << endl;

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
  TString inputFileName = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/aug18/reco-aug18.root"; // inputFile_Data_Aug2018_Be6cm
  TString label = "reclev";
  // TString inputFileName = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/aug18/reco-mupmum.root"; // inputFile_MC_Aug2018_Be6cm
  // TString label = "MCreclev";

  double zEndTarget = 10.*(457.9+3.-84.6);   // [mm] - dataset: AUGUST 2018    Be target 6 cm
  //double zEndTarget = 10.*(460.93+3.-82.78); // [mm] - dataset: SEPTEMBER 2018 Be target 6 cm and C target 6cm
  //double zEndTarget = 10.*(460.93+1.-82.78); // [mm] - dataset: SEPTEMBER 2018 C  target 2 cm 

  // ----------------------------------


  cout << label << endl;
  bool isMC = false;                        // for Data
  if(label.Contains("MC") ){ isMC = true; } // for MC
  
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
  Double_t vtx_x;
  Double_t vtx_z;
  Double_t vtx_chi2;
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
  inputTree->SetBranchAddress("vtx_x",          &vtx_x);
  inputTree->SetBranchAddress("vtx_z",          &vtx_z);
  inputTree->SetBranchAddress("vtx_chi2",       &vtx_chi2);
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
  TH1F* hist_emittanceValue_mup = new TH1F("hist_emittanceValue_mup","Emittance Value mup; nm #times rad; entries",30,1500.,4500.);
  TH1F* hist_emittanceValue_mum = new TH1F("hist_emittanceValue_mum","Emittance Value mum; nm #times rad; entries",30,1500.,4500.);



  // -9999. -> extrapolate_track_x, to be used to process MC as if data and data
  // 0.     -> track points @ 30 and 31, FOR TEST PURPOSES
  // >0.    -> Geant4 points @ 30 and 31 smeared by sigma_x, FOR TEST PURPOSES, smeares by Gaus(1.,sigma_x/1000.)
  Double_t sigma_x=-9999.;



  // ---------------------------
  // --- get tree entries
  Long64_t entries = inputTree->GetEntries();
  cout<<" Reading input file ..."<<endl;


  // estimate is done on a reference plane 
  Double_t z_ref  = zEndTarget; // [mm] z ref (end of the target)


  // repeat 300 times the procedure
  for(Int_t j=0; j<300; j++){

    if(j%10 == 0) { cout<<" Running iteration n "<<j<<endl; }


    // def vectors for emittance estimate
    vector<double> vec_emittance_x_mup; 
    vec_emittance_x_mup.clear();
    vector<double> vec_emittance_xprime_mup;
    vec_emittance_xprime_mup.clear();

    vector<double> vec_emittance_x_mum;
    vec_emittance_x_mum.clear();
    vector<double> vec_emittance_xprime_mum;
    vec_emittance_xprime_mum.clear();


  
    // --- define Trandom variable
    TRandom3* qwerty = new TRandom3(0); 


    // --- while cycle from 0 to nRecoEventsTot forming the subsample
    Int_t nEvents = 0;
    while(nEvents < nRecoEventsTot){

      inputTree->GetEntry(entries*qwerty->Rndm());


      if( p_mup > 0. && p_mum > 0. ) {


	// --- if data or MC as if data require a valid e+/mu+/mu- vertex constraint
	// if( label.Contains("reclev") && vtx_x<-9990. ) continue; // backward compatibility studies
	if( label.Contains("reclev") ){
	  if( vtx_chi2>=500. ) continue; // 9999. 
	  if( chi2Si5MuP>500. ) continue;
	  if( chi2Si5MuM>500. ) continue;
	}


        
        // -----------------------------
        // --- emittance 


	// --- mu+
        Double_t x_det30_atZref_mup=0,pTot_genLev_mup=0,x_prime_ondet30_mup=0;
        if( label == "MC" ){
          // x  (extrapolation on reference plane): x_det30_atZref_mup = x_onDet30 - (z_onDet30 - z_ref)* px_onDet30 / pz_onDet30
          x_det30_atZref_mup = gen_pos_mup[0] - (gen_pos_mup[2] - z_ref)*(gen_pos_mup[3]/gen_pos_mup[5]); 
          // x' (extrapolation on reference plane): the direction remain the same as on det30 or at vtx
          pTot_genLev_mup = sqrt(gen_pos_mup[3]*gen_pos_mup[3] + gen_pos_mup[4]*gen_pos_mup[4] + gen_pos_mup[5]*gen_pos_mup[5]);
          x_prime_ondet30_mup = gen_pos_mup[3] / pTot_genLev_mup;
	}
        if( label.Contains("reclev") ){
          // use full track and extrapolate_track_x
          Double_t chi2_mup=9999;
	  Double_t x_ext_mup=-9999,x_ext_err_mup=-9999;
	  Double_t dx_on_dz_ext_mup=-9999,dx_on_dz_ext_err_mup=-9999;
          chi2_mup = extrapolate_track_x(inputFileName,z_ref, x_pos_mup, x_pos_mup_err, z_x_pos_mup, 
					 x_ext_mup, x_ext_err_mup, dx_on_dz_ext_mup, dx_on_dz_ext_err_mup);
          x_det30_atZref_mup = x_ext_mup;
          pTot_genLev_mup = p_mup;
          x_prime_ondet30_mup = dx_on_dz_ext_mup;
          // use measured track points on 30 and 31 for x_det30_atZref_mum
          if( sigma_x==0 ){
            Double_t DeltaX=(x_pos_mup[3]-x_pos_mup[4]); Double_t DeltaZ=(z_x_pos_mup[3]-z_x_pos_mup[4]);
            x_det30_atZref_mup = (x_pos_mup[4] + (z_ref-z_x_pos_mup[4])*DeltaX/DeltaZ);
	  }
          // use smeared Geant4 points on 30 and 31 for x_det30_atZref_mum
          if( sigma_x>0 ){
            TRandom3* prandom = new TRandom3(0);
            Double_t r30=prandom->Gaus(0.,sigma_x/1000.);
            Double_t r31=prandom->Gaus(0.,sigma_x/1000.);
            // cout << r30 << " " << r31 << endl;
            r31+=gen_pos_mup[6]; r30+=gen_pos_mup[0];
            Double_t DeltaX=(r31-r30); Double_t DeltaZ=(gen_pos_mup[8]-gen_pos_mup[2]);
            x_det30_atZref_mup = (r30 + (z_ref-gen_pos_mup[2])*DeltaX/DeltaZ);
            delete prandom;
	  }
	}
        
        // emittance of mu+
        Double_t x_emittance_mup       = x_det30_atZref_mup;
        Double_t x_prime_emittance_mup = x_prime_ondet30_mup;
        // --- fill vectors
        vec_emittance_x_mup     .push_back(x_emittance_mup);
        vec_emittance_xprime_mup.push_back(x_prime_emittance_mup);



	// --- mu-
        Double_t x_det30_atZref_mum=0,pTot_genLev_mum=0,x_prime_ondet30_mum=0;
        if( label == "MC" ){
          // x  (extrapolation on reference plane): x_det30_atZref_mum = x_onDet30 - (z_onDet30 - z_ref)* px_onDet30 / pz_onDet30
          x_det30_atZref_mum = gen_pos_mum[0] - (gen_pos_mum[2] - z_ref)*(gen_pos_mum[3]/gen_pos_mum[5]); 
          // x' (extrapolation on reference plane): the direction remain the same as on det30 or at vtx
          pTot_genLev_mum = sqrt(gen_pos_mum[3]*gen_pos_mum[3] + gen_pos_mum[4]*gen_pos_mum[4] + gen_pos_mum[5]*gen_pos_mum[5]);
          x_prime_ondet30_mum = gen_pos_mum[3] / pTot_genLev_mum; 
	}
        if( label.Contains("reclev") ){
          // use full track and extrapolate_track_x
          Double_t chi2_mum=9999;
          Double_t x_ext_mum=-9999,x_ext_err_mum=-9999;
          Double_t dx_on_dz_ext_mum=-9999,dx_on_dz_ext_err_mum=-9999;
          chi2_mum = extrapolate_track_x(inputFileName,z_ref, x_pos_mum, x_pos_mum_err, z_x_pos_mum,
                                         x_ext_mum, x_ext_err_mum, dx_on_dz_ext_mum, dx_on_dz_ext_err_mum);
          x_det30_atZref_mum = x_ext_mum;
          pTot_genLev_mum = p_mum;
          x_prime_ondet30_mum = dx_on_dz_ext_mum;
          // use measured track points on 30 and 31 for x_det30_atZref_mum
          if( sigma_x==0 ){
            Double_t DeltaX=(x_pos_mum[3]-x_pos_mum[4]); Double_t DeltaZ=(z_x_pos_mum[3]-z_x_pos_mum[4]);
            x_det30_atZref_mum = (x_pos_mum[4] + (z_ref-z_x_pos_mum[4])*DeltaX/DeltaZ);
	  }
          // use smeared Geant4 points on 30 and 31 for x_det30_atZref_mum
          if( sigma_x>0 ){
            TRandom3* prandom = new TRandom3(0);
            Double_t r30=prandom->Gaus(0.,sigma_x/1000.);
            Double_t r31=prandom->Gaus(0.,sigma_x/1000.);
            // cout << r30 << " " << r31 << endl;
            r31+=gen_pos_mum[6]; r30+=gen_pos_mum[0]; 
            Double_t DeltaX=(r31-r30); Double_t DeltaZ=(gen_pos_mum[8]-gen_pos_mum[2]);
            x_det30_atZref_mum = (r30 + (z_ref-gen_pos_mum[2])*DeltaX/DeltaZ);
            delete prandom;
	  }
	}


        // emittance of mu-
        Double_t x_emittance_mum       = x_det30_atZref_mum;
        Double_t x_prime_emittance_mum = x_prime_ondet30_mum;
        // --- fill vectors
        vec_emittance_x_mum     .push_back(x_emittance_mum);
        vec_emittance_xprime_mum.push_back(x_prime_emittance_mum);


	nEvents++;

      } // end if (p_mup > 0. && p_mum > 0.)

    }//end for while cycle form 0 to nRecoEventsTot forming the subsample

  
  
    // ---------------------------------
    //  compute emittance from vectors
    // ---------------------------------

    hist_emittanceValue_mup->Fill(getemittance(vec_emittance_x_mup, vec_emittance_xprime_mup)); 
    hist_emittanceValue_mum->Fill(getemittance(vec_emittance_x_mum, vec_emittance_xprime_mum)); 

    cout<<"muP "<<getemittance(vec_emittance_x_mup, vec_emittance_xprime_mup)<<endl;
    cout<<"muM "<<getemittance(vec_emittance_x_mum, vec_emittance_xprime_mum)<<endl;


  }//end for cycle from 0 to 300 repeating the procedure 300 times


  
  TCanvas* c_emittanceValue_mup = new TCanvas("c_emittanceValue_mup","c_emittanceValue_mup");
  c_emittanceValue_mup->cd();
  hist_emittanceValue_mup->Draw();
  if(isMC){
    c_emittanceValue_mup->SaveAs("emittanceUnc_mup_MC.png");
  }else{
    c_emittanceValue_mup->SaveAs("emittanceUnc_mup.png");
  }


  TCanvas* c_emittanceValue_mum = new TCanvas("c_emittanceValue_mum","c_emittanceValue_mum");
  c_emittanceValue_mum->cd();
  hist_emittanceValue_mum->Draw();
  if(isMC){
    c_emittanceValue_mum->SaveAs("emittanceUnc_mum_MC.png");
  }else{
    c_emittanceValue_mum->SaveAs("emittanceUnc_mum.png");
  }


}
