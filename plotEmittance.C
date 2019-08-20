// *****************************************************
//
// author: Alessandra Cappati
//         14/03/2019
// 
// usage: specify the input files (Data and MC), 
//        the output directory, and other options 
//        at the end of the script
//
// run with:
//        root -l -b -q plotEmittance.C++
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





// doTheHistos function: read root file and do histos 
void doTheHistos(TString inputFileName, TString label, double zEndTarget, double zPosDet20, TString plotOutputPath){

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

  inputTree->SetBranchAddress("chi2Si5MuM",	&chi2Si5MuM);
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


  // def vectors for emittance estimate

  // DATA vectors  
  vector<double> vec_PositronBeamEmittance_x_1eplus_Data;
  vector<double> vec_PositronBeamEmittance_xprime_1eplus_Data;
  vector<double> vec_PositronBeamEmittance_x_moreThan1eplus_Data;
  vector<double> vec_PositronBeamEmittance_xprime_moreThan1eplus_Data;

  // MC vectors
  vector<double> vec_emittanceControl_emittance_x_positron_MC;
  vector<double> vec_emittanceControl_emittance_xprime_positron_MC;

  vector<double> vec_emittance_x_mup_MC;
  vector<double> vec_emittance_xprime_mup_MC;
  vector<double> vec_emittanceControl_emittance_x_mup_MC;
  vector<double> vec_emittanceControl_emittance_xprime_mup_MC;

  vector<double> vec_emittance_x_mum_MC;
  vector<double> vec_emittance_xprime_mum_MC;
  vector<double> vec_emittanceControl_emittance_x_mum_MC;
  vector<double> vec_emittanceControl_emittance_xprime_mum_MC;



  // def histos limits

  int    h_n_bins_rawEmitt     =  100;
  double h_min_x_rawEmitt      = -30.;   // [mm]
  double h_max_x_rawEmitt      =  30.;   // [mm]
  double h_min_xprime_rawEmitt = -0.002; // [rad]
  double h_max_xprime_rawEmitt =  0.002; // [rad]
  if( !isMC ){
    h_n_bins_rawEmitt = 25;
    h_min_x_rawEmitt*=2.; 
    h_max_x_rawEmitt*=2.;
    h_min_xprime_rawEmitt*=2.; 
    h_max_xprime_rawEmitt*=2.;
  }


  double h_min_x_emitt      = -0.3;   // [mm]
  double h_max_x_emitt      =  0.3;   // [mm]
  double h_min_xprime_emitt = -0.002; // [rad]
  double h_max_xprime_emitt =  0.002; // [rad]



  // def histos 
 
  // --- DATA HISTOS
  TH1F* hist_npos_Data = new TH1F("hist_npos_Data", "N positrons", 11, -0.5, 10.5);                                           
  TH1F* hist_xbe_positrons_Data = new TH1F("hist_xbe_positrons_Data", "Positron: Be exit point (mm)", h_n_bins_rawEmitt, h_min_x_rawEmitt,      h_max_x_rawEmitt);      
  TH1F* hist_the_positrons_Data = new TH1F("hist_the_positrons_Data", "Positron: theta exit (rad)",   h_n_bins_rawEmitt, h_min_xprime_rawEmitt, h_max_xprime_rawEmitt);    

  // on det20
  TH1F* hist_npos_Data_det20          = new TH1F("hist_npos_Data_det20", "N positrons on det20", 11, -0.5, 10.5);                                           
  TH1F* hist_xbe_positrons_Data_det20 = new TH1F("hist_xbe_positrons_Data_det20", "Positron: Be exit point on det20 (mm)", 500, -30., 30.);      
  TH1F* hist_the_positrons_Data_det20 = new TH1F("hist_the_positrons_Data_det20", "Positron: theta exit on det20 (rad)", 500, -2e-3, 2e-3);    


  // --- raw emittance
  TH1F* hist1D_PositronBeamEmittance_x______1eplus_Data = new TH1F("hist1D_PositronBeamEmittance_x______1eplus_Data","hist1D_PositronBeamEmittance_x______1eplus_Data",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt); 
  TH1F* hist1D_PositronBeamEmittance_xprime_1eplus_Data = new TH1F("hist1D_PositronBeamEmittance_xprime_1eplus_Data","hist1D_PositronBeamEmittance_xprime_1eplus_Data",h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
  TH2F* hist2D_PositronBeamEmittance_emitt__1eplus_Data = new TH2F("hist2D_PositronBeamEmittance_emitt__1eplus_Data","hist2D_PositronBeamEmittance_emitt__1eplus_Data",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt,h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
 
  TH1F* hist1D_PositronBeamEmittance_x______moreThan1eplus_Data = new TH1F("hist1D_PositronBeamEmittance_x______moreThan1eplus_Data","hist1D_PositronBeamEmittance_x______moreThan1eplus_Data",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt);
  TH1F* hist1D_PositronBeamEmittance_xprime_moreThan1eplus_Data = new TH1F("hist1D_PositronBeamEmittance_xprime_moreThan1eplus_Data","hist1D_PositronBeamEmittance_xprime_moreThan1eplus_Data",h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
  TH2F* hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data = new TH2F("hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data","hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt,100,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);

  

  // --- MC HISTOS
  //2D emittance plots
  TH2D* hist2D_emittance_x_mup_MC = new TH2D("hist2D_emittance_x_mup_MC","hist2D_emittance_x_mup_MC",100, h_min_x_emitt, h_max_x_emitt, 100, h_min_xprime_emitt, h_max_xprime_emitt);
  TH2D* hist2D_emittance_x_mum_MC = new TH2D("hist2D_emittance_x_mum_MC","hist2D_emittance_x_mum_MC",100, h_min_x_emitt, h_max_x_emitt, 100, h_min_xprime_emitt, h_max_xprime_emitt);
  //1D emittance plots
  TH1D* hist1D_emittance_x_mup_MC       = new TH1D("hist1D_emittance_x_mup_MC",      "hist1D_emittance_x_mup_MC",      100, h_min_x_emitt,      h_max_x_emitt);
  TH1D* hist1D_emittance_x_prime_mup_MC = new TH1D("hist1D_emittance_x_prime_mup_MC","hist1D_emittance_x_prime_mup_MC",100, h_min_xprime_emitt, h_max_xprime_emitt);
  TH1D* hist1D_emittance_x_mum_MC       = new TH1D("hist1D_emittance_x_mum_MC",      "hist1D_emittance_x_mum_MC",      100, h_min_x_emitt,      h_max_x_emitt);
  TH1D* hist1D_emittance_x_prime_mum_MC = new TH1D("hist1D_emittance_x_prime_mum_MC","hist1D_emittance_x_prime_mum_MC",100, h_min_xprime_emitt, h_max_xprime_emitt);

  // --- raw emittance
  // 1D emittance control plots 
  TH1D* hist1D_emittanceControl_x_atZref_eplus_MC       = new TH1D("hist1D_emittanceControl_x_atZref_eplus_MC","hist1D_emittanceControl_x_atZref_eplus_MC",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt);
  TH1D* hist1D_emittanceControl_x_prime_atZref_eplus_MC = new TH1D("hist1D_emittanceControl_x_prime_atZref_eplus_MC","hist1D_emittanceControl_x_prime_atZref_eplus_MC",h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);

  TH1D* hist1D_emittanceControl_x_onDet30_mup_MC = new TH1D("hist1D_emittanceControl_x_onDet30_mup_MC","hist1D_emittanceControl_x_onDet30_mup_MC",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt);
  TH1D* hist1D_emittanceControl_x_atZref_mup_MC  = new TH1D("hist1D_emittanceControl_x_atZref_mup_MC" ,"hist1D_emittanceControl_x_atZref_mup_MC" ,h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt);
  TH1D* hist1D_emittanceControl_xprime_mup_MC    = new TH1D("hist1D_emittanceControl_xprime_mup_MC"   ,"hist1D_emittanceControl_xprime_mup_MC"   ,h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
  TH1D* hist1D_emittanceControl_px_mup_MC        = new TH1D("hist1D_emittanceControl_px_mup_MC"       ,"hist1D_emittanceControl_px_mup_MC"       ,100,-60.,60.);
  TH1D* hist1D_emittanceControl_py_mup_MC        = new TH1D("hist1D_emittanceControl_py_mup_MC"       ,"hist1D_emittanceControl_py_mup_MC"       ,100,-60.,60.);
  TH1D* hist1D_emittanceControl_pz_mup_MC        = new TH1D("hist1D_emittanceControl_pz_mup_MC"       ,"hist1D_emittanceControl_pz_mup_MC"       ,100,10000.,35000.);
  TH1D* hist1D_emittanceControl_pTot_mup_MC      = new TH1D("hist1D_emittanceControl_pTot_mup_MC"     ,"hist1D_emittanceControl_pTot_mup_MC"     ,100,10000.,35000.);

  TH1D* hist1D_emittanceControl_x_onDet30_mum_MC = new TH1D("hist1D_emittanceControl_x_onDet30_mum_MC","hist1D_emittanceControl_x_onDet30_mum_MC",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt);
  TH1D* hist1D_emittanceControl_x_atZref_mum_MC  = new TH1D("hist1D_emittanceControl_x_atZref_mum_MC" ,"hist1D_emittanceControl_x_atZref_mum_MC" ,h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt);
  TH1D* hist1D_emittanceControl_xprime_mum_MC    = new TH1D("hist1D_emittanceControl_xprime_mum_MC"   ,"hist1D_emittanceControl_xprime_mum_MC"   ,h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
  TH1D* hist1D_emittanceControl_px_mum_MC        = new TH1D("hist1D_emittanceControl_px_mum_MC"       ,"hist1D_emittanceControl_px_mum_MC"       ,100,-60.,60.);
  TH1D* hist1D_emittanceControl_py_mum_MC        = new TH1D("hist1D_emittanceControl_py_mum_MC"       ,"hist1D_emittanceControl_py_mum_MC"       ,100,-60.,60.);
  TH1D* hist1D_emittanceControl_pz_mum_MC        = new TH1D("hist1D_emittanceControl_pz_mum_MC"       ,"hist1D_emittanceControl_pz_mum_MC"       ,100,10000.,35000.);
  TH1D* hist1D_emittanceControl_pTot_mum_MC      = new TH1D("hist1D_emittanceControl_pTot_mum_MC"     ,"hist1D_emittanceControl_pTot_mum_MC"     ,100,10000.,35000.);
  // 2D emittance control plots
  TH2D* hist2D_emittanceControl_emittance_positron_MC = new TH2D("hist2D_emittanceControl_emittance_positron_MC","hist2D_emittanceControl_emittance_positron_MC",h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt,h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
  TH2D* hist2D_emittanceControl_emittance_mup_MC      = new TH2D("hist2D_emittanceControl_emittance_mup_MC"     ,"hist2D_emittanceControl_emittance_mup_MC"     ,h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt,h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);
  TH2D* hist2D_emittanceControl_emittance_mum_MC      = new TH2D("hist2D_emittanceControl_emittance_mum_MC"     ,"hist2D_emittanceControl_emittance_mum_MC"     ,h_n_bins_rawEmitt,h_min_x_rawEmitt,h_max_x_rawEmitt,h_n_bins_rawEmitt,h_min_xprime_rawEmitt,h_max_xprime_rawEmitt);

  // -9999. -> extrapolate_track_x, to be used to process MC as if data and data
  // 0.     -> track points @ 30 and 31, FOR TEST PURPOSES
  // >0.    -> Geant4 points @ 30 and 31 smeared by sigma_x, FOR TEST PURPOSES, smeares by Gaus(1.,sigma_x/1000.)
  Double_t sigma_x=-9999.;

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
      hist_npos_Data->Fill(nxbe);

      for(Int_t j=0; j<nxbe; j++){
        hist_xbe_positrons_Data->Fill(xpatz[j]); //position [mm]
        hist_the_positrons_Data->Fill(thatz[j]); //angle [rad]
  
        if(j==1){
          hist1D_PositronBeamEmittance_x______1eplus_Data->Fill(xpatz[j]);
          hist1D_PositronBeamEmittance_xprime_1eplus_Data->Fill(thatz[j]);
          hist2D_PositronBeamEmittance_emitt__1eplus_Data->Fill(xpatz[j],thatz[j]);

          vec_PositronBeamEmittance_x_1eplus_Data     .push_back(xpatz[j]);
	  vec_PositronBeamEmittance_xprime_1eplus_Data.push_back(thatz[j]);

        }else if(j>1){

          hist1D_PositronBeamEmittance_x______moreThan1eplus_Data->Fill(xpatz[j]);
          hist1D_PositronBeamEmittance_xprime_moreThan1eplus_Data->Fill(thatz[j]);
          hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->Fill(xpatz[j],thatz[j]);

          vec_PositronBeamEmittance_x_moreThan1eplus_Data     .push_back(xpatz[j]);
          vec_PositronBeamEmittance_xprime_moreThan1eplus_Data.push_back(thatz[j]);
        }
      }

      // --- if z_pos is on det20
      Double_t xpatz_det20[10], thatz_det20[10];
      Double_t zpos_det20(zPosDet20); //Z position of det20
      Int_t nxbe_det20 = ReturnBeamInfo(subdet, xh, zh, nhits, zpos_det20, xpatz_det20, thatz_det20);  // number of positrons
      hist_npos_Data_det20->Fill(nxbe_det20);

      for(Int_t j=0; j<nxbe_det20; j++){
        hist_xbe_positrons_Data_det20->Fill(xpatz_det20[j]); //position [mm]
        hist_the_positrons_Data_det20->Fill(thatz_det20[j]); //angle [rad]
      }


    }//end if !isMC
    // ----------------------------------------------------------------------


    
    // --- condition for candidate events
    if( p_mup > 0. && p_mum > 0. ) {

      // --- if data or MC as if data require a valid e+/mu+/mu- vertex constraint
      // if( label.Contains("reclev") && vtx_x<-9990. ) continue; // backward compatibility studies
      if( label.Contains("reclev") ){
	if( vtx_chi2>=500. ) continue; // 9999. 
        if( chi2Si5MuP>500. ) continue;
        if( chi2Si5MuM>500. ) continue;
      }

      // -----------------------------
      // --- emittance in x 1D and 2D histos
      // -----------------------------
      // x  = x  @det30 - x  incoming e+
      // x' = x' @det30 - x' incoming e+
      // estimate is done on a reference plane                                                                                          
      Double_t z_ref  = zEndTarget; // [mm] z ref (end of the target)

      Double_t x_atZref_eplus=0.;
      Double_t x_prime_atZref_eplus=0.;
      if(isMC){
        
        // --- e+ incoming
        // px of e+ = Cx mu- * En of mu- + Cx mu+ * En of mu+ = px of mu- + px of mu+
        Double_t px_eplus = gen_vtx_mum[3]*gen_vtx_mum[6] + gen_vtx_mup[3]*gen_vtx_mup[6];
        Double_t py_eplus = gen_vtx_mum[4]*gen_vtx_mum[6] + gen_vtx_mup[4]*gen_vtx_mup[6];
        Double_t pz_eplus = gen_vtx_mum[5]*gen_vtx_mum[6] + gen_vtx_mup[5]*gen_vtx_mup[6];
        // x  of e+ (extrapolation on reference plane)  = x_vtx - (z_vtx - z_ref)*(px_e+/pz_e+)
        // N.B. gen_vtx_mup[0] = gen_vtx_mum[0] and gen_vtx_mup[2] = gen_vtx_mum[2] 
        x_atZref_eplus = gen_vtx_mup[0] - (gen_vtx_mup[2] - z_ref)*(px_eplus/pz_eplus);
        // x' of e+ (extrapolation on reference plane)  = px / p  (the direction remain the same as at vtx)
        x_prime_atZref_eplus = px_eplus / sqrt(px_eplus*px_eplus + py_eplus*py_eplus + pz_eplus*pz_eplus);
      }else{

        // --- e+ incoming
        // use the positron giving the best vertex chi2 fit, these hits are flagged
        Double_t x_10=0.; Double_t z_10=0.; Double_t x_20=0.; Double_t z_20=0.;
        for(Int_t i=0;i<nhits;i++){
	  if( subdet[i]==10 && itrack[i]==-11 ){
            if( xh[i]>-9990. ) x_10=xh[i];
            if( zh[i]>-9990. ) z_10=zh[i];
	  }
          if( subdet[i]==20 && itrack[i]==-11 ){
            if( xh[i]>-9990. ) x_20=xh[i];
            if( zh[i]>-9990. ) z_20=zh[i];
          }
	}
        x_prime_atZref_eplus = ((x_20-x_10)/(z_20-z_10)); // DeltaX/DeltaZ
        x_atZref_eplus = (x_20 + (z_ref-z_20)*x_prime_atZref_eplus);
      }

      // emittance e+ control plots
      hist1D_emittanceControl_x_atZref_eplus_MC->Fill(x_atZref_eplus);
      hist1D_emittanceControl_x_prime_atZref_eplus_MC->Fill(x_prime_atZref_eplus);
      hist2D_emittanceControl_emittance_positron_MC->Fill(x_atZref_eplus,x_prime_atZref_eplus);

      vec_emittanceControl_emittance_x_positron_MC     .push_back(x_atZref_eplus);
      vec_emittanceControl_emittance_xprime_positron_MC.push_back(x_prime_atZref_eplus);
 
      if( true ){

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

        // emittance mu+ control plots
        hist1D_emittanceControl_x_onDet30_mup_MC->Fill(gen_pos_mup[0]);
        hist1D_emittanceControl_x_atZref_mup_MC->Fill(x_det30_atZref_mup);
        hist1D_emittanceControl_xprime_mup_MC->Fill(x_prime_ondet30_mup);
        if( label == "MC" ){
          hist1D_emittanceControl_px_mup_MC->Fill(gen_pos_mup[3]);
          hist1D_emittanceControl_py_mup_MC->Fill(gen_pos_mup[4]);
          hist1D_emittanceControl_pz_mup_MC->Fill(gen_pos_mup[5]);
	}
        hist1D_emittanceControl_pTot_mup_MC->Fill(pTot_genLev_mup);
	hist2D_emittanceControl_emittance_mup_MC->Fill(x_det30_atZref_mup,x_prime_ondet30_mup); 
 
        vec_emittanceControl_emittance_x_mup_MC     .push_back(x_det30_atZref_mup);
        vec_emittanceControl_emittance_xprime_mup_MC.push_back(x_prime_ondet30_mup);
        
        // emittance of mu+
        Double_t x_emittance_mup = x_det30_atZref_mup - x_atZref_eplus;
        Double_t x_prime_emittance_mup = x_prime_ondet30_mup - x_prime_atZref_eplus;
        hist2D_emittance_x_mup_MC->Fill(x_emittance_mup, x_prime_emittance_mup);
        hist1D_emittance_x_mup_MC->Fill(x_emittance_mup);
        hist1D_emittance_x_prime_mup_MC->Fill(x_prime_emittance_mup);

        vec_emittance_x_mup_MC     .push_back(x_emittance_mup);
        vec_emittance_xprime_mup_MC.push_back(x_prime_emittance_mup);



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

        // emittance mu- control plots
        hist1D_emittanceControl_x_onDet30_mum_MC->Fill(gen_pos_mum[0]);
        hist1D_emittanceControl_x_atZref_mum_MC->Fill(x_det30_atZref_mum);
        hist1D_emittanceControl_xprime_mum_MC->Fill(x_prime_ondet30_mum);
        if( label == "MC" ){
          hist1D_emittanceControl_px_mum_MC->Fill(gen_pos_mum[3]);
          hist1D_emittanceControl_py_mum_MC->Fill(gen_pos_mum[4]);
          hist1D_emittanceControl_pz_mum_MC->Fill(gen_pos_mum[5]);
	}
        hist1D_emittanceControl_pTot_mum_MC->Fill(pTot_genLev_mum);
       	hist2D_emittanceControl_emittance_mum_MC->Fill(x_det30_atZref_mum,x_prime_ondet30_mum);

        vec_emittanceControl_emittance_x_mum_MC     .push_back(x_det30_atZref_mum);
        vec_emittanceControl_emittance_xprime_mum_MC.push_back(x_prime_ondet30_mum);


        // emittance of mu-
        Double_t x_emittance_mum = x_det30_atZref_mum - x_atZref_eplus;
        Double_t x_prime_emittance_mum = x_prime_ondet30_mum - x_prime_atZref_eplus;
        hist2D_emittance_x_mum_MC->Fill(x_emittance_mum, x_prime_emittance_mum);
        hist1D_emittance_x_mum_MC->Fill(x_emittance_mum);
        hist1D_emittance_x_prime_mum_MC->Fill(x_prime_emittance_mum);

        vec_emittance_x_mum_MC     .push_back(x_emittance_mum);
        vec_emittance_xprime_mum_MC.push_back(x_prime_emittance_mum);


        //cout<<x_prime_atZref_eplus<<"    "<<x_prime_ondet30_mup<<"    "<<x_prime_ondet30_mum<<endl;

      } // true

    } // end if (p_mup > 0. && p_mum > 0.)

  }//end over tree entries 



  // ------------------
  //  plot histos
  // ------------------

  gStyle->SetOptStat(1);
 
  

  // beam info 
  // DATA ONLY
  if(!isMC){
  
    TCanvas * c_pos = new TCanvas("c_pos","c_pos",800, 1200);
    c_pos->Divide(1,3);
    c_pos->cd(1); hist_npos_Data->Draw();
    c_pos->cd(2); hist_xbe_positrons_Data->Draw();
    c_pos->cd(3); hist_the_positrons_Data->Draw();
    c_pos->SaveAs((plotOutputPath + "/" + c_pos->GetName() + ".png"));

    // on det20
    TCanvas * c_pos_det20 = new TCanvas("c_pos_det20","c_pos_det20",800, 1200);
    c_pos_det20->Divide(1,3);
    c_pos_det20->cd(1); hist_npos_Data_det20->Draw();
    c_pos_det20->cd(2); 
    hist_xbe_positrons_Data_det20->Draw();
    c_pos_det20->cd(3);
    gStyle->SetOptStat(1); 
    TF1* g1 = new TF1("g1","gaus", -0.0004, 0.0004);
    g1->SetParameter(0,4500.);
    g1->SetParameter(1,0.);
    g1->SetParameter(2,0.0004);
    g1->SetLineColor(kRed);
    hist_the_positrons_Data_det20->Fit("g1","R");
    hist_the_positrons_Data_det20->Draw();
    g1->Draw("samel");
    c_pos_det20->SaveAs((plotOutputPath + "/" + c_pos_det20->GetName() + ".png"));
    gStyle->SetOptStat(0);
    
    TCanvas* c_PositronBeamEmittance_x_1eplus = new TCanvas("c_PositronBeamEmittance_x_1eplus","c_PositronBeamEmittance_x_1eplus");
    c_PositronBeamEmittance_x_1eplus->cd();
    hist1D_PositronBeamEmittance_x______1eplus_Data->SetLineColor(kViolet);
    hist1D_PositronBeamEmittance_x______1eplus_Data->Draw();
    c_PositronBeamEmittance_x_1eplus->SaveAs((plotOutputPath + "/" + c_PositronBeamEmittance_x_1eplus->GetName() + ".png"));
    
    TCanvas* c_PositronBeamEmittance_xprime_1eplus = new TCanvas("c_PositronBeamEmittance_xprime_1eplus","c_PositronBeamEmittance_xprime_1eplus");
    c_PositronBeamEmittance_xprime_1eplus->cd();
    hist1D_PositronBeamEmittance_xprime_1eplus_Data->SetLineColor(kViolet);
    hist1D_PositronBeamEmittance_xprime_1eplus_Data->Draw();
    c_PositronBeamEmittance_xprime_1eplus->SaveAs((plotOutputPath + "/" + c_PositronBeamEmittance_xprime_1eplus->GetName() + ".png"));
    
    TCanvas* c_PositronBeamEmittance_emitt_1eplus = new TCanvas("c_PositronBeamEmittance_emitt_1eplus","c_PositronBeamEmittance_emitt_1eplus");
    c_PositronBeamEmittance_emitt_1eplus->cd();
    hist2D_PositronBeamEmittance_emitt__1eplus_Data->SetTitle("emittance: incoming beam with 1 e^{+}");
    hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetXaxis()->SetTitle("x [mm]");
    hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetYaxis()->SetTitle("x' [rad]");
    hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kCool);
    hist2D_PositronBeamEmittance_emitt__1eplus_Data->Draw("COLZ");
    // --- Binned formula
    //float emittanceValue_positronBeam_1eplus_Data = sqrt(hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetCovariance(1,1)*hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetCovariance(2,2) - hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetCovariance(2,1)*hist2D_PositronBeamEmittance_emitt__1eplus_Data->GetCovariance(1,2)); 
    // --- Unbinned formula
    float emittanceValue_positronBeam_1eplus_Data = getemittance(vec_PositronBeamEmittance_x_1eplus_Data, vec_PositronBeamEmittance_xprime_1eplus_Data); 
    TPaveText* pv_positronBeam_1eplus = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_positronBeam_1eplus->AddText(Form("#epsilon = %.2f",emittanceValue_positronBeam_1eplus_Data));
    pv_positronBeam_1eplus->AddText("nm #times rad");
    pv_positronBeam_1eplus->SetFillColor(kWhite);
    pv_positronBeam_1eplus->SetBorderSize(0);
    pv_positronBeam_1eplus->SetTextFont(40);
    pv_positronBeam_1eplus->SetTextSize(0.05);
    pv_positronBeam_1eplus->SetTextFont(42);
    pv_positronBeam_1eplus->SetTextAlign(22); //text centering
    pv_positronBeam_1eplus->Draw();
    c_PositronBeamEmittance_emitt_1eplus->Update();
    c_PositronBeamEmittance_emitt_1eplus->SaveAs((plotOutputPath + "/" + c_PositronBeamEmittance_emitt_1eplus->GetName() + ".png"));
    
    
    TCanvas* c_PositronBeamEmittance_x_moreThan1eplus = new TCanvas("c_PositronBeamEmittance_x_moreThan1eplus","c_PositronBeamEmittance_x_moreThan1eplus");
    c_PositronBeamEmittance_x_moreThan1eplus->cd();
    hist1D_PositronBeamEmittance_x______moreThan1eplus_Data->SetLineColor(kViolet);
    hist1D_PositronBeamEmittance_x______moreThan1eplus_Data->Draw();
    c_PositronBeamEmittance_x_moreThan1eplus->SaveAs((plotOutputPath + "/" + c_PositronBeamEmittance_x_moreThan1eplus->GetName() + ".png"));
    
    TCanvas* c_PositronBeamEmittance_xprime_moreThan1eplus = new TCanvas("c_PositronBeamEmittance_xprime_moreThan1eplus","c_PositronBeamEmittance_xprime_moreThan1eplus");
    c_PositronBeamEmittance_xprime_moreThan1eplus->cd();
    hist1D_PositronBeamEmittance_xprime_moreThan1eplus_Data->SetLineColor(kViolet);
    hist1D_PositronBeamEmittance_xprime_moreThan1eplus_Data->Draw();
    c_PositronBeamEmittance_xprime_moreThan1eplus->SaveAs((plotOutputPath + "/" + c_PositronBeamEmittance_xprime_moreThan1eplus->GetName() + ".png"));
    
    TCanvas* c_PositronBeamEmittance_emitt_moreThan1eplus = new TCanvas("c_PositronBeamEmittance_emitt_moreThan1eplus","c_PositronBeamEmittance_emitt_moreThan1eplus");
    c_PositronBeamEmittance_emitt_moreThan1eplus->cd();
    hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->SetTitle("emittance: incoming beam with more than 1 e^{+}");
    hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetXaxis()->SetTitle("x [mm]");
    hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetYaxis()->SetTitle("x' [rad]");
    hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kCool);
    hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->Draw("COLZ");
    // --- Binned formula
    // float emittanceValue_positronBeam_moreThan1eplus_Data = sqrt(hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetCovariance(1,1)*hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetCovariance(2,2) - hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetCovariance(2,1)*hist2D_PositronBeamEmittance_emitt__moreThan1eplus_Data->GetCovariance(1,2));
    float emittanceValue_positronBeam_moreThan1eplus_Data = getemittance(vec_PositronBeamEmittance_x_moreThan1eplus_Data, vec_PositronBeamEmittance_xprime_moreThan1eplus_Data);
    TPaveText* pv_positronBeam_moreThan1eplus = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_positronBeam_moreThan1eplus->AddText(Form("#epsilon = %.2f",emittanceValue_positronBeam_moreThan1eplus_Data));
    pv_positronBeam_moreThan1eplus->AddText("nm #times rad");
    pv_positronBeam_moreThan1eplus->SetFillColor(kWhite);
    pv_positronBeam_moreThan1eplus->SetBorderSize(0);
    pv_positronBeam_moreThan1eplus->SetTextFont(40);
    pv_positronBeam_moreThan1eplus->SetTextSize(0.05);
    pv_positronBeam_moreThan1eplus->SetTextFont(42);
    pv_positronBeam_moreThan1eplus->SetTextAlign(22); //text centering
    pv_positronBeam_moreThan1eplus->Draw();
    c_PositronBeamEmittance_emitt_moreThan1eplus->Update();
    c_PositronBeamEmittance_emitt_moreThan1eplus->SaveAs((plotOutputPath + "/" + c_PositronBeamEmittance_emitt_moreThan1eplus->GetName() + ".png"));

  }
  // end beam info (DATA ONLY)


  // MC HISTOS
  if( true ){

    // 2D emittance plots
    TCanvas* c_x_emittance_mup = new TCanvas("c_x_emittance_mup","c_x_emittance_mup"); 
    c_x_emittance_mup->cd();
    hist2D_emittance_x_mup_MC->SetTitle("emittance #mu^{+}");
    hist2D_emittance_x_mup_MC->GetXaxis()->SetTitle("x [mm]");
    hist2D_emittance_x_mup_MC->GetYaxis()->SetTitle("x' [rad]");
    hist2D_emittance_x_mup_MC->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kBird);
    hist2D_emittance_x_mup_MC->Draw("COLZ");
    // --- Binned formula
    // float emittanceValue_mup = sqrt(hist2D_emittance_x_mup_MC->GetCovariance(1,1)*hist2D_emittance_x_mup_MC->GetCovariance(2,2) - hist2D_emittance_x_mup_MC->GetCovariance(2,1)*hist2D_emittance_x_mup_MC->GetCovariance(1,2));
    // --- Unbinned formula
    if( false ){
      cout << vec_emittance_x_mup_MC.size() << endl;
      cout << vec_emittance_xprime_mup_MC.size() << endl;
      for(UInt_t i_dbg=0;i_dbg<vec_emittance_x_mup_MC.size();i_dbg++){
        cout << vec_emittance_x_mup_MC[i_dbg] << " " << vec_emittance_xprime_mup_MC[i_dbg] << endl;
      }
    }
    float emittanceValue_mup = getemittance(vec_emittance_x_mup_MC, vec_emittance_xprime_mup_MC);
    cout << "emittanceValue_mup = " << emittanceValue_mup << endl;
    TPaveText* pv_x_emittance_mup = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_x_emittance_mup->AddText(Form("#epsilon = %.2f",emittanceValue_mup));
    pv_x_emittance_mup->AddText("nm #times rad");
    pv_x_emittance_mup->SetFillColor(kWhite);
    pv_x_emittance_mup->SetBorderSize(0);
    pv_x_emittance_mup->SetTextFont(40);
    pv_x_emittance_mup->SetTextSize(0.05);
    pv_x_emittance_mup->SetTextFont(42);
    pv_x_emittance_mup->SetTextAlign(22); //text centering
    pv_x_emittance_mup->Draw();
    c_x_emittance_mup->Update();
    c_x_emittance_mup->SaveAs((plotOutputPath + "/" + c_x_emittance_mup->GetName() + ".png"));
    
    TCanvas* c_x_emittance_mum = new TCanvas("c_x_emittance_mum","c_x_emittance_mum"); 
    c_x_emittance_mum->cd();
    hist2D_emittance_x_mum_MC->SetTitle("emittance #mu^{-}");
    hist2D_emittance_x_mum_MC->GetXaxis()->SetTitle("x [mm]");
    hist2D_emittance_x_mum_MC->GetYaxis()->SetTitle("x' [rad]");
    hist2D_emittance_x_mum_MC->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kBird);
    hist2D_emittance_x_mum_MC->Draw("COLZ");
    // -- Binned formula
    // float emittanceValue_mum = sqrt(hist2D_emittance_x_mum_MC->GetCovariance(1,1)*hist2D_emittance_x_mum_MC->GetCovariance(2,2) - hist2D_emittance_x_mum_MC->GetCovariance(2,1)*hist2D_emittance_x_mum_MC->GetCovariance(1,2));
    // --- Unbinned formula
    float emittanceValue_mum = getemittance(vec_emittance_x_mum_MC, vec_emittance_xprime_mum_MC);
    cout << "emittanceValue_mum = " << emittanceValue_mum << endl;
    TPaveText* pv_x_emittance_mum = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_x_emittance_mum->AddText(Form("#epsilon = %.2f",emittanceValue_mum));
    pv_x_emittance_mum->AddText("nm #times rad");
    pv_x_emittance_mum->SetFillColor(kWhite);
    pv_x_emittance_mum->SetBorderSize(0);
    pv_x_emittance_mum->SetTextFont(40);
    pv_x_emittance_mum->SetTextSize(0.05);
    pv_x_emittance_mum->SetTextFont(42);
    pv_x_emittance_mum->SetTextAlign(22); //text centering
    pv_x_emittance_mum->Draw();
    c_x_emittance_mum->Update();
    c_x_emittance_mum->SaveAs((plotOutputPath + "/" + c_x_emittance_mum->GetName() + ".png"));
    
    //1D emittance plots  
    TCanvas* c_x_emittance_1D_mup = new TCanvas("c_x_emittance_1D_mup","c_x_emittance_1D_mup");
    c_x_emittance_1D_mup->cd();
    hist1D_emittance_x_mup_MC->Draw();
    c_x_emittance_1D_mup->SaveAs((plotOutputPath + "/" + c_x_emittance_1D_mup->GetName() + ".png"));
    
    TCanvas* c_x_prime_emittance_1D_mup = new TCanvas("c_x_prime_emittance_1D_mup","c_x_prime_emittance_1D_mup");
    c_x_prime_emittance_1D_mup->cd();
    hist1D_emittance_x_prime_mup_MC->Draw();
    c_x_prime_emittance_1D_mup->SaveAs((plotOutputPath + "/" + c_x_prime_emittance_1D_mup->GetName() + ".png"));
    
    TCanvas* c_x_emittance_1D_mum = new TCanvas("c_x_emittance_1D_mum","c_x_emittance_1D_mum");
    c_x_emittance_1D_mum->cd();
    hist1D_emittance_x_mum_MC->Draw();
    c_x_emittance_1D_mum->SaveAs((plotOutputPath + "/" + c_x_emittance_1D_mum->GetName() + ".png"));
    
    TCanvas* c_x_prime_emittance_1D_mum = new TCanvas("c_x_prime_emittance_1D_mum","c_x_prime_emittance_1D_mum");
    c_x_prime_emittance_1D_mum->cd();
    hist1D_emittance_x_prime_mum_MC->Draw();
    c_x_prime_emittance_1D_mum->SaveAs((plotOutputPath + "/" + c_x_prime_emittance_1D_mum->GetName() + ".png"));

    // -------------------------------------- 
    // save 1D emittance plots in root file
    TFile* fOut_1DemittanceHistos = new TFile("fout_1DemittancePlots_Sep.root","recreate");
    fOut_1DemittanceHistos->cd();
    hist1D_emittance_x_mup_MC      ->Write(hist1D_emittance_x_mup_MC      ->GetName());
    hist1D_emittance_x_prime_mup_MC->Write(hist1D_emittance_x_prime_mup_MC->GetName());
    hist1D_emittance_x_mum_MC      ->Write(hist1D_emittance_x_mum_MC      ->GetName());
    hist1D_emittance_x_prime_mum_MC->Write(hist1D_emittance_x_prime_mum_MC->GetName());
    fOut_1DemittanceHistos->Close();
    // --------------------------------------
    
    
    //1D emittance control plots 
    TCanvas* c_emittanceControl_x_atZref_eplus = new TCanvas("c_emittanceControl_x_atZref_eplus","c_emittanceControl_x_atZref_eplus");
    c_emittanceControl_x_atZref_eplus->cd();
    hist1D_emittanceControl_x_atZref_eplus_MC->SetLineColor(kMagenta);
    hist1D_emittanceControl_x_atZref_eplus_MC->Draw();
    c_emittanceControl_x_atZref_eplus->SaveAs((plotOutputPath + "/" + c_emittanceControl_x_atZref_eplus->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_x_prime_atZref_eplus = new TCanvas("c_emittanceControl_x_prime_atZref_eplus","c_emittanceControl_x_prime_atZref_eplus");
    c_emittanceControl_x_prime_atZref_eplus->cd();
    hist1D_emittanceControl_x_prime_atZref_eplus_MC->SetLineColor(kMagenta);
    hist1D_emittanceControl_x_prime_atZref_eplus_MC->Draw();  
    c_emittanceControl_x_prime_atZref_eplus->SaveAs((plotOutputPath + "/" + c_emittanceControl_x_prime_atZref_eplus->GetName() + ".png"));
    

    TCanvas* c_emittanceControl_x_onDet30_mup = new TCanvas("c_emittanceControl_x_onDet30_mup","c_emittanceControl_x_onDet30_mup");
    c_emittanceControl_x_onDet30_mup->cd();
    hist1D_emittanceControl_x_onDet30_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_x_onDet30_mup_MC->Draw();
    c_emittanceControl_x_onDet30_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_x_onDet30_mup->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_x_atZref_mup = new TCanvas("c_emittanceControl_x_atZref_mup","c_emittanceControl_x_atZref_mup");
    c_emittanceControl_x_atZref_mup->cd();
    hist1D_emittanceControl_x_atZref_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_x_atZref_mup_MC->Draw();
    c_emittanceControl_x_atZref_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_x_atZref_mup->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_xprime_mup = new TCanvas("c_emittanceControl_xprime_mup","c_emittanceControl_xprime_mup");
    c_emittanceControl_xprime_mup->cd();
    hist1D_emittanceControl_xprime_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_xprime_mup_MC->Draw();
    c_emittanceControl_xprime_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_xprime_mup->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_px_mup = new TCanvas("c_emittanceControl_px_mup","c_emittanceControl_px_mup");
    c_emittanceControl_px_mup->cd();
    hist1D_emittanceControl_px_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_px_mup_MC->Draw();
    c_emittanceControl_px_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_px_mup->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_py_mup = new TCanvas("c_emittanceControl_py_mup","c_emittanceControl_py_mup");
    c_emittanceControl_py_mup->cd();
    hist1D_emittanceControl_py_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_py_mup_MC->Draw();
    c_emittanceControl_py_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_py_mup->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_pz_mup = new TCanvas("c_emittanceControl_pz_mup","c_emittanceControl_pz_mup");
    c_emittanceControl_pz_mup->cd();
    hist1D_emittanceControl_pz_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_pz_mup_MC->Draw();
    c_emittanceControl_pz_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_pz_mup->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_pTot_mup = new TCanvas("c_emittanceControl_pTot_mup","c_emittanceControl_pTot_mup");
    c_emittanceControl_pTot_mup->cd();
    hist1D_emittanceControl_pTot_mup_MC->SetLineColor(kRed);
    hist1D_emittanceControl_pTot_mup_MC->Draw();
    c_emittanceControl_pTot_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_pTot_mup->GetName() + ".png"));
    

    TCanvas* c_emittanceControl_x_onDet30_mum = new TCanvas("c_emittanceControl_x_onDet30_mum","c_emittanceControl_x_onDet30_mum");
    c_emittanceControl_x_onDet30_mum->cd();
    hist1D_emittanceControl_x_onDet30_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_x_onDet30_mum_MC->Draw();
    c_emittanceControl_x_onDet30_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_x_onDet30_mum->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_x_atZref_mum = new TCanvas("c_emittanceControl_x_atZref_mum","c_emittanceControl_x_atZref_mum");
    c_emittanceControl_x_atZref_mum->cd();
    hist1D_emittanceControl_x_atZref_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_x_atZref_mum_MC->Draw();
    c_emittanceControl_x_atZref_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_x_atZref_mum->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_xprime_mum = new TCanvas("c_emittanceControl_xprime_mum","c_emittanceControl_xprime_mum");
    c_emittanceControl_xprime_mum->cd();
    hist1D_emittanceControl_xprime_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_xprime_mum_MC->Draw();
    c_emittanceControl_xprime_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_xprime_mum->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_px_mum = new TCanvas("c_emittanceControl_px_mum","c_emittanceControl_px_mum");
    c_emittanceControl_px_mum->cd();
    hist1D_emittanceControl_px_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_px_mum_MC->Draw();
    c_emittanceControl_px_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_px_mum->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_py_mum = new TCanvas("c_emittanceControl_py_mum","c_emittanceControl_py_mum");
    c_emittanceControl_py_mum->cd();
    hist1D_emittanceControl_py_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_py_mum_MC->Draw();
    c_emittanceControl_py_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_py_mum->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_pz_mum = new TCanvas("c_emittanceControl_pz_mum","c_emittanceControl_pz_mum");
    c_emittanceControl_pz_mum->cd();
    hist1D_emittanceControl_pz_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_pz_mum_MC->Draw();
    c_emittanceControl_pz_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_pz_mum->GetName() + ".png"));
    
    TCanvas* c_emittanceControl_pTot_mum = new TCanvas("c_emittanceControl_pTot_mum","c_emittanceControl_pTot_mum");
    c_emittanceControl_pTot_mum->cd();
    hist1D_emittanceControl_pTot_mum_MC->SetLineColor(kBlue);
    hist1D_emittanceControl_pTot_mum_MC->Draw();
    c_emittanceControl_pTot_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_pTot_mum->GetName() + ".png"));
    
    // 2D emittance control plots 
    TCanvas* c_emittanceControl_emittance_positron = new TCanvas("c_emittanceControl_emittance_positron","c_emittanceControl_emittance_positron");
    c_emittanceControl_emittance_positron->cd();
    hist2D_emittanceControl_emittance_positron_MC->SetTitle("raw emittance e^{+}");
    hist2D_emittanceControl_emittance_positron_MC->GetXaxis()->SetTitle("x [mm]");
    hist2D_emittanceControl_emittance_positron_MC->GetYaxis()->SetTitle("x' [rad]");
    hist2D_emittanceControl_emittance_positron_MC->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kCool);
    hist2D_emittanceControl_emittance_positron_MC->Draw("COLZ");
    // --- Binned formula
    // float emittanceValue_raw_positron = sqrt(hist2D_emittanceControl_emittance_positron_MC->GetCovariance(1,1)*hist2D_emittanceControl_emittance_positron_MC->GetCovariance(2,2) - hist2D_emittanceControl_emittance_positron_MC->GetCovariance(2,1)*hist2D_emittanceControl_emittance_positron_MC->GetCovariance(1,2));
    // --- Unbinned formula
    float emittanceValue_raw_positron_MC = getemittance(vec_emittanceControl_emittance_x_positron_MC, vec_emittanceControl_emittance_xprime_positron_MC);
    TPaveText* pv_emittanceControl_emittance_positron = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_emittanceControl_emittance_positron->AddText(Form("#epsilon = %.2f",emittanceValue_raw_positron_MC));
    pv_emittanceControl_emittance_positron->AddText("nm #times rad");
    pv_emittanceControl_emittance_positron->SetFillColor(kWhite);
    pv_emittanceControl_emittance_positron->SetBorderSize(0);
    pv_emittanceControl_emittance_positron->SetTextFont(40);
    pv_emittanceControl_emittance_positron->SetTextSize(0.05);
    pv_emittanceControl_emittance_positron->SetTextFont(42);
    pv_emittanceControl_emittance_positron->SetTextAlign(22); //text centering
    pv_emittanceControl_emittance_positron->Draw();
    c_emittanceControl_emittance_positron->Update();
    c_emittanceControl_emittance_positron->SaveAs((plotOutputPath + "/" + c_emittanceControl_emittance_positron->GetName() + ".png"));
    
    
    TCanvas* c_emittanceControl_emittance_mup = new TCanvas("c_emittanceControl_emittance_mup","c_emittanceControl_emittance_mup");
    c_emittanceControl_emittance_mup->cd();
    hist2D_emittanceControl_emittance_mup_MC->SetTitle("raw emittance #mu^{+}");
    hist2D_emittanceControl_emittance_mup_MC->GetXaxis()->SetTitle("x [mm]");
    hist2D_emittanceControl_emittance_mup_MC->GetYaxis()->SetTitle("x' [rad]");
    hist2D_emittanceControl_emittance_mup_MC->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kAvocado);
    hist2D_emittanceControl_emittance_mup_MC->Draw("COLZ");
    // --- Binned formula
    //float emittanceValue_raw_mup = sqrt(hist2D_emittanceControl_emittance_mup_MC->GetCovariance(1,1)*hist2D_emittanceControl_emittance_mup_MC->GetCovariance(2,2) - hist2D_emittanceControl_emittance_mup_MC->GetCovariance(2,1)*hist2D_emittanceControl_emittance_mup_MC->GetCovariance(1,2));
    // --- Unbinned formula
    float emittanceValue_raw_mup = getemittance(vec_emittanceControl_emittance_x_mup_MC, vec_emittanceControl_emittance_xprime_mup_MC);
    TPaveText* pv_emittanceControl_emittance_mup = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_emittanceControl_emittance_mup->AddText(Form("#epsilon = %.2f",emittanceValue_raw_mup));
    pv_emittanceControl_emittance_mup->AddText("nm #times rad");
    pv_emittanceControl_emittance_mup->SetFillColor(kWhite);
    pv_emittanceControl_emittance_mup->SetBorderSize(0);
    pv_emittanceControl_emittance_mup->SetTextFont(40);
    pv_emittanceControl_emittance_mup->SetTextSize(0.05);
    pv_emittanceControl_emittance_mup->SetTextFont(42);
    pv_emittanceControl_emittance_mup->SetTextAlign(22); //text centering
    pv_emittanceControl_emittance_mup->Draw();
    c_emittanceControl_emittance_mup->Update();
    c_emittanceControl_emittance_mup->SaveAs((plotOutputPath + "/" + c_emittanceControl_emittance_mup->GetName() + ".png"));
    
    
    TCanvas* c_emittanceControl_emittance_mum = new TCanvas("c_emittanceControl_emittance_mum","c_emittanceControl_emittance_mum");
    c_emittanceControl_emittance_mum->cd();
    hist2D_emittanceControl_emittance_mum_MC->SetTitle("raw emittance #mu^{-}");
    hist2D_emittanceControl_emittance_mum_MC->GetXaxis()->SetTitle("x [mm]");
    hist2D_emittanceControl_emittance_mum_MC->GetYaxis()->SetTitle("x' [rad]");
    hist2D_emittanceControl_emittance_mum_MC->GetYaxis()->SetTitleOffset(1.4);
    gStyle->SetPalette(kAvocado);
    hist2D_emittanceControl_emittance_mum_MC->Draw("COLZ");
    // --- Binned formula
    // float emittanceValue_raw_mum = sqrt(hist2D_emittanceControl_emittance_mum_MC->GetCovariance(1,1)*hist2D_emittanceControl_emittance_mum_MC->GetCovariance(2,2) - hist2D_emittanceControl_emittance_mum_MC->GetCovariance(2,1)*hist2D_emittanceControl_emittance_mum_MC->GetCovariance(1,2));
    // --- Unbinned formula
    float emittanceValue_raw_mum = getemittance(vec_emittanceControl_emittance_x_mum_MC, vec_emittanceControl_emittance_xprime_mum_MC);
    TPaveText* pv_emittanceControl_emittance_mum = new TPaveText(0.15,0.75,0.35,0.85,"brNDC");
    pv_emittanceControl_emittance_mum->AddText(Form("#epsilon = %.2f",emittanceValue_raw_mum));
    pv_emittanceControl_emittance_mum->AddText("nm #times rad");
    pv_emittanceControl_emittance_mum->SetFillColor(kWhite);
    pv_emittanceControl_emittance_mum->SetBorderSize(0);
    pv_emittanceControl_emittance_mum->SetTextFont(40);
    pv_emittanceControl_emittance_mum->SetTextSize(0.05);
    pv_emittanceControl_emittance_mum->SetTextFont(42);
    pv_emittanceControl_emittance_mum->SetTextAlign(22); //text centering
    pv_emittanceControl_emittance_mum->Draw();
    c_emittanceControl_emittance_mum->Update();
    c_emittanceControl_emittance_mum->SaveAs((plotOutputPath + "/" + c_emittanceControl_emittance_mum->GetName() + ".png"));
  
  } // true

  cout<<" Plots done! =) "<<endl; 


}





// main function 
void plotEmittance(){

  // define input files 
  TString inputFile_Data_Aug2018_Be6cm = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/aug18/reco-aug18.root";
  TString inputFile_MC_Aug2018_Be6cm   = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/aug18/reco-mupmum.root";
  TString inputFile_MC_Sep2018_Be6cm   = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/with_vtx_fit/reco-mupmum-Be6cm.root";
  TString inputFile_MC_Sep2018_C6cm    = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/reco-mupmum-C6cm.root";
  TString inputFile_MC_Sep2018_C2cm    = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/reco-mupmum-C2cm.root";
  TString inputFile_MC_Sep2018_Be6cm_new = "/afs/cern.ch/user/a/abertoli/public/lemma/reco/sep18/with_vtx_fit/reco-mupmum-Be6cm-GausGaus.root"; 
  
  // define output path and make output directory
  //TString plotOutputPath = "190327_Emittance_August2018_targetBe6cm_DATA";
  //TString plotOutputPath = "Emittance_August2018_targetBe6cm_MC";
  //TString plotOutputPath = "190327_Emittance_September2018_targetBe6cm_MC";
  //TString plotOutputPath = "190327_Emittance_September2018_targetC6cm_MC";
  //TString plotOutputPath = "190327_Emittance_September2018_targetC2cm_MC";
  //TString plotOutputPath = "190516_Emittance_Sep18_Be6cm_GausGaus";
  TString plotOutputPath = "test_plotEmittance";
  gSystem->Exec(("mkdir -p "+plotOutputPath));

  // choose type of target
  double zEndTarget = 10.*(457.9+3.-84.6);   // [mm] - dataset: AUGUST 2018    Be target 6 cm
  //double zEndTarget = 10.*(460.93+3.-82.78); // [mm] - dataset: SEPTEMBER 2018 Be target 6 cm and C target 6cm
  //double zEndTarget = 10.*(460.93+1.-82.78); // [mm] - dataset: SEPTEMBER 2018 C  target 2 cm 

  // choose configuration
  double zPosDet20 = 10.*(444.5-84.6); // [mm] - AUGUST 2018
  // double zPosDet20 = 10.*(441.63-82.78); // [mm] - SEPTEMBER 2018
 
  // --- call do the histos function
  // arguments: input file, label for data or MC

  doTheHistos(inputFile_Data_Aug2018_Be6cm, "reclev",   zEndTarget, zPosDet20, plotOutputPath);
  //doTheHistos(inputFile_MC_Aug2018_Be6cm,   "MC",       zEndTarget, zPosDet20, plotOutputPath);
  //doTheHistos(inputFile_MC_Aug2018_Be6cm,   "MCreclev", zEndTarget, zPosDet20, plotOutputPath);
  //doTheHistos(inputFile_MC_Sep2018_Be6cm,   "MC",       zEndTarget, zPosDet20, plotOutputPath);
  //doTheHistos(inputFile_MC_Sep2018_C6cm,    "MC",       zEndTarget, zPosDet20, plotOutputPath); 
  //doTheHistos(inputFile_MC_Sep2018_C2cm,    "MC",       zEndTarget, zPosDet20, plotOutputPath);
  // gauss profile of input beam
  //doTheHistos(inputFile_MC_Sep2018_Be6cm_new, "MC",     zEndTarget, zPosDet20, plotOutputPath);

}
