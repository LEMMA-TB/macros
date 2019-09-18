
// root -b -q xprimeoverlap.C++
// tested with ROOT 6.18/00

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


using namespace std;


void xprimeoverlap(){

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2(true);


  // DATA file
  TFile *inFile_Aug18_data = TFile::Open("fout_1DemittancePlots_Aug18data.root");
 
  TH1D* x_prime_mup_data = (TH1D*)inFile_Aug18_data->Get("hist1D_emittance_x_prime_mup_MC");
  TH1D* x_prime_mum_data = (TH1D*)inFile_Aug18_data->Get("hist1D_emittance_x_prime_mum_MC");

  // MC file
  TFile *inFile_Aug18_MC = TFile::Open("fout_1DemittancePlots_Aug18MC.root");

  TH1D* x_prime_mup_MC = (TH1D*)inFile_Aug18_MC->Get("hist1D_emittance_x_prime_mup_MC");
  TH1D* x_prime_mum_MC = (TH1D*)inFile_Aug18_MC->Get("hist1D_emittance_x_prime_mum_MC");

  
  // add histos
  x_prime_mup_data ->Add(x_prime_mum_data);
  x_prime_mup_MC   ->Add(x_prime_mum_MC);


  // normalize hists
  x_prime_mup_data->Scale(1./x_prime_mup_data->GetEntries());
  x_prime_mup_MC->Scale(1./x_prime_mup_MC->GetEntries());

  cout << "number of bins along X: " << x_prime_mup_data->GetXaxis()->GetNbins() << endl;
  Int_t rebin_factor=4;
  x_prime_mup_data->Rebin(rebin_factor);
  x_prime_mup_MC  ->Rebin(rebin_factor);

  // cosmetics
  x_prime_mup_data->SetTitle("uncorrected x'");
  x_prime_mup_data->SetMaximum(0.15);
  x_prime_mup_data->GetXaxis()->SetTitle("x' [rad]");
  x_prime_mup_data->GetYaxis()->SetTitle("a.u.");
  x_prime_mup_data->SetMarkerStyle(20);
  x_prime_mup_data->SetMarkerColor(kBlack);
  x_prime_mup_data->SetLineColor(kBlack);

  x_prime_mup_MC->SetLineColor(kGreen+2);
  x_prime_mup_MC->SetFillColor(kGreen-9);

  TCanvas *c_out = new TCanvas();
  c_out->cd();
  x_prime_mup_data->Draw("pe"); 
  x_prime_mup_MC->Draw("samehisto");
  x_prime_mup_data->Draw("samepe");

  TLegend* l = new TLegend(0.80,0.84,0.98,0.97);
  l->AddEntry(x_prime_mup_data,"Data","pl");
  l->AddEntry(x_prime_mup_MC,  "MC",  "f");
  l->SetFillColor(kWhite);
  l->SetLineColor(kBlack);
  l->SetTextFont(43);
  l->SetTextSize(20);
  l->Draw();

  c_out->Update();
  c_out->SaveAs("xprimeoverlap.png");

}
