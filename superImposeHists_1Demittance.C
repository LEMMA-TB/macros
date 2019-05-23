// ******************************************************************
// this script reads x and x' histos from 2 files (Aug18 and Sep18)
// and superimposes them (both for mu+ and mu-)
// -- input files have to be produced with plotEmittance.C script
//
//
// run with: 
//   root -l -b -q superImposeHists_1Demittance.C++
//
//*******************************************************************

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


void superImposeHists_1Demittance(){

  TFile *inFile_Aug18 = TFile::Open("fout_1DemittancePlots_Aug.root");
  
  TH1F* hist1D_emittance_x_mup_MC_Aug       = (TH1F*)inFile_Aug18->Get("hist1D_emittance_x_mup_MC");
  TH1F* hist1D_emittance_x_prime_mup_MC_Aug = (TH1F*)inFile_Aug18->Get("hist1D_emittance_x_prime_mup_MC");
  TH1F* hist1D_emittance_x_mum_MC_Aug       = (TH1F*)inFile_Aug18->Get("hist1D_emittance_x_mum_MC");
  TH1F* hist1D_emittance_x_prime_mum_MC_Aug = (TH1F*)inFile_Aug18->Get("hist1D_emittance_x_prime_mum_MC");

  
  TFile *inFile_Sep18 = TFile::Open("fout_1DemittancePlots_Sep.root");
  
  TH1F* hist1D_emittance_x_mup_MC_Sep       = (TH1F*)inFile_Sep18->Get("hist1D_emittance_x_mup_MC");
  TH1F* hist1D_emittance_x_prime_mup_MC_Sep = (TH1F*)inFile_Sep18->Get("hist1D_emittance_x_prime_mup_MC");
  TH1F* hist1D_emittance_x_mum_MC_Sep       = (TH1F*)inFile_Sep18->Get("hist1D_emittance_x_mum_MC");
  TH1F* hist1D_emittance_x_prime_mum_MC_Sep = (TH1F*)inFile_Sep18->Get("hist1D_emittance_x_prime_mum_MC");



  TString plotOutputPath = "190523_superImposeHist_Aug18_Sep18_Be6cm";
  gSystem->Exec(("mkdir -p "+plotOutputPath));


  gStyle->SetOptStat(0);
  

  TCanvas* c_emittance_x_mup = new TCanvas("c_emittance_x_mup","c_emittance_x_mup");
  c_emittance_x_mup->cd();
  hist1D_emittance_x_mup_MC_Aug->SetLineColor(kRed);
  hist1D_emittance_x_mup_MC_Aug->SetLineWidth(2);
  hist1D_emittance_x_mup_MC_Aug->SetMaximum(1.1 * max(hist1D_emittance_x_mup_MC_Aug->GetMaximum(),hist1D_emittance_x_mup_MC_Sep->GetMaximum()));
  hist1D_emittance_x_mup_MC_Aug->SetTitle("x #mu^{+}");
  hist1D_emittance_x_mup_MC_Aug->GetXaxis()->SetTitle("x [mm]");
  hist1D_emittance_x_mup_MC_Aug->GetYaxis()->SetTitle("counts");
  hist1D_emittance_x_mup_MC_Aug->Draw("hist");
  hist1D_emittance_x_mup_MC_Sep->SetLineColor(kBlue);
  hist1D_emittance_x_mup_MC_Sep->SetLineWidth(1);
  hist1D_emittance_x_mup_MC_Sep->Draw("same hist");
  TLegend* l_x_mup = new TLegend(0.76,0.55,0.98,0.96);
  l_x_mup->AddEntry(hist1D_emittance_x_mup_MC_Aug,"August 2018","f");
  l_x_mup->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_mup_MC_Aug->GetEntries()),"");
  l_x_mup->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_mup_MC_Aug->GetMean()),"");  
  l_x_mup->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_mup_MC_Aug->GetRMS()),"");  
  l_x_mup->AddEntry(hist1D_emittance_x_mup_MC_Sep, "September 2018", "f");
  l_x_mup->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_mup_MC_Sep->GetEntries()),"");
  l_x_mup->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_mup_MC_Sep->GetMean()),"");
  l_x_mup->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_mup_MC_Sep->GetRMS()),"");
  l_x_mup->SetFillColor(kWhite);
  l_x_mup->SetLineColor(kBlack);
  l_x_mup->SetTextFont(43);
  l_x_mup->SetTextSize(14);
  l_x_mup->Draw();
  c_emittance_x_mup->Update();
  c_emittance_x_mup->SaveAs((plotOutputPath + "/c_emittance_x_mup.png"));


  TCanvas* c_emittance_x_prime_mup = new TCanvas("c_emittance_x_prime_mup","c_emittance_x_prime_mup");
  c_emittance_x_prime_mup->cd();
  hist1D_emittance_x_prime_mup_MC_Aug->SetLineColor(kRed);
  hist1D_emittance_x_prime_mup_MC_Aug->SetLineWidth(2);
  hist1D_emittance_x_prime_mup_MC_Aug->SetMaximum(1.1 * max(hist1D_emittance_x_prime_mup_MC_Aug->GetMaximum(),hist1D_emittance_x_prime_mup_MC_Sep->GetMaximum()));
  hist1D_emittance_x_prime_mup_MC_Aug->SetTitle("x' #mu^{+}");
  hist1D_emittance_x_prime_mup_MC_Aug->GetXaxis()->SetTitle("x' [rad]");
  hist1D_emittance_x_prime_mup_MC_Aug->GetYaxis()->SetTitle("counts");
  hist1D_emittance_x_prime_mup_MC_Aug->Draw("hist");
  hist1D_emittance_x_prime_mup_MC_Sep->SetLineColor(kBlue);
  hist1D_emittance_x_prime_mup_MC_Sep->SetLineWidth(1);
  hist1D_emittance_x_prime_mup_MC_Sep->Draw("same hist");
  TLegend* l_x_prime_mup = new TLegend(0.76,0.55,0.98,0.96);
  l_x_prime_mup->AddEntry(hist1D_emittance_x_prime_mup_MC_Aug,"August 2018","f");
  l_x_prime_mup->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_prime_mup_MC_Aug->GetEntries()),"");
  l_x_prime_mup->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_prime_mup_MC_Aug->GetMean()),"");  
  l_x_prime_mup->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_prime_mup_MC_Aug->GetRMS()),"");  
  l_x_prime_mup->AddEntry(hist1D_emittance_x_prime_mup_MC_Sep, "September 2018", "f");
  l_x_prime_mup->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_prime_mup_MC_Sep->GetEntries()),"");
  l_x_prime_mup->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_prime_mup_MC_Sep->GetMean()),"");
  l_x_prime_mup->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_prime_mup_MC_Sep->GetRMS()),"");
  l_x_prime_mup->SetFillColor(kWhite);
  l_x_prime_mup->SetLineColor(kBlack);
  l_x_prime_mup->SetTextFont(43);
  l_x_prime_mup->SetTextSize(14);
  l_x_prime_mup->Draw();
  c_emittance_x_prime_mup->Update();
  c_emittance_x_prime_mup->SaveAs((plotOutputPath + "/c_emittance_x_prime_mup.png"));



  TCanvas* c_emittance_x_mum = new TCanvas("c_emittance_x_mum","c_emittance_x_mum");
  c_emittance_x_mum->cd();
  hist1D_emittance_x_mum_MC_Aug->SetLineColor(kRed);
  hist1D_emittance_x_mum_MC_Aug->SetLineWidth(2);
  hist1D_emittance_x_mum_MC_Aug->SetMaximum(1.1 * max(hist1D_emittance_x_mum_MC_Aug->GetMaximum(),hist1D_emittance_x_mum_MC_Sep->GetMaximum()));
  hist1D_emittance_x_mum_MC_Aug->SetTitle("x #mu^{-}");
  hist1D_emittance_x_mum_MC_Aug->GetXaxis()->SetTitle("x [mm]");
  hist1D_emittance_x_mum_MC_Aug->GetYaxis()->SetTitle("counts");
  hist1D_emittance_x_mum_MC_Aug->Draw("hist");
  hist1D_emittance_x_mum_MC_Sep->SetLineColor(kBlue);
  hist1D_emittance_x_mum_MC_Sep->SetLineWidth(1);
  hist1D_emittance_x_mum_MC_Sep->Draw("same hist");
  TLegend* l_x_mum = new TLegend(0.76,0.55,0.98,0.96);
  l_x_mum->AddEntry(hist1D_emittance_x_mum_MC_Aug,"August 2018","f");
  l_x_mum->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_mum_MC_Aug->GetEntries()),"");
  l_x_mum->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_mum_MC_Aug->GetMean()),"");  
  l_x_mum->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_mum_MC_Aug->GetRMS()),"");  
  l_x_mum->AddEntry(hist1D_emittance_x_mum_MC_Sep, "September 2018", "f");
  l_x_mum->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_mum_MC_Sep->GetEntries()),"");
  l_x_mum->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_mum_MC_Sep->GetMean()),"");
  l_x_mum->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_mum_MC_Sep->GetRMS()),"");
  l_x_mum->SetFillColor(kWhite);
  l_x_mum->SetLineColor(kBlack);
  l_x_mum->SetTextFont(43);
  l_x_mum->SetTextSize(14);
  l_x_mum->Draw();
  c_emittance_x_mum->Update();
  c_emittance_x_mum->SaveAs((plotOutputPath + "/c_emittance_x_mum.png"));


  TCanvas* c_emittance_x_prime_mum = new TCanvas("c_emittance_x_prime_mum","c_emittance_x_prime_mum");
  c_emittance_x_prime_mum->cd();
  hist1D_emittance_x_prime_mum_MC_Aug->SetLineColor(kRed);
  hist1D_emittance_x_prime_mum_MC_Aug->SetLineWidth(2);
  hist1D_emittance_x_prime_mum_MC_Aug->SetMaximum(1.1 * max(hist1D_emittance_x_prime_mum_MC_Aug->GetMaximum(),hist1D_emittance_x_prime_mum_MC_Sep->GetMaximum()));
  hist1D_emittance_x_prime_mum_MC_Aug->SetTitle("x' #mu^{-}");
  hist1D_emittance_x_prime_mum_MC_Aug->GetXaxis()->SetTitle("x' [rad]");
  hist1D_emittance_x_prime_mum_MC_Aug->GetYaxis()->SetTitle("counts");
  hist1D_emittance_x_prime_mum_MC_Aug->Draw("hist");
  hist1D_emittance_x_prime_mum_MC_Sep->SetLineColor(kBlue);
  hist1D_emittance_x_prime_mum_MC_Sep->SetLineWidth(1);
  hist1D_emittance_x_prime_mum_MC_Sep->Draw("same hist");
  TLegend* l_x_prime_mum = new TLegend(0.76,0.55,0.98,0.96);
  l_x_prime_mum->AddEntry(hist1D_emittance_x_prime_mum_MC_Aug,"August 2018","f");
  l_x_prime_mum->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_prime_mum_MC_Aug->GetEntries()),"");
  l_x_prime_mum->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_prime_mum_MC_Aug->GetMean()),"");  
  l_x_prime_mum->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_prime_mum_MC_Aug->GetRMS()),"");  
  l_x_prime_mum->AddEntry(hist1D_emittance_x_prime_mum_MC_Sep, "September 2018", "f");
  l_x_prime_mum->AddEntry((TObject*)0,Form("entries: %.2f",hist1D_emittance_x_prime_mum_MC_Sep->GetEntries()),"");
  l_x_prime_mum->AddEntry((TObject*)0,Form("mean: %.6f",   hist1D_emittance_x_prime_mum_MC_Sep->GetMean()),"");
  l_x_prime_mum->AddEntry((TObject*)0,Form("stDev: %.6f",  hist1D_emittance_x_prime_mum_MC_Sep->GetRMS()),"");
  l_x_prime_mum->SetFillColor(kWhite);
  l_x_prime_mum->SetLineColor(kBlack);
  l_x_prime_mum->SetTextFont(43);
  l_x_prime_mum->SetTextSize(14);
  l_x_prime_mum->Draw();
  c_emittance_x_prime_mum->Update();
  c_emittance_x_prime_mum->SaveAs((plotOutputPath + "/c_emittance_x_prime_mum.png"));
  


}
