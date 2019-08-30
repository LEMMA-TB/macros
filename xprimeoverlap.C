
// root -b -q xprimeoverlap.C
// tested with ROOT 6.18/00

void xprimeoverlap(){

  gStyle->SetOptStat(0);

  TString data_file="1DemittancePlots_Aug18data.root";
  TString MC_file="1DemittancePlots_Aug18MC.root";

  if( !gSystem->AccessPathName( data_file ) ||
      !gSystem->AccessPathName( MC_file ) ){
    cout << "missing input files ..." << endl;
    exit;
  }

  TFile *p_data_file = new TFile(data_file);
  TFile *p_MC_file = new TFile(MC_file);

  TH1F *phxp_mup_data = (TH1F*)p_data_file->Get("hist1D_emittance_x_prime_mup_MC");
  TH1F *phxp_mum_data = (TH1F*)p_data_file->Get("hist1D_emittance_x_prime_mum_MC");

  TH1F *phxp_mup_MC = (TH1F*)p_MC_file->Get("hist1D_emittance_x_prime_mup_MC");
  TH1F *phxp_mum_MC = (TH1F*)p_MC_file->Get("hist1D_emittance_x_prime_mum_MC");

  TH1F *phxp_data = (TH1F*)phxp_mup_data->Clone("emittance_x_prime_data"); phxp_data->Reset();
  TH1F *phxp_MC = (TH1F*)phxp_mup_MC->Clone("emittance_x_prime_MC"); phxp_MC->Reset();

  phxp_data->Add(phxp_mup_data,phxp_mum_data);
  phxp_MC->Add(phxp_mup_MC,phxp_mum_MC);

  phxp_data->Scale(1./phxp_data->GetEntries());
  phxp_MC->Scale(1./phxp_MC->GetEntries());

  cout << "number of bins along X: " << phxp_data->GetXaxis()->GetNbins() << endl;
  Int_t rebin_factor=4;
  phxp_data->Rebin(rebin_factor);
  phxp_MC->Rebin(rebin_factor);

  // cosmetics
  phxp_data->SetTitle("uncorrected x'");
  phxp_data->SetMaximum(0.15);
  phxp_data->GetXaxis()->SetTitle("x' [rad]");
  phxp_data->GetYaxis()->SetTitle("a.u.");

  TCanvas *c_out = new TCanvas();
  c_out->cd();
  phxp_data->Draw("pe"); phxp_MC->Draw("samehisto");
  c_out->SaveAs("xprimeoverlap.png");

}
