
TRandom *r=new TRandom();

void Emitest(){

  // ROOT 6.02.01

  Double_t x[5000];
  Double_t xprime[5000];

  Int_t n_points=5000;

  // units: mm rad
  TH2F *phxxprime = new TH2F("hxxprime","emittance",100,-0.5,0.5,100,-0.5,0.5);
  TH1F *phx = new TH1F("hx","x",100,-0.5,0.5);
  TH1F *phxprime = new TH1F("hxprime","xprime",100,-0.5,0.5);

  // fill the arrays in case of an alternative definition of emittance
  for(Int_t i=0;i<n_points;i++){
    // cout << i << endl;
    x[i]=r->Gaus(0,150e-3);
    xprime[i]=r->Gaus(0,150e-3);
  }

  // fill histos
  for(Int_t j=0;j<n_points;j++){
    // cout << j << endl;
    phx->Fill(x[j]);
    phxprime->Fill(xprime[j]);
    phxxprime->Fill(x[j],xprime[j]);
  }

  TCanvas* cEmittance = new TCanvas("cEmittance"); cEmittance->Divide(2,2);
  cEmittance->cd(1); phx->Draw();
  cEmittance->cd(2); phxprime->Draw();
  cEmittance->cd(3); phxxprime->Draw("BOX");

  // compute emittance from histo
  Double_t emittance = sqrt(phxxprime->GetCovariance(1,1)*phxxprime->GetCovariance(2,2)-
                            phxxprime->GetCovariance(2,1)*phxxprime->GetCovariance(1,2) );
  cout << emittance << endl;

}
