
TRandom *r=new TRandom();

Double_t getemittance(Int_t n_points, Double_t x[5000], Double_t xp[5000]){

  Double_t emittance=0.;

  Double_t x2=0.;
  Double_t xp2=0.;
  Double_t xxp=0.;
  for(Int_t i=0;i<n_points;i++){
    x2+=(x[i]*x[i]);
    xp2+=(xp[i]*xp[i]);
    xxp+=(x[i]*xp[i]);
  }
  x2*=1./float(n_points);
  xp2*=1./float(n_points);
  xxp*=1./float(n_points);

  emittance=TMath::Sqrt(x2*xp2-xxp*xxp);

  return emittance;

}

void Emitest(){

  // ROOT 6.02.01

  Double_t x[5000];
  Double_t xprime[5000];

  Int_t n_points=5000;

  Int_t n_bins=50; Double_t fatness=150.;
  Double_t x_range=5.*fatness*1e-3; // mm
  Double_t xprime_range=5.*fatness*1e-6; // rad
  // units: mm rad
  TH2F *phxxprime = new TH2F("hxxprime","emittance",n_bins,-1.*x_range,x_range,
   n_bins,-1.*xprime_range,xprime_range);
  TH1F *phx = new TH1F("hx","x",n_bins,-1.*x_range,x_range);
  TH1F *phxprime = new TH1F("hxprime","xprime",n_bins,-1.*xprime_range,xprime_range);

  // fill the arrays in case of an alternative definition of emittance
  for(Int_t i=0;i<n_points;i++){
    // cout << i << endl;
    x[i]=r->Gaus(0,fatness*1e-3);
    xprime[i]=r->Gaus(0,fatness*1e-6);
  }
  Double_t emittance_unbin=getemittance(n_points,x,xprime);
  cout << emittance_unbin << endl;

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
