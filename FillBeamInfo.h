
// *****************************************************
//
// author: Marcello Rotondo
//         13/11/2018
// 
// usage: useful to easily access information from 
//        positrons in the beam
//
// *****************************************************

#include <math.h>

// Function prototype: implemented below 
void fit_line(
         double x[], 
         double y[], 
         int ndata, 
         double sig[], 
         double &a, 
         double &b, 
         double &siga, 
         double &sigb, 
         double &chi2);

Int_t ReturnBeamInfo(/*Inputs*/
                     const Int_t subdet[], 
                     const Double_t xh[], 
                     const Double_t zh[], 
                     const Int_t nhits,
                     Double_t zpos,    /* z position where tracks are extrapolated */
                     /*Outputs*/
                     Double_t xpatz[], /* x position @ z */
                     Double_t thatz[]  /* exit angel */)   
{
// ---------------------------------------------------------------------
// Return number of positron in the beam
//  Fill vectors 
//    xpatz: x coordinate at zpos position and angle
//    thatz: angle at position zpos (this actually does not depend on z
// ---------------------------------------------------------------------
    Double_t SiHitError(0.015); // 15um for Si 10 and Si 20
    Double_t BeamSpread(0.000400); // 400urad
 
    Int_t npatz=0;
    for(Int_t i=0; i<10; i++){
      xpatz[i]=-.999;
      thatz[i]=-.999;
    }

    // --------------------------------------------------------
    // Loop over hits, fill temporary vector with coordinates
    // in Si10 and Si20
    // --------------------------------------------------------
    Int_t nx10, nx20;
    Double_t x10[10], x20[10], z10[10], z20[10];
    nx10=0; nx20=0;
      for(Int_t i=0; i<10; i++){ 
        x10[i]=-999.; x20[i]=-999.; z10[i]=-999.; z20[i]=-999.;
    }
 
    for(Int_t jhit=0; jhit<nhits; jhit++) {
      if(xh[jhit]>-999){
        if(subdet[jhit]==10) {
          x10[nx10]=xh[jhit];
          z10[nx10]=zh[jhit];
          nx10++;
        }//if(subdet[jhit]) 
        if(subdet[jhit]==20) {
          x20[nx20]=xh[jhit];
          z20[nx20]=zh[jhit];
          nx20++;
        }//if(subdet[jhit]... 
      }//if(xh[...
    }//for(Int_t...

    // -----------------------------------------------------------
    // Loop over all hits in Si10, combine with any hit in Si20
    //      store only of theta_track is < 220urad (beam spread)
    // ----------------------------------------------------------- 
    int np=2; // 2 hits
    double x[np], z[np], ex[np];

    for(Int_t i=0; i < nx10 ; i++){
      for(Int_t j=0; j < nx20 ; j++){
        // only 10 tracks
        if( npatz > 9 ) continue;
        x[0] = x10[i]; x[1] = x20[j];
        z[0] = z10[i]; z[1] = z20[j];
        ex[0] = SiHitError; ex[1]=SiHitError; //Si10 & Si20 single hit uncertainty is 30um
        // Track:   X = a + b*Z
        double a,b,ae,be,chi2;
        fit_line(z,x,np,ex,a,b,ae,be,chi2);
        double theta = atan(b);
        //cout << " A: " << a << " " << "B: " << b <<endl;
        if( abs(theta) < BeamSpread ){
          xpatz[npatz] = a + b*zpos;
          thatz[npatz] = theta;
          npatz++;
        } // if consistent with beam 
      } // si 20
    } // si 10 
      
    return npatz;
}


void fit_line(
         double x[], 
         double y[], 
         int ndata, 
         double sig[], 
         double &a, 
         double &b, 
         double &siga, 
         double &sigb, 
         double &chi2)
/*
 * From Numerical Recipes in C 
 * Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations
sig[1..ndata], fit them to a straight line y = a + bx by minimizing χ2. Returned are
a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and the
goodness-of-fit probability q (that the fit would have chi2 this large or larger). 
*/

{
  int i;
  double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  b=0.0;
  a=0.0;

  ss=0.0;
  for (i=0;i<ndata;i++) { 
    wt=1.0/pow(sig[i],2);
    ss += wt;
    sx += x[i]*wt;
    sy += y[i]*wt;
  }

  sxoss=sx/ss;
  for (i=0;i<ndata;i++) {
    t=(x[i]-sxoss)/sig[i];
    st2 += t*t;
    b += t*y[i]/sig[i];
  }

  b /= st2; 
  a=(sy-sx*b)/ss;
  siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  sigb=sqrt(1.0/st2);
  chi2=0.0; 
  
  for (i=0;i<ndata;i++) chi2 += pow(((y[i]-a-b*x[i])/sig[i]),2);
} 
