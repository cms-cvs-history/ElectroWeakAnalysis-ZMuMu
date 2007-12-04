#include "TMinuit.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TStyle.h"
#include <iostream>
#include <cmath>
using namespace std;
static const int fitXMin = 74, fitXMax = 106;
static const int nRebinZmt = 20;
static const int nRebinZms = 40;

static size_t fullBins = 0;

double BW (double *x, double * p) {
  double m = p[0];
  double g = p[1];
  double m2 = p[0]*p[0];
  double g2 = p[1]*p[1];
  double g2m2 = g2/m2;
  double s = x[0] * x[0] ;
  double deltaS = s - m2;
  double Zline = 0, Int = 0;
  if ( fabs(deltaS/m2) < 16 ) {
    double prop = ( deltaS * deltaS + s * s * g2m2 );
   Zline =  2 /(TMath::Pi()) * g * s / prop ;
   Int =  m * 5 * ( s - m2 ) / prop;  
  }
  double Gline = p[2]/s;
  return ( (1-fabs(p[3])) * Zline + Gline + p[3] * Int );
}

double BWZmm (double *x, double * p) {
  // x[0] = E ; sp = E'^2 ; y = E'; H(E) = /int{BW(E') * g(E-E') dE'}
  //double Ni = p[4];
  //double Nf = p[3];  
  //double m = p[1];
  //double g = p[2];
  // double * p = { &m, &g, &Nf, &Ni};
  double sigma = p[7];
  double sigma2 = sigma*sigma;
  double mean = p[8];
  const double ng = 1.0 / sqrt(2*TMath::Pi());
  static const int bins = 100;
  double xMin = x[0] - 3. * sigma, xMax = x[0] + 3. * sigma, DeltaX = xMax - xMin;
  double dx =  DeltaX  / bins;
  double f = 0; double fbw =0; double gau =0;
  for( int n = 0; n < bins; ++n ) {
    /// down 
    double y = x[0] + (n+.5) * dx;
    fbw = BW(&y, &p[1]) ;
    gau = exp( -.5 * (x[0] - mean - y)*(x[0] -mean - y)  / sigma2 ) * ng/sigma    ;
    f +=  fbw * gau ;
    //up
    double yy = x[0] - (n+.5) * dx ;
    fbw = BW(&yy, &p[1]) ;
    gau =  exp( -.5 * (x[0] - mean - yy)* (x[0] - mean - yy) / sigma2 )* ng/sigma ;
    f +=  fbw * gau ;
   
        }   
 double BWG =    f * DeltaX / bins  ;
  double l = p[18];
  double a = p[19];
  double b = p[20];
  double N1 = exp (- l * fitXMax) / (l* l * l) * (-l * (l * fitXMax * a + a + l) - b * (l * l * fitXMax * fitXMax +  2 * l * fitXMax + 2 ));
  
  double N2 = exp (- l * fitXMin) / (l* l * l) * (-l * (l * fitXMin * a + a + l) - b * (l * l * fitXMin * fitXMin + 2 * l  * fitXMin + 2 )); 
 return  BWG * p[0] /(10)  * p[5] * p[5] * p[6] * p[6] +  p[17] / (N1-N2) * exp( -l * x[0] )  * (1 + x[0] * ( a + x[0] * b));

  
}

double BWZmt (double *x, double * p) {
  // x[0] = E ; sp = E'^2 ; y = E'; H(E) = /int{BW(E') * g(E-E') dE'}
  //double Ni = p[4];
  //double Nf = p[3];  
  //double m = p[1];
  //double g = p[2];
  // double * p = { &m, &g, &Nf, &Ni};
  double sigma = p[9];
  double sigma2 = sigma*sigma;
  double mean = p[10];
  const double ng = 1.0 / sqrt(2*TMath::Pi());
  static const int bins = 100;
  double xMin = x[0] - 3. * sigma, xMax = x[0] + 3. * sigma, DeltaX = xMax - xMin;
  double dx =  DeltaX  / bins;
  double f = 0; double fbw =0; double gau =0;
  for( int n = 0; n < bins; ++n ) {
    /// down 
    double y = x[0] + (n+.5) * dx;
    fbw = BW(&y, &p[1]) ; 
    gau = exp( -.5 * (x[0] - mean - y)*(x[0] -mean - y)  / sigma2 ) * ng/sigma    ;
    f +=  fbw * gau ;
    //up
    double yy = x[0] - (n+.5) * dx ;
    fbw = BW(&yy, &p[1]) ;
    gau =  exp( -.5 * (x[0] - mean - yy)* (x[0] - mean - yy) / sigma2 )* ng/sigma ;
    f +=  fbw * gau ;
   
        }   
  double BWZ =   f * DeltaX / bins  ;
  double l = p[14];
  double a = p[15];
  double b = p[16];
  
  double N1 = exp (- l * fitXMax) / (l* l * l) * (-l * (l * fitXMax * a + a + l) - b * (l * l * fitXMax * fitXMax +  2 * l * fitXMax + 2 ));
  
  double N2 = exp (- l * fitXMin) / (l* l * l) * (-l * (l * fitXMin * a + a + l) - b * (l * l * fitXMin * fitXMin + 2 * l  * fitXMin + 2 )); 
  
  
  return  2 * p[0] * BWZ *2* (p[5]* p[5] * p[6]* (1.0-p[6])) +  p[13]/(N1-N2) * exp( -l * x[0] )  * (1 + x[0] * ( a + x[0] * b));
  
}

double BWZms (double *x, double * p) {
  // x[0] = E ; sp = E'^2 ; y = E'; H(E) = /int{BW(E') * g(E-E') dE'}
  //double Ni = p[4];
  //double Nf = p[3];  
  //double m = p[1];
  //double g = p[2];
  // double * p = { &m, &g, &Nf, &Ni};
  double sigma = p[11];
  double sigma2 = sigma*sigma;
  double mean = p[12];
  const double ng = 1.0 / sqrt(2*TMath::Pi());
  static const int bins = 100;
  double xMin = x[0] - 3. * sigma, xMax = x[0] + 3. * sigma, DeltaX = xMax - xMin;
  double dx =  DeltaX  / bins;
  double f = 0; double fbw =0; double gau =0;
  for( int n = 0; n < bins; ++n ) {
    /// down 
    double y = x[0] + (n+.5) * dx;
    fbw = BW(&y, &p[1]) ;
    gau = exp( -.5 * (x[0] - mean - y)*(x[0] -mean - y)  / sigma2 ) * ng/sigma    ;
    f +=  fbw * gau ;
    //up
    double yy = x[0] - (n+.5) * dx ;
    fbw = BW(&yy, &p[1]) ;
    gau =  exp( -.5 * (x[0] - mean - yy)* (x[0] - mean - yy) / sigma2 )* ng/sigma ;
    f +=  fbw * gau ;
    
  }   
  double BWZ =   f * DeltaX / bins  ;
  
  double l = p[22];
  double a = p[23];
  double b = p[24];
  
  double N1 = exp (- l * fitXMax) / (l* l * l) * (-l * (l * fitXMax * a + a + l) - b * (l * l * fitXMax * fitXMax +  2 * l * fitXMax + 2 ));
  
  double N2 = exp (- l * fitXMin) / (l* l * l) * (-l * (l * fitXMin * a + a + l) - b * (l * l * fitXMin * fitXMin + 2 * l  * fitXMin + 2 ));
  return   4*  p[0] * BWZ *2* (p[6]* p[6] * p[5]* (1.0-p[5])) 
+  p[21]/(N1-N2) * exp( -l * x[0] )  * (1 + x[0] * ( a + x[0] * b)) ;
  
  
  
}


double pol2(double *x,double *p){
  
  double l = p[1];
  double a = p[2];
  double b = p[3];

  double N1 = exp (- l * fitXMax) / (l* l * l) * (-l * (l * fitXMax * a + a + l) - b * (l * l * fitXMax * fitXMax +  2 * l * fitXMax + 2 ));
  
  double N2 = exp (- l * fitXMin) / (l* l * l) * (-l * (l * fitXMin * a + a + l) - b * (l * l * fitXMin * fitXMin + 2 * l  * fitXMin + 2 )); 
  
  return  p[0]/(N1-N2) * exp(-l* x[0]) * (1+ x[0]*( a + x[0] * b));
}



void fcn( int &, double *, double & f, double * par, int ) {
  static bool firstTime = true;
  static const size_t nBins = 2000; 
  static const size_t n = nBins, nn = nBins/nRebinZmt, nnn = nBins/nRebinZms;
  static double h[n] , dh[n];
  static double hh[nn], dhh[nn];
  static double hhh[nnn], dhhh[nnn];
  static double xMin, xMax, deltaX, ddeltaX, dddeltaX;
  if ( firstTime ) {
    firstTime = false;
    //    gROOT->ProcessLine(".L fit_function.C+");
 

    //Zmm senza back
 TFile * ZToLL_file1 = new TFile("test/EffTrackEffStandAloneHisto/ZMM_ZMT_ZMSA_HistonoIso_Spring07.root","read");
    TH1D * zToMuMu = (TH1D*) ZToLL_file1->Get("zToMM");
    TH1D * zToSingleTrackMu = (TH1D*) ZToLL_file1->Get("zToSingleTrackMu");
 TH1D * zToSingleStandAloneMu = (TH1D*) ZToLL_file1->Get("zToSingleStandAloneMu");
 zToSingleTrackMu->Sumw2();
 zToSingleTrackMu->Rebin(2);
 zToSingleStandAloneMu->Sumw2();
 zToSingleStandAloneMu->Rebin(4);
 
  
    xMin = zToMuMu->GetXaxis()->GetXmin();
    xMax = zToMuMu->GetXaxis()->GetXmax();
    double delta = xMax - xMin;
    deltaX = delta / n;
    ddeltaX = delta / nn;
    dddeltaX = delta / nnn;

      for( size_t i = 0; i < n; ++i ) {
      h[i] = zToMuMu->GetBinContent(i+1);
      dh[i] = zToMuMu->GetBinError(i+1);
      double x = xMin + ( i +.5 ) * deltaX;
      if( x > fitXMin && x < fitXMax && (h[i]>0) ) fullBins =fullBins++; 
    }
    
    for( size_t i = 0; i < nn; ++i ) {
      hh[i] = zToSingleTrackMu->GetBinContent(i+1);
      dhh[i] = zToSingleTrackMu->GetBinError(i+1);
      double x = xMin + ( i +.5 ) * ddeltaX;
      if( x > fitXMin && x < fitXMax &&  (hh[i]>0)) fullBins =fullBins++; 
    }

    
    for( size_t i = 0; i < nnn; ++i ) {
      hhh[i] = zToSingleStandAloneMu->GetBinContent(i+1);
      dhhh[i] = zToSingleStandAloneMu->GetBinError(i+1);
      double x = xMin + ( i  + .5) * dddeltaX;
      if( x > fitXMin && x < fitXMax &&  (hhh[i]>0)) fullBins =fullBins++; 
    }
    
  }
  static double var[1];
  double f1 =0, f2=0, f3=0;
  f=0; 
  cout << " f = " << f;
  for( size_t i = 0; i < n; ++ i ) {
    double x = xMin + ( i + .5 ) * deltaX;
    if ( x > fitXMin && x < fitXMax && h[i] > 0    ) {
      var[0] = x; 
      double l = BWZmm( var, par );
      double delta = l- h[i];
      // f += (delta * delta) / (l) ;    
      f1 += (delta * delta) / (dh[i]* dh[i]) ;
      
    }
  }
  cout << " --> " << f1;
  for( size_t i = 0; i < nn; ++ i ) {
    double x = xMin + ( i +.5 ) * ddeltaX;
    if ( x > fitXMin && x < fitXMax && hh[i]>0   ) {
      var[0] = x; 
      double ll = BWZmt( var, par ) ;
      double ddelta = ll - hh[i];
      
      f2 += (ddelta * ddelta) / ( 2 * dhh[i] * dhh[i]) ;
     
    }
  }
  cout << " --> " << f2;  
  for( size_t i = 0; i < nnn; ++ i ) {
    double x = xMin + ( i +.5) * dddeltaX;
    if ( x > fitXMin && x < fitXMax &&  hhh[i]>0 ) {
      var[0] = x; 
      double lll = BWZms( var, par );
      double dddelta = lll- hhh[i];
      f3 += (dddelta * dddelta) / ( 4 * dhhh[i] * dhhh[i] )  ;     
    }
  }
  cout << " --> " << f3 << endl;
  f = f1 + f2 + f3;
  cout << "f  === " << f << endl;
  
}
 
void zMinuit3Independentfit_conv_Spring07noIso() {
  
  cout << ">>> starting fit program" << endl;


  cout << ">>> fit init" << endl;
  
  unsigned int nPars = 25;
  TMinuit minuit( nPars );
  minuit.SetFCN( fcn );
  int ierflg = 0;
  double arglist[ 10 ];
  arglist[0] = -1;
  minuit.mnexcm( "SET PRINT", arglist, 1, ierflg );
  if ( ierflg != 0 ) cerr << "set print: AARGH!!" << endl;
  arglist[0] = 1;
  minuit.mnexcm( "SET ERR", arglist, 1, ierflg );
  if ( ierflg != 0 ) cerr << "set err: AARGH!!" << endl;
  // set parameters
  string names[] = { "P0", "mZmm", "GZmm", 
		     "f_gamma", "f_int", 
		     "eff_track","eff_standalone", "sigmaZmm", "meanZmm","sigmaZmt", "meanZmt", "sigmaZms", "meanZms","f_backZmt", "l_backZmt", "a_backZmt", "b_backZmt" , "f_backZmm", "l_backZmm", "a_backZmm", "b_backZmm"
		     , "f_backZms", "l_backZms", "a_backZms", "b_backZms"
};
  

  cout << ">>> set fit parms" << endl;

  minuit.mnparm( 0, names[0].c_str(),70000  ,1000 , 1000, 1000000, ierflg ); 
  
  //  minuit.FixParameter( 0 );
  minuit.mnparm( 1, names[1].c_str(), 91.3, .1, 70., 110, ierflg );
  // minuit.FixParameter( 1 );
  minuit.mnparm( 2, names[2].c_str(), 2.5, 1, 1, 20, ierflg ); 
  //minuit.FixParameter( 2 ); 
  minuit.mnparm( 3, names[3].c_str(), 0, 10, -100, 1000000, ierflg ); 
  minuit.FixParameter( 3 );
  minuit.mnparm( 4, names[4].c_str(),0, .01, -100, 100, ierflg ); 
  // minuit.FixParameter( 4 );
  minuit.mnparm( 5, names[5].c_str(), 0.98, .01, 0.,1, ierflg ); 
  //  minuit.FixParameter( 5 );
  minuit.mnparm( 6, names[6].c_str(), 0.97, .01, 0., 1,ierflg ); 
  //minuit.FixParameter( 6 );
  minuit.mnparm( 7, names[7].c_str(), 1., .1, -5, 5, ierflg ); 
  // minuit.FixParameter( 7);
  minuit.mnparm( 8, names[8].c_str(), 0, .001, -.5, .5, ierflg ); 
  minuit.FixParameter( 8 );


 minuit.mnparm( 9, names[9].c_str(), 1, .1, -3, 3, ierflg ); 
 //minuit.FixParameter( 9);
 minuit.mnparm( 10, names[10].c_str(), 1., .1, -10, 10, ierflg ); 
 //minuit.FixParameter( 10);
 minuit.mnparm( 11, names[11].c_str(), 5, .1, -15, 15, ierflg ); 
 minuit.mnparm( 12, names[12].c_str(), 0., .1, -1, 1, ierflg ); 
  minuit.mnparm( 13, names[13].c_str(), 8000, 100, 0, 100000, ierflg ); 
  minuit.mnparm( 14, names[14].c_str(),- 0.01,  .001, -1, 1, ierflg ); 
  minuit.mnparm( 15, names[15].c_str(), 1000, 1, 0., 10000, ierflg ); 
  minuit.mnparm( 16, names[16].c_str(), -5, 1, -1000, 1000, ierflg ); 

  minuit.mnparm( 17, names[17].c_str(), 500, 1, 30, 2000, ierflg ); 
  minuit.mnparm( 18, names[18].c_str(), 0.002,  .00001, 0.001, 1, ierflg ); 
  minuit.mnparm( 19, names[19].c_str(), 100, 1, 1, 10000, ierflg ); 
  minuit.mnparm( 20, names[20].c_str(), 0, .0001, -1, 1, ierflg ); 
  minuit.FixParameter( 20);

  minuit.mnparm( 21, names[21].c_str(), 100, 1, 0, 400, ierflg ); 
  // minuit.FixParameter( 21 );
  minuit.mnparm( 22, names[22].c_str(), 0.001,  .00001, 0, 1, ierflg ); 
  minuit.mnparm( 23, names[23].c_str(), 100, .1, -10000, 10000, ierflg ); 
  minuit.mnparm( 24, names[24].c_str(), 0, .001, -1, 1, ierflg );
  minuit.FixParameter( 24);
 
  
  cout << ">>> run fit" << endl;
  
  arglist[0] = 5000;
  arglist[1] = .01;
  minuit.mnexcm("MIGRAD",arglist, 2, ierflg );
  if ( ierflg != 0 ) cerr << "migrad: AARGH!!" << endl;
  // minuit.Migrad();
  // minuit.mnmnos();
  minuit.mnmatu(1);
  
  cout << ">>> fit completed" << endl;

  double amin;
  double edm, errdef;
  int nvpar, nparx;
  minuit.mnstat( amin, edm, errdef, nvpar, nparx, ierflg );  
  minuit.mnprin( 3, amin );
  unsigned int ndof = fullBins - minuit.GetNumFreePars();
  cout << "Chi-2 = " << amin << "/" << ndof << " = " << amin/ndof 
       << "; prob: " << TMath::Prob( amin, ndof )
       << endl;
  //minuit.mnmatu(1);
 
   
   TF1 * fitZmm = new TF1("fitZmm",BWZmm,fitXMin,fitXMax,nPars); 
   TF1 * fitZmt = new TF1("fitZmt",BWZmt,fitXMin,fitXMax,nPars); 
   TF1 * fitZms = new TF1("fitZms",BWZms,fitXMin,fitXMax,nPars); 
   fitZmm->SetLineColor(kRed);
   fitZmt->SetLineColor(kRed);
   fitZms->SetLineColor(kRed);
   TF1 * bkgZmm = new TF1("bkgZmm",pol2,fitXMin,fitXMax,4);
   TF1 * bkgZms = new TF1("bkgZms",pol2,fitXMin,fitXMax,4);
    TF1 * bkgZmt = new TF1("bkgZmt",pol2,fitXMin,fitXMax,4);
   
   
   gStyle->SetOptFit(1111);
   // gStyle->SetOptLogy();
   
   for( size_t i = 0; i < nPars; ++ i ) {
     double p, dp;
     minuit.GetParameter( i, p, dp );
     cout << names[i] << " = " << p << "+/-" << dp << endl;
   }
   
   for( size_t i = 0; i < nPars; ++ i ) {   
     double p, dp;
     minuit.GetParameter( i, p, dp );
     fitZmm->SetParameter(i,p);
     //bkgZmm->SetParameter(j,0.);
   }
   
   for( size_t i = 0; i < nPars; ++ i ) {   
     double p, dp;
     minuit.GetParameter( i, p, dp );
     fitZmt->SetParameter(i,p);
     //bkgZmt->SetParameter(i,0.);
   }
   
   for( size_t i = 0; i < nPars; ++ i ) {   
     double p, dp;
     minuit.GetParameter( i, p, dp );
     fitZms->SetParameter(i,p);
   }
    
  for( size_t i = 13; i < 17; ++ i ) {   
  double p, dp;
  minuit.GetParameter( i, p, dp );
  bkgZmt->SetParameter(i-13,p);
}
for( size_t i = 17; i < 21; ++ i ) {   
  double p, dp;
  minuit.GetParameter( i, p, dp );
  bkgZmm->SetParameter(i-17,p);
}

for( size_t i = 21; i < 25; ++ i ) {   
  double p, dp;
  minuit.GetParameter( i, p, dp );
  bkgZms->SetParameter(i-21,p);
}

TFile * ZToLL_file1 = new TFile("test/EffTrackEffStandAloneHisto/ZMM_ZMT_ZMSA_HistonoIso_Spring07.root","read");
    TH1D * zToMuMu = (TH1D*) ZToLL_file1->Get("zToMM");
    TH1D * zToSingleTrackMu = (TH1D*) ZToLL_file1->Get("zToSingleTrackMu");

  TH1D * zToSingleStandAloneMu = (TH1D*) ZToLL_file1->Get("zToSingleStandAloneMu");
  zToSingleTrackMu->Rebin(2);
  
  zToSingleStandAloneMu->Rebin(4);
   gStyle->SetOptLogy();  
   TCanvas c;
  
    zToMuMu->Draw("e");
    //zToMuMu->GetXaxis()->SetRangeUser(85,97);

    zToMuMu->GetXaxis()->SetRangeUser(fitXMin-10, fitXMax+10);
    zToMuMu->SetXTitle("#mu #mu mass (GeV/c^{2})");
    zToMuMu->SetTitle("fit of #mu #mu mass without isolation cut");
    zToMuMu->SetYTitle("Events / 0.1 GeV/c^{2}");


   /* 
  TH1D * h = new TH1D("h","h", 2000 , 80, 100);
   for (int i = 0; i < h->GetNbinsX(); i ++) {
     h->SetBinContent(i,zToMuMu->GetBinContent(8000 + i) );
     h->SetBinError(i,zToMuMu->GetBinError(8000 + i) ) ; 
      }
      // bkgZmm->Draw("same");
  
   //   fitZmm->GetXaxis()->SetLimits(fitXMin,fitXMax);
   h->Draw("e");
   */
   fitZmm->Draw("same");
   bkgZmm->Draw("same");
   c.SaveAs("zmm_noIso.eps");

 
   TCanvas c1;
   
   zToSingleTrackMu->Draw("e");

zToSingleTrackMu->GetXaxis()->SetRangeUser(fitXMin-20, fitXMax+20);
   zToSingleTrackMu->SetXTitle("#mu + (unmatched) track mass (GeV/c^{2})");
    zToSingleTrackMu->SetTitle("fit of #mu + (unmatched) track  mass without isolation cut");
    zToSingleTrackMu->SetYTitle("Events /  2 GeV/c^{2}");
   //  bkgZmt->Draw("same");
   fitZmt->Draw("same");
    bkgZmt->Draw("same");
  c1.SaveAs( "zmt_noIso.eps");
   
   TCanvas c2;
   
   zToSingleStandAloneMu->Draw("e");
zToSingleStandAloneMu->GetXaxis()->SetRangeUser(fitXMin-20, fitXMax+20);
    zToSingleStandAloneMu->SetXTitle("#mu + (unmatched) standalone mass (GeV/c^{2})");
    zToSingleStandAloneMu->SetTitle("fit of #mu + (unmatched) standalone  mass without isolation cut");
    zToSingleStandAloneMu->SetYTitle("Events / 4 GeV/c^{2}");
    bkgZms->Draw("same");
   fitZms->Draw("same");
   c2.SaveAs( "zms_noIso.eps");
   

   double Int = fitZmm->Integral(fitXMin,fitXMax);
   double IntHisto = zToMuMu->Integral((10*fitXMin)+1, (10*fitXMax)+1);

   cout <<"Integral of Zmm == "<<Int<<"vs Histo Integral "<<IntHisto<<endl;
   
   //   fitZmm->SetParameter(4,0);
   //double p4, dp4;
   //minuit.GetParameter( 4, p4, dp4 );
   
    double p0, dp0;
   minuit.GetParameter( 0, p0, dp0 );
   double p5, dp5;
   minuit.GetParameter( 5, p5, dp5 );
   double p6, dp6;
   minuit.GetParameter( 6, p6, dp6 );
   double p17, dp17;
   minuit.GetParameter( 17, p17, dp17 );

   // double Norm  =   Int/ p0;
   //  double Norm = 
   // cout <<"Norm == "<<Norm <<endl;
    
   fitZmm->SetParameter(17,0);
   double IntS = fitZmm->Integral(fitXMin ,fitXMax);
   cout <<"Integral of Zmm Signal == "<<IntS<<endl;

   double Norm  = (p0 * p5 * p5* p6* p6 + p17)/ IntHisto;
 
   cout <<"Norm == "<<Norm <<endl;
   
   double N0 = p0 /Norm ;
   double dN0p0 = dp0 ;
   double dN0p17 = dp17 ;

   double dN0p5 = dp5/ (p5*p5);
   double dN0p6 = dp6/ (p6*p6);

   double dNorm= TMath::Sqrt(dN0p0 *dN0p0 + dN0p17 * dN0p17 + dN0p5 * dN0p5 + dN0p6 * dN0p6); 
     
   double dN0 = TMath::Sqrt((dp0*dp0) + (dNorm * dNorm) / (Norm* Norm* Norm *Norm));
   cout <<"Stima di N0 == "<<N0<<" p0 "<< dN0  <<endl;
   TFile * ZToLL_file2 = new TFile("test/zSkimHisto/zSpring07withZMCHisto_utx_PtandEtaandMassCut.root","read");
    TH1D * ZMC = (TH1D*) ZToLL_file2->Get("ZAnalyzer/Z0MCmass"); 
    cout<< "---------vs MC value "<<ZMC->Integral((2 * fitXMin) +1 , (2 * fitXMax) + 1)<<endl;; 
    cout<<"ZMC integral (40,200) == "<<ZMC->Integral((2*40) +1 , (2 * 200) + 1)<<endl; 

}


int main() {
  //  gROOT->Reset();
  zMinuit3Independentfit_conv_Spring07noIso();
  return 0;
}



