#include "PhysicsTools/Utilities/interface/Gaussian.h"
#include "PhysicsTools/Utilities/interface/HistoChiSquare.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "PhysicsTools/Utilities/interface/Parameter.h"
#include "PhysicsTools/Utilities/interface/Product.h"
#include "PhysicsTools/Utilities/interface/Constant.h"
#include "PhysicsTools/Utilities/interface/RootFunctionAdapter.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

using namespace std;
using namespace boost;

int main() { 
  TFile * ZToLL_file1 = new TFile("histoZ_zmumu_1.root","read");
  TDirectory *Histos = (TDirectory*) ZToLL_file1->GetDirectory("ZHisto");
  TDirectory *MCHistos = (TDirectory*) Histos->GetDirectory("ZMCHisto");
  TH1D * zMass = (TH1D*) MCHistos->Get("ZMCMass");
  gROOT->SetStyle("Plain");
  function::Parameter yield("Yield", 1000);
  function::Parameter mean("Mean", 91.2);
  function::Parameter sigma("Sigma", 2.5);
  function::Parameter dyield("Yield Error", 0);
  function::Parameter dmean("Mean Error", 0);
  function::Parameter dsigma("Sigma Error", 0); 
  cout << "set pars: " << endl;
  cout << yield << " ; " << dyield << endl; 
  cout << mean << " ; " << dmean << endl; 
  cout << sigma << " ; " << dsigma << endl; 
  function::Gaussian gaus(mean, sigma);
  function::Constant c(yield);
  typedef function::Product<function::Constant, function::Gaussian> FitFunction;
  FitFunction f = c * gaus;
  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
  ChiSquared chi2(f, zMass, 80, 120);
  int fullBins = chi2.degreesOfFreedom();
  cout << "N. deg. of freedom: " << fullBins << endl;
  fit::RootMinuit<ChiSquared> minuit(3, chi2, true);
  minuit.setParameter(0, yield, 10, 100, 100000);
  minuit.setParameter(1, mean, 0.001, 80, 100);
  minuit.setParameter(2, sigma, 0.1, -5., 5.);
  double amin = minuit.minimize();
  cout << "fullBins = " << fullBins 
       << "; free pars = " << minuit.getNumberOfFreeParameters() 
       << endl;
  unsigned int ndof = fullBins - minuit.getNumberOfFreeParameters();
  cout << "Chi^2 = " << amin << "/" << ndof << " = " << amin/ndof 
       << "; prob: " << TMath::Prob( amin, ndof )
       << endl;
  dyield = minuit.getParameterError(0);
  cout << yield << " ; " << dyield << endl;
  dmean = minuit.getParameterError(1);
  cout << mean << " ; " << dmean << endl;
  dsigma = minuit.getParameterError(2);
  cout << sigma << " ; " << dsigma << endl;
  TF1 fun = root::tf1("fun", f, 0, 200, yield, mean, sigma);
  fun.SetParameters(yield.value(), mean.value(), sigma.value());
  fun.SetParNames(yield.name().c_str(), mean.name().c_str(), sigma.name().c_str());
  fun.SetLineColor(kRed);
  TCanvas canvas;
  zMass->Draw("e");
  fun.Draw("same");
  canvas.SaveAs("ZMassMCGaussHistoFunFit_zmumu_1.eps");
  return 0;
}
