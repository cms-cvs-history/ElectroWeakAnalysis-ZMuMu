#include "PhysicsTools/Utilities/interface/Gaussian.h"
#include "PhysicsTools/Utilities/interface/HistoChiSquare.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "PhysicsTools/Utilities/interface/Parameter.h"
#include "PhysicsTools/Utilities/interface/Sum.h"
#include "PhysicsTools/Utilities/interface/Product.h"
#include "PhysicsTools/Utilities/interface/Constant.h"
#include "PhysicsTools/Utilities/interface/RootFunctionAdapter.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

using namespace std;
using namespace boost;

int main() { 
  TFile * ZToLL_file1 = new TFile("zMassRes.root","read");
  TH1D * zMass = (TH1D*) ZToLL_file1->Get("zMassRes");
  gROOT->SetStyle("Plain");
  function::Parameter yield1("Yield 1", 1000);
  function::Parameter yield2("Yield 2", 1000);
  function::Parameter mean1("Mean 1", 0);
  function::Parameter sigma1("Sigma 1", 1.);
  function::Parameter mean2("Mean 2", -1.);
  function::Parameter sigma2("Sigma 2", 1.);
  function::Parameter dyield1("Yield 1 Error", 0);
  function::Parameter dyield2("Yield 2 Error", 0);
  function::Parameter dmean1("Mean 1 Error", 0);
  function::Parameter dsigma1("Sigma 1 Error", 0); 
  function::Parameter dmean2("Mean 2 Error", 0);
  function::Parameter dsigma2("Sigma 2 Error", 0); 
  cout << "set pars: " << endl;
  cout << yield1 << " ; " << dyield1 << endl;
  cout << yield2 << " ; " << dyield2 << endl;
  cout << mean1 << " ; " << dmean1 << endl; 
  cout << sigma1 << " ; " << dsigma1 << endl; 
  cout << mean2 << " ; " << dmean2 << endl; 
  cout << sigma2 << " ; " << dsigma2 << endl; 
  function::Gaussian gaus1(mean1, sigma1);
  function::Gaussian gaus2(mean2, sigma2);
  function::Constant c1(yield1);
  function::Constant c2(yield2);
  typedef function::Product<function::Constant, function::Gaussian> ConstGaussian;
  ConstGaussian cog1 = c1 * gaus1;
  ConstGaussian cog2 = c2 * gaus2;
  typedef function::Sum<ConstGaussian, ConstGaussian> FitFunction;
  FitFunction f = cog1 + cog2;
  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
  ChiSquared chi2(f, zMass, -20, 20);
  int fullBins = chi2.degreesOfFreedom();
  cout << "N. deg. of freedom: " << fullBins << endl;
  fit::RootMinuit<ChiSquared> minuit(6, chi2, true);
  minuit.setParameter(0, yield1, 10, 100, 100000);
  minuit.setParameter(1, yield2, 10, 100, 100000);
  minuit.setParameter(2, mean1, 0.001, -10, 10);
  minuit.setParameter(3, sigma1, 0.1, -5., 5.);
  minuit.setParameter(4, mean2, 0.001, -11, 9.);
  minuit.setParameter(5, sigma2, 0.1, -5., 5.);
  double amin = minuit.minimize();
  cout << "fullBins = " << fullBins 
       << "; free pars = " << minuit.getNumberOfFreeParameters() 
       << endl;
  unsigned int ndof = fullBins - minuit.getNumberOfFreeParameters();
  cout << "Chi^2 = " << amin << "/" << ndof << " = " << amin/ndof 
       << "; prob: " << TMath::Prob( amin, ndof )
       << endl;
  dyield1 = minuit.getParameterError(0);
  cout << yield1 << " ; " << dyield1 << endl;
  dyield2 = minuit.getParameterError(1);
  cout << yield2 << " ; " << dyield2 << endl;
  dmean1 = minuit.getParameterError(2);
  cout << mean1 << " ; " << dmean1 << endl;
  dsigma1 = minuit.getParameterError(3);
  cout << sigma1 << " ; " << dsigma1 << endl;
  dmean2 = minuit.getParameterError(4);
  cout << mean2 << " ; " << dmean2 << endl;
  dsigma2 = minuit.getParameterError(5);
  cout << sigma2 << " ; " << dsigma2 << endl;
  vector<shared_ptr<double> > pars;
  pars.push_back(yield1.ptr());
  pars.push_back(yield2.ptr());
  pars.push_back(mean1.ptr());
  pars.push_back(sigma1.ptr());
  pars.push_back(mean2.ptr());
  pars.push_back(sigma2.ptr());
  TF1 fun = root::tf1("fun", f, -60, 60, pars);
  fun.SetParNames(yield1.name().c_str(), yield2.name().c_str(), 
		  mean1.name().c_str(), sigma1.name().c_str(), 
		  mean2.name().c_str(), sigma2.name().c_str());
  fun.SetLineColor(kRed);
  TCanvas canvas;
  zMass->Draw("e");
  fun.Draw("same");
  canvas.SaveAs("ZMassResGGHistoFunFit_zmumu_1.eps");
  return 0;
}
