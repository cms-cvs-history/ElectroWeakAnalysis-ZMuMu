#include "PhysicsTools/Utilities/interface/Constant.h"
#include "PhysicsTools/Utilities/interface/Parameter.h"
#include "PhysicsTools/Utilities/interface/Product.h"
#include "PhysicsTools/Utilities/interface/Sum.h"
#include "PhysicsTools/Utilities/interface/Convolution.h"
#include "PhysicsTools/Utilities/interface/BreitWigner.h"
#include "PhysicsTools/Utilities/interface/Gaussian.h"
#include "PhysicsTools/Utilities/interface/HistoChiSquare.h"
#include "PhysicsTools/Utilities/interface/RootFunctionAdapter.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost;

int main() { 
  TFile * ZToLL_file1 = new TFile("out.root","read");
  TDirectory *Histos = (TDirectory*) ZToLL_file1->GetDirectory("zMuMuAnalyzer");
  TH1D * zMass = (TH1D*) Histos->Get("ZMuMumass");
  gROOT->SetStyle("Plain");
  cout << ">>> histo loaded" << endl;
  function::Parameter yield1("Yield 1", 600);
  function::Parameter yield2("Yield 2", 220);
  function::Parameter mass("Mass", 91.2);
  function::Parameter gamma("Gamma", 2.5);
  function::Parameter mean1("Mean 1", 0.08);
  function::Parameter sigma1("Sigma 1", 1.11);
  function::Parameter mean2("Mean 2", -1.);
  function::Parameter sigma2("Sigma 2", 3.9);
  cout << ">>> set pars: " << endl;
  cout << yield1 << endl; 
  cout << yield2 << endl; 
  cout << mass  << endl; 
  cout << gamma << endl; 
  cout << mean1 << endl; 
  cout << sigma1 << endl; 
  cout << mean2  << endl; 
  cout << sigma2 << endl; 
  function::Constant c1(yield1);
  function::Constant c2(yield2);
  function::BreitWigner bw(mass, gamma);
  function::Gaussian gaus1(mean1, sigma1);
  function::Gaussian gaus2(mean2, sigma2);
  double range = 3 * max(sigma1.value(), sigma2.value());
  typedef function::Convolution<function::BreitWigner, function::Gaussian> ConvBwG;
  ConvBwG cbg1(bw, gaus1, -range , range, 200);
  ConvBwG cbg2(bw, gaus2, -range , range, 200);
  typedef function::Product<function::Constant, ConvBwG> ConstConvBwG;
  ConstConvBwG cocobg1 = c1 * cbg1;
  ConstConvBwG cocobg2 = c2 * cbg2;
  typedef function::Sum<ConstConvBwG, ConstConvBwG> FitFunction;
  FitFunction f = cocobg1 + cocobg2;
  cout << "set functions" << endl;
  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
  ChiSquared chi2(f, zMass, 80, 120);
  int fullBins = chi2.degreesOfFreedom();
  cout << "N. deg. of freedom: " << fullBins << endl;
  fit::RootMinuit<ChiSquared> minuit(chi2, true);
  minuit.addParameter(yield1, 10, 100, 100000);
  minuit.addParameter(yield2, 10, 100, 100000);
  minuit.addParameter(mass, .1, 70., 110);
  minuit.addParameter(gamma, 1, 1, 10);
  minuit.addParameter(mean1, 0.001, -0.5, 0.5);
  minuit.fixParameter(mean1.name());
  minuit.addParameter(sigma1, 0.1, -5., 5.);
  minuit.fixParameter(sigma1.name());
  minuit.addParameter(mean2, 0.001, -5., 5.);
  minuit.fixParameter(mean2.name());
  minuit.addParameter(sigma2, 0.1, -5., 5.);
  minuit.fixParameter(sigma2.name());
  minuit.minimize();
  minuit.printFitResults();
  vector<shared_ptr<double> > pars;
  pars.push_back(yield1.ptr());
  pars.push_back(yield2.ptr());
  pars.push_back(mass.ptr());
  pars.push_back(gamma.ptr());
  pars.push_back(mean1.ptr());
  pars.push_back(sigma1.ptr());
  pars.push_back(mean2.ptr());
  pars.push_back(sigma2.ptr());
  TF1 fun = root::tf1("fun", f, 0, 200, pars);
  fun.SetParNames(yield1.name().c_str(), yield2.name().c_str(), 
		  mass.name().c_str(), gamma.name().c_str(), 
 		  mean1.name().c_str(), sigma1.name().c_str(), 
		  mean2.name().c_str(), sigma2.name().c_str());
  fun.SetLineColor(kRed);
  TCanvas canvas;
  zMass->Draw("e");
  fun.Draw("same");
  canvas.SaveAs("ZMassRecoCoBwGGHistoFunFit_out.ps");
  return 0;
}
