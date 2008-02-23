#include "PhysicsTools/Utilities/interface/ZLineShape.h"
#include "PhysicsTools/Utilities/interface/HistoChiSquare.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "PhysicsTools/Utilities/interface/Parameter.h"
#include "PhysicsTools/Utilities/interface/Product.h"
#include "PhysicsTools/Utilities/interface/Constant.h"
#include "PhysicsTools/Utilities/interface/RootFunctionAdapter.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost;

int main() { 
  TFile * ZToLL_file1 = new TFile("histoZMass.root","read");
  TDirectory *Histos = (TDirectory*) ZToLL_file1->GetDirectory("ZMassHisto");
  TH1D * zMass = (TH1D*) Histos->Get("mZ");
  gROOT->SetStyle("Plain");
  cout << "histo loaded" << endl;
  function::Parameter yield("Yield", 500);
  function::Parameter mass("Mass", 91.2);
  function::Parameter gamma("Gamma", 2.5);
  function::Parameter f_gamma("Photon factor", 0);
  function::Parameter f_int("Interference factor", 0.001);
  function::Parameter dyield("Yield Error", 0);
  function::Parameter dmass("Mass Error", 0);
  function::Parameter dgamma("Gamma Error", 0);
  function::Parameter df_gamma("Photon factor Error", 0);
  function::Parameter df_int("Interference factor Error", 0);
  cout << "set pars: " << endl;
  cout << yield << " ; " << dyield << endl; 
  cout << mass << " ; " << dmass << endl; 
  cout << gamma << " ; " << dgamma << endl; 
  cout << f_gamma << " ; " << df_gamma << endl; 
  cout << f_int << " ; " << df_int << endl; 
  function::ZLineShape zls(mass, gamma, f_gamma, f_int);
  function::Constant c(yield);
  typedef function::Product<function::Constant, function::ZLineShape> FitFunction;
  FitFunction f = c * zls;
  cout << "set functions" << endl;
  vector<shared_ptr<double> > pars;
  pars.push_back(yield.ptr());
  pars.push_back(mass.ptr());
  pars.push_back(gamma.ptr());
  pars.push_back(f_gamma.ptr());
  pars.push_back(f_int.ptr());
  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
  ChiSquared chi2(f, zMass, 80, 120);
  int fullBins = chi2.degreesOfFreedom();
  cout << "N. deg. of freedom: " << fullBins << endl;
  fit::RootMinuit<ChiSquared> minuit(5, chi2, true);
  minuit.setParameter(0, yield, 10, 100, 100000);
  minuit.setParameter(1, mass, .1, 70., 110);
  minuit.setParameter(2, gamma, 1, 1, 10);
  minuit.setParameter(3, f_gamma, 0.1, -100, 1000);
  minuit.fixParameter(3);
  minuit.setParameter(4, f_int, .0001, -1000000, 1000000);
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
  dmass = minuit.getParameterError(1);
  cout << mass << " ; " << dmass << endl;
  dgamma = minuit.getParameterError(2);
  cout << gamma << " ; " << dgamma << endl;
  df_gamma = minuit.getParameterError(3);
  cout << f_gamma << " ; " << df_gamma << endl;
  df_int = minuit.getParameterError(4);
  cout << f_int << " ; " << df_int << endl;
  TF1 fun = root::tf1("fun", f, 0, 200, pars);
  fun.SetParNames(yield.name().c_str(), mass.name().c_str(), gamma.name().c_str(), 
		  f_gamma.name().c_str(), f_int.name().c_str());
  fun.SetLineColor(kRed);
  TCanvas canvas;
  zMass->Draw("e");
  fun.Draw("same");
  canvas.SaveAs("ZMassRecoZLSHistoFunFit.eps");
  return 0;
}
