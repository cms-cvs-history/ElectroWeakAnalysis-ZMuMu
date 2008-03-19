#include "PhysicsTools/Utilities/interface/Parameter.h"
#include "PhysicsTools/Utilities/interface/Constant.h"
#include "PhysicsTools/Utilities/interface/Number.h"
#include "PhysicsTools/Utilities/interface/Product.h"
#include "PhysicsTools/Utilities/interface/Sum.h"
#include "PhysicsTools/Utilities/interface/Difference.h"
#include "PhysicsTools/Utilities/interface/Gaussian.h"
#include "PhysicsTools/Utilities/interface/HistoChiSquare.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "PhysicsTools/Utilities/interface/RootFunctionAdapter.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include <boost/shared_ptr.hpp>
#include <boost/program_options.hpp>
using namespace boost;
namespace po = boost::program_options;

#include <iostream>
#include <algorithm> 
#include <iterator>
#include <string>
#include <vector>
using namespace std;
using namespace ::function;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " ")); 
    return os;
}

int main(int ac, char *av[]) { 
  try {
      double fMin, fMax;
      string ext;
      po::options_description desc("Allowed options");
      desc.add_options()
	("help", "produce help message")
	("include-path,I", po::value< vector<string> >(), 
	 "include path")
	("input-file", po::value< vector<string> >(), "input file")
	("min,m", po::value<double>(&fMin)->default_value(-20), "minimum value for fit range")
	("max,M", po::value<double>(&fMax)->default_value(20), "maximum value for fit range")
	("gauss", "fit to a gaussian")
	("2gauss", "fit to a linear combination of two gaussians")
	("3gauss", "fit to a linear combination of three gaussians")
	("output-file,O", po::value<string>(&ext)->default_value(".ps"), 
	 "output file format")
	;
      
      po::positional_options_description p;
      p.add("input-file", -1);
      
      po::variables_map vm;
      po::store(po::command_line_parser(ac, av).
	    options(desc).positional(p).run(), vm);
      po::notify(vm);
           
      if (vm.count("help")) {
	cout << "Usage: options_description [options]\n";
	cout << desc;
	return 0;
      }
      
      if (vm.count("include-path")) {
	cout << "Include paths are: " 
	     << vm["include-path"].as< vector<string> >() << "\n";
      }
      
      vector<string> v_file;
      vector<TH1D*> v_ZMassResHistos;
      vector<string> v_eps;
      
      if (vm.count("input-file")) {
	  cout << "Input files are: " 
	       << vm["input-file"].as< vector<string> >() << "\n";
	  v_file = vm["input-file"].as< vector<string> >();
	  for(vector<string>::const_iterator it = v_file.begin(); 
	      it != v_file.end(); ++it) { 
	    TFile * root_file = new TFile(it->c_str(),"read");
	    TDirectory *Histos = (TDirectory*) root_file->GetDirectory("ZHisto");
	    TDirectory *RecoHistos = (TDirectory*) Histos->GetDirectory("ZResolutionHisto");
	    TH1D * zMass = (TH1D*) RecoHistos->Get("ZToMuMuRecoMassResolution");
	    zMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})"); 
	    v_ZMassResHistos.push_back(zMass);
	    gROOT->SetStyle("Plain");
	    string f_string = *it;
	    replace(f_string.begin(), f_string.end(), '.', '_');
	    string eps_string = f_string + ext;
	    v_eps.push_back(eps_string);
	    cout << ">>> histogram loaded\n";
	  }
	  cout << v_file.size() << ", " << v_ZMassResHistos.size() << ", " << v_eps.size() << endl;
	  cout <<">>> Input files loaded\n";
	}
      //Parameters for fit
      Parameter yield("Yield", 10000);
      Parameter alpha("Alpha", 0);
      Parameter beta("Beta", 0);
      Parameter mean("Mean", 0);
      Parameter mean2("Mean 2", 0);
      Parameter mean3("Mean 3", 0);
      Parameter sigma1("Sigma 1", 1.);
      Parameter sigma2("Sigma 2", 2.5);
      Parameter sigma3("Sigma 3", 9.);
      Parameter dyield("Yield Error", 0);
      Parameter dalpha("Alpha Error", 0);
      Parameter dbeta("Beta Error", 0);
      Parameter dmean("Mean Error", 0);
      Parameter dmean2("Mean 2 Error", 0);
      Parameter dmean3("Mean 3 Error", 0);
      Parameter dsigma1("Sigma 1 Error", 0); 
      Parameter dsigma2("Sigma 2 Error", 0); 
      Parameter dsigma3("Sigma 3 Error", 0); 
      
      if (vm.count("gauss")) {
	cout << "Fitting histograms in input files to a Gaussian\n"; 
	cout << ">>> set pars: " << endl;
	cout << yield << " ; " << dyield << endl; 
	cout << mean << " ; " << dmean << endl; 
	cout << sigma1 << " ; " << dsigma1 << endl;
	for(size_t i = 0; i < v_ZMassResHistos.size(); ++i) { 
	  TH1D * zMass = v_ZMassResHistos[i]; 
	  zMass->Rebin(4); //remember...
	  Gaussian gaus(mean, sigma1);
	  Constant c(yield);
	  typedef Product<Constant, Gaussian> FitFunction;
	  FitFunction f = c * gaus;
	  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	  ChiSquared chi2(f, zMass, fMin, fMax);
	  int fullBins = chi2.degreesOfFreedom();
	  cout << "N. deg. of freedom: " << fullBins << endl;
	  fit::RootMinuit<ChiSquared> minuit(3, chi2, true);
	  minuit.setParameter(0, yield, 10, 100, 10000000);
	  minuit.setParameter(1, mean, 0.001, -1., 1.);
	  //minuit.fixParameter(1);
	  minuit.setParameter(2, sigma1, 0.1, -5., 5.);
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
	  dsigma1 = minuit.getParameterError(2);
	  cout << sigma1 << " ; " << dsigma1 << endl;
	  TF1 fun = root::tf1("fun", f, fMin, fMax, yield, mean, sigma1);
	  fun.SetParNames(yield.name().c_str(), mean.name().c_str(), sigma1.name().c_str());
	  fun.SetLineColor(kRed);
	  fun.SetNpx(1000);
	  TCanvas *canvas = new TCanvas("canvas");
	  zMass->Draw("e");
	  fun.Draw("same");	
	  string epsFilename = "ZMassResFitG_" + v_eps[i];
	  canvas->SaveAs(epsFilename.c_str());
	  canvas->SetLogy();
	  string epsLogFilename = "ZMassResFitG_Log_" + v_eps[i];
	  canvas->SaveAs(epsLogFilename.c_str());
	}
      }
      
      if (vm.count("2gauss")) {
	cout << "Fitting histograms in input files to a linear combination of two Gaussians\n"; 
	cout << "set pars: " << endl;
	cout << yield << " ; " << dyield << endl;
	cout << alpha << " ; " << dalpha << endl;
	cout << mean << " ; " << dmean << endl; 
	cout << mean2 << " ; " << dmean2 << endl; 
	cout << sigma1 << " ; " << dsigma1 << endl; 
	cout << sigma2 << " ; " << dsigma2 << endl; 
	for(size_t i = 0; i < v_ZMassResHistos.size(); ++i) { 
	  TH1D * zMass = v_ZMassResHistos[i]; 
	  zMass->Rebin(4); //remember...
	  Gaussian gaus1(mean, sigma1);
	  Gaussian gaus2(mean, sigma2);
	  typedef Product<Constant, Gaussian> G1;
	  typedef Product<Difference<Number,Constant>, Gaussian> G2;
	  typedef Product<Constant,Sum<G1, G2> > FitFunction;
	  Number _1(1);
	  FitFunction f = yield*(alpha*gaus1 + (_1 - alpha)*gaus2);
	  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	  ChiSquared chi2(f, zMass, fMin, fMax);
	  int fullBins = chi2.degreesOfFreedom();
	  cout << "N. deg. of freedom: " << fullBins << endl;
	  fit::RootMinuit<ChiSquared> minuit(5, chi2, true);
	  minuit.setParameter(0, yield, 10, 100, 10000000);
	  minuit.setParameter(1, alpha, 0.1, -1., 1.);
	  minuit.setParameter(2, mean, 0.001, -1., 1.);
	  //minuit.fixParameter(2);
	  minuit.setParameter(3, sigma1, 0.1, -5., 5.);
	  minuit.setParameter(4, sigma2, 0.1, -5., 5.);
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
	  dalpha = minuit.getParameterError(1);
	  cout << alpha << " ; " << dalpha << endl;
	  dmean = minuit.getParameterError(2);
	  cout << mean << " ; " << dmean << endl;
	  dsigma1 = minuit.getParameterError(3);
	  cout << sigma1 << " ; " << dsigma1 << endl;
	  dsigma2 = minuit.getParameterError(4);
	  cout << sigma2 << " ; " << dsigma2 << endl;
	  vector<shared_ptr<double> > pars;
	  pars.push_back(yield.ptr());
	  pars.push_back(alpha.ptr());
	  pars.push_back(mean.ptr());
	  pars.push_back(sigma1.ptr());
	  pars.push_back(sigma2.ptr());
	  TF1 fun = root::tf1("fun", f, fMin, fMax, pars);
	  fun.SetParNames(yield.name().c_str(), alpha.name().c_str(),
			  mean.name().c_str(), sigma1.name().c_str(), sigma2.name().c_str());
	  fun.SetLineColor(kRed);
	  fun.SetNpx(1000);
	  TCanvas *canvas = new TCanvas("canvas");
	  zMass->Draw("e");
	  fun.Draw("same");	
	  string epsFilename = "ZMassResFitGG_" + v_eps[i];
	  canvas->SaveAs(epsFilename.c_str());
	  canvas->SetLogy();
	  string epsLogFilename = "ZMassResFitGG_Log_" + v_eps[i];
	  canvas->SaveAs(epsLogFilename.c_str());
	}
      }
      
      if (vm.count("3gauss")) {
	cout << "Fitting histograms in input files to a linear combination of three Gaussians\n";
	cout << "set pars: " << endl;
	cout << yield << " ; " << dyield << endl;
	cout << alpha << " ; " << dalpha << endl;
	cout << beta << " ; " << dbeta << endl;
	cout << mean << " ; " << dmean << endl; 
	cout << mean2 << " ; " << dmean << endl; 
	cout << mean3 << " ; " << dmean3 << endl; 
	cout << sigma1 << " ; " << dsigma1 << endl; 
	cout << sigma2 << " ; " << dsigma2 << endl; 
	cout << sigma3 << " ; " << dsigma3 << endl; 
	for(size_t i = 0; i < v_ZMassResHistos.size(); ++i) { 
	  TH1D * zMass = v_ZMassResHistos[i]; 
	  zMass->Rebin(4); //remember...
	  Gaussian gaus1(mean, sigma1);
	  Gaussian gaus2(mean2, sigma2);
	  Gaussian gaus3(mean3, sigma3);
	  Constant a(alpha), b(beta), c(yield);
	  Number _1(1);
	  typedef Product<Constant, Gaussian> G1;
	  typedef Sum<G1, G1> SumG1;
	  typedef Sum<Constant, Constant> ConstSum;
	  typedef Product<Difference<Number,ConstSum>, Gaussian> G2;
	  typedef Product<Constant,Sum<SumG1, G2> > FitFunction;
	  FitFunction f = c*(a*gaus1 + b*gaus2 + (_1 - (a+b))*gaus3);
	  typedef fit::HistoChiSquare<FitFunction> ChiSquared;
	  ChiSquared chi2(f, zMass, fMin, fMax);
	  int fullBins = chi2.degreesOfFreedom();
	  cout << "N. deg. of freedom: " << fullBins << endl;
	  fit::RootMinuit<ChiSquared> minuit(9, chi2, true);
	  minuit.setParameter(0, yield, 10, 100, 10000000);
	  minuit.setParameter(1, alpha, 0.1, -1., 1.);
	  minuit.setParameter(2, beta, 0.1, -1., 1.);
	  minuit.setParameter(3, mean, 0.001, -1., 1.);
	  //minuit.fixParameter(3);
	  minuit.setParameter(4, mean2, 0.001, -1., 1.);
	  //minuit.fixParameter(4);
	  minuit.setParameter(5, mean3, 0.001, -1., 1.);
	  //minuit.fixParameter(5);
	  minuit.setParameter(6, sigma1, 0.1, -5., 5.);
	  minuit.setParameter(7, sigma2, 0.1, -5., 5.);
	  minuit.setParameter(8, sigma3, 0.1, -20., 20.);
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
	  dalpha = minuit.getParameterError(1);
	  cout << alpha << " ; " << dalpha << endl;
	  dbeta = minuit.getParameterError(2);
	  cout << beta << " ; " << dbeta << endl;
	  dmean = minuit.getParameterError(3);
	  cout << mean << " ; " << dmean << endl;
	  dmean2 = minuit.getParameterError(4);
	  cout << mean2 << " ; " << dmean2 << endl;
	  dmean3 = minuit.getParameterError(5);
	  cout << mean3 << " ; " << dmean3 << endl;
	  dsigma1 = minuit.getParameterError(6);
	  cout << sigma1 << " ; " << dsigma1 << endl;
	  dsigma2 = minuit.getParameterError(7);
	  cout << sigma2 << " ; " << dsigma2 << endl;
	  dsigma3 = minuit.getParameterError(8);
	  cout << sigma3 << " ; " << dsigma3 << endl;
	  vector<shared_ptr<double> > pars;
	  pars.push_back(yield.ptr());
	  pars.push_back(alpha.ptr());
	  pars.push_back(beta.ptr());
	  pars.push_back(mean.ptr());
	  pars.push_back(mean2.ptr());
	  pars.push_back(mean3.ptr());
	  pars.push_back(sigma1.ptr());
	  pars.push_back(sigma2.ptr());
	  pars.push_back(sigma3.ptr());
	  TF1 fun = root::tf1("fun", f, fMin, fMax, pars);
	  fun.SetParNames(yield.name().c_str(), alpha.name().c_str(), beta.name().c_str(), 
			  mean.name().c_str(), mean2.name().c_str(), mean3.name().c_str(), 
			  sigma1.name().c_str(), sigma2.name().c_str(), sigma3.name().c_str());
	  fun.SetLineColor(kRed);
	  //fun.SetNpx(100000);
	  TCanvas *canvas = new TCanvas("canvas");
	  zMass->Draw("e");
	  fun.Draw("same");	
	  string epsFilename = "ZMassResFitGGG_" + v_eps[i];
	  canvas->SaveAs(epsFilename.c_str());
	  canvas->SetLogy();
	  string epsLogFilename = "ZMassResFitGGG_Log_" + v_eps[i];
	  canvas->SaveAs(epsLogFilename.c_str());
	}
      }
      cout << "It works!\n";
  }
  catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }
  return 0;  
}