#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooChi2Var.h"
#include "RooGlobalFunc.h" 
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "RooMCStudy.h"
#include "RooRandom.h"

using namespace RooFit;

// A function that get histogram, restricting it on fit range  and sets contents to 0.01 if entries are too small


TH1F * getHisto(TFile * file, const char * name, double fMin, double fMax, unsigned int rebin) {
  TObject * h = file->Get(name);
  if(h == 0)
    // throw edm::Exception(edm::errors::Configuration) 
    cout  << "Can't find object " << name << "\n";
  TH1F * histo = dynamic_cast<TH1F*>(h);
  // histo->Sumw2();
  //throw edm::Exception(edm::errors::Configuration) 
  if(histo == 0)
    cout << "Object " << name << " is of type " << h->ClassName() << ", not TH1\n";
  
  TH1F * new_histo = new TH1F(name, name, (int) (fMax-fMin), fMin, fMax);
  //new_histo->Sumw2();
  int bin_num=0;
  for (int i = (int)fMin; i <= (int)fMax; ++i ) {
    bin_num= (i - (int)fMin + 1);  
    new_histo->SetBinContent( bin_num, histo->GetBinContent(i) );				    
  } 
  delete histo;
  new_histo->Rebin(rebin);
  

   for(int i = 1; i <= new_histo->GetNbinsX(); ++i) {
if(new_histo->GetBinContent(i) == 0.00) {
  cout<< " WARNING: histo " << name << " has 0 enter in bin number " << i << endl;   
   }
   }

//   for(int i = 1; i <= new_histo->GetNbinsX(); ++i) {
//     if(new_histo->GetBinContent(i) < 0.0001) {
//       new_histo->SetBinContent(i, 0.0001);
//      new_histo->SetBinError(i, 0.0001);
//     }
//   }


return new_histo;
}
 
// a function that create fromm a model pdf a toy RooDataHist
 
RooDataHist * genHistFromModelPdf(const char * name, RooAbsPdf *model,  RooRealVar *var,  unsigned int seed,  double ScaleLumi,  int rebin ) {
  //var->setBins(60 / rebin); 
  // adjust rebin to range.....
  RooRandom::randomGenerator()->SetSeed(34 + seed); 
  double genEvents =  ScaleLumi * model->expectedEvents(*var);
  RooDataSet * data = model->generate(*var ,  (int) genEvents ,    Extended(kTRUE) );
  cout<< " data->numEntries() for name " << name << " == " << data->numEntries()<< endl;
 
   double integ = model->expectedEvents(*var);
   cout<< " expected events for " << name << " = "<< integ << endl; 
  RooAbsData *binned_data = data->binnedClone();

  TH1 * toy_hist = binned_data->createHistogram( name, *var, Binning(rebin)  );
 RooDataHist * toy_rooHist = new RooDataHist(name, name , RooArgList(*var), toy_hist );
 
 
 return toy_rooHist; 
}


// funtion to create the pdf used for the fit, need the histo model, should be zmm except for zmusta case..... 
RooHistPdf * createHistPdf( const char * name, TH1F *model,  RooRealVar *var,  int rebin){
  TH1F * model_clone = new TH1F(*model);
  model_clone-> Rebin( rebin ); 
  RooDataHist * model_dataHist = new RooDataHist(name, name , RooArgList(*var), model_clone ); 
  RooHistPdf * HistPdf = new RooHistPdf(name, name, *var, *model_dataHist,  0);
  //delete  model_clone;
  //delete  model_dataHist;
  return HistPdf;
}







int main() {
  TCanvas  c1;
  c1.SetLogy(kTRUE);
  // non funge!!!!!
  c1.Divide(2,3);
  
  gROOT->SetBatch(kTRUE);
   gROOT->SetStyle("Plain");
  
  
  double fMin = 60.;
  double fMax = 120.;
  
  
  RooRealVar * mass = new RooRealVar("mass", "mass (GeV)", fMin, fMax );
  //  mass->setRange("range", fMin,fMax);
  
  TFile * root_file = new TFile("Analisi_45pb.root", "read");
  TFile * out_root_file = new TFile("Fit_Result.root", "recreate");
  
  int rebinZMuMu =1, rebinZMuSa = 12/*(int) (fMax-fMin)*/,   rebinZMuTk = 1,  rebinZMuMuNoIso=1,  rebinZMuMuHlt =  (int)(fMax - fMin); 
  int range = (int)(fMax - fMin);

  int numberOfBins = range/rebinZMuSa  + range/rebinZMuTk + range/rebinZMuMuNoIso + 2* range/rebinZMuMuHlt ;

  // zmm histograms used for pdf 

  TH1F * zmm = getHisto(root_file, "goodZToMuMuPlots/zMass",fMin, fMax, rebinZMuMu );  
  
  // zmsta used for pdf
  TH1F *zmsta = getHisto(root_file, "zmumuSaMassHistogram/zMass",fMin, fMax, 1);
 
  // histogramms to fit.....
 TH1F *zmmsta = getHisto(root_file, "goodZToMuMuOneStandAloneMuonPlots/zMass",fMin, fMax,  rebinZMuSa/rebinZMuMu);  

  TH1F *zmt = getHisto(root_file, "goodZToMuMuOneTrackPlots/zMass", fMin, fMax, rebinZMuTk / rebinZMuMu);  
  TH1F *zmmNotIso = getHisto(root_file,"nonIsolatedZToMuMuPlots/zMass",  fMin, fMax, rebinZMuMuNoIso / rebinZMuMu); 
  TH1F *zmm1hlt = getHisto(root_file, "goodZToMuMu1HLTPlots/zMass", fMin, fMax, rebinZMuMuHlt / rebinZMuMu);   
  TH1F *zmm2hlt = getHisto(root_file, "goodZToMuMu2HLTPlots/zMass", fMin, fMax, rebinZMuMuHlt / rebinZMuMu);   
  

// creating a pdf for Zmt
  
  /*TH1F *zmm_forZMuTk = new TH1F(*zmm);
  zmm_forZMuTk-> Rebin( rebinZMuTk/ rebinZMuMu ); 
  RooDataHist * zMass_forZMuTk = new RooDataHist("zMass_forZMuTk", "good Z mass _forZMuTk" , RooArgList(*mass), zmm_forZMuTk ); 
  RooHistPdf * ZmtPdf = new RooHistPdf("ZmtPdf", "Zmt pdf from histogram", *mass, *zMass_forZMuTk,  0); */
  
  RooHistPdf * ZmtPdf = createHistPdf( "ZmtPdf", zmm, mass, rebinZMuTk/ rebinZMuMu  );

  // creating a pdf for Zmm not iso
   RooHistPdf * ZmmNoIsoPdf = createHistPdf( "ZmmNoIsoPdf", zmm, mass, rebinZMuMuNoIso/ rebinZMuMu  );
  
// creating a pdf for Zms from zmsta!!!
  RooHistPdf * ZmsPdf = createHistPdf( "ZmsPdf", zmsta, mass, rebinZMuSa/ rebinZMuMu  );

  // creating a pdf for Zmmhlt
  RooHistPdf * ZmmHltPdf = createHistPdf( "ZmmHltPdf", zmm, mass, rebinZMuMuHlt/ rebinZMuMu  );



 // fit parameters

RooRealVar *a0 = new RooRealVar("a0","coefficient 0", 1,0.,2.) ;
//a0->setConstant(kTRUE) ;
RooRealVar *a1 = new RooRealVar("a1","coefficient 1", -0.001,-100,100.) ;
//a1->setConstant(kTRUE) ;
RooRealVar *a2 = new RooRealVar("a2","coefficient 2", 0.0,-100.,100.) ;
//a2->setConstant(kTRUE) ;

RooRealVar *alpha = new RooRealVar("alpha","coefficient alpha", -0.01 , -1., 0.) ;
//alpha->setConstant(kTRUE) ;
RooRealVar Yield("Yield","Yield",15000.,0.,10000000.) ;
//Yield.setConstant(kTRUE);

RooRealVar nbkg_mutrk("nbkg_mutrk","background _mutrk fraction",500,0.,100000.) ;

RooRealVar *b0 = new RooRealVar("b0","coefficient 0", 1,0.,2.) ;
//b0->setConstant(kTRUE) ;
RooRealVar *b1 = new RooRealVar("b1","coefficient 1", -0.001,-100,100.) ;
//b1->setConstant(kTRUE) ;
RooRealVar *b2 = new RooRealVar("b2","coefficient 2", 0.0,-100.,100.) ;
//b2->setConstant(kTRUE) ;

 RooRealVar *beta = new RooRealVar("beta","coefficient beta", -0.01,-1.0 , 0. ) ;
 //beta->setConstant(kTRUE) ;

 RooRealVar nbkg_mumuNotIso("nbkg_mumuNotIso","background fraction",500,0.,100000.) ;
 RooRealVar eff_tk("eff_tk","signal _mutrk fraction",.99,0.8,1.0) ;
 //eff_tk.setConstant(kTRUE) ;   
RooRealVar eff_iso("eff_iso","eff mumuNotIso",.99,0.8,1.0) ;
//eff_iso.setConstant(kTRUE) ;  
RooRealVar eff_sa("eff_sa","eff musta",0.99,0.8,1.0) ;
//eff_sa.setConstant(kTRUE) ; 
RooRealVar eff_hlt("eff_hlt","eff 1hlt",0.99, 0.8,1.0) ;
//eff_hlt.setConstant(kTRUE) ;
 


 //efficiency term

 //zMuMuEff1HLTTerm = _2 * (effTk ^ _2) *  (effSa ^ _2) * (effIso ^ _2) * effHLT * (_1 - effHLT); 
 RooFormulaVar zMuMu1HLTEffTerm("zMuMu1HLTEffTerm", " @4 * (2.* (@0)^2 * (@1)^2 * (@2)^2 * @3 *(1.- @3))", RooArgList(eff_tk, eff_sa,eff_iso, eff_hlt, Yield));
 
//zMuMuEff2HLTTerm = (effTk ^ _2) *  (effSa ^ _2) * (effIso ^ _2) * (effHLT ^ _2) ; 
 RooFormulaVar zMuMu2HLTEffTerm("zMuMu2HLTEffTerm", " @4 * ((@0)^2 * (@1)^2 * (@2)^2 * (@3)^2)", RooArgList(eff_tk, eff_sa,eff_iso, eff_hlt, Yield));

 //zMuMuNoIsoEffTerm = (effTk ^ _2) * (effSa ^ _2) * (_1 - (effIso ^ _2)) * (_1 - ((_1 - effHLT)^_2));
 RooFormulaVar zMuMuNoIsoEffTerm("zMuMuNoIsoEffTerm", " @4 * ((@0)^2 * (@1)^2 * (1.- (@2)^2) * (1.- ((1.-@3)^2)))", RooArgList(eff_tk, eff_sa,eff_iso, eff_hlt, Yield));

 //zMuTkEffTerm = _2 * (effTk ^ _2) * effSa * (_1 - effSa) * (effIso ^ _2) * effHLT;
 RooFormulaVar zMuTkEffTerm("zMuTkEffTerm", " @4 * (2. *(@0)^2 * @1 * (1.-@1)* (@2)^2 * @3)", RooArgList(eff_tk, eff_sa,eff_iso, eff_hlt, Yield));

//zMuSaEffTerm = _2 * (effSa ^ _2) * effTk * (_1 - effTk) * (effIso ^ _2) * effHLT;
 RooFormulaVar zMuSaEffTerm("zMuSaEffTerm", "@4 * (2. *(@1)^2 * @0 * (1.-@0)* (@2)^2 * @3)", RooArgList(eff_tk, eff_sa,eff_iso, eff_hlt, Yield));


 // creating model for the  fit 
// z mu track


  RooPlot * massFrame_mutrk = mass->frame() ;
 RooDataHist * zmtMass = new RooDataHist("zmtMass", "good z mu track" , RooArgList(*mass), zmt );
  RooGenericPdf *bkg_mutrk = new RooGenericPdf("bkg_mutrk","zmt bkg_model", "exp(alpha*mass) * ( a0 + a1 * mass + a2 * mass^2 )", RooArgSet( *mass, *alpha, *a0, *a1,*a2) );
 RooFormulaVar fracSigMutrk("fracSigMutrk", "@0 / (@0 + @1)", RooArgList(zMuTkEffTerm, nbkg_mutrk ));
 RooAddPdf * model_mutrk= new RooAddPdf("model_mutrk","model_mutrk",RooArgList(*ZmtPdf,*bkg_mutrk),RooArgList(zMuTkEffTerm , nbkg_mutrk)) ;

 


// z mu mu not Iso
 RooPlot * massFrame_mumuNotIso = mass->frame() ;
 RooDataHist * zmmNotIsoMass = new RooDataHist("ZmmNotIso", "good z mu mu not isolated" , RooArgList(*mass), zmmNotIso );
 // creating background pdf for zmu mu not Iso
  RooGenericPdf *bkg_mumuNotIso = new RooGenericPdf("bkg_mumuNotIso","zmumuNotIso bkg_model", "exp(beta * mass) * (b0 + b1 * mass + b2 * mass^2)", RooArgSet( *mass, *beta, *b0, *b1,*b2) );
RooFormulaVar fracSigMuMuNoIso("fracSigMuMuNoIso", "@0 / (@0 + @1)", RooArgList(zMuMuNoIsoEffTerm, nbkg_mumuNotIso ));
 RooAddPdf * model_mumuNotIso= new RooAddPdf("model_mumuNotIso","model_mumuNotIso",RooArgList(*ZmmNoIsoPdf,*bkg_mumuNotIso), RooArgList(zMuMuNoIsoEffTerm, nbkg_mumuNotIso)); 

// z mu sta
 RooPlot * massFrame_musta = mass->frame() ;
 RooDataHist * zmsMass = new RooDataHist("zmsMass", "good z mu sta mass" , RooArgList(*mass), zmmsta );
 RooGenericPdf model_musta("model_musta",  " ZmsPdf * zMuSaEffTerm ", RooArgSet( *ZmsPdf, zMuSaEffTerm)) ;
 RooExtendPdf * eZmsSig= new RooExtendPdf("eZmsSig","extended signal p.d.f for zms ",*ZmsPdf,  zMuSaEffTerm ) ;
 
 // z mu mu HLT

RooPlot * massFrame_mumu1hlt = mass->frame() ;
RooPlot * massFrame_mumu2hlt = mass->frame() ;
 // count ZMuMu Yield
 double nZMuMu = 0.;
 double nZMuMu1HLT = 0.;
 double  nZMuMu2HLT = 0.;
 unsigned int nBins = zmm->GetNbinsX();
 double xMin = zmm->GetXaxis()->GetXmin();
 double xMax = zmm->GetXaxis()->GetXmax();
 double deltaX =(xMax - xMin) / nBins;
 for(size_t i = 0; i < nBins; ++i) { 
   
   double x = xMin + (i +.5) * deltaX;
   if(x > fMin && x < fMax){
    nZMuMu += zmm->GetBinContent(i+1);
  
    nZMuMu1HLT += zmm1hlt->GetBinContent(i+1);
  
    nZMuMu2HLT += zmm2hlt->GetBinContent(i+1);
   }
 }

cout << ">>> count of ZMuMu yield in the range [" << fMin << ", " << fMax << "]: " << nZMuMu << endl;
cout << ">>> count of ZMuMu (1HLT) yield in the range [" << fMin << ", " << fMax << "]: "  << nZMuMu1HLT << endl;
cout << ">>> count of ZMuMu (2HLT) yield in the range [" << fMin << ", " << fMax << "]: "  << nZMuMu2HLT << endl;
// we set eff_hlt

 eff_hlt.setVal( 1. / (1. + (nZMuMu1HLT/ (2 * nZMuMu2HLT))) ) ;
 //eff_hlt.setConstant(kTRUE);
 // starting point for the corrected Yield near ZMuMu Yield 

 Yield.setVal(nZMuMu);

 // creating the pdf for z mu mu 1hlt
 RooDataHist * HistoZCount1HLT = new RooDataHist("HistoZCount1HLT", "good Z mu mu 1hlt" , RooArgList(*mass),   zmm1hlt ); 
 RooExtendPdf * eZmm1hltSig= new RooExtendPdf("eZmm1hltSig","extended signal p.d.f for zmm 1hlt",*ZmmHltPdf,  zMuMu1HLTEffTerm ) ;



// z mu mu 2HLT
 // creating the pdf for z mu mu 2hlt 
 RooDataHist * HistoZCount2HLT = new RooDataHist("HistoZCount2HLT", "good Z mu mu 2hlt" , RooArgList(*mass),   zmm2hlt ); 
RooExtendPdf * eZmm2hltSig= new RooExtendPdf("eZmm2hltSig","extended signal p.d.f for zmm 2hlt",*ZmmHltPdf,  zMuMu2HLTEffTerm ) ;
 

 // creting the chi2s 
RooChi2Var *chi2_mutrk = new RooChi2Var("chi2_mutrk","chi2_mutrk",*model_mutrk,*zmtMass, Extended(kTRUE)) ;

 RooChi2Var *chi2_mumuNotIso = new RooChi2Var("chi2_mumuNotIso","chi2_mumuNotIso",*model_mumuNotIso,*zmmNotIsoMass,Extended(kTRUE)) ;

 RooChi2Var *chi2_musta = new RooChi2Var("chi2_musta","chi2_musta",*eZmsSig, *zmsMass,  Extended(kTRUE) ) ;

 RooChi2Var *chi2_mu1hlt = new RooChi2Var("chi2_mu1hlt","chi2_mu1hlt", *eZmm1hltSig, *HistoZCount1HLT,  Extended(kTRUE) ) ;
 RooChi2Var *chi2_mu2hlt = new RooChi2Var("chi2_mu2hlt","chi2_mu2hlt", *eZmm2hltSig, *HistoZCount2HLT,  Extended(kTRUE) ) ;


 // adding the chi2
 RooAddition totChi2("totChi2","chi2_mutrk + chi2_mumuNotIso + chi2_musta   + chi2_mu1hlt  +  chi2_mu2hlt " , RooArgSet(*chi2_mutrk   ,*chi2_mumuNotIso  , *chi2_musta  ,*chi2_mu1hlt  ,  *chi2_mu2hlt  )) ;

 RooMinuit m_tot(totChi2) ;
 m_tot.migrad();
 // m_tot.hesse();
 RooFitResult* r_chi2_tot = m_tot.save() ;
 cout << "==> Chi2 Fit results " << endl ;
 r_chi2_tot->Print("v") ;
 r_chi2_tot->floatParsFinal().Print("s") ;
 int NumberOfFreeParameters =  r_chi2_tot->floatParsFinal().getSize() ;
 cout<<"Normalized Chi2  = " << totChi2.getVal()/ (numberOfBins - NumberOfFreeParameters)<<endl; 
fracSigMutrk.Print();
fracSigMuMuNoIso.Print();
 
 
 //plotting
c1.cd(1);
 
 zmtMass->plotOn(massFrame_mutrk, LineColor(kBlue)) ;
 model_mutrk->plotOn(massFrame_mutrk,LineColor(kRed)) ;
 model_mutrk->plotOn(massFrame_mutrk,Components(*bkg_mutrk),LineColor(kGreen)) ;
 massFrame_mutrk->SetTitle("Z -> #mu track");
 // massFrame_mutrk->GetYaxis()->SetLogScale();
 massFrame_mutrk->Draw(); 
 
 c1.cd(2);
 
 zmmNotIsoMass->plotOn(massFrame_mumuNotIso, LineColor(kBlue)) ;
 model_mumuNotIso->plotOn(massFrame_mumuNotIso,LineColor(kRed)) ;
 model_mumuNotIso->plotOn(massFrame_mumuNotIso,Components(*bkg_mumuNotIso), LineColor(kGreen)) ;
 massFrame_mumuNotIso->SetTitle("Z -> #mu #mu not isolated");
 massFrame_mumuNotIso->Draw(); 
 
 
 c1.cd(3);
 
 zmsMass->plotOn(massFrame_musta, LineColor(kBlue)) ;
 eZmsSig->plotOn(massFrame_musta,LineColor(kRed)) ;
 massFrame_musta->SetTitle("Z -> #mu sta");
 massFrame_musta->Draw(); 
 
 
 
 c1.cd(4);
 
 HistoZCount1HLT->plotOn(massFrame_mumu1hlt, LineColor(kBlue)) ;
 eZmm1hltSig->plotOn(massFrame_mumu1hlt,LineColor(kRed)) ;
 massFrame_mumu1hlt->SetTitle("Z -> #mu #mu 1hlt");
 massFrame_mumu1hlt->Draw(); 
 
 c1.cd(5);
 
 HistoZCount2HLT->plotOn(massFrame_mumu2hlt, LineColor(kBlue)) ;
 eZmm2hltSig->plotOn(massFrame_mumu2hlt,LineColor(kRed)) ;
 massFrame_mumu2hlt->SetTitle("Z -> #mu #mu 2hlt");
 massFrame_mumu2hlt->Draw(); 
 
 
 c1.SaveAs("mass.eps");



 // test on toy....
 
   int n_experiments = 2;

 for (int i =1; i<= n_experiments; ++i){  

   double Lumi = 45.;
   // Lumi should be intented as pb-1 and passed from outside
   double scaleLumi=  Lumi / 45.0 ; // 45 is the current lumi correspondent to the histogram.....
 




 RooPlot * toy_massFrame_mutrk = mass->frame();
 RooPlot * toy_massFrame_mumuNotIso = mass->frame();
 RooPlot * toy_massFrame_musta = mass->frame();
 RooPlot * toy_massFrame_mumu1hlt = mass->frame();
 RooPlot * toy_massFrame_mumu2hlt = mass->frame();







 Yield.setVal(scaleLumi* 17000);
// ostream s; 
//s << r_chi2_tot->floatParsFinal().find("eff_hlt")->Print();
//  eff_hlt.setVal(0.916175);
//  eff_hlt.setConstant( kTRUE) ;
//  eff_tk.setVal(0.998111);
//  eff_tk.setConstant( kTRUE) ;
//  eff_sa.setVal(0.998111);
//  eff_sa.setConstant( 0.989971) ;
//  eff_iso.setVal(0.982936);
//  eff_iso.setConstant( kTRUE) ;
 //eff_tk.setVal(0.998);
 //eff_iso.setVal(0.983);
 //eff_sa.setVal(0.989);
 //a0->setVal(0.1);
 //  a1->setVal(-0.001);
 a2->setVal(0);
 a2->setConstant(kTRUE);
 //alpha->setVal(-0.015);
 // b0->setVal(0.1);
 // b1->setVal(-0.01);
 b2->setVal(0);
   b2->setConstant(kTRUE);
  //beta->setVal(-0.006);
//  nbkg_mutrk.setVal(770);
 nbkg_mutrk.setVal(scaleLumi* 830);
  nbkg_mutrk.setConstant(kTRUE);
 nbkg_mumuNotIso.setVal(scaleLumi* 1236);
 nbkg_mumuNotIso.setConstant(kTRUE);




 RooDataHist *toy_zmtMass = genHistFromModelPdf("toy_zmtMass", model_mutrk, mass, i, scaleLumi, range/rebinZMuTk);
 RooDataHist *toy_zmmNotIsoMass = genHistFromModelPdf("toy_zmmNotIsoMass", model_mumuNotIso, mass, i, scaleLumi,  range/rebinZMuMuNoIso);
 RooDataHist *toy_zmsMass = genHistFromModelPdf("toy_zmsMass", eZmsSig, mass, i, scaleLumi, range/rebinZMuSa );
 RooDataHist *toy_zmm1hltMass = genHistFromModelPdf("toy_zmm1hltMass", eZmm1hltSig, mass, i, scaleLumi,  range/rebinZMuMuHlt );
 RooDataHist *toy_zmm2hltMass = genHistFromModelPdf("toy_zmm2hltMass", eZmm2hltSig, mass,   i, scaleLumi, range/rebinZMuMuHlt);


RooChi2Var *toy_chi2_mutrk = new RooChi2Var("toy_chi2_mutrk","toy_chi2_mutrk",*model_mutrk,*toy_zmtMass, Extended(kTRUE)) ;

 RooChi2Var *toy_chi2_mumuNotIso = new RooChi2Var("toy_chi2_mumuNotIso","toy_chi2_mumuNotIso",*model_mumuNotIso,*toy_zmmNotIsoMass,Extended(kTRUE)) ;

 RooChi2Var *toy_chi2_musta = new RooChi2Var("toy_chi2_musta","toy_chi2_musta",*eZmsSig, *toy_zmsMass,  Extended(kTRUE) ) ;

 RooChi2Var *toy_chi2_mu1hlt = new RooChi2Var("toy_chi2_mu1hlt","toy_chi2_mu1hlt", *eZmm1hltSig, *toy_zmm1hltMass,  Extended(kTRUE) ) ;


 RooChi2Var *toy_chi2_mu2hlt = new RooChi2Var("toy_chi2_mu2hlt","toy_chi2_mu2hlt", *eZmm2hltSig, *toy_zmm2hltMass,  Extended(kTRUE) ) ;


 // adding the chi2
 RooAddition toy_totChi2("toy_totChi2","toy_chi2_mutrk + toy_chi2_mumuNotIso + toy_chi2_musta   + toy_chi2_mu1hlt  +  toy_chi2_mu2hlt " , RooArgSet(*toy_chi2_mutrk   ,*toy_chi2_mumuNotIso  , *toy_chi2_musta  ,*toy_chi2_mu1hlt  ,  *toy_chi2_mu2hlt  )) ;
 
 TCanvas cTk;
 toy_zmtMass->plotOn(toy_massFrame_mutrk, LineColor(kBlue)) ;
 toy_massFrame_mutrk->Draw();
 cTk.SaveAs("toyTk.eps");


 RooMinuit toy_m_tot(toy_totChi2) ;


 toy_m_tot.migrad();
 // m_tot.hesse();
 RooFitResult* r_toy_chi2_tot = toy_m_tot.save() ;
 cout << "==> Toy_Chi2 Fit results " << endl ;
 r_toy_chi2_tot->Print( " v ") ;
 r_toy_chi2_tot->floatParsFinal().Print("s");
 r_toy_chi2_tot->Write( ) ;
 // const int NumberOfFreeParameters = ;
 // cout<<"Number of bins = "<<numberOfBins<<" ***  getNPar = " << NumberOfFreeParameters <<endl;   
 NumberOfFreeParameters = r_toy_chi2_tot->floatParsFinal().getSize();
 cout<<"Normalized Toy_Chi2 OF TOY NUMBER " << i << "  = " << toy_totChi2.getVal()/ (numberOfBins - NumberOfFreeParameters)<<endl; 




 TCanvas cT;

 cT.Divide(2,3);
cT.cd(1);
 toy_zmtMass->plotOn(toy_massFrame_mutrk, LineColor(kBlue)) ;
 model_mutrk->plotOn(toy_massFrame_mutrk,LineColor(kRed)) ;
 model_mutrk->plotOn(toy_massFrame_mutrk,Components(*bkg_mutrk),LineColor(kGreen)) ;
 toy_massFrame_mutrk->SetTitle("Z -> #mu track");
 toy_massFrame_mutrk->Draw(); 

 cT.cd(2);
 
 toy_zmmNotIsoMass->plotOn(toy_massFrame_mumuNotIso, LineColor(kBlue)) ;
 model_mumuNotIso->plotOn(toy_massFrame_mumuNotIso,LineColor(kRed)) ;
 model_mumuNotIso->plotOn(toy_massFrame_mumuNotIso,Components(*bkg_mumuNotIso), LineColor(kGreen)) ;
 toy_massFrame_mumuNotIso->SetTitle("Z -> #mu #mu not isolated");
 toy_massFrame_mumuNotIso->Draw(); 
 
 
 cT.cd(3);
 
 toy_zmsMass->plotOn(toy_massFrame_musta, LineColor(kBlue)) ;
 eZmsSig->plotOn(toy_massFrame_musta,LineColor(kRed)) ;
 toy_massFrame_musta->SetTitle("Z -> #mu sta");
 toy_massFrame_musta->Draw(); 
 
 
 
 cT.cd(4);
 
 toy_zmm1hltMass->plotOn(toy_massFrame_mumu1hlt, LineColor(kBlue)) ;
 eZmm1hltSig->plotOn(toy_massFrame_mumu1hlt,LineColor(kRed)) ;
 toy_massFrame_mumu1hlt->SetTitle("Z -> #mu #mu 1hlt");
 toy_massFrame_mumu1hlt->Draw(); 
 
 cT.cd(5);
 
 toy_zmm2hltMass->plotOn(toy_massFrame_mumu2hlt, LineColor(kBlue)) ;
 eZmm2hltSig->plotOn(toy_massFrame_mumu2hlt,LineColor(kRed)) ;
 toy_massFrame_mumu2hlt->SetTitle("Z -> #mu #mu 2hlt");
 toy_massFrame_mumu2hlt->Draw(); 
 cT.SaveAs("toy_mass.eps");   







 delete toy_zmtMass;
delete toy_zmmNotIsoMass;
 delete toy_zmsMass;
 delete toy_zmm1hltMass;
 delete toy_zmm2hltMass;
 delete toy_chi2_mutrk;
 delete toy_chi2_mumuNotIso;
 delete toy_chi2_musta;
delete  toy_chi2_mu1hlt;
delete  toy_chi2_mu2hlt;
 delete toy_massFrame_mutrk;
 delete  toy_massFrame_mumuNotIso;
 delete  toy_massFrame_musta;
 delete  toy_massFrame_mumu1hlt;
delete  toy_massFrame_mumu2hlt;



 }


 /* how to read the fit result in root 
TH1D h_Yield("h_Yield", "h_Yield", 100, 10000, 30000)
for (int i =0: i < 100; i++){ 
RooFitResult* r = gDirectory->Get(Form("toy_totChi2;%d)",i)
//r->floatParsFinal().Print("s");
// without s return a list,  can we get the number?
 RooFitResult* r = gDirectory->Get("toy_totChi2;1")
// chi2
r->minNll();
//distamce form chi2.....
//r->edm();
// yield
r->floatParsFinal()[0]->Print();
//RooAbsReal * l = r->floatParsFinal()->first()
RooAbsReal * y = r->floatParsFinal()->find("Yield");
h_Yield->Fill(y->getVal());
}

  */




 delete root_file;
 delete out_root_file;
 return 0;
 
}
