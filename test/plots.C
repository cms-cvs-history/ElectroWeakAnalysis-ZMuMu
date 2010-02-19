{
  TFile *file = TFile::Open("NtupleLoose_test.root");
  TTree * Events = dynamic_cast< TTree *> (file->Get("Events"));
  TFile * output_file = TFile::Open("histo.root", "RECREATE");
  
  // zGolden plots
  TCut cut_zGolden("zGoldenMass>20 && zGoldenDau1Pt> 20 && zGoldenDau2Pt>20 && zGoldenDau1Iso< 3.0 && zGoldenDau2Iso < 3.0");
  TDirectory * dir = output_file->mkdir("goodZToMuMuPlots");
  dir->cd();
  TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  //  Events->Draw("zGoldenMass");
  Events->Project("zMass", "zGoldenMass", cut_zGolden );
  zMass->Write();
  output_file->cd("/");

  // zGolden1HLT plots
  // TCut cut_zGolden1HLT("zGolden1HLTMass>20 && zGoldenDau1Pt> 20 && zGolden1HLTDau2Pt>20 && zGolden1HLTDau1Iso< 3.0 && zGolden1HLTDau2Iso < 3.0");
  // to be equivalent to:
  TCut cut2_zGolden1HLT("zGoldenMass>20 && zGoldenDau1Pt> 20 && zGoldenDau2Pt>20 && zGoldenDau1Iso< 3.0 && zGoldenDau2Iso < 3.0 && ((zGoldenDau1HLTBit==1 && zGoldenDau2HLTBit==0) || (zGoldenDau1HLTBit==0 && zGoldenDau2HLTBit==1))");
  TDirectory * dir = output_file->mkdir("goodZToMuMu1HLTPlots");
  dir->cd();
  TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  // TH1F * zMass2 = new TH1F("zMass2", "zMass2", 200, 0, 200);
  //  Events->Draw("zGoldenMass");
  //Events->Project("zMass", "zGolden1HLTMass", cut_zGolden1HLT );
   Events->Project("zMass", "zGoldenMass", cut2_zGolden1HLT );
 zMass->Write();
 //zMass2->Write();

  output_file->cd("/");


 // zGolden2HLT plots 
  // TCut cut_zGolden2HLT("zGolden2HLTMass>20 && zGoldenDau1Pt> 20 && zGolden2HLTDau2Pt>20 && zGolden2HLTDau1Iso< 3.0 && zGolden2HLTDau2Iso < 3.0");
  TCut cut2_zGolden2HLT("zGoldenMass>20 && zGoldenDau1Pt> 20 && zGoldenDau2Pt>20 && zGoldenDau1Iso< 3.0 && zGoldenDau2Iso < 3.0 && ((zGoldenDau1HLTBit==1 && zGoldenDau2HLTBit==1) )");  
  TDirectory * dir = output_file->mkdir("goodZToMuMu2HLTPlots");
  dir->cd();
  //TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  //  Events->Draw("zGoldenMass");
  // Events->Project("zMass", "zGolden2HLTMass", cut_zGolden2HLT );
  Events->Project("zMass", "zGoldenMass", cut2_zGolden2HLT );

  zMass->Write();
  //zMass2->Write();

  output_file->cd("/");

 // zGoldenNotIso plots
  TCut cut_zGoldenNotIso("zGoldenMass>20 && zGoldenDau1Pt> 20 && zGoldenDau2Pt>20 && ( zGoldenDau1Iso> 3.0 || zGoldenDau2Iso > 3.0 )");
  TDirectory * dir = output_file->mkdir("nonIsolatedZToMuMuPlots");
  dir->cd();
  TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  //  Events->Draw("zGoldenMass");
  Events->Project("zMass", "zGoldenMass", cut_zGoldenNotIso );
  zMass->Write();
  output_file->cd("/");

// zMuTrk plots
  TCut cut_zMuTrk("zMuTrkMass>20 && zMuTrkDau1Pt> 20 && zMuTrkDau2Pt>20 && zMuTrkDau1Iso< 3.0 && zMuTrkDau2Iso < 3.0");
  TDirectory * dir = output_file->mkdir("goodZToMuMuOneTrackPlots");
  dir->cd();
  TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  //  Events->Draw("zGoldenMass");
  Events->Project("zMass", "zMuTrkMass", cut_zMuTrk );
  zMass->Write();
  output_file->cd("/");


// zMuSta plots
  TCut cut_zMuSta("zMuStaMass>20 && zMuStaDau1Pt> 20 && zMuStaDau2Pt>20 && zMuStaDau1Iso< 3.0 && zMuStaDau2Iso < 3.0");
  TDirectory * dir = output_file->mkdir("goodZToMuMuOneStandAloneMuonPlots");
  dir->cd();
  TH1F * zMass = new TH1F("zMass", "zMass", 200, 0, 200);
  //  Events->Draw("zGoldenMass");
  Events->Project("zMass", "zMuStaMass", cut_zMuSta );
  zMass->Write();
  output_file->cd("/");


  //(mu1.phi -mu2.phi)
  TH1F * deltaPhi = new TH1F("deltaPhi", "deltaPhi", 120, 0, 6.);
  TH1F * h2 = new TH1F("h2", "h2", 120, 0, 6. );
  TH1F * h3 = new TH1F("h3", "h3", 120, 0, 6. );
  //  Events->Draw("zGoldenMass");

  /*    result = phi1 - phi2;
040     while (result > M_PI) result -= 2*M_PI;
041     while (result <= -M_PI) result += 2*M_PI;
042     return result;
  */ 



   Events->Project("deltaPhi", "abs(zGoldenDau1Phi - zGoldenDau2Phi)", "-TMath::Pi() < (zGoldenDau1Phi - zGoldenDau2Phi) < TMath::Pi()" );
   Events->Project("h2", "abs(zGoldenDau1Phi - zGoldenDau2Phi - 2 * TMath::Pi())", "(zGoldenDau1Phi - zGoldenDau2Phi) > TMath::Pi()" );
   Events->Project("h3", "abs(zGoldenDau1Phi - zGoldenDau2Phi + 2 * TMath::Pi())", "(zGoldenDau1Phi - zGoldenDau2Phi) <=  -TMath::Pi()" );

   // mu1PhiMinusMu2Phi ->Add(  h2 , h3 );
   //mu1PhiMinusMu2Phi->Draw();
   // output_file->cd("/");
   deltaPhi->Add(h2, h3);
   deltaPhi->Write();
  
   delete h2;
   delete h3;



   // mu1.eta -mu2.eta
 TH1F * deltaEta = new TH1F("deltaEta", "deltaEta", 120, 0, 6.);
   Events->Project("deltaEta", "abs(zGoldenDau1Eta - zGoldenDau2Eta)" );
   deltaEta->Write();

   // mu1.eta -mu2.eta
 TH1F * zGoldenPt = new TH1F("zGoldenPt", "zGoldenPt", 200, 0, 200);
   Events->Project("zGoldenPt", "zGoldenPt" );
   zGoldenPt->Write();

   // isolations...

 TH1F * TrkIso = new TH1F("TrkIso", "TrkIso", 1000, 0, 100);
 TH1F * h2 = new TH1F("h2", "h2", 1000, 0, 100);
   Events->Project("TkIso", "zGoldenDau1TrkIso" );
   Events->Project("h2", "zGoldenDau2TrkIso" );
   TrkIso->Add(h2);

   TrkIso->Write();
   delete h2;

 TH1F * EcalIso = new TH1F("EcalIso", "EcalIso", 1000, 0, 100);
 TH1F * h2 = new TH1F("h2", "h2", 1000, 0, 100);
   Events->Project("TkIso", "zGoldenDau1EcalIso" );
   Events->Project("h2", "zGoldenDau2EcalIso" );
   EcalIso->Add(h2);

   EcalIso->Write();
   delete h2;

 TH1F * HcalIso = new TH1F("HcalIso", "HcalIso", 1000, 0, 100);
 TH1F * h2 = new TH1F("h2", "h2", 1000, 0, 100);
   Events->Project("TkIso", "zGoldenDau1HcalIso" );
   Events->Project("h2", "zGoldenDau2HcalIso" );
   HcalIso->Add(h2);

   HcalIso->Write();
   delete h2;

   // eta of muon not triggered

   //   Events->Draw("zGoldenDau2Eta", "zGoldenDau1HLTBit==0");
 TH1F * muNotTriggeredEta = new TH1F("muNotTriggeredEta", "muNotTriggeredEta", 240, -6, 6.);
 TH1F * h2 = new TH1F("h2", "h2", 240, -6, 6.);

 Events->Project("muNotTriggeredEta","zGoldenDau1Eta", "zGoldenDau1HLTBit==0");
 Events->Project("h2","zGoldenDau2Eta", "zGoldenDau2HLTBit==0");

 muNotTriggeredEta->Add(h2);
 muNotTriggeredEta->Write();
 delete h2;
    output_file->Close();



 
 
}
