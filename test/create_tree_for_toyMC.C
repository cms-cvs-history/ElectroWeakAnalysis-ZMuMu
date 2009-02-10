/*************************/
/*                       */ 
/* author: Pasquale Noli */
/* INFN Naples           */
/* Create TTree from     */
/* fit on Toy Montecarlo */
/*                       */
/*************************/
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "TObjArray.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include <iomanip>


void create_tree_for_toyMC()()
{
  gROOT->Reset();

  Double_t Y;
  Double_t dY;
  Double_t Tk;
  Double_t dTk;
  Double_t Sa;
  Double_t dSa;
  Double_t Iso;
  Double_t dIso;
  Double_t Hlt;
  Double_t dHlt;
  Double_t chi2;

  TFile *f;
  TTree *tree;
  
  f = new TFile("fitResult.root","RECREATE");
  tree = new TTree("tree"," C data from ASCII file");

  tree->Branch("Y",&Y,"Y/D");
  tree->Branch("dY",&dY,"dY/D");
  tree->Branch("Tk",&Tk," Tk/D");
  tree->Branch("dTk",&dTk," dTk/D");
  tree->Branch("Sa",&Sa," Sa/D");
  tree->Branch("dSa",&dSa," dSa/D");
  tree->Branch("Iso",&Iso," Iso/D");
  tree->Branch("dIso",&dIso," dIso/D");
  tree->Branch("Hlt",&Hlt," Hlt/D");
  tree->Branch("dHlt",&dHlt," dHlt/D");
  tree->Branch("chi2",&chi2," chi2/D");

  ifstream fin;
  fin.open("fitResult.txt");
  while(!(fin.eof())){
    fin >> Y >> dY >> Tk >>  dTk >>  Sa >>
      dSa >>  Iso >> dIso >>  Hlt >>  dHlt >>chi2;
    tree->Fill();
  }

  tree->Print();
  f->Write();
  f->Close();


}
