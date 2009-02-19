/*************************/
/*                       */
/* author: Pasquale Noli */
/* INFN Naples           */
/* macro to save the eps */ 
/* of pulls              */
/*                       */
/*************************/

{
  TFile *f = TFile::Open("fitResult.root"); 
  TH1D frameTrk("frameTrk", "track eff.", 100, -10, 10);
  TH1D frameSa("frameSa", "stand-alone eff.", 100, -10, 10);
  TH1D frameIso("frameIso", "isolation eff.", 100, -10, 10);
  TH1D frameHlt("frameHlt", "HLT eff.", 100, -10, 10);
  tree->Project("frameTrk","(Tk-Tk_true)/dTk", "abs((Tk-Tk_true)/dTk)>5.e-3");
  tree->Project("frameSa","(Sa-Sa_true)/dSa", "abs((Sa-Sa_true)/dSa)>5.e-3"); 
  tree->Project("frameIso", "(Iso-Iso_true)/dIso", "abs((Iso-Iso_true)/dIso)>5.e-3");
  tree->Project("frameHlt", "(Iso-Iso_true)/dIso", "abs((Iso-Iso_true)/dIso)>5.e-3");
  frameTrk.Fit("gaus");
  frameSa.Fit("gaus");
  frameIso.Fit("gaus");
  frameHlt.Fit("gaus");
  TCanvas *c1 = new TCanvas("c1","pulls",10,10,900,900);
  gStyle->SetOptStat(111111);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetFitFormat("5.3g");
  c1->Divide (2,2);
  c1->cd(1);
  frameTrk.Draw();
  c1->cd(2);
  frameSa.Draw();
  c1->cd(3);
  frameIso.Draw();
  c1->cd(4);
  frameHlt.Draw();
  c1->Draw();
  c1->SaveAs("pulls.eps");
}

