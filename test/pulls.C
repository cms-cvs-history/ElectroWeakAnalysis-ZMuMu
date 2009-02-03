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
  TCanvas *c1 = new TCanvas("c1","pulls",10,10,900,900);
  c1->Divide (2,2);
  c1->cd(1);
  tree->Fit("gaus","(Tk-0.998364)/dTk");
  c1->cd(2);
  tree->Fit("gaus","(Sa-0.998364)/dSa"); 
  c1->cd(3);
  tree->Fit("gaus","(Iso-0.998364)/dIso");
  c1->cd(4);
  tree->Fit("gaus","(Hlt-0.998364)/dHlt");
  c1->Draw();
  c1->SaveAs("pulls.eps");

}

