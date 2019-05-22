#ifndef __CINT__
#include "RooGlobalFunc.h"
#else
// Refer to a class implemented in libRooFit to force its loading
// via the autoloader.
class Roo2DKeysPdf;
#endif

#include <time.h>
#include <TROOT.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TTimeStamp.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <iostream>
#include <TMath.h> 
#include "TH1D.h"
#include "TString.h"
#include "RooAddPdf.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "RooPrintable.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCmdArg.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooKeysPdf.h"
#include "RooMCStudy.h"
#include "RooFitResult.h"
#include <sstream>
#include <TLatex.h>
#include "ConfigParams.h"
#include <TLegend.h>
#include <TLine.h>
#include <RooExponential.h>
#include <RooHist.h>
#include <TStyle.h>
#include "TApplication.h"
#include "RooStats/Heaviside.h"
#include "TFormula.h"
#include "RooIntegrator1D.h"
#include "RooRealIntegral.h"
#include "RooTFnBinding.h" 
#include "Math/DistFunc.h"
#include "RooKeysPdf.h"
#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction3Binding.h"
#endif
using namespace std;
using namespace RooFit;
using namespace RooStats;
//RooMsgService::instance().deleteStream(0);
//RooMsgService::instance().deleteStream(0);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//        _______  ______ _______ __   _ ______ 
// |      |______ |  ____ |______ | \  | |     \
// |_____ |______ |_____| |______ |  \_| |_____/
//                                              
/////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

TLegend* getLegend(const char* signalText, ConfigParams* myconfigs) {
  TLegend* legend = new TLegend(myconfigs->m_legendPos[0],myconfigs->m_legendPos[1],myconfigs->m_legendPos[2],myconfigs->m_legendPos[3]);
  TLine* l1 = new TLine();
  l1->SetLineColor(kRed);
  l1->SetLineWidth(5.);

  legend->AddEntry(l1,signalText,"L"); 
  legend->SetTextFont(myconfigs->m_font);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetMargin(myconfigs->m_legendMargin);
  legend->SetTextSize(myconfigs->m_legendTextSize);

  return legend;
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
// _______  _____   _____                 _____ 
// |_____| |_____] |     | |      |      |     |
// |     | |       |_____| |_____ |_____ |_____|
//                                              
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//
 Double_t apollo(Double_t *x,Double_t *par){
   Double_t xi = (x[0]-par[0])/par[1];
   Double_t xi2 = xi*xi;
   if (xi > -par[3]) {return exp(-par[2]*sqrt(1+ xi2));}
   Double_t a2 = par[3]*par[3];
  
   Double_t B =  -1.*(par[3] - par[4]*sqrt(1+a2)/(par[3]*par[2]));
   Double_t A = exp(-par[2]*sqrt(1+a2)) *pow(B+par[3],par[4]);
   
   return A*pow(B-xi,-par[4]);
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
// _______ _____ _______
// |______   |      |   
// |       __|__    |   
//                      
////////////////////////////////////////////////////
////////////////////////////////////////////////////

void fit(ConfigParams* myconfigs, int forsig=0, int usesigd = 1, int usesigs = 1, int usetrig = 1, int usecssample = 0, int usetightsample = 1, double mmin=5350,double mmax=5900){

    TFile* f;
    if (usetightsample) {
     RooRealVar MLb_DTF("MLb_DTF","MLb_DTF",mmin,mmax);
     RooRealVar mean("mean","",5620.4,5600,5640);
     f=new TFile("/sum/lb/Lb2LcDs.root");
    }
  
  TTree* tree = (TTree*) f->Get("");
  double mlb;
  tree->SetBranchStatus("*",1);
  tree->SetBranchAddress("MLb_DTF",&mlb);

  Int_t numentries = tree->GetEntries();

  RooDataSet *data = new RooDataSet("data","",RooArgSet(MLb_DTF));

  for (int i =0;i<numentries;i++) {
    tree->GetEntry(i);
    if(mlb < mmin || mlb > mmax) continue;
    MLb_DTF.setVal(mlb);
    data->add(RooArgSet(MLb_DTF));
  }  
  data->Print();

  // fit variables to all the PDFs
  RooRealVar nsig("nsig","",5000,3000,8500);// B->dds
  RooRealVar ncomb("ncomb","",4800,3000,17400);
  RooRealVar nlcdsst("nlcdsst","",6500,3400,15400);

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//                                                  
//    /\     //
//   /  \    //
//  / /\ \   //
// / /  \ \  //////  
///_/    \_\ //  //
//           //////
///////////////////////////////////////////////////
///////////////////////////////////////////////////
// Define signal, Apollo Func
  RooPlot* debugframe=MLb_DTF.frame();
//Declare Func/PDF Parameters
Double_t m0(5269),sigma(5),b(10),a(1),n(10);

RooRealVar msigma("msigma","",25.42802,1,30);
RooRealVar m_b("b","",.61728,.05,5);
RooRealVar m_a("a","",2.92796,0.1,4);
RooRealVar m_n("n","",2.48844,.1,4);
mean.setConstant(kTRUE);
//msigma.setConstant(kTRUE);
//m_n.setConstant(kTRUE);
//m_a.setConstant(kTRUE);
//m_b.setConstant(kTRUE);//*/


RooFormulaVar masssig("masssig","","@0/@1",RooArgList(msigma,m_b));
//Create bound function 
  TF1 *apollofunc= new TF1("apollofunc",apollo,mmin,mmax,5);

  apollofunc->SetParameters(m0,sigma,b,a,n);
  RooAbsReal *L_bfunc=bindFunction(apollofunc,MLb_DTF,RooArgList(mean,masssig,m_b,m_a,m_n));
//Create 0 bkg signal
  RooRealVar unit("unit","",0);
unit.setConstant(kTRUE);
RooChebychev fake("fake","",MLb_DTF,RooArgList(unit));
//Create Func PDF(B_0)
RooAbsPdf *Lb_pdf= new RooRealSumPdf("Lb_pdf","",fake,*L_bfunc,unit);
///////////////////////////////////////////////////
TFile *lowBFile = new TFile("LcDsst.KEYS.root");
  RooAbsPdf* pdf_PARTRECO = (RooKeysPdf*)lowBFile->Get("pdf_PARTRECO");

  RooRealVar a0("a0", "a0", -.003, -.15, 0);  
  RooExponential comb_expo("comb_expo","",MLb_DTF,a0); //Combinatorics PDFs 
  



  RooRealVar  widesiglikeback1width("widesiglikeback1width","widesiglikeback1width",33.6);
  RooGaussian widesiglikeback1("widesiglikeback1","widesiglikeback1",MLb_DTF,mean,widesiglikeback1width);

  RooRealVar charmlessfrac("charmlessfrac","charmlessfrac",myconfigs->fixedChlessReflForDsDs,0.0,0.1);
  RooFormulaVar  ncharmless("ncharmless","@0*@1",RooArgList(nsig,charmlessfrac));
// Keep charmless pdf mass window shape constant
charmlessfrac.setVal(myconfigs->fixedChlessReflForLcDs);
charmlessfrac.setConstant(kTRUE);

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//
// Define totalPDF and fit to dataset
//
////////////////////////////////////////////////////
////////////////////////////////////////////////////

  RooAbsPdf* pdf = new RooAddPdf("pdf","",RooArgList(*Lb_pdf,comb_expo,*pdf_PARTRECO,widesiglikeback1), RooArgList(nsig, ncomb,nlcdsst,ncharmless));
cout<<"Welcome to: "<<endl;
 /* data->plotOn(debugframe);
  pdf->plotOn(debugframe);
  debugframe->Draw();return 0;// */
  const char time[100]=TDatime().AsString();
  cout<<time<<endl;
  pdf->fitTo(*data,PrintLevel(-1),Extended());

////////////////////////////////////////////////////
////////////////////////////////////////////////////
//  _____          _____  _______
// |_____] |      |     |    |   
// |       |_____ |_____|    |   
//                               
////////////////////////////////////////////////////
////////////////////////////////////////////////////

//Initialize canvas and RooPlot frame, plot data on frame
  TCanvas *c = new TCanvas("c", "c",myconfigs->m_canvWidth,myconfigs->m_canvHeight); 
  c->cd();
  RooPlot* plot = MLb_DTF.frame(5350,5900,110); 
 // plot->SetTitle("LHCb");
  //plot->GetXaxis()->SetNDivisions(4);
  data->plotOn(plot,Binning(myconfigs->m_nBinsData),MarkerStyle(myconfigs->m_markerStyle));
 
  //frame Settings
  plot->GetXaxis()->SetTitle("#Lambda_{c}^{-} D^{+}_{s} Mass [MeV/c^{2}]");
  plot->GetYaxis()->SetTitle("Candidates / (5 MeV/c^{2})");
  plot->GetYaxis()->SetTitleOffset(1.1);
  plot->GetYaxis()->SetNdivisions(505);
  plot->GetYaxis()->SetTitleFont(myconfigs->m_axisfont);
  plot->GetYaxis()->SetTitleSize(myconfigs->m_titlesize);
  plot->GetXaxis()->SetTitleSize(myconfigs->m_titlesize);
  plot->GetXaxis()->SetTitleFont(myconfigs->m_axisfont);
  plot->GetXaxis()->SetLabelSize(.06);
  plot->GetYaxis()->SetLabelSize(.06);

  pdf->plotOn(plot,Components("Lb_pdf,comb_expo,pdf_PARTRECO,widesiglikeback1"),LineStyle(kSolid),LineColor(kBlack),DrawOption("F"),FillColor(kWhite)); 
  RooHist* pull_hist = (*plot).pullHist(); 
  pull_hist->GetYaxis()->SetTitle("Fit pull /(5 MeV/c^{2})");
    pdf->plotOn(plot,Components("widesiglikeback1,comb_expo,pdf_PARTRECO,"),LineColor(15),DrawOption("F"),FillColor(15));
  pdf->plotOn(plot,Components("comb_expo,pdf_PARTRECO"),LineColor(kBlack),DrawOption("F"),FillColor(kCyan));

   pdf->plotOn(plot,Components("comb_expo"),LineColor(kBlack),DrawOption("F"),FillColor(kCyan-5));
     pdf->plotOn(plot,Components("Lb_pdf,comb_expo,pdf_PARTRECO"),LineColor(kBlack));

  pdf->plotOn(plot,Components("Lb_pdf"),LineColor(kRed)); 

  TLegend* legend = getLegend("#Lambda_{b}^{0}#rightarrow #Lambda_{c}^{#bf{-}} D^{#bf{+}}_{s}",myconfigs); 



  TH1F* lcdssthistoforleg = new TH1F("lcdssthisto","",10,0,1);  
  lcdssthistoforleg->SetFillColor(kCyan);
  legend->AddEntry(lcdssthistoforleg,"#Lambda_{b}^{0}#rightarrow #Lambda_{c}^{#bf{-}} D_{s}#kern[-0.75]{#lower[-0.1]{*}}^{#bf{+}}","F"); 

  TH1F* combhistoforleg = new TH1F("combhisto","",10,0,1);//Combinatorics
  combhistoforleg->SetFillColor(kCyan-5);
  legend->AddEntry(combhistoforleg,"Combinatorics","F"); 
 
  TH1F* widesiglikeback1histoforleg = new TH1F("widesiglikeback1histoforleg","widesiglikeback1histoforleg",10,0,1);
  widesiglikeback1histoforleg->SetFillColor(15);
  legend->AddEntry(widesiglikeback1histoforleg,"#Lambda_{b}^{0}#rightarrow #Lambda_{c}^{#bf{-}} K^{#bf{+}}K^{#bf{-}}#pi^{#bf{+}}","F");


  data->plotOn(plot,Binning(myconfigs->m_nBinsData),MarkerStyle(myconfigs->m_markerStyle));
  plot->Draw("E");
  legend->Draw("same");

  ////TLatex


  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42);
  myLatex->SetTextColor(1);
  myLatex->SetTextSize(0.04);
  myLatex->SetTextAlign(12);
  myLatex->SetNDC(kTRUE);
  myLatex->SetTextSize(0.035);

  //Pull plots
  //

  TCanvas* pull = new TCanvas("pull","pull",myconfigs->m_canvWidth,myconfigs->m_canvHeight);
  pull->cd();  
  RooPlot* pull_plot = MLb_DTF.frame(5350,5900,myconfigs->m_nBinsData);
  pull_plot->GetXaxis()->SetTitle("#Lambda_{c}^{-} D^{+}_{s} Mass [MeV/c^{2}]");
  pull_plot->GetYaxis()->SetTitle("Fit Pull");
  pull_plot->GetYaxis()->SetTitleFont(myconfigs->m_axisfont);
  pull_plot->GetXaxis()->SetTitleFont(myconfigs->m_axisfont);
  pull_plot->GetYaxis()->SetTitleSize(myconfigs->m_titlesize);
  pull_plot->GetXaxis()->SetTitleSize(myconfigs->m_titlesize);

  pull_plot->addPlotable(pull_hist,"P");
  pull_plot->Draw(); 
  
    c->Print("/FitResults_ForPaper/Lb2LcDsFit"+myconfigs->m_suffix+".pdf");
    c->SetLogy();
    plot->GetYaxis()->SetRangeUser(5,6000);    
    c->Print("/FitResults_ForPaper/Lb2LcDsFitLog"+myconfigs->m_suffix+".pdf");
    pull->Print("/FitResults_ForPaper/Lb2LcDsFitPull"+myconfigs->m_suffix+".pdf");
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//  _____   ______ _____ __   _ _______
// |_____] |_____/   |   | \  |    |   
// |       |    \_ __|__ |  \_|    |   
//                                     
////////////////////////////////////////////////////
////////////////////////////////////////////////////
  
a0.Print();  
nsig.Print();
ncomb.Print(); 
nlcdsst.Print();
ncharmless.Print();
mean.Print();
m_b.Print();
m_n.Print();
m_a.Print();
msigma.Print();
masssig.Print();
cout<<"START: "<<time<<endl;
cout<<"END: "<<TDatime().AsString()<<endl;

//cout<<""<<endl<<""<<;>texresults.txt


}
int main(int argc, char* argv[],ConfigParams* myconfigs){
TApplication* rootapp= new TApplication("why",&argc,argv); //Required to work
ConfigParams* figs=new ConfigParams;//ConfigParams* myconfigs
fit(figs);
return 0;
}

