#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "TH1D.h"
#include "TRandom.h"
#include "roc_header.h"

using namespace RooFit;

int main(int arg0, int argc, char* argv[]){
int a;
TApplication* rootapp= new TApplication("why",&argc,argv); //Required to work

TFile *f= new TFile("/sum/roc/chainedtree_noDDs_windowcuts.root");
TTree *t= (TTree*)f->Get("T");

///////////////////////////////////////////// Defining PDF variables
RooRealVar MD("MD","",1820,1920);
RooRealVar MDs("MDs","",1920,2020);
RooRealVar MB_DTF("MB_DTF","",5220,5400);

RooRealVar MDmean("MDmean","",1870,1820,1920);
RooRealVar MDsmean("MDsmean","",1970,1920,2020);
RooRealVar sig("sig","",10,0,15);
RooRealVar alpha("alpha","",.0005,0,1);

///////////////////////////////////////////// Defining PDFs, 0->D, 1->Ds
if (arg0==0){
RooGaussian gaussfit("gauss","",MD,MDmean,sig);
RooExponential expfit("exp","",MD,alpha);
RooRealVar gfrac("gfrac","",.5,0,1);

RooAddPdf totalfit("totalfit","",RooArgList(gaussfit,expfit),gfrac);
ofstream g("C:\\sum\\roc\\D.opt");
double bdt_min=t.GetMinimum("d_bdt");
double bdt_max=t.GetMaximum("d_bdt");
cout<<bdt_min<<" "<< bdt_max<< endl;
MD.Print();
MDmean.Print();
gaussfit.Print();
expfit.Print();
}

if (arg0==1){
RooGaussian gaussfit("gauss","",MDs,MDsmean,sig);
RooExponential expfit("exp","",MDs,alpha);
RooRealVar gfrac("gfrac","",.5,0,1);

RooAddPdf totalfit("totalfit","",RooArgList(gaussfit,expfit),gfrac);
ofstream g("C:\\sum\\roc\\Ds.opt");
double bdt_min=t.GetMinimum("ds_bdt");
double bdt_max=t.GetMaximum("ds_bdt");
}

///////////////////////////////////////////// Fitting for every bdt value
double n=.15;
for(n; n<bdt_max; n+=.01){


    if(arg0==0){
RooRealVar d_bdt("d_bdt","",n,1);
RooDataSet dataset("dataset","",RooArgSet(MD,d_bdt),Import(*t));
    dataset.Print();}
    if(arg0==1){
RooRealVar ds_bdt("ds_bdt","",n,1);
RooDataSet dataset("dataset","",RooArgSet(MDs,ds_bdt),Import(*t));}

RooRealVar nsig("nsig","signal events",0,100000); 
RooRealVar nbkg("nbkg","signal background events",0,100000); 
RooAddPdf model("model","",RooArgList(gaussfit,expfit),RooArgList(nsig,nbkg)); 
model.fitTo(dataset,PrintLevel(-1), Extended()); 
//nsig.Print();
//nbkg.Print();
if(n==bdt_min){a =nsig.getVal();
cout <<"will not crash" << endl;
g<<endl;}//creates a first blank line, and defined eff parameter a
float eff=nsig.getVal()/a;
/*if(n>.24 && n<.25){
RooPlot* frame = MDs.frame(Title("MD1 Eff Fit"));
dataset.plotOn(frame);
model.plotOn(frame,Components(RooArgSet(expfit)),LineColor(kRed));
model.plotOn(frame,LineColor(kBlue));
TCanvas* c= new TCanvas("c","");
c->Divide(2,2);
c->cd(1);
frame->Draw();
}
if(n>.25 && n<.26){
RooPlot* frame2 = MDs.frame(Title("MD2 Eff Fit"));
dataset.plotOn(frame2);
model.plotOn(frame2,Components(RooArgSet(expfit)),LineColor(kRed));
model.plotOn(frame2,LineColor(kBlue));
c->cd(2);
frame2->Draw();
}
if(n>.26 && n<.27){
RooPlot* frame3 = MDs.frame(Title("MD3 Eff Fit"));
dataset.plotOn(frame3);
model.plotOn(frame3,Components(RooArgSet(expfit)),LineColor(kRed));
model.plotOn(frame3,LineColor(kBlue));
c->cd(3);
frame3->Draw();
}
if(n>.34 && n<.35){
RooPlot* frame4 = MDs.frame(Title("MD4 Eff Fit"));
dataset.plotOn(frame4);
model.plotOn(frame4,Components(RooArgSet(expfit)),LineColor(kRed));
model.plotOn(frame4,LineColor(kBlue));
c->cd(4);
frame4->Draw();
}
*/

g<<n<<"\t"<< nsig.getVal() <<"\t"<< nbkg.getVal()<<"\t"<< eff<<endl;
//totalfit.fitTo(dataset);
//alpha.Print();
//mean.Print();
//sig.Print();
}
g.close();
/*
////////////////////////////////////////////// Plotting

RooPlot* frame1 = MD.frame(Title("Binned MB_DTF")) ;

dataset.plotOn(frame1) ;
totalfit.plotOn(frame1,Components(RooArgSet(expfit)),LineColor(kRed));
totalfit.plotOn(frame1,LineColor(kBlue));
cout << "TotalEvents= "<< totalfit.getVal() <<" "<< gaussfit.getVal() <<endl;

////////////////////////////////////////////// Drawing
TCanvas* c= new TCanvas("c","");

frame1->Draw();
*/
return cumulative_bdt_out(arg0);
}

/*
for(all different bdt)

{
Fit(only ones bdt value)
write line to file of parameters
}
*/

