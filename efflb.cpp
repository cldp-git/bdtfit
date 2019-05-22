///lb2lcDS setup
void efflb(){
const char histname[100]="Lc_NC_sig";
const char hist2name[100]="Ds_NC_sig";
const char histfilename[100]="/sum/eggs/dfromb_effs.root";
const char mcdataname[100]="/sum/eggs/lb2lcds.mc.root";
const char Dmass[100]="lab5_MM";//L_c
const char Dsmass[100]="lab1_MM";
const char Dchi2[100]="lab5_FDCHI2_OWNPV";
const char Dschi2[100]="lab1_FDCHI2_OWNPV";
const char DPX[100]="lab5_PX";
const char DsPX[100]="lab1_PX";
const char DPY[100]="lab5_PY";
const char DsPY[100]="lab1_PY";
const char outfilename1[100]="C:\\sum\\eggs\\Lb_lc_mc.opt";
const char outfilename2[100]="C:\\sum\\eggs\\Lb_ds_mc.opt";

ofstream g1(outfilename1);
ofstream g2(outfilename2);
Double_t d_total_eff,ds_total_eff;
TFile *histfile = new TFile(histfilename);
//Declare hist pointer & variables to load tree values
TH3F *h,*hh;
histfile->GetObject(histname,h);
histfile->GetObject(hist2name,hh);

TFile *mcdatafile= new TFile(mcdataname);
// Declare Tree, variables
TTree* mcdata=(TTree*) mcdatafile->Get("T");

Double_t md,mds,d_chi2,ds_chi2,d_px,d_py,ds_px,ds_py,d_pt,ds_pt,d_bdt,ds_bdt,d_total,d_partotal,ds_total,ds_partotal,d_eff,ds_eff;
int dcounts,dscounts;
bool d_preventdoublecount,ds_preventdoublecount;
//
//set input/output branch addresses
mcdata->SetBranchStatus("*",1);
mcdata->SetBranchAddress(Dmass,&md);//L_c
mcdata->SetBranchAddress(Dsmass,&mds);
mcdata->SetBranchAddress(Dchi2,&d_chi2);
mcdata->SetBranchAddress(Dschi2,&ds_chi2);
mcdata->SetBranchAddress(DPX,&d_px);
mcdata->SetBranchAddress(DsPX,&ds_px);
mcdata->SetBranchAddress(DPY,&d_py);
mcdata->SetBranchAddress(DsPY,&ds_py);
double d_bdt_max=h->Project3D("z")->GetXaxis()->GetXmax();

double d_bdt_min=h->Project3D("z")->GetXaxis()->GetXmin();

double ds_bdt_max=hh->Project3D("z")->GetXaxis()->GetXmax();

double ds_bdt_min=hh->Project3D("z")->GetXaxis()->GetXmin();
int n=mcdata->GetEntries();
int binID,binID2,binID3;
/*cout<<d_bdt_max<<" "<<d_bdt_min<<" "<<ds_bdt_min<<" "<<ds_bdt_max<<endl;
h->Project3D("z")->Draw();
return 0; // */


for(double mcout=-.28;mcout<.5;mcout+=.005){
d_total_eff=0;
ds_total_eff=0;
    if (mcout==-.28){g1<<endl;g2<<endl;}//adds empty line to file

///////////////////////////////////////////////
////////////////////////////////////////////////

for(int i=0;i<n;i++){//reinitialize variables to 0
    d_total=0; d_partotal=0; ds_total=0; ds_partotal=0; d_eff=0; ds_eff=0; d_preventdoublecount=0;ds_preventdoublecount=0;
    mcdata->GetEntry(i);
    d_pt=sqrt(d_px*d_px+d_py*d_py);
    ds_pt=sqrt(ds_px*ds_px+ds_py*ds_py);

//////////////////////////////////////////v
    for(double d_bdtval=d_bdt_min;
            d_bdtval<d_bdt_max+.005;
            d_bdtval+=.01)
    {
        binID=h->FindFixBin(d_pt,d_chi2,d_bdtval);
        if(binID==binID3)
            continue;
        binID3=binID;

        d_total+=h->GetBinContent(binID);
        if(mcout<d_bdtval)
            d_partotal+=h->GetBinContent(binID);


        if(fabs(d_bdtval-d_bdt_max)<.005 && !d_preventdoublecount){
            if (d_total>0) d_eff=d_partotal/d_total;
            else {d_eff=-1;cout<<"bug?"<<endl;};
if (d_eff>1) d_eff=1;
            d_total_eff+=d_eff;
        //    cout<<d_total_eff<<endl;
            d_preventdoublecount=1;
            //dcounts++;
        }
    }
//////////////////////////////////////////v
    for(double ds_bdtval=ds_bdt_min;ds_bdtval<ds_bdt_max+.005;ds_bdtval+=.01){

//cout<<ds_bdtval<<endl;
        binID=hh->FindFixBin(ds_pt,ds_chi2,ds_bdtval);
        if(binID==binID2)
            continue;
        binID2=binID;
 
        ds_total+=hh->GetBinContent(binID);
        if(mcout<ds_bdtval)
            ds_partotal+=hh->GetBinContent(binID);
        

        if(fabs(ds_bdtval-ds_bdt_max)<.005 && !ds_preventdoublecount){
            if (ds_total>0) ds_eff=ds_partotal/ds_total;
            else {ds_eff=-1;cout<<"bug?"<<endl;}
            if (ds_eff>1) ds_eff=1;
            ds_total_eff+=ds_eff;
            ds_preventdoublecount=1;
            //dscounts++;
        }
 

    }


}
double d_out=d_total_eff/n;
double ds_out=ds_total_eff/n;
g1<<mcout<<"\t"<<d_out<<endl;
g2<<mcout<<"\t"<<ds_out<<endl;
//fout.cd();
//tout.Write();
}


g1.close();
g2.close();

//Here is some code to look at the projections. X:PT,Y:CHI2,Z:BDT
//cout<<h->Print("all")<<endl;>readme.txt
//
//
/*TCanvas *c=new TCanvas("c","");
h->Project3D("z")->Draw();
c->SetTitle("BDT Projection");
c->Print("/sum/opt/bdtprojection.pdf");

h->Project3D("y")->Draw();
c->SetTitle("chi2 Projection");
c->Print("/sum/opt/chi2projection.pdf");

h->Project3D("x")->Draw();
c->SetTitle("PT Projection");
c->Print("/sum/opt/ptprojection.pdf");*/

}
