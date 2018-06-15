#include "TH2F.h"
#include "RooProdPdf.h"
void D0_deltam_sim_plot(){

using namespace RooFit;
    gStyle->SetPalette(kBird);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

// Getting real data Delta Mass and D0 Mass
TFile *_file1 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/simultaneous_plot/DaVinciTuples_S24r1_D0M_deltam_Dataset.root");

// Defining dataset and variable deltam
RooDataSet* data = (RooDataSet*) _file1->Get("data");
const RooArgSet* deltam1 = data->get(0);
RooRealVar* deltam = (RooRealVar*)&(*deltam1)["deltam"];
RooRealVar* D0_M = (RooRealVar*)&(*deltam1)["D0_M"];


//                          SIGNAL                             //


// Creating guassians and bifurcated guassian signal pdfs for deltam 
RooRealVar mean_del("mean_del","mean_del",145,140,170) ; 
RooRealVar sigma_del("sigma_del","sigma_del",1,0.,10.) ; 
RooGaussian gauss_del("gauss_del","gauss_del",*deltam,mean_del,sigma_del) ; 

RooRealVar mean_del2("mean_del2","mean_del2",145,140,170) ; 
RooRealVar sigma_del2("sigma_del2","sigma_del2",1,0.,10.) ; 
RooGaussian gauss_del2("gauss_del2","gauss_del2",*deltam,mean_del2,sigma_del2) ;

RooRealVar mean_bi("mean_bi","mean_bi",146,140,150) ;
RooRealVar sigma_L_bi("sigma_L_bi","sigma_L_bi",2.5,0.01,5) ;
RooRealVar sigma_R_bi("sigma_R_bi","sigma_R_bi",2.5,0.01,5) ;
RooBifurGauss gauss_del_bi("gauss_del_bi","gauss_del_bi",*deltam,mean_bi,sigma_L_bi, sigma_R_bi) ;

// Total delta, signal (two gaussians + bifurcated gaussian)
RooRealVar cdel_1("cdel_1", "", 0.5, 0, 1) ;
RooRealVar cdel_2("cdel_2", "", 0.5, 0, 1) ;
RooAddPdf del_sig("del_sig","", RooArgList(gauss_del,gauss_del2,gauss_del_bi), RooArgList( cdel_1, cdel_2), true);

// Creating guassians signal pdf for D0
RooRealVar mean_D0M("mean_D0M","mean_D0M",1850,1710,2000) ; 
RooRealVar sigma_D0M("sigma_D0M","sigma_D0M",5,0.,80.) ; 
RooGaussian gauss_D0M("gauss_D0M","gauss_D0M",*D0_M,mean_D0M,sigma_D0M) ; 

//RooRealVar mean_D0M2("mean_D0M2","mean_D0M2",1856,1850,) ; 
RooRealVar sigma_D0M2("sigma_D0M2","sigma_D0M2",5,0.,80.) ; 
RooGaussian gauss_D0M2("gauss_D0M2","gauss_D0M2",*D0_M,mean_D0M,sigma_D0M2) ;

// Total D0 signal (two gaussians)
RooRealVar cD0_1("cD0_1", "", 0.5, 0, 1) ;
RooAddPdf D0_sig("D0_sig","", RooArgList(gauss_D0M,gauss_D0M2), cD0_1);

//                    BACKGROUND                  // 


// D0 combinatorial background (straight line)
RooRealVar a0("a0","a0",0.5,-1,1) ;
RooChebychev D0_com_bkg("D0_com_bkg","D0_com_bkg",*D0_M,a0) ;

// Deltam combinatorial background (modified exponential)
RooRealVar M_th_del("M_th_del","M_th_del",139,170) ;
M_th_del.setVal(139.57);   // pi plus/minus mass as threshold value
M_th_del.setConstant(); 
RooRealVar a_del("a_del","a_del",20,0,60) ;
RooRealVar b_del("b_del","b_del", -20,-100,100) ;
RooArgList(*deltam,M_th_del,a_del,b_del);
RooGenericPdf bkg_del("bkg_del","Combinatorial Background plot","(1-exp(-(deltam-M_th_del)/a_del))*(deltam/M_th_del)^b_del",RooArgList(*deltam,M_th_del,a_del,b_del)) ;


//                  COMBINING           //


// Total signal (deltam: 2 gauss + bifurcated , D0: 2 gauss)
RooProdPdf sig_pdf("sig_pdf","sig_pdf",RooArgSet(del_sig,D0_sig));

// Total Combinatorial background (deltam: modified exponential , D0: straight line)
RooProdPdf com_bkgpdf("com_bkgpdf","com_bkgpdf",RooArgSet(bkg_del,D0_com_bkg));

// Total Random pi background (deltam: modified exponential, D0: same gaussian as signal
RooProdPdf ranpi_bkgpdf("ranpi_bkgpdf","ranpi_bkgpdf",RooArgSet(bkg_del,D0_sig));

//Here the bkg_del pdf is used for both the random pi and combinatorial background

//                   FINAL  PDF               //


RooRealVar nsignal("nsignal", "", data->numEntries()*0.3, 0, data->numEntries()*1.5) ;
RooRealVar n_com_bkg("n_com_bkg", "", data->numEntries()*0.35, 0, data->numEntries()*1.5) ;
RooRealVar n_ranpi_bkg("n_ranpi_bkg", "", data->numEntries()*0.35, 0, data->numEntries()*1.5) ;

// Sum of signal, random pi and combinatorial background components
RooAddPdf totalpdf("totalpdf", "", RooArgList(sig_pdf, com_bkgpdf, ranpi_bkgpdf), RooArgList(nsignal, n_com_bkg,n_ranpi_bkg)) ;



//                  FITTING AND PLOTTING              //


//Fitting and printing plots
 totalpdf.fitTo(*data);

auto c = new TCanvas();
RooPlot* D0_frame = D0_M->frame();
data->plotOn(D0_frame);
totalpdf.plotOn(D0_frame);
totalpdf.plotOn(D0_frame, Components("sig_pdf"), LineColor(kBlue), LineStyle(kDashed)) ;
totalpdf.plotOn(D0_frame, Components("gauss_D0M"), LineColor(kGreen), LineStyle(kDashed)) ;
totalpdf.plotOn(D0_frame, Components("gauss_D0M2"), LineColor(kOrange), LineStyle(kDashed)) ;
totalpdf.plotOn(D0_frame, Components("com_bkgpdf"), LineColor(kMagenta), LineStyle(kDashed)) ;
totalpdf.plotOn(D0_frame, Components("ranpi_bkgpdf"), LineColor(kRed), LineStyle(kDashed)) ;
D0_frame->Draw();

auto c2 = new TCanvas();
RooPlot* del_frame = deltam->frame();
data->plotOn(del_frame);
totalpdf.plotOn(del_frame);

 totalpdf.Print("t");
 totalpdf.plotOn(del_frame, Components("sig_pdf"), LineColor(kBlue), LineStyle(kDashed)) ;
 totalpdf.plotOn(del_frame, Components("gauss_del_bi"), LineColor(kCyan), LineStyle(kDashed)) ;
totalpdf.plotOn(del_frame, Components("gauss_del"), LineColor(kGreen), LineStyle(kDashed)) ;
totalpdf.plotOn(del_frame, Components("gauss_del2"), LineColor(kOrange), LineStyle(kDashed)) ;
totalpdf.plotOn(del_frame, Components("com_bkgpdf"), LineColor(kMagenta), LineStyle(kDashed)) ;
totalpdf.plotOn(del_frame, Components("ranpi_bkgpdf"), LineColor(kRed), LineStyle(kDashed)) ;
del_frame->Draw();
}
