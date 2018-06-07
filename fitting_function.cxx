void fitting_function(){

// Getting dataset deltam 

using namespace RooFit;
TFile *_file0 = TFile::Open("DaVinciTuples_S24r1_part_Kpipi0_Dataset.root");

// Defining dataset and variable deltam
RooDataSet* data = (RooDataSet*) _file0->Get("data");
const RooArgSet* deltam1 = data->get(0);
RooRealVar* deltam = (RooRealVar*)&(*deltam1)["deltam"];

// Plotting dataset
RooPlot* frame = deltam->frame();
//data->plotOn(frame);
//frame->Draw();

// Setting up gaussian 1 and 2 for signal
RooRealVar mean("mean","mean",0,140,170) ;
RooRealVar sigma("sigma","sigma",3,0.1,10) ; //arbitrary atm
RooGaussian gauss1("gauss1","signal component 1",*deltam,mean,sigma) ;

RooRealVar mean2("mean2","mean2",0,140,170) ;
RooRealVar sigma2("sigma2","sigma2",3,0.1,4) ; //arbitrary atm
RooGaussian gauss2("gauss2","signal component 2",*deltam,mean2,sigma2) ;

// Sum signal components to single PDF
RooRealVar sigfrac("sigfrac", "composite signal fraction", 0.5, 0. , 1.);  //signal fraction ; 
RooAddPdf sig("sig","Signal", RooArgList(gauss1,gauss2),sigfrac);
 
//setting up bifurcated gaussian for signal
/* 
  RooRealVar mean_bi("mean_bi","mean_bi",0,140,170) ;
  RooRealVar sigma_L_bi("sigma_L_bi","sigma_L_bi",3,0.1,10) ;
  RooRealVar sigma_R_bi("sigma_R_bi","sigma_R_bi",3,0.1,10) ;
  RooBifurGauss gauss_bi("gauss_bi","gauss_bi",*deltam,mean_bi,sigma_L_bi, sigma_R_bi) ;
  gauss_bi.fitTo(data) ;
  gauss_bi.plotOn(frame) ; 
*/

// Build modified exponential background PDF
RooRealVar M_th("M_th","M_th",139,170) ;
M_th.setVal(139.57);   // pi plus/minus mass as threshold value
M_th.setConstant(); 
RooRealVar a("a","a",0.5,0,500) ;
RooRealVar b("b","b",1,-100,100) ;
// a.setConstant(kTRUE); // For tinkering with only one variable

// Creating generic PDF for modified exponential background
 RooArgList(*deltam,M_th,a,b);
 RooGenericPdf background("background","Modified-Exponential Background plot","(1-exp(-(deltam-M_th)/a))*(deltam/M_th)^b",RooArgList(*deltam,M_th,a,b)) ;

 // Sum composite signal and background components

RooRealVar fsig("fsig", "total signal fraction", 0.5, 0. , 1.);  //signal fraction ; 

// model(x) = fsig*sig(x) +(1 - fsig)*background(x) 
RooAddPdf model("model","model", RooArgList(sig, background),fsig);  
//also include bifuracted gaussian

 model.Print("t");
 /* 
 model.fitTo(*data);
 model.plotOn(frame);
 frame->Draw();
 */
}
