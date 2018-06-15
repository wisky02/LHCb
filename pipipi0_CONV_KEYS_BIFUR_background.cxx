void pipipi0_CONV_KEYS_BIFUR_background(){

using namespace RooFit;

// Getting Monte-Carlo data set
TFile *_file2 = TFile::Open("./DaVinciTuples_MC_S28_Matched_pipipi0_Dataset.root");

// Defining dataset and variable deltam
RooDataSet* data2 = (RooDataSet*) _file2->Get("data");
const RooArgSet* deltam3 = data2->get(0);
RooRealVar* deltam2 = (RooRealVar*)&(*deltam3)["deltam"];

// Defining KEYS pdf
RooKeysPdf kest2("kest2","kest2",*deltam2,*data2,RooKeysPdf::NoMirror);

// Setting up bifurcated gaussian for signal
RooRealVar mean_bi("mean_bi","mean_bi",145.4,140,170) ;
RooRealVar sigma_L_bi("sigma_L_bi","sigma_L_bi",3,0.1,15) ;
RooRealVar sigma_R_bi("sigma_R_bi","sigma_R_bi",3,0.1,15) ;
RooBifurGauss gauss_bi("gauss_bi","gauss_bi",*deltam2,mean_bi,sigma_L_bi, sigma_R_bi) ;

//Getting real data set
TFile *_file0 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/pipipi0_DaVinciDev_v44r3/DaVinciTuples_S24r1_part_pipipi0_Dataset.root");

// Defining dataset and variable deltam
RooDataSet* data = (RooDataSet*) _file0->Get("data");
const RooArgSet* deltam1 = data->get(0);
RooRealVar* deltam = (RooRealVar*)&(*deltam1)["deltam"];

// Plotting dataset
RooPlot* frame = deltam->frame();
data->plotOn(frame);
frame->Draw();

// Build modified exponential background PDF
RooRealVar M_th("M_th","M_th",139,170) ;
M_th.setVal(139.57);   // pi plus/minus mass as threshold value
M_th.setConstant(); 
RooRealVar a("a","a",20,0,60) ;
RooRealVar b("b","b", -20,-100,100) ;

// Creating generic PDF for modified exponential background
RooArgList(*deltam,M_th,a,b);
RooGenericPdf bkgpdf("bkgpdf","Modified-Exponential Background plot","(1-exp(-(deltam-M_th)/a))*(deltam/M_th)^b",RooArgList(*deltam,M_th,a,b)) ;
 RooRealVar nsignal("nsignal", "", data->numEntries()*0.3, 0, data->numEntries()*1.5) ; 
 RooRealVar nbkg("nbkg", "", data->numEntries()*0.7, 0, data->numEntries()*1.5) ;

// Creating convoluting guassing
RooRealVar mean("mean","mean",0) ;
mean.setConstant();
RooRealVar sigma("sigma","sigma",0.1,0,2) ;
RooGaussian gauss_conv("gauss_conv","gauss_conv",*deltam2,mean,sigma) ;

// Construct convolution pdf
deltam->setBins(10000,"cache");
RooFFTConvPdf conv_pdf("conv_pdf","KEYS (X) gauss",*deltam2,kest2,gauss_conv) ;

// Signal pdf (MC + convoluting gaussian + bifurcated gaussian)
RooRealVar cgaus1("cgaus1", "", 0.5, 0, 1) ;
RooAddPdf signal_pdf("signal_pdf","", RooArgList(conv_pdf,gauss_bi), cgaus1, true);

// Sum signal and background components 
// model(x) = fsig*sig(x) +(1 - fsig)*background(x) 
RooAddPdf totalpdf("totalpdf", "", RooArgList(signal_pdf, bkgpdf), RooArgList(nsignal, nbkg)) ;

// Fitting and printing plots
totalpdf.fitTo(*data);
totalpdf.plotOn(frame);
totalpdf.Print("t");

totalpdf.plotOn(frame, Components("kest2"), LineColor(kRed), LineStyle(kDashed)) ;
totalpdf.plotOn(frame, Components("conv_pdf"), LineColor(kGreen), LineStyle(kDashed)) ;
totalpdf.plotOn(frame, Components("gauss_bi"), LineColor(kCyan), LineStyle(kDashed)) ;
totalpdf.plotOn(frame, Components("signal_pdf"), LineColor(kOrange), LineStyle(kDashed)) ;
totalpdf.plotOn(frame, Components("bkgpdf"), LineColor(kCyan), LineStyle(kDashed)) ;
frame->Draw();
 
}
