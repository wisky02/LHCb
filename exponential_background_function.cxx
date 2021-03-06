/*// fitting and plotting to real data 
TFile file("DaVinciTuples_S24r1_part_pipipi0.root") 
using namespace RooFit ;
DecayTree->Draw("lab0_DTF_D0Mass_M >> h(100,2005,2030)")
TH1* hh = (TH1*) gDirectory->Get("h") 
RooRealVar lab0_DTF_D0Mass_M("lab0_DTF_D0Mass_M","lab0_DTF_D0Mass_M",2005,2030) ;
RooDataHist data("data","dataset with lab0_DTF_D0Mass_M",lab0_DTF_D0Mass_M,hh) ; 
RooPlot* frame = lab0_DTF_D0Mass_M.frame()
data.plotOn(frame)
*/
/////// getting dataset deltam /////////
using namespace RooFit;

  void import_data_set() {
const RooArgSet* deltam1 = data->get(0);
(*deltam1)["deltam"];
RooRealVar* deltam = (RooRealVar*)&(*deltam1)["deltam"];
RooPlot* frame = deltam->frame();
}
/*
data->plotOn(frame);
frame->Draw();
*/

///////////// Signal fitting and plotting /////////////////

// setting up gaussian 1 for signal
/*  RooRealVar mean("mean","mean",0,2005,2030) ;
  RooRealVar sigma("sigma","sigma",3,0.1,10) ;
  RooGaussian gauss("gauss","gauss",lab0_DTF_D0Mass_M,mean,sigma) ;
  gauss.fitTo(data) ;
  gauss.plotOn(frame) ;
*/
// setting up gaussian 2 for signal
/* RooRealVar mean2("mean2","mean2",0,2005,2030) ;
  RooRealVar sigma2("sigma2","sigma2",3,0.1,10) ;
  RooGaussian gauss2("gauss2","gauss2",x,mean2,sigma2) ;
  gauss2.fitTo(data) ;
  gauss2.plotOn(frame) ;
*/
//setting up bifurcated gaussian for signal
/*  RooRealVar mean_bi("mean_bi","mean_bi",0,2005,2030) ;
  RooRealVar sigma_L_bi("sigma_L_bi","sigma_L_bi",3,0.1,10) ;
  RooRealVar sigma_R_bi("sigma_R_bi","sigma_R_bi",3,0.1,10) ;
  RooBifurGauss gauss_bi("gauss_bi","gauss_bi",x,mean_bi,sigma_L_bi, sigma_R_bi) ;
  gauss_bi.fitTo(data) ;
  gauss_bi.plotOn(frame) ; 
*/
//Fitting modified exponential to background
void background_fit() {
RooRealVar M_th("M_th","M_th",139,170) ;
M_th.setVal(139.57);   // pi plus of minus mass as threshold value
M_th.setConstant(); 
RooRealVar a("a","a",2,0,500) ;
RooRealVar b("b","b",-20,-100,100) ;

///////////// Background fitting and plotting /////////////////

RooGenericPdf background("background","Modified-Exponential Background plot","(1-exp(-(deltam-M_th)/a))*(deltam/M_th)^b",RooArgList(*deltam,M_th,a,b)) 

 background.fitTo(data) ;
 background.plotOn(frame) 
}
//plotting
frame->Draw()
