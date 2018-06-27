//Input data from decay times (distiguished for D0 and D0Bar
//Input data from D0_M and deltam for both types
//Produce Pdfs to find the shape parameters for the D0_M and deltam fits by plotting simultaneously
//Fix these paramaters and leave the yield variables as floats
//Put the decay time data into bins width 0.5 originally
//Loop over these bins with the shape parameters to find the yield parameters that fit the data. The yeild parameters are plotted to the two mass variables but obviously only the ones that are in the histoframe decay range
//Plot two histograms (D0 and D0bar) containing the yield values for each decay time and their uncertainties
//Use this information to find the assymetry from: (D0-D0bar)/(D0+D0bar)
#include <iostream>
#include <vector>
#include "RooDataSet.h" 
#include "RooRealVar.h"
#include "RooAbsArg.h"
#include "math.h"
void Time_dependence_Kpipi0(){

  using namespace RooFit;
  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  using namespace std;

  //                        Data Input                          //


  // Getting D0 decay time data
  TFile *_file1 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/Kpipi0_DaVinciDev_v44r3/DaVinci_Kpipi0_D0_decaytime_Dataset.root ");

  // Defining dataset and D0 decay time
  RooDataSet* data = (RooDataSet*) _file1->Get("data");
  const RooArgSet* decay_D0 = data->get(0);
  RooRealVar* D0_decay = (RooRealVar*)&(*decay_D0)["ctau"];

  // Getting D0_BAR decay time data
  TFile *_file2 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/Kpipi0_DaVinciDev_v44r3/DaVinci_Kpipi0_D0BAR_decaytime_Dataset.root ");

  // Defining dataset and D0_BAR decay time
  RooDataSet* data_BAR = (RooDataSet*) _file2->Get("data");
  const RooArgSet* decay_BAR = data_BAR->get(0);
  RooRealVar* D0_BAR_decay = (RooRealVar*)&(*decay_BAR)["ctau_bar"];

  // Getting real data Delta Mass and D0 Mass
  TFile *_file3 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/Kpipi0_DaVinciDev_v44r3/DaVinci_Kpipi0_newBDT_Dataset.root ");

  // Defining dataset and variables D0_M and deltam
  RooDataSet* data_mass = (RooDataSet*) _file3->Get("data");
  const RooArgSet* deltam1 = data_mass->get(0);
  RooRealVar* deltam = (RooRealVar*)&(*deltam1)["deltam"];
  RooRealVar* D0_M = (RooRealVar*)&(*deltam1)["D0_M"];


  //                          SIGNAL: SHAPE PARAMETER                             //


  // Creating guassians and bifurcated guassian signal pdfs for deltam 
  RooRealVar mean_del("mean_del","mean_del",145,140,150) ; 
  RooRealVar sigma_del("sigma_del","sigma_del",1,0.,10.) ; 
  RooGaussian gauss_del("gauss_del","gauss_del",*deltam,mean_del,sigma_del) ; 

  RooRealVar mean_del2("mean_del2","mean_del2",145,140,150) ; 
  RooRealVar sigma_del2("sigma_del2","sigma_del2",1,0.,10.) ; 
  RooGaussian gauss_del2("gauss_del2","gauss_del2",*deltam,mean_del2,sigma_del2) ;

  RooRealVar mean_bi("mean_bi","mean_bi",146,140,150) ;
  RooRealVar sigma_L_bi("sigma_L_bi","sigma_L_bi",2.5,0.01,5) ;
  RooRealVar sigma_R_bi("sigma_R_bi","sigma_R_bi",2.5,0.01,5) ;
  RooBifurGauss gauss_del_bi("gauss_del_bi","gauss_del_bi",*deltam,mean_bi,sigma_L_bi, sigma_R_bi) ;

  // Total delta signal (two gaussians + bifurcated gaussian)
  RooRealVar cdel_1("cdel_1", "", 0.5, 0, 1) ;
  RooRealVar cdel_2("cdel_2", "", 0.5, 0, 1) ;
  RooAddPdf del_sig("del_sig","", RooArgList(gauss_del,gauss_del2,gauss_del_bi), RooArgList( cdel_1, cdel_2), true);

  // Creating guassians signal pdf for D0
  RooRealVar mean_D0M("mean_D0M","mean_D0M",1850,1710,2000) ; 
  RooRealVar sigma_D0M("sigma_D0M","sigma_D0M",5,0.,80.) ; 
  RooGaussian gauss_D0M("gauss_D0M","gauss_D0M",*D0_M,mean_D0M,sigma_D0M) ; 

  //RooRealVar mean_D0M2("mean_D0M2","mean_D0M2",1856,1850,) ; 
  RooRealVar sigma_D0M2("sigma_D0M2","sigma_D0M2",5,0.,80.) ; 
  RooGaussian gauss_D0M2("gauss_D0M2","gauss_D0M2",*D0_M,mean_D0M,sigma_D0M2) ; //using same mean as gauss_D0M

  // Total D0 signal (two gaussians)
  RooRealVar cD0_1("cD0_1", "", 0.5, 0, 1) ;
  RooAddPdf D0_sig("D0_sig","", RooArgList(gauss_D0M,gauss_D0M2), cD0_1);

  //                    BACKGROUND: SHAPE PARAMETER                  // 


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


  //                  COMBINING: SHAPE PARAMETER           //


  // Total signal (deltam: 2 gauss + bifurcated , D0: 2 gauss)
  RooProdPdf sig_pdf("sig_pdf","sig_pdf",RooArgSet(del_sig,D0_sig));

  // Total Combinatorial background (deltam: modified exponential , D0: straight line)
  RooProdPdf com_bkgpdf("com_bkgpdf","com_bkgpdf",RooArgSet(bkg_del,D0_com_bkg));

  // Total Random pi background (deltam: modified exponential, D0: same gaussian as signal
  RooProdPdf ranpi_bkgpdf("ranpi_bkgpdf","ranpi_bkgpdf",RooArgSet(bkg_del,D0_sig));

  //Here the bkg_del pdf is used for both the random pi and combinatorial background

  //                   FINAL  PDF: SHAPE PARAMETER               //


  RooRealVar nsignal("nsignal", "", data_mass->numEntries()*0.3, 0, data_mass->numEntries()*1.5) ;
  RooRealVar n_com_bkg("n_com_bkg", "", data_mass->numEntries()*0.35, 0, data_mass->numEntries()*1.5) ;
  RooRealVar n_ranpi_bkg("n_ranpi_bkg", "", data_mass->numEntries()*0.35, 0, data_mass->numEntries()*1.5) ;

  // Sum of signal, random pi and combinatorial background components
  RooAddPdf totalpdf("totalpdf", "", RooArgList(sig_pdf, com_bkgpdf, ranpi_bkgpdf), RooArgList(nsignal, n_com_bkg,n_ranpi_bkg)) ;


  //                  FITTING AND PLOTTING: SHAPE PARAMETER              //


  //Fitting and printing plots
  totalpdf.fitTo(*data_mass);

 
  // auto c = new TCanvas();
  // RooPlot* D0_frame = D0_M->frame();
  // data_mass->plotOn(D0_frame);
  // totalpdf.plotOn(D0_frame);
  // totalpdf.plotOn(D0_frame, Components("sig_pdf"), LineColor(kBlue), LineStyle(kDashed)) ;
  // totalpdf.plotOn(D0_frame, Components("gauss_D0M"), LineColor(kGreen), LineStyle(kDashed)) ;
  // totalpdf.plotOn(D0_frame, Components("gauss_D0M2"), LineColor(kOrange), LineStyle(kDashed)) ;
  // totalpdf.plotOn(D0_frame, Components("com_bkgpdf"), LineColor(kMagenta), LineStyle(kDashed)) ;
  // totalpdf.plotOn(D0_frame, Components("ranpi_bkgpdf"), LineColor(kRed), LineStyle(kDashed)) ;
  // D0_frame->Draw();

  // auto c2 = new TCanvas();
  // RooPlot* del_frame = deltam->frame();
  // data_mass->plotOn(del_frame);
  // totalpdf.plotOn(del_frame);

  // totalpdf.Print("t");
  // totalpdf.plotOn(del_frame, Components("sig_pdf"), LineColor(kBlue), LineStyle(kDashed)) ;
  // totalpdf.plotOn(del_frame, Components("gauss_del_bi"), LineColor(kCyan), LineStyle(kDashed)) ;
  // totalpdf.plotOn(del_frame, Components("gauss_del"), LineColor(kGreen), LineStyle(kDashed)) ;
  // totalpdf.plotOn(del_frame, Components("gauss_del2"), LineColor(kOrange), LineStyle(kDashed)) ;
  // totalpdf.plotOn(del_frame, Components("com_bkgpdf"), LineColor(kMagenta), LineStyle(kDashed)) ;
  // totalpdf.plotOn(del_frame, Components("ranpi_bkgpdf"), LineColor(kRed), LineStyle(kDashed)) ;
  // del_frame->Draw();


  // auto c3 = new TCanvas();
  // RooPlot* D0_BAR_decay_frame = D0_BAR_decay->frame();
  // data_BAR->plotOn(D0_BAR_decay_frame);
  // D0_BAR_decay_frame->Draw();

  // auto c4 = new TCanvas();
  // RooPlot* D0_decay_frame = D0_decay->frame();
  // data->plotOn(D0_decay_frame);
  // D0_decay_frame->Draw();


  //                    GETTING COMBINED DATA                            //


  // Getting combined data for D0
  TFile *_file4 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/Kpipi0_DaVinciDev_v44r3/DaVinci_Kpipi0_D0_FULL_Dataset.root");

  // Defining dataset and variables D0_M, deltam and ctau for D0
  RooDataSet* data_D0_comb = (RooDataSet*) _file4->Get("data");
  const RooArgSet* variables_D0 = data_D0_comb->get(0);
  RooRealVar* ctau_D0 = (RooRealVar*)&(*variables_D0)["ctau"];
  RooRealVar* deltam_D0 = (RooRealVar*)&(*variables_D0)["deltam"];
  RooRealVar* D0_M_D0 = (RooRealVar*)&(*variables_D0)["D0_M"];

  // Getting combined data for D0_bar
  TFile *_file5 = TFile::Open("/home/ppe/l/ldickson/documents/lhcb_project/small_sample/analysis/Kpipi0_DaVinciDev_v44r3/DaVinci_Kpipi0_D0BAR_FULL_Dataset.root");

  // Defining dataset and variables D0_M, deltam and ctau for D0_bar
  RooDataSet* data_D0BAR_comb = (RooDataSet*) _file5->Get("data");
  const RooArgSet* variables_D0BAR = data_D0BAR_comb->get(0);
  RooRealVar* ctau_D0BAR = (RooRealVar*)&(*variables_D0BAR)["ctau_bar"];
  RooRealVar* deltam_D0BAR = (RooRealVar*)&(*variables_D0BAR)["deltam"];
  RooRealVar* D0_M_D0BAR = (RooRealVar*)&(*variables_D0BAR)["D0_M"];


 //                    REDUCING DATASETS AND CREATING VECTORS POINTING TO THESE FOR D0 AND D0BAR                            //
  //Declaring data subsets D0
  RooDataSet* d1 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>= 0 && ctau < 0.5");
  RooDataSet* d2 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=0.5 && ctau<1.");
  RooDataSet* d3 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=1. && ctau <1.5");
  RooDataSet* d4 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=1.5 && ctau <2.");
  RooDataSet* d5 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=2. && ctau <2.5");
  RooDataSet* d6 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=2.5 && ctau<3.");
  RooDataSet* d7 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=3. && ctau<3.5");
  RooDataSet* d8 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=3.5 && ctau <4.");
  RooDataSet* d9 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=4. && ctau <4.5");
  RooDataSet* d10 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *D0_M_D0, *ctau_D0),"ctau>=4.5 && ctau <5.");

  //Creating data subsets vector D0
  vector<RooDataSet*> dataset_d0_vect;
  dataset_d0_vect.push_back(d1);
  dataset_d0_vect.push_back(d2);
  dataset_d0_vect.push_back(d3);
  dataset_d0_vect.push_back(d4);
  dataset_d0_vect.push_back(d5);
  dataset_d0_vect.push_back(d6);
  dataset_d0_vect.push_back(d7);
  dataset_d0_vect.push_back(d8);
  dataset_d0_vect.push_back(d9);
  dataset_d0_vect.push_back(d10);
 
  // int y;  //for checking that the vector was created properly
  //  cout << "d0 dataset contains:";
  //  for (y=0; y<dataset_d0_vect.size(); y++)
  //    cout << " " << dataset_d0_vect.at(y);

  //  cout << endl;

  //Declaring data subsets D0_bar
  RooDataSet* dbar1 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=0 && ctau_bar<0.5");
  RooDataSet* dbar2 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=0.5 && ctau_bar<1.");
  RooDataSet* dbar3 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=1. && ctau_bar<1.5");
  RooDataSet* dbar4 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=1.5 && ctau_bar<2.");
  RooDataSet* dbar5 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=2. && ctau_bar<2.5");
  RooDataSet* dbar6 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=2.5 && ctau_bar<3.");
  RooDataSet* dbar7 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=3. && ctau_bar<3.5");
  RooDataSet* dbar8 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=3.5 && ctau_bar<4.");
  RooDataSet* dbar9 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=4. && ctau_bar<4.5");
  RooDataSet* dbar10 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *D0_M_D0BAR, *ctau_D0BAR),"ctau_bar>=4.5 && ctau_bar<5.");
 
  //Creating data subsets vector D0_bar
  vector<RooDataSet*> dataset_d0bar_vect;
  dataset_d0bar_vect.push_back(dbar1);
  dataset_d0bar_vect.push_back(dbar2);
  dataset_d0bar_vect.push_back(dbar3);
  dataset_d0bar_vect.push_back(dbar4);
  dataset_d0bar_vect.push_back(dbar5);
  dataset_d0bar_vect.push_back(dbar6);
  dataset_d0bar_vect.push_back(dbar7);
  dataset_d0bar_vect.push_back(dbar8);
  dataset_d0bar_vect.push_back(dbar9);
  dataset_d0bar_vect.push_back(dbar10);
 
  //for checking that the vector was created properly
  // int z;  
  //   cout << "d0 bar contains:";
  //   for (z=0; z<dataset_d0bar_vect.size(); z++)
  //     cout << " " << dataset_d0bar_vect.at(z);

  //   cout << endl;




  //                    SETTING SHAPE PARAMETERS CONSTANT                          //
  //looping over paramaters setting them constant other than the yeild values for D0 and D0bar
  RooArgSet* vars = totalpdf.getParameters(*data_mass);
  RooFIter iter=vars->fwdIterator() ;
  // vars->Print("v");
  RooAbsArg* arg;
  RooRealVar *rrvar;
  vector<RooRealVar*> d0_vect;
  int i =0;
  while((arg=iter.next())) {
    rrvar = dynamic_cast<RooRealVar*>(arg);
    // arg->Print("v");
    if(rrvar){
      switch(i) {
      case 11:{    //n_com_bkg
	break;}
      case 12:{    //n_ranpi_bkg 
	break;}
      case 13: {  //nsignal
	break;}
      default: {
	rrvar->setConstant(); 
	d0_vect.push_back(rrvar); }         
      }
    }
    else { cout << "issue with rrvar"; }
    i++;
    
  }

  //printing contents of d0_vect (pointers)
  // int f;
  // int d0_size = dataset_d0_vect.size();
  // cout << "myvector contains:";
  // for (f=0; f<=d0_size; f++)
  //   cout << " " << d0_vect.at(f);

  // cout << endl;

  //                    FITTING NEW TIME-DEPENDANT SIGNAL VALUES AND PUTTING INTO VECTOR                            //
  // D0 decay time itteration over subdatasets from vector
  int k;
  int data_d0_size=dataset_d0_vect.size(); 
  RooRealVar *nsignal_vals_D0 ;
  vector<double> nsignal_D0_vect_double;
  vector<double> nsignal_D0_error_vect_double;
  RooAbsArg* coefficient;
  double nsignal_D0_double;
  double nsignal_D0_error_double;
  for (k=0; k<data_d0_size; k++) {
    totalpdf.fitTo(*dataset_d0_vect[k]);
    coefficient = totalpdf.coefList().at(0);  //nsignal(0), n_com_bkg(1), n_randpi_bkg(2) after removing shape parameters
    nsignal_vals_D0 =dynamic_cast<RooRealVar*>(coefficient);
    nsignal_D0_double = nsignal_vals_D0->getVal();
    nsignal_D0_error_double =  nsignal_vals_D0->getError();

    if(nsignal_vals_D0){
      nsignal_D0_vect_double.push_back(nsignal_D0_double);
      nsignal_D0_error_vect_double.push_back(nsignal_D0_error_double);
    }
    else { cout << "ISSUE WITH CREATING SIGNAL VECTOR- NSIGNAL_D0_VECT"; }
  }

  //D0_bar decay time itteration over subdatasets from vector
  int p; 
  int data_D0bar_size = dataset_d0bar_vect.size();
  RooRealVar *nsignal_vals_D0bar ;
  vector<double> nsignal_D0bar_vect_double;
  vector<double> nsignal_D0bar_error_vect_double;
  double nsignal_D0bar_double;
  double nsignal_D0bar_error_double;
  for (p=0; p<data_D0bar_size; p++) {
    totalpdf.fitTo(*dataset_d0bar_vect[p]);
    coefficient = totalpdf.coefList().at(0);  //nsignal(0), n_com_bkg(1), n_randpi_bkg(2) after removing shape parameters
    nsignal_vals_D0bar =dynamic_cast<RooRealVar*>(coefficient);
    nsignal_D0bar_double = nsignal_vals_D0bar->getVal();
    nsignal_D0bar_error_double =  nsignal_vals_D0bar->getError();

    if(nsignal_vals_D0bar){
      nsignal_D0bar_vect_double.push_back(nsignal_D0bar_double);
      nsignal_D0bar_error_vect_double.push_back(nsignal_D0bar_error_double);
    }
    else { cout << "ISSUE WITH CREATING SIGNAL VECTOR- NSIGNAL_D0_VECT"; }
  }

  // CHECKING IF SIGNAL VECTOR IS PRODUCED CORRECTLY 
  int o;
  cout << "myvector contains:";
  for (o=0; o<10; o++)
    cout << " " <<nsignal_D0bar_error_vect_double.at(o);

  cout << endl;


  //                 PLOTTING DATA ON HIST              //
 
  //  D0  and D0bar  //
  const Int_t NBINS = 10;
  Double_t edges[NBINS + 1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};

  TH1D* D0_sig_hist = new TH1D("D0_sig_hist","Hist with variable bin width", NBINS, edges); 

  TH1D* D0bar_sig_hist = new TH1D("D0bar_sig_hist","Hist with variable bin width", NBINS, edges);
  int q;
  D0_sig_hist->Sumw2();
  D0bar_sig_hist->Sumw2();
  for (q=0; q<NBINS; q++) {
    // D0 //
    nsignal_D0_double = nsignal_D0_vect_double[q];
    //   nsignal_D0_error_vect_double = nsignal_D0_error_vect_double[q];
    D0_sig_hist->SetBinContent(q+1,nsignal_D0_double);
    //  D0_sig_hist->SetError(nsignal_D0_error_vect_double);

    // D0 BAR //
    nsignal_D0bar_double = nsignal_D0bar_vect_double[q];
    //nsignal_D0bar_error_vect_double = nsignal_D0bar_error_vect_double[q];
    D0bar_sig_hist->SetBinContent(q+1,nsignal_D0bar_double);
    //D0bar_sig_hist->SetBinError(q+1 )			     
  }  
  auto c6 =  new TCanvas();
  D0_sig_hist->Draw();

  auto c7 = new TCanvas();
  D0bar_sig_hist->Draw();

  TH1D *h3 = (TH1D*)D0_sig_hist->Clone("h3");
  TH1D *h4 = (TH1D*)D0_sig_hist->Clone("h4");
  h3->Add(D0bar_sig_hist,-1);
  h4->Add(D0bar_sig_hist);
  h3->Divide(h4);

  //     CALCULATING ERRORS     //
  //err = sqrt( errD0^2*(dD0(f))^2 + errD0bar^2*(dD0bar(f))^2 )

  int u;  
  double dD0_denom;
  double dD0_numer;
  double dD0_differential;
  double dD0_term1;
  double dD0bar_denom;
  double dD0bar_numer;
  double dD0bar_differential;
  double dD0bar_term2;
  double error_square_sum;
  vector<double> error_vect;
  for (u=0; u<data_d0_size; u++) {
    // D0 partial differntial squared  
    dD0_numer = 4*pow(nsignal_D0bar_vect_double[u],2);
    dD0_denom = pow(nsignal_D0_vect_double[u]+ nsignal_D0bar_vect_double[u],4);
    dD0_differential = dD0_numer/dD0_denom;
    dD0_term1 = dD0_differential*pow(nsignal_D0_error_vect_double[u],2);

    //D0bar partial differential squared
    dD0bar_numer = 4*pow(nsignal_D0_vect_double[u],2);
    dD0bar_denom = pow(nsignal_D0_vect_double[u]+ nsignal_D0bar_vect_double[u],4);
    dD0bar_differential = dD0bar_numer/dD0bar_denom;
    dD0bar_term2 =  dD0bar_differential*pow(nsignal_D0bar_error_vect_double[u],2);

    // Total error
    error_square_sum = sqrt(dD0_term1+dD0bar_term2);
    error_vect.push_back(error_square_sum);
    //  h3->SetBinError(u , error_square_sum);
 }

 int r;
  cout << "myvector contains:";
  for (r=0; r<10; r++)
    cout << " " <<error_vect.at(r);

  cout << endl;


  auto c8 = new TCanvas();
  h3->Draw();


  // FOR SEEING THE REDUCED DATA PLOTS 
  // RooPlot* example_frame = ctau_D0->frame();
  // d1->plotOn(example_frame);
  // totalpdf.fitTo(*d1);
  // totalpdf.plotOn(example_frame);
  // example_frame->Draw();




}
