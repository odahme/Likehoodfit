
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooMinuitMCMC.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include <TApplication.h>
using namespace RooFit ;





int main(int argc, char **argv) {

  //create a window
  TApplication theApp("demoApplication",&argc,argv);
  //create a canvas

  TCanvas c1("c1","c1",1,1,1920,1080);
  //TCanvas c2("c2","c2",1,1,1920,1080);
  // TCanvas c3("c3","c3",1,1,1920,1080);


    // Observable
  RooRealVar x("x","x",-20,20) ;

  // Model (intentional strong correlations)
  RooRealVar mean_g1("mean1","mean of g2",3.0,-5.0,5.0) ;
  RooRealVar sigma_g1("sigma1","width of g1",1.0,0.0,3) ;
  RooGaussian g1("g1","g1",x,mean_g1,sigma_g1) ;

  RooRealVar mean_g2("mean2","mean of g1",-3.0,-5.0,5.0) ;
  RooRealVar sigma_g2("sigma2","width of g2",1.5,0.0,4.0) ;
  RooGaussian g2("g2","g2",x,mean_g2,sigma_g2) ;

  RooRealVar frac("frac","frac",0.5,0.0,1.0) ;
  RooAddPdf model("model","model",RooArgList(g1,g2),frac) ;

  // Generate 1000 events
  RooDataSet* data = model.generate(x,1000) ;

  // Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = model.createNLL(*data);

  RooMinuit mi(*nll);
  mi.minos();


  RooMinuitMCMC m(*nll);
  m.mcmc(2000,100);
  RooPlot* frame = x.frame();
  data->plotOn(frame);
  model.plotOn(frame);
  c1.cd();
  frame->Draw();
  c1.SaveAs("twogausfit.png");
  m.saveCornerPlotAs("twogauscorner.pdf");


  //turns off the program with mous clic
  theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
  //starts the canvas
  theApp.Run();

  return 1;
}
