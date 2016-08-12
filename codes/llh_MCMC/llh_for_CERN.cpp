/*
Likehood fit for CERN
*/

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TApplication.h>
#include "TMath.h"
#include <omp.h>
#include "TFile.h"
#include <fstream>
#include <string>
#include <TVectorD.h>
#include "TMatrix.h"
#include "TDecompChol.h"
#include "assert.h"
#include "TGraph.h"
#include <unistd.h>
using namespace std;

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooMinuitMCMC.h"
#include "RooFitResult.h"
#include "TH1.h"
using namespace RooFit ;

double f_gllh(double mu, double sigma, TVectorD gauss_points)
  {
    double calc = 0;
    for (int i = 0; i < gauss_points.GetNrows(); i++) {
      calc += log(TMath::Gaus(gauss_points[i],mu,sigma,true));
    }

    return -calc;
  }


void rf101_basics() {
  // S e t u p   m o d e l
  // ---------------------

  // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
  RooRealVar x("x","x",-10,10) ;
  RooRealVar mean("mean","mean of gaussian",1,-10,10) ;
  RooRealVar sigma("sigma","width of gaussian",1,0.1,10) ;

  // Build gaussian p.d.f in terms of x,mean and sigma
  RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma) ;

  // Construct plot frame in 'x'
  RooPlot* xframe = x.frame(Title("Gaussian p.d.f.")) ;


  // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
  // ---------------------------------------------------------------------------

  // Plot gauss in frame (i.e. in x)
  gauss.plotOn(xframe) ;

  // Change the value of sigma to 3
  sigma.setVal(3) ;

  // Plot gauss in frame (i.e. in x) and draw frame on canvas
  gauss.plotOn(xframe,LineColor(kRed)) ;


  // G e n e r a t e   e v e n t s
  // -----------------------------

  // Generate a dataset of 1000 events in x from gauss
  RooDataSet* data = gauss.generate(x,10000) ;

  // Make a second plot frame in x and draw both the
  // data and the p.d.f in the frame
  RooPlot* xframe2 = x.frame(Title("Gaussian p.d.f. with data")) ;
  data->plotOn(xframe2) ;
  gauss.plotOn(xframe2) ;


  // F i t   m o d e l   t o   d a t a
  // -----------------------------

  // Fit pdf to data
  gauss.fitTo(*data) ;

  // Print values of mean and sigma (that now reflect fitted values and errors)
  mean.Print() ;
  sigma.Print() ;

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("rf101_basics","rf101_basics",800,400) ;
  c->Divide(2) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.6) ; xframe->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; xframe2->GetYaxis()->SetTitleOffset(1.6) ; xframe2->Draw() ;


}



int main(int argc, char **argv) {

TRandom3 *rnd = new TRandom3(35);

//create a window
TApplication theApp("demoApplication",&argc,argv);
//create a canvas

TCanvas c1("c1","c1",1,1,1920,1080);
TCanvas c2("c2","c2",1,1,1920,1080);
//TCanvas c3("c3","c3",1,1,1024,768);

TH1F * hist = new TH1F("hist", "hist", 100, -2, 8);


//simulate gauss distributed points

unsigned long gp = 1000; //number of points
TVectorD gauss_points(gp);
double real_sigma = 1.5;
double real_mu = 3;
for (int j = 0; j < gauss_points.GetNrows(); j++) {
  gauss_points[j] = rnd->Gaus(real_mu,real_sigma);
  hist->Fill(gauss_points[j]);
}

//initialise

unsigned int nparams = 2; //number of parameters
unsigned int nstat = 1000;//number of tries
double maxstep = 0.01; //maximum step size
double ubmu = 4; //upper bound for mu
double lbmu = 1.5; //lower bound for mu
double ubsigma = 3; //upper bound for sigma
double lbsigma = 0.5; //lower bound for sigma
double alphastar = 0.234; //forced acceptance rate







ofstream rejectstream;
ofstream acceptstream;

Bool_t accepted;
unsigned int ntested = 0;
unsigned int naccepted = 0;

TVectorD last(nparams);
TVectorD curr(nparams);

unsigned int nplotp = 50; //number of points in plot
TVectorD llh_p(nstat);
TVectorD mu_p(nstat);
double minllh = 1e32;
double bestmu = 42;

unsigned int nlast = 200;
std::vector<bool> lastaccepted;

TMatrixDSym identity(nparams);
identity.UnitMatrix();

TMatrixDSym S1(nparams);
S1.Zero();
S1 = identity;

TMatrixDSym SNminusone(nparams);
SNminusone = S1;
TMatrixDSym SN(nparams);
SN.Zero();

double first_mu = rnd->Uniform(1,4);
double first_simga = rnd->Uniform(1,2);

last[0] = first_mu;
last[1] = first_simga;

//MCMC


rejectstream.open("list_of_rejected_candidates.txt");
acceptstream.open("List_of_accepted_candidates.txt");
for (unsigned int i = 0; i < nstat; i++) {
  curr = last; //use value of last for current then vary

  TVectorD WN(nparams);
  for (int j = 0; j < WN.GetNrows() ; j++) {
    WN[j] = rnd->Gaus(0.0, maxstep);
  }
  TVectorD SW(SNminusone*WN);
  curr += SW;

  if (curr[0] < lbmu) {
    curr[0] = lbmu;
  }
  if (curr[0] > ubmu) {
    curr[0] = ubmu;
  }
  if (curr[1] < lbsigma) {
    curr[1] = lbsigma;
  }
  if (curr[1] > ubsigma) {
  curr[1] = ubsigma;
  }




// std::cout << "sigma = "<< curr[1] << std::endl;

  //get llh for current and last
  double llh_last = f_gllh(last[0],last[1],gauss_points);
  double llh_curr = f_gllh(curr[0],curr[1],gauss_points);

if (llh_curr < minllh) {
  minllh = llh_curr;
  bestmu = curr[0];
}

  llh_p[i] = llh_curr;
  mu_p[i] = curr[0];

//plotting everz 50th point of the llh
/*
if (i % (nstat/nplotp) == 0) {
  double llhplotp[inplotp+1];
  double muplotp[inplotp+1];
  std::cout << "mu = "<< curr[0] << std::endl;
  std::cout << "llh = " << llh_curr << std::endl;
  llh_p[inplotp] = llh_curr;
  mu_p[inplotp] = curr[0];
  for (unsigned int j = 0; j < inplotp; j++) {
    llhplotp[j] = llh_p[j];
    muplotp[j] = mu_p[j];
  }

  TGraph * gr = new TGraph(inplotp,muplotp,llhplotp);
  c1.cd();
  gr->Draw("A*");
  c1.Update();
  inplotp++;
  usleep(1000000);
}
*/
//evolution of the llh
  //
  // if (i > nplotp) {
  //   double llhplotp[nplotp];
  //   double muplotp[nplotp];
  //   for (unsigned int inplotp = 0; inplotp < nplotp; inplotp++) {
  //     llhplotp[inplotp] = llh_p[i-inplotp];
  //     muplotp[inplotp] = mu_p[i-inplotp];
  //   }
  //
  //
  //   TGraph * gr = new TGraph(nplotp,muplotp,llhplotp);
  //   c2.cd();
  //   gr->Draw("A*");
  //   c2.Update();
  //   usleep(100000);
  // }



  double alpha = std::min(1.0, exp(llh_last - llh_curr));
//  std::cout << "alpha = " << alpha << std::endl;
  double r = rnd->Uniform(0,1);
  accepted = false;
  double acceptrate = double(naccepted)/double(ntested);


  if (r < alpha) {
    accepted = true;
    naccepted++;
    //success, candidate filled below
    last = curr;
  } else {
    //std::cout << "candidate found" << std::endl;
    accepted = false;
    //first save rejected candidate to file
    for (int line = 0; line < curr.GetNrows(); line++) {
      rejectstream << curr[line] << "\n";
      //std::cout << curr[line] << std::endl;
    }
    //then reset to last candidate
    curr = last;
  }

  for (int line = 0; line < curr.GetNrows(); line++) {
    acceptstream << curr[line] << "\n";
    //std::cout << curr[line] << std::endl;
  }



  lastaccepted.insert(lastaccepted.begin(), accepted);
  double lastratio = 0.0;
  if (lastaccepted.size() > nlast)
  {
  //lastaccepted.push_back(accepted);
  unsigned int nlastaccepted = 0;
  for (unsigned int j=0; j<lastaccepted.size(); j++)
  if (lastaccepted.at(j))
    nlastaccepted++;
  lastratio = double(nlastaccepted)/double(lastaccepted.size());
  lastaccepted.pop_back();
  }

  // if (accepted){
  //   std::cout << "set " << i << " accepted, rate " << std::fixed << std::setprecision(3) << acceptrate << " last " << nlast << " " << lastratio << std::endl;
  // }
  // else{
  //   std::cout << "set " << i << " rejected, rate " << std::fixed << std::setprecision(3) << acceptrate << " last " << nlast << " " << lastratio << std::endl;
  // }

  ntested++;

  //update S matrix
  TMatrixDSym SNminusoneT(SNminusone);
  SNminusoneT.T();
  double etan = std::min(1.0, nparams*pow(double(i), -2.0/3.0));
  TMatrixDSym WNWNT(nparams);
  WNWNT.Zero();
  for (int row = 0; row < WNWNT.GetNrows(); row++) {
    for (int col = 0; col < WNWNT.GetNcols(); col++) {
      WNWNT[row][col] = WN[row]*WN[col]/WN.Norm2Sqr();
    }
  }
  TMatrixDSym SNSNT(identity + WNWNT*etan*(alpha-alphastar));
  SNSNT = SNSNT.Similarity(SNminusone);

  //SNSNT = (SNminusone*identity*SNminusoneT);
  TDecompChol chol(SNSNT);
  bool success = chol.Decompose();
  assert(success);
  TMatrixD SNT = chol.GetU();
  TMatrixD SN(SNT);
  SN.T();
  for (int row = 0; row < SN.GetNrows(); row++) {
    for (int col = 0; col < SN.GetNcols(); col++) {
      SNminusone[row][col] = SN[row][col];
    }
  }
}

rejectstream.close();
acceptstream.close();
std::cout << ntested << " points were tested" << std::endl;
std::cout << naccepted << " points were accepted" << std::endl;
std::cout << "The accept fraction is " << double(naccepted)/double(ntested) << std::endl;

//checking error at 68% confidence

unsigned int npininterval = 0;
unsigned int trials = 0;
TVectorD mup_for_error(int(nstat*0.68));
while (double(npininterval)/double(nstat) < 0.68) {
  TVectorD imup_for_error(int(nstat*0.68));
  unsigned int npinintervalcurr = 0;
  double interval = double(trials)/(nstat*100.0);
  for (unsigned int j = 0; j < nstat; j++) {
    if (mu_p[j] < bestmu+interval && mu_p[j] > bestmu-interval) {
      imup_for_error[npinintervalcurr] = mu_p[j];
      npinintervalcurr++;
    }
  }
  trials++;
  npininterval = npinintervalcurr;
  mup_for_error = imup_for_error;
}
 double mue[2];
 mue[0] = -1e32;
 mue[1] = +1e32;
for (unsigned int i = 0; i < nstat*0.68; i++) {
  if (mup_for_error[i] > mue[0]) {
    mue[0] = mup_for_error[i];
  }
  if (mup_for_error[i] < mue[1]) {
    mue[1] = mup_for_error[i];
  }
}


// for (unsigned int i = 0; i < nstat; i++) {
//   llh_p[i] -= minllh;
// }
// double halfmin = 0.5;
// double halfmax = 0.5;
// Bool_t not2mu = true;
// while (not2mu) {
//   unsigned int imu = 0;
//   for (unsigned int i = 0; i < nstat; i++) {
//     if (llh_p[i] < halfmax && llh_p[i] > halfmin) {
//       imu++;
//     }
//   }
//   if (imu == 2) {
//     not2mu = false;
//   } else {
//     halfmin -= 0.001;
//     halfmax += 0.001;
//   }
//   std::cout << "imu = "<< imu << std::endl;
//   std::cout << "halfmax ="<< halfmax << std::endl;
// }
//

// unsigned int imu = 0;
// for (unsigned int i = 0; i < nstat; i++) {
//   if (llh_p[i] < halfmax && llh_p[i] > halfmin) {
//     mue[imu] = llh_p[i];
//     imu++;
//   }
// }



// TGraph *gr = new TGraph(mu_p,llh_p);
// c1.cd();
// gr->Draw("A*");
// c1.SaveAs("lastllh.png");
//c2.cd();
//hist->Fit("gaus");
//hist->Draw();
//c2.SaveAs("gauss_points.png");

//root fit stuff start
RooRealVar x("x","x",-10,10) ;
//
RooRealVar mean("mean","mean of gaussian",first_mu,lbmu,ubmu) ;
RooRealVar sigma("sigma","width of gaussian",first_simga,lbsigma,ubsigma) ;
RooGaussian gauss("gauss","gaussian PDF",x,mean,sigma) ;
RooDataSet data("data", "data",RooArgSet(x));

 for (unsigned int i = 0; i < gp; i++) {
   x = gauss_points[i];
   data.add(RooArgSet(x));
}

// Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = gauss.createNLL(data) ;

// Create MINUIT interface object
  RooMinuitMCMC m(*nll) ;

nll->Print();
m.mcmc(1000,100);


TGraph profile = m.getProfile("sigma",kTRUE);
TGraph stepProfil = m.getStepProfile("mean",kTRUE);
//data.Print("v");
//cout << endl ;
//gauss.fitTo(data);
//mean.Print() ;
//std::cout <<"mu(MCMC) = "<< bestmu << "+-" << (mue[0] - mue[1]) << std::endl;

//data->plotOn(xframe) ;
//model.plotOn(xframe) ;
//nll->plotOn(xframe);
//xframe->SetAxisRange(-10,10);

//root fit stuff end
c1.cd();
profile.Draw();
c2.cd();
stepProfil.Draw();

// TVectorD testrange(5400);
// TVectorD testval(5400);
//
// int n = 0;
// for (int i = 800; i < 5200; i++) {
//   testrange[n] = i/1000.0;
//   mean.setVal(i/1000.0);
//   sigma.setVal(real_sigma);
//   RooArgList testset(mean,sigma);
//   testval[n] = nll->getVal(testset);
// //  std::cout << testval[n] << std::endl;
//   n++;
//   // testset.at(0)->Print();
//   // testset.at(1)->Print();
// }
// TGraph *testg = new TGraph(testrange,testval);
// c1.cd();
// testg->Draw("A*");


//turns off the program with mous clic
theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
//starts the canvas
theApp.Run();

  return 1;
}
