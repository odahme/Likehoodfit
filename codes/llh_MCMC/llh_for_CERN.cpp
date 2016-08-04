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


// double f_gllh(double mu, double sigma, TVectorD gauss_points)
// {
//   double calc = gauss_points.GetNrows() * log(sqrt(2*TMath::Pi()*pow(sigma,2)));
//   for (unsigned int i = 0; i < gauss_points.GetNrows(); i++) {
//     calc += pow(gauss_points[i]-mu,2)/(2*pow(sigma,2));
//   }
//   return -calc;
// }
double f_gllh(double mu, double sigma, TVectorD gauss_points)
{
  double calc = 0;
  for (unsigned int i = 0; i < gauss_points.GetNrows(); i++) {
    calc += log(TMath::Gaus(gauss_points[i],mu,sigma,true));
  }

  return -calc;
}


int main(int argc, char **argv) {
unsigned long jobID = 0;
TRandom3 *rnd = new TRandom3(35);

//create a window
TApplication theApp("demoApplication",&argc,argv);
//create a canvas

TCanvas c1("c1","c1",1,1,1920,1080);
//TCanvas c2("c2","c2",1,1,1024,768);
//TCanvas c3("c3","c3",1,1,1024,768);

TH1F * hist = new TH1F("hist", "hist", 100, -2, 8);


//simulate gauss distributed points

unsigned long gp = 1000; //number of points
TVectorD gauss_points(gp);
double real_sigma = 1.5;
double real_mu = 3;
for (unsigned int j = 0; j < gauss_points.GetNrows(); j++) {
  gauss_points[j] = rnd->Gaus(real_mu,real_sigma);
  hist->Fill(gauss_points[j]);
}

//initialise

unsigned int nparams = 2; //number of parameters
unsigned int nstat = 1000;//number of tries
double maxstep = 0.01; //maximum step size
double ubmu = 5; //upper bound for mu
double lbmu = 1; //lower bound for mu
double ubsigma = 3; //upper bound for sigma
double lbsigma = 0.5; //lower bound for sigma
double alphastar = 0.234; //forced acceptance rate

ofstream rejectstream;
ofstream acceptstream;

Bool_t accepted;
Bool_t inchain;
unsigned int ntested = 0;
unsigned int naccepted = 0;

TVectorD last(nparams);
TVectorD curr(nparams);

unsigned int nplotp = 50; //number of points in plot
double llh_p[nstat];
double mu_p[nstat];
double minllh = 1e32;

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

double first_mu = rnd->Uniform(1,5);
double first_simga = rnd->Uniform(1,2);

last[0] = first_mu;
last[1] = first_simga;

//MCMC


rejectstream.open("list_of_rejected_candidates.txt");
acceptstream.open("List_of_accepted_candidates.txt");
for (unsigned int i = 0; i < nstat; i++) {
  curr = last; //use value of last for current then vary

  TVectorD WN(nparams);
  for (unsigned int j = 0; j < WN.GetNrows() ; j++) {
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
}*/

//evolution of the llh

  if (i > nplotp) {
    double llhplotp[nplotp];
    double muplotp[nplotp];
    for (unsigned int inplotp = 0; inplotp < nplotp; inplotp++) {
      llhplotp[inplotp] = llh_p[i-inplotp];
      muplotp[inplotp] = mu_p[i-inplotp];
    }


    TGraph * gr = new TGraph(nplotp,muplotp,llhplotp);
    c1.cd();
    gr->Draw("A*");
    c1.Update();
    //usleep(500);
  }



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
    inchain = false;
    //first save rejected candidate to file
    for (unsigned int line = 0; line < curr.GetNrows(); line++) {
      rejectstream << curr[line] << "\n";
      //std::cout << curr[line] << std::endl;
    }
    //then reset to last candidate
    curr = last;
  }

  inchain = true;

  for (unsigned int line = 0; line < curr.GetNrows(); line++) {
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

  if (accepted){
    std::cout << "set " << i << " accepted, rate " << std::fixed << std::setprecision(3) << acceptrate << " last " << nlast << " " << lastratio << std::endl;
  }
  else{
    std::cout << "set " << i << " rejected, rate " << std::fixed << std::setprecision(3) << acceptrate << " last " << nlast << " " << lastratio << std::endl;
  }

  ntested++;

  //update S matrix
  TMatrixDSym SNminusoneT(SNminusone);
  SNminusoneT.T();
  double etan = std::min(1.0, nparams*pow(double(i), -2.0/3.0));
  TMatrixDSym WNWNT(nparams);
  WNWNT.Zero();
  for (unsigned int row = 0; row < WNWNT.GetNrows(); row++) {
    for (unsigned int col = 0; col < WNWNT.GetNcols(); col++) {
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
  for (unsigned int row = 0; row < SN.GetNrows(); row++) {
    for (unsigned int col = 0; col < SN.GetNcols(); col++) {
      SNminusone[row][col] = SN[row][col];
    }
  }
}

rejectstream.close();
acceptstream.close();
std::cout << ntested << " points were tested" << std::endl;
std::cout << naccepted << " points were accepted" << std::endl;
std::cout << "The accept fraction is " << double(naccepted)/double(ntested) << std::endl;

double mue[2];
unsigned int imu = 0;
for (unsigned int i = 0; i < nstat; i++) {
  llh_p[i] -= minllh;
  if (llh_p[i] < 0.51 && llh_p[i] > 0.49) {
    std::cout << "mu = "<< mu_p[i] << std::endl;
    mue[imu] = mu_p[i];
    imu++;
  }
}

std::cout << "error at 68p = " << abs(mue[0] - mue[1]) << std::endl;


TGraph * gr = new TGraph(nstat,mu_p,llh_p);
c1.cd();
//gr->Draw("A*");
c1.SaveAs("lastllh.png");
//c2.cd();
//hist->Fit("gaus");
//hist->Draw();
//c2.SaveAs("gauss_points.png");

//turns off the program with mous clic
theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
//starts the canvas
theApp.Run();

  return 1;
}
