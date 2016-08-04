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

  return calc;
}


int main(int argc, char **argv) {
unsigned long jobID = 0;
TRandom3 *rnd = new TRandom3(35);

//create a window
TApplication theApp("demoApplication",&argc,argv);
//create a canvas

TCanvas c1("c1","c1",1,1,1024,768);
//TCanvas c2("c2","c2",1,1,1024,768);
//TCanvas c3("c3","c3",1,1,1024,768);



//simulate gauss distributed points

unsigned long gp = 100; //number of points
TVectorD gauss_points(gp);
double real_sigma = 1.5;
double real_mu = 3;
for (unsigned int j = 0; j < gauss_points.GetNrows(); j++) {
  gauss_points[j] = rnd->Gaus(real_mu,real_sigma);
}

//initialise

unsigned int nparams = 2; //number of parameters
unsigned int nstat = 30;//number of tries

ofstream mystream;

Bool_t accepted;
Bool_t inchain;
unsigned int ntested = 0;
unsigned int naccepted = 0;

TVectorD last(nparams);
TVectorD curr(nparams);

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

double alphastar = 0.234; //forced acceptance rate

double llh_p[nstat];
double mu_p[nstat];

for (unsigned int i = 0; i < nstat; i++) {
  curr = last; //use value of last for current then vary

  TVectorD WN(nparams);
  for (unsigned int j = 0; j < WN.GetNrows() ; j++) {
    WN[j] = rnd->Gaus(0.0, 1.0);
  }
  TVectorD SW(SNminusone*WN);
  curr += SW;

std::cout << "mu = "<< curr[0] << std::endl;
std::cout << "sigma = "<< curr[1] << std::endl;

  //get llh for current and last
  double llh_last = f_gllh(last[0],last[1],gauss_points);
  double llh_curr = f_gllh(curr[0],curr[1],gauss_points);

  std::cout << "llh = " << llh_curr << std::endl;
  llh_p[i] = llh_curr;
  mu_p[i] = curr[0];

  double alpha = std::min(1.0, exp(abs(llh_curr - llh_last)));
  std::cout << "alpha = " << alpha << std::endl;
  double r = rnd->Uniform(0,1);
  accepted = false;
  double acceptrate = double(naccepted)/double(ntested);

  if (r < alpha) {
    accepted = true;
    naccepted++;
    //success, candidate filled below
    last = curr;
  } else {
    accepted = false;
    inchain = false;
    //first save rejected candidate to file
    mystream.open("list_of_candidates.txt");
    for (unsigned int line = 0; line < curr.GetNrows(); line++) {
      mystream << curr[i] << "\n";
    }
    mystream.close();
    //then reset to last candidate
    curr = last;
  }


  inchain = true;





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

std::cout << ntested << " points were tested" << std::endl;
std::cout << naccepted << " points were accepted" << std::endl;
std::cout << "The accept fraction is " << double(naccepted)/double(ntested) << std::endl;


TGraph * gr = new TGraph(nstat,mu_p,llh_p);
c1.cd();
gr->Draw("A*");

//turns off the program with mous clic
theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
//starts the canvas
theApp.Run();

  return 1;
}
