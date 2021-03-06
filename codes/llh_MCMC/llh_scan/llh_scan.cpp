
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


double f_gllh(double mu, double sigma, TVectorD gauss_points)
{
  double calc = 0;
  for (unsigned int i = 0; i < gauss_points.GetNrows(); i++) {
    calc += log(TMath::Gaus(gauss_points[i],mu,sigma,true));
  }

  return -calc;
}

int main(int argc, char **argv) {

TApplication theApp("demoApplication",&argc,argv);

TCanvas c1("c1","c1",1,1,1024,768);
TCanvas c2("c2","c2",1,1,1024,768);


TRandom3 *rnd = new TRandom3(35);

unsigned long gp = 1000; //number of points
TVectorD gauss_points(gp);
const double real_sigma = 1.5;
const double real_mu = 3.0;
for (unsigned int j = 0; j < gauss_points.GetNrows(); j++) {
  gauss_points[j] = rnd->Gaus(real_mu,real_sigma);
}


unsigned int count = 10000;
double mu_start = 2.9;
double mu_end = 3.1;
TVectorD llh_p(count);
TVectorD mu_p(count);

double sigma = 1.5;
double mu = mu_start;

double llhmin = 1e32;
double calc = 0;
double x = abs(mu_end - mu_start)/count;
for (unsigned int i = 0; i < count; i++) {
  calc = f_gllh(mu,sigma,gauss_points);
  mu_p[i] = mu;
  // std::cout << "mu = "<< mu << std::endl;
  llh_p[i] = calc;
  if (calc < llhmin) {
    llhmin = calc;
  }
  mu += x;
}
for (size_t i = 0; i < count; i++) {
  llh_p[i] -= llhmin;
  //std::cout << "mu"<<i<<" = "<< mu_p[i] << std::endl;
}


mu = 3;
unsigned int sigma_count = 300;
double sigma_p[sigma_count];
double llh_p2[sigma_count];

std::cout << "OBEN" << std::endl;
// for (unsigned int i = 20; i < sigma_count+20; i++) {
//   std::cout << "i = "<<i << std::endl;
//   sigma = i/100.0;
//   calc = f_gllh(mu,sigma,gauss_points);
//   sigma_p[i] = sigma;
//   // std::cout << "sigma = "<< sigma << std::endl;
//   // std::cout << "llhs = "<< calc << std::endl;
//   llh_p2[i] = calc;
// }

std::cout << "UNTEN" << std::endl;

std::cout << "llhmin = "<< llhmin << std::endl;

TGraph *gr = new TGraph(mu_p,llh_p);
TGraph *gr2 = new TGraph(sigma_count,sigma_p,llh_p2);
c1.cd();
gr->SetTitle("negtive Log-Likelihood");
gr->GetXaxis()->SetTitle("parameter of intrest");
gr->GetYaxis()->SetTitle("nll value");
gr->Draw("A*");
c1.SaveAs("lastmu.png");
c2.cd();
gr2->Draw("A*");
c2.SaveAs("lastsigma.png");

//turns off the program with mous clic
theApp.Connect("TCanvas","Closed()","TApplication",&theApp,"Terminate()");
//starts the canvas
theApp.Run();

return 1;

}
