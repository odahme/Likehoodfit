 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Copyright (c) 2000-2005, Regents of the University of California          * 
  *                          and Stanford University. All rights reserved.    * 
  *                                                                           * 
  * Redistribution and use in source and binary forms,                        * 
  * with or without modification, are permitted according to the terms        * 
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
  *****************************************************************************/ 

 // -- CLASS DESCRIPTION [PDF] -- 
 // Your description goes here... 

 #include <iostream> 

// class declaration include file below retrieved from workspace code storage
#include "RooKstMMPdf_Unfold.h"
#include "RooAbsReal.h" 
#include "TMath.h"

#include "RooKstMMPdf_My_Unfold.h"



 ClassImp(RooKstMMPdf_Unfold) 

 
 RooKstMMPdf_Unfold::RooKstMMPdf_Unfold()
 {

   RooRealVar cosThetaL("cosThetaL","cos thetaL", -1.0, 1.0);
   RooRealVar cosThetaK("cosThetaK","cos thetaK", -1.0, 1.0);
   RooRealVar phi("phi","phi", -1.*TMath::Pi(), TMath::Pi());
   RooRealVar FL("FL","FL", 0.0, 1.0);
   RooRealVar P4("P4","P4",-1.0, 1.0);
   RooRealVar P5("P5","P5",-1.0, 1.0);
   RooRealVar P6("P6","P6",-1.0, 1.0);
   RooRealVar P8("P8","P8",-1.0, 1.0);
   RooRealVar AT2("AT2","AT2", -1.0, 1.0);
   RooRealVar ATRe("ATRe","ATRe", -1.0, 1.0);
   RooRealVar ATIm("ATIm","ATIm", -1.0, 1.0);
   RooRealVar As("As","As", -1.0, 1.0);
   RooRealVar A4s("A4s","A4s", -1.0, 1.0);
   RooRealVar A5s("A5s","A5s", -1.0, 1.0);
   RooRealVar A7s("A7s","A7s", -1.0, 1.0);
   RooRealVar A8s("A8s","A8s", -1.0, 1.0);
   RooRealVar Fs("Fs","Fs", 0.0, 1.0);
   RooRealVar q2("q2","q2", 0.0, 1000.0);

     RooKstMMPdf_Unfold("RooKstMMPdf_Unfold", "RooKstMMPdf_Unfold", cosThetaL, cosThetaK, phi, FL, P4, P5, P6, P8, AT2, ATRe, ATIm, As, A4s, A5s, A7s, A8s, Fs, q2);
   
}



RooKstMMPdf_Unfold::RooKstMMPdf_Unfold(const char *name, const char *title, 
	     RooAbsReal& _cosThetaL,
	     RooAbsReal& _cosThetaK,
	     RooAbsReal& _phi,
	     RooAbsReal& _FL,
	     RooAbsReal& _P4,
	     RooAbsReal& _P5,
	     RooAbsReal& _P6,
	     RooAbsReal& _P8,
	     RooAbsReal& _AT2,
	     RooAbsReal& _ATRe,
	     RooAbsReal& _ATIm,
	     RooAbsReal& _As,
	     RooAbsReal& _A4s,
	     RooAbsReal& _A5s,
	     RooAbsReal& _A7s,
	     RooAbsReal& _A8s,
             RooAbsReal& _Fs,
	     RooAbsReal& _q2			       
	     ) :
  RooAbsPdf(name,title), 
  cosThetaL("cosThetaL","cos thetaL", this, _cosThetaL),
  cosThetaK("cosThetaK","cos thetaK", this, _cosThetaK),
  phi("phi","phi", this, _phi),
  FL("FL","FL", this, _FL),
  P4("P4","P4",this,_P4),
  P5("P5","P5",this,_P5),
  P6("P6","P6",this,_P6),  
  P8("P8","P8",this,_P8),  
  AT2("AT2","AT2",this,_AT2),
  ATRe("ATRe","ATRe",this,_ATRe),
  ATIm("ATIm","ATIm",this,_ATIm),
  As("As","As",this,_As),
  A4s("A4s","A4s",this,_A4s),
  A5s("A5s","A5s",this,_A5s),
  A7s("A7s","A7s",this,_A7s),
  A8s("A8s","A8s",this,_A8s),
  Fs("Fs","Fs",this,_Fs),  
  q2("q2","q2",this,_q2){ 

} 

     
RooKstMMPdf_Unfold::RooKstMMPdf_Unfold(const RooKstMMPdf_Unfold& other, const char* name) :  
  RooAbsPdf(other,name),
  cosThetaL("cosThetaL", this, other.cosThetaL),
  cosThetaK("cosThetaK", this, other.cosThetaK),
  phi("phi", this,other.phi),
  FL("FL", this,other.FL),
  P4("P4",this,other.P4),
  P5("P5",this,other.P5),
  P6("P6",this,other.P6),
  P8("P8",this,other.P8),
  AT2("AT2",this,other.AT2),
  ATRe("ATRe",this,other.ATRe),
  ATIm("ATIm",this,other.ATIm),
  As("As",this,other.As),
  A4s("A4s",this,other.A4s),
  A5s("A5s",this,other.A5s),
  A7s("A7s",this,other.A7s),
  A8s("A8s",this,other.A8s),
  Fs("Fs",this,other.Fs),
  q2("q2",this,other.q2){ 
} 
     

Double_t RooKstMMPdf_Unfold::evaluate() const 
{ 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   
  double sinThetaK = sqrt( 1 - cosThetaK*cosThetaK );
  double sinThetaL = sqrt( 1 - cosThetaL*cosThetaL );

  double sinTheta2K =  1 - cosThetaK*cosThetaK; 
  double sinTheta2L =  1 - cosThetaL*cosThetaL; 

  double sin2ThetaK = 2*sinThetaK*cosThetaK;
  //double cos2ThetaK = 2*cosThetaK*cosThetaK - 1;

  double sin2ThetaL = 2*sinThetaL*cosThetaL;
  double cos2ThetaL = 2*cosThetaL*cosThetaL - 1;
  

  double value  = 3.0/4*( 1 - FL )* sinTheta2K;
  value+= FL * cosThetaK * cosThetaK;
  value+= 1.0/4*(1 - FL)* sinTheta2K * cos2ThetaL;
  value+= (-1)*FL * cosThetaK*cosThetaK * cos2ThetaL;
  value+= 1/2.0*(1 - FL)* AT2 * sinTheta2K * sinTheta2L * TMath::Cos( 2*phi ); 
  value+= (1- FL)*ATRe*sinTheta2K*cosThetaL;
  value+=(1./2.)*(1- FL)*ATIm*sinTheta2K*sinTheta2L*TMath::Sin( 2*phi );
  value+= P4*sqrt(FL*(1-FL))* sin2ThetaK * sin2ThetaL * TMath::Cos( phi );
  value+= P5*sqrt(FL*(1-FL))* sin2ThetaK * sinThetaL * TMath::Cos( phi );
  value+= P6*sqrt(FL*(1-FL))* sin2ThetaK * sinThetaL * TMath::Sin( phi );
  value+= P8*sqrt(FL*(1-FL))* sin2ThetaK * sin2ThetaL * TMath::Sin( phi );
  
    
  //// Adding the Swave 
  value*=(1-Fs);
  value+=(2./3.)*Fs*sinThetaL*sinThetaL;
  value+=(4./3.)*As*sinThetaL*sinThetaL*cosThetaK;
  value+=A4s*sinThetaK*sin2ThetaL*TMath::Cos( phi );
  value+=A5s*sinThetaK*sinThetaL*TMath::Cos( phi );
  value+=A7s*sinThetaK*sinThetaL*TMath::Sin( phi );
  value+=A8s*sinThetaK*sin2ThetaL*TMath::Sin( phi );

   value*=RooKstMMPdf_My_Unfold::eff(q2, cosThetaL, cosThetaK,phi);

  return value ;
  
  
} 





