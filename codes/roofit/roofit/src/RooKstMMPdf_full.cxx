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
#include "RooKstMMPdf_full.h"
 #include "RooAbsReal.h" 
#include "TMath.h"


 ClassImp(RooKstMMPdf_full) 

 
 RooKstMMPdf_full::RooKstMMPdf_full()
 {

   RooRealVar cosThetaL("cosThetaL","cos thetaL", -1.0, 1.0);
   RooRealVar cosThetaK("cosThetaK","cos thetaK", -1.0, 1.0);
   RooRealVar phi("phi","phi", -1.*TMath::Pi(), TMath::Pi());
   RooRealVar S1s("S1s","S1s",-1.0, 1.0);
   RooRealVar S1c("S1c","S1c",-1.0, 1.0);
   RooRealVar S2s("S2s","S2s",-1.0, 1.0);
   RooRealVar S2c("S2c","S2c",-1.0, 1.0);

   RooRealVar S3("S3","S3",-1.0, 1.0);
   RooRealVar S4("S4","S4",-1.0, 1.0);
   RooRealVar S5("S5","S5",-1.0, 1.0);
   RooRealVar S6("S6","S6",-1.0, 1.0);
   RooRealVar S7("S7","S7",-1.0, 1.0);  
   RooRealVar S8("S8","S8",-1.0, 1.0);
   RooRealVar S9("S9","S9",-1.0, 1.0); 
   RooRealVar As("As","As", -1.0, 1.0);
   RooRealVar A4s("A4s","A4s", -1.0, 1.0);
   RooRealVar A5s("A5s","A5s", -1.0, 1.0);
   RooRealVar A7s("A7s","A7s", -1.0, 1.0);
   RooRealVar A8s("A8s","A8s", -1.0, 1.0);
   RooRealVar Fs("Fs","Fs", 0.0, 1.0);

     RooKstMMPdf_full("RooKstMMPdf_full", "RooKstMMPdf_full", cosThetaL, cosThetaK, phi, S1s,S1c,S2s,S2c, S3, S4, S5, S6, S7, S8, S9, As, A4s, A5s, A7s, A8s, Fs);
   
}



RooKstMMPdf_full::RooKstMMPdf_full(const char *name, const char *title, 
			       RooAbsReal& _cosThetaL,
			       RooAbsReal& _cosThetaK,
			       RooAbsReal& _phi,
			       RooAbsReal& _S1s,
			       RooAbsReal& _S1c,
			       RooAbsReal& _S2s,
			       RooAbsReal& _S2c,
			       RooAbsReal& _S3,
			       RooAbsReal& _S4,
			       RooAbsReal& _S5,
			       RooAbsReal& _S6,
			       RooAbsReal& _S7,
			       RooAbsReal& _S8,
			       RooAbsReal& _S9,
			       RooAbsReal& _As,
			       RooAbsReal& _A4s,
			       RooAbsReal& _A5s,
			       RooAbsReal& _A7s,
			       RooAbsReal& _A8s,
			       RooAbsReal& _Fs
	     ) :
  RooAbsPdf(name,title), 
  cosThetaL("cosThetaL","cos thetaL", this, _cosThetaL),
  cosThetaK("cosThetaK","cos thetaK", this, _cosThetaK),
  phi("phi","phi", this, _phi),
  S1s("S1s","S1s",this,_S1s),
  S1c("S1c","S1c",this,_S1c),
  S2s("S2s","S2s",this,_S2s),
  S2c("S2c","S2c",this,_S2c),
  S3("S3","S3",this,_S3),
  S4("S4","S4",this,_S4),
  S5("S5","S5",this,_S5),
  S6("S6","S6",this,_S6),  
  S7("S7","S7",this,_S7),
  S8("S8","S8",this,_S8),
  S9("S9","S9",this,_S9),
  As("As","As",this,_As),
  A4s("A4s","A4s",this,_A4s),
  A5s("A5s","A5s",this,_A5s),
  A7s("A7s","A7s",this,_A7s),
  A8s("A8s","A8s",this,_A8s),
  Fs("Fs","Fs",this,_Fs)  
{ 

} 

     
RooKstMMPdf_full::RooKstMMPdf_full(const RooKstMMPdf_full& other, const char* name) :  
  RooAbsPdf(other,name),
  cosThetaL("cosThetaL", this, other.cosThetaL),
  cosThetaK("cosThetaK", this, other.cosThetaK),
  phi("phi", this,other.phi),
  S1s("S1s",this,other.S1s),
  S1c("S1c",this,other.S1c),
  S2s("S2s",this,other.S2s),
  S2c("S2c",this,other.S2c),
  S3("S3",this,other.S3),
  S4("S4",this,other.S4),
  S5("S5",this,other.S5),
  S6("S6",this,other.S6),
  S7("S7",this,other.S7),
  S8("S8",this,other.S8),
  S9("S9",this,other.S9),
  As("As",this,other.As),
  A4s("A4s",this,other.A4s),
  A5s("A5s",this,other.A5s),
  A7s("A7s",this,other.A7s),
  A8s("A8s",this,other.A8s),
  Fs("Fs",this,other.Fs){ 
} 
     

Double_t RooKstMMPdf_full::evaluate() const 
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
  

  double value  = S1s* sinTheta2K;
  value+= S1c  * cosThetaK * cosThetaK;
  value+= S2s * sinTheta2K * cos2ThetaL;
  value+= S2c * cosThetaK*cosThetaK * cos2ThetaL;
  value+= S3*sinTheta2K * sinTheta2L * TMath::Cos( 2*phi );  
  value+= S4* sin2ThetaK * sin2ThetaL * TMath::Cos( phi );
  value+= S5* sin2ThetaK * sinThetaL * TMath::Cos( phi );
  value+= S6*sinTheta2K*cosThetaL;
  value+= S7* sin2ThetaK * sinThetaL * TMath::Sin( phi );
  value+= S8* sin2ThetaK * sin2ThetaL * TMath::Sin( phi );
  value+= S9*sinTheta2K*sinTheta2L*TMath::Sin( 2*phi );
  //// Adding the Swave 
  value*=(1-Fs);
  value+=(2./3.)*Fs*sinThetaL*sinThetaL;
  value+=(4./3.)*As*sinThetaL*sinThetaL*cosThetaK;
  value+=A4s*sinThetaK*sin2ThetaL*TMath::Cos( phi );
  value+=A5s*sinThetaK*sinThetaL*TMath::Cos( phi );
  value+=A7s*sinThetaK*sinThetaL*TMath::Sin( phi );
  value+=A8s*sinThetaK*sin2ThetaL*TMath::Sin( phi );

  
  return value ;
  
  
} 





