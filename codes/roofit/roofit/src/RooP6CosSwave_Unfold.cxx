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
#include "RooP6CosSwave_Unfold.h"
 #include "RooAbsReal.h" 
#include "TMath.h"

#include "RooKstMMPdf_My_Unfold.h"


 ClassImp(RooP6CosSwave_Unfold) 

 
 RooP6CosSwave_Unfold::RooP6CosSwave_Unfold()
 {

   RooRealVar cosThetaL("cosThetaL","cos thetaL", -1.0, 1.0);
   RooRealVar cosThetaK("cosThetaK","cos thetaK", -1.0, 1.0);
   RooRealVar phi("phi","phi", 0.0, TMath::Pi());
   RooRealVar FL("FL","FL", 0.0, 1.0);
   RooRealVar P6("P6","P6",0.0, 1.0);
   RooRealVar AT2("AT2","AT2", 0.0, 1.0);
   RooRealVar As("As","As", -1.0, 1.0);
   RooRealVar A6s("A6s","A6s", -1.0, 1.0);
   RooRealVar Fs("Fs","Fs", 0.0, 1.0);
   RooRealVar q2("q2","q", 0.0, 1000.0);

     RooP6CosSwave_Unfold("RooP6CosSwave_Unfold", "RooP6CosSwave_Unfold", cosThetaL, cosThetaK, phi, FL, P6, AT2, As, A6s, Fs, q2);

  
   
}



RooP6CosSwave_Unfold::RooP6CosSwave_Unfold(const char *name, const char *title, 
	     RooAbsReal& _cosThetaL,
	     RooAbsReal& _cosThetaK,
	     RooAbsReal& _phi,
	     RooAbsReal& _FL,
	     RooAbsReal& _P6,
	     RooAbsReal& _AT2,
	     RooAbsReal& _As,
	     RooAbsReal& _A6s,
	     RooAbsReal& _Fs,
	     RooAbsReal& _q2  ) :

  RooAbsPdf(name,title), 
  cosThetaL("cosThetaL","cos thetaL", this, _cosThetaL),
  cosThetaK("cosThetaK","cos thetaK", this, _cosThetaK),
  phi("phi","phi", this, _phi),
  FL("FL","FL", this, _FL),
  P6("P6","P6",this,_P6),
  AT2("AT2","AT2",this,_AT2),
  As("As","As",this,_As),
  A6s("A6s","A6s",this,_A6s),
  Fs("Fs","Fs",this,_Fs),
  q2("q2","q2",this,_q2)

  
{ 

} 

     
RooP6CosSwave_Unfold::RooP6CosSwave_Unfold(const RooP6CosSwave_Unfold& other, const char* name) :  
  RooAbsPdf(other,name),
  cosThetaL("cosThetaL", this, other.cosThetaL),
  cosThetaK("cosThetaK", this, other.cosThetaK),
  phi("phi", this,other.phi),
  FL("FL", this,other.FL),
  P6("P6",this,other.P6),
  AT2("AT2",this,other.AT2),
  As("As",this,other.As),
  A6s("A6s",this,other.A6s),
  Fs("Fs",this,other.Fs),
  q2("q2",this,other.q2){ 
} 
     


Double_t RooP6CosSwave_Unfold::evaluate() const 
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
  value+= P6*sqrt(FL*(1-FL))* sin2ThetaK * sinThetaL * TMath::Sin( phi );
  //// Adding the Swave 
  value*=(1-Fs);
  value+=(2./3.)*Fs*sinThetaL*sinThetaL;
  value+=(4./3.)*As*sinThetaL*sinThetaL*cosThetaK;
  value+=A6s*sinThetaK*sinThetaL*TMath::Sin( phi );

  value*=RooKstMMPdf_My_Unfold::eff(q2, cosThetaL, cosThetaK,phi);
    
  return value ;
  
  
} 





