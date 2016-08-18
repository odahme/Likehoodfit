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
#include "RooMoriondPdf.h"
 #include "RooAbsReal.h" 
#include "TMath.h"


 ClassImp(RooMoriondPdf) 

 
 RooMoriondPdf::RooMoriondPdf()
 {

   RooRealVar cosThetaL("cosThetaL","cos thetaL", -1.0, 1.0);
   RooRealVar cosThetaK("cosThetaK","cos thetaK", -1.0, 1.0);
   RooRealVar phi("phi","phi", 0.0, TMath::Pi());
   RooRealVar FL("FL","FL", 0.0, 1.0);
   RooRealVar AT2("AT2","AT2", 0.0, 1.0);
   RooRealVar ATRe("ATRe","ATRe", -1.0, 1.0);
   RooRealVar ATIm("ATIm","ATIm", -1.0, 1.0);
   RooRealVar As("As","As", -1.0, 1.0);
   RooRealVar Fs("Fs","Fs", 0.0, 1.0);

   RooMoriondPdf("RooMoriondPdf", "RooMoriondPdf", cosThetaL, cosThetaK, phi, FL, ATRe, ATIm, AT2, As, Fs);

  
   
}



RooMoriondPdf::RooMoriondPdf(const char *name, const char *title, 
	     RooAbsReal& _cosThetaL,
	     RooAbsReal& _cosThetaK,
	     RooAbsReal& _phi,
	     RooAbsReal& _FL,
	     RooAbsReal& _ATRe,
	     RooAbsReal& _ATIm,
	     RooAbsReal& _AT2,
	     RooAbsReal& _As,
	     RooAbsReal& _Fs
	     ) :
  RooAbsPdf(name,title), 
  cosThetaL("cosThetaL","cos thetaL", this, _cosThetaL),
  cosThetaK("cosThetaK","cos thetaK", this, _cosThetaK),
  phi("phi","phi", this, _phi),
  FL("FL","FL", this, _FL),
  ATRe("ATRe","ATRe",this,_ATRe),
  ATIm("ATIm","ATIm",this,_ATIm),
  AT2("AT2","AT2",this,_AT2),
  As("As","As",this,_As),
  Fs("Fs","Fs",this,_Fs)

  
{ 

} 

     
RooMoriondPdf::RooMoriondPdf(const RooMoriondPdf& other, const char* name) :  
  RooAbsPdf(other,name),
  cosThetaL("cosThetaL", this, other.cosThetaL),
  cosThetaK("cosThetaK", this, other.cosThetaK),
  phi("phi", this,other.phi),
  FL("FL", this,other.FL),
  ATRe("ATRe",this,other.ATRe),
  ATIm("ATIm",this,other.ATIm),
  AT2("AT2",this,other.AT2),
  As("As",this,other.As),
  Fs("Fs",this,other.Fs){ 
} 
     


Double_t RooMoriondPdf::evaluate() const 
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
  //// Adding the Swave 
  value*=(1-Fs);
  value+=(2./3.)*Fs*sinThetaL*sinThetaL;
  value+=As*sinThetaL*sinThetaL*cosThetaK;
    
  return value ;
  
  
} 





