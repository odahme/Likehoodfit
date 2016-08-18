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

#ifndef ROOS4COSSWAVEU
#define ROOS4COSSWAVEU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS4CosSwave_Unfold : public RooAbsPdf {
public:
  
  RooS4CosSwave_Unfold();
  
  RooS4CosSwave_Unfold(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S4,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A4s,
	   RooAbsReal& _Fs,
	   RooAbsReal& _q2
	   );
  RooS4CosSwave_Unfold(const RooS4CosSwave_Unfold& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS4CosSwave_Unfold(*this,newname); }
  inline virtual ~RooS4CosSwave_Unfold() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy S4 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A4s ;
  RooRealProxy Fs ;  
  RooRealProxy q2 ;

  
  Double_t evaluate() const ;

private:

  ClassDef(RooS4CosSwave_Unfold,1) // Your description goes here...
};
 
#endif
