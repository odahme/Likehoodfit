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

#ifndef ROOS4COSSWAVE
#define ROOS4COSSWAVE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS4CosSwave : public RooAbsPdf {
public:
  
  RooS4CosSwave();
  
  RooS4CosSwave(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S4,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A4s,
	   RooAbsReal& _Fs
	   );
  RooS4CosSwave(const RooS4CosSwave& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS4CosSwave(*this,newname); }
  inline virtual ~RooS4CosSwave() { }


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
  
  Double_t evaluate() const ;

private:

  ClassDef(RooS4CosSwave,1) // Your description goes here...
};
 
#endif
