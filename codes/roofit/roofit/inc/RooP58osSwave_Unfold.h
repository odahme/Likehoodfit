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

#ifndef ROOP8COSSWAVE
#define ROOP8COSSWAVE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooP8CosSwave : public RooAbsPdf {
public:
  
  RooP8CosSwave();
  
  RooP8CosSwave(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _P8,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A7s,
	   RooAbsReal& _Fs
	   );
  RooP8CosSwave(const RooP8CosSwave& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooP8CosSwave(*this,newname); }
  inline virtual ~RooP8CosSwave() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy P8 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A7s ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooP8CosSwave,1) // Your description goes here...
};
 
#endif
