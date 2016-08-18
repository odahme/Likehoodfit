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

#ifndef ROOP6COSSWAVE
#define ROOP6COSSWAVE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooP6CosSwave : public RooAbsPdf {
public:
  
  RooP6CosSwave();
  
  RooP6CosSwave(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _P6,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A6s,
	   RooAbsReal& _Fs
	   );
  RooP6CosSwave(const RooP6CosSwave& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooP6CosSwave(*this,newname); }
  inline virtual ~RooP6CosSwave() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy P6 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A6s ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooP6CosSwave,1) // Your description goes here...
};
 
#endif
