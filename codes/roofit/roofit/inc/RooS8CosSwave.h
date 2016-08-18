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

#ifndef ROOS8COSSWAVE
#define ROOS8COSSWAVE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS8CosSwave : public RooAbsPdf {
public:
  
  RooS8CosSwave();
  
  RooS8CosSwave(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S8,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A7s,
	   RooAbsReal& _Fs
	   );
  RooS8CosSwave(const RooS8CosSwave& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS8CosSwave(*this,newname); }
  inline virtual ~RooS8CosSwave() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy S8 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A7s ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooS8CosSwave,1) // Your description goes here...
};
 
#endif
