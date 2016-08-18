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

#ifndef ROOS7COSSWAVE
#define ROOS7COSSWAVE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS7CosSwave : public RooAbsPdf {
public:
  
  RooS7CosSwave();
  
  RooS7CosSwave(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S7,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A7s,
	   RooAbsReal& _Fs
	   );
  RooS7CosSwave(const RooS7CosSwave& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS7CosSwave(*this,newname); }
  inline virtual ~RooS7CosSwave() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy S7 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A7s ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooS7CosSwave,1) // Your description goes here...
};
 
#endif
