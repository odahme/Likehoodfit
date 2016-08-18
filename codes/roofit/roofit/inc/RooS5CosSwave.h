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

#ifndef ROOS5COSSWAVE
#define ROOS5COSSWAVE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS5CosSwave : public RooAbsPdf {
public:
  
  RooS5CosSwave();
  
  RooS5CosSwave(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S5,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A5s,
	   RooAbsReal& _Fs
	   );
  RooS5CosSwave(const RooS5CosSwave& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS5CosSwave(*this,newname); }
  inline virtual ~RooS5CosSwave() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy S5 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A5s ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooS5CosSwave,1) // Your description goes here...
};
 
#endif
