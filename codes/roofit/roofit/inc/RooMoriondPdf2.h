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

#ifndef ROOMORIONDPDF2
#define ROOMORIONDPDF2

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooMoriondPdf2 : public RooAbsPdf {
public:
  
  RooMoriondPdf2();
  
  RooMoriondPdf2(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _AFB,
	   RooAbsReal& _ATIm,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _Fs
	   );
  RooMoriondPdf2(const RooMoriondPdf2& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMoriondPdf2(*this,newname); }
  inline virtual ~RooMoriondPdf2() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy AFB ;
  RooRealProxy ATIm ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooMoriondPdf2,1) // Your description goes here...
};
 
#endif
