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

#ifndef ROOMORIONDPDF
#define ROOMORIONDPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooMoriondPdf : public RooAbsPdf {
public:
  
  RooMoriondPdf();
  
  RooMoriondPdf(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _ATRe,
	   RooAbsReal& _ATIm,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _Fs
	   );
  RooMoriondPdf(const RooMoriondPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooMoriondPdf(*this,newname); }
  inline virtual ~RooMoriondPdf() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy ATRe ;
  RooRealProxy ATIm ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy Fs ;  
  
  Double_t evaluate() const ;

private:

  ClassDef(RooMoriondPdf,1) // Your description goes here...
};
 
#endif
