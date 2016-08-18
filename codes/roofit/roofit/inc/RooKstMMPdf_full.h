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

#ifndef ROOKSTMMPDFFULL
#define ROOKSTMMPDFFULL

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooKstMMPdf_full : public RooAbsPdf {
public:
  
  RooKstMMPdf_full();
  
  RooKstMMPdf_full(const char *name, const char *title,
		 RooAbsReal& _cosThetaL,
		 RooAbsReal& _cosThetaK,
		 RooAbsReal& _phi,
		 RooAbsReal& _S1s,
		 RooAbsReal& _S1c,
		 RooAbsReal& _S2s,
		 RooAbsReal& _S2c,
		 RooAbsReal& _S3,
		 RooAbsReal& _S4,
		 RooAbsReal& _S5,
		 RooAbsReal& _S6,
		 RooAbsReal& _S7,
		 RooAbsReal& _S8,
		 RooAbsReal& _S9,
		 RooAbsReal& _As,
		 RooAbsReal& _A4s,
		 RooAbsReal& _A5s,
		 RooAbsReal& _A7s,
		 RooAbsReal& _A8s,
		 RooAbsReal& _Fs);

  RooKstMMPdf_full(const RooKstMMPdf_full& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooKstMMPdf_full(*this,newname); }
  inline virtual ~RooKstMMPdf_full() { }

  Double_t evaluate() const ;


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy S1s ;
  RooRealProxy S1c ;
  RooRealProxy S2s ;
  RooRealProxy S2c ;
  RooRealProxy S3 ;
  RooRealProxy S4 ;
  RooRealProxy S5 ;
  RooRealProxy S6 ;
  RooRealProxy S7 ;
  RooRealProxy S8 ;
  RooRealProxy S9 ;

  RooRealProxy As ;
  RooRealProxy A4s ;
  RooRealProxy A5s ;
  RooRealProxy A7s ;
  RooRealProxy A8s ;
  RooRealProxy Fs ;  
  
  //Double_t evaluate() const ;

private:

  ClassDef(RooKstMMPdf_full,1) // Your description goes here...
};
 
#endif
