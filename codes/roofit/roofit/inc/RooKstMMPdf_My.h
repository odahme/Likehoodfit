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

#ifndef ROOKSTMMPDFMY
#define ROOKSTMMPDFMY

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooKstMMPdf_My : public RooAbsPdf {
public:
  
  RooKstMMPdf_My();
  
  RooKstMMPdf_My(const char *name, const char *title,
		 RooAbsReal& _cosThetaL,
		 RooAbsReal& _cosThetaK,
		 RooAbsReal& _phi,
		 RooAbsReal& _FL,
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

  RooKstMMPdf_My(const RooKstMMPdf_My& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooKstMMPdf_My(*this,newname); }
  inline virtual ~RooKstMMPdf_My() { }

  Double_t evaluate() const ;


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL;
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

  ClassDef(RooKstMMPdf_My,1) // Your description goes here...
};
 
#endif
