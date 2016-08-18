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

#ifndef ROOKSTMMPDFMYU
#define ROOKSTMMPDFMYU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooKstMMPdf_My_Unfold : public RooAbsPdf {
public:
  
  RooKstMMPdf_My_Unfold();
  
  RooKstMMPdf_My_Unfold(const char *name, const char *title,
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
		 RooAbsReal& _Fs,
		 RooAbsReal& _q2);

  RooKstMMPdf_My_Unfold(const RooKstMMPdf_My_Unfold& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooKstMMPdf_My_Unfold(*this,newname); }
  inline virtual ~RooKstMMPdf_My_Unfold() { }
  static  double eff(double, double, double, double);
  Double_t evaluate() const ;


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
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
  RooRealProxy q2 ;

  
  //Double_t evaluate() const ;

private:

  ClassDef(RooKstMMPdf_My_Unfold,1) // Your description goes here...
};
 
#endif
