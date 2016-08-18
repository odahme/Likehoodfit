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

#ifndef ROOKSTMMPDFU
#define ROOKSTMMPDFU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooKstMMPdf_Unfold : public RooAbsPdf {
public:
  
  RooKstMMPdf_Unfold();
  
  RooKstMMPdf_Unfold(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _P4,
	   RooAbsReal& _P5,
	   RooAbsReal& _P6,
	   RooAbsReal& _P8,
	   RooAbsReal& _AT2,
	   RooAbsReal& _ATRe,
	   RooAbsReal& _ATIm,
	   RooAbsReal& _As,
	   RooAbsReal& _A4s,
	   RooAbsReal& _A5s,
	   RooAbsReal& _A7s,
	   RooAbsReal& _A8s,
	   RooAbsReal& _Fs,
	   RooAbsReal& _q2 );

  RooKstMMPdf_Unfold(const RooKstMMPdf_Unfold& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooKstMMPdf_Unfold(*this,newname); }
  inline virtual ~RooKstMMPdf_Unfold() { }

  Double_t evaluate() const ;


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy P4 ;
  RooRealProxy P5 ;
  RooRealProxy P6 ;
  RooRealProxy P8 ;
  RooRealProxy AT2 ;
  RooRealProxy ATRe ;
  RooRealProxy ATIm ;
  RooRealProxy As ;
  RooRealProxy A4s ;
  RooRealProxy A5s ;
  RooRealProxy A7s ;
  RooRealProxy A8s ;
  RooRealProxy Fs ;  
  RooRealProxy q2 ;
  
  //Double_t evaluate() const ;

private:

  ClassDef(RooKstMMPdf_Unfold,1) // Your description goes here...
};
 
#endif
