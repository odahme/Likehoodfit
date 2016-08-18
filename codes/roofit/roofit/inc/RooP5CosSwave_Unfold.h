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

#ifndef ROOP5COSSWAVEU
#define ROOP5COSSWAVEU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooP5CosSwave_Unfold : public RooAbsPdf {
public:
  
  RooP5CosSwave_Unfold();
  
  RooP5CosSwave_Unfold(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _P5,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A5s,
	   RooAbsReal& _Fs,
	   RooAbsReal& q2
	   );
  RooP5CosSwave_Unfold(const RooP5CosSwave_Unfold& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooP5CosSwave_Unfold(*this,newname); }
  inline virtual ~RooP5CosSwave_Unfold() { }


protected:
  RooRealProxy cosThetaL ;
  RooRealProxy cosThetaK ;
  RooRealProxy phi ;
  RooRealProxy FL ;
  RooRealProxy P5 ;
  RooRealProxy AT2 ;
  RooRealProxy As ;
  RooRealProxy A5s ;
  RooRealProxy Fs ;  
  RooRealProxy q2;
  Double_t evaluate() const ;

private:

  ClassDef(RooP5CosSwave_Unfold,1) // Your description goes here...
};
 
#endif
