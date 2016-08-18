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

#ifndef ROOS8COSSWAVEU
#define ROOS8COSSWAVEU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS8CosSwave_Unfold : public RooAbsPdf {
public:
  
  RooS8CosSwave_Unfold();
  
  RooS8CosSwave_Unfold(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S8,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A7s,
		       RooAbsReal& _Fs,
		       RooAbsReal& _q2
	   );
  RooS8CosSwave_Unfold(const RooS8CosSwave_Unfold& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS8CosSwave_Unfold(*this,newname); }
  inline virtual ~RooS8CosSwave_Unfold() { }


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
  RooRealProxy q2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooS8CosSwave_Unfold,1) // Your description goes here...
};
 
#endif
