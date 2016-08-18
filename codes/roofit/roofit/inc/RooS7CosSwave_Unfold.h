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

#ifndef ROOS7COSSWAVEU
#define ROOS7COSSWAVEU

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h" 

class RooS7CosSwave_Unfold : public RooAbsPdf {
public:
  
  RooS7CosSwave_Unfold();
  
  RooS7CosSwave_Unfold(const char *name, const char *title,
	   RooAbsReal& _cosThetaL,
	   RooAbsReal& _cosThetaK,
	   RooAbsReal& _phi,
	   RooAbsReal& _FL,
	   RooAbsReal& _S7,
	   RooAbsReal& _AT2,
	   RooAbsReal& _As,
	   RooAbsReal& _A7s,
		       RooAbsReal& _Fs,
		       RooAbsReal& _q2
	   );
  RooS7CosSwave_Unfold(const RooS7CosSwave_Unfold& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooS7CosSwave_Unfold(*this,newname); }
  inline virtual ~RooS7CosSwave_Unfold() { }


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
  RooRealProxy q2;

  Double_t evaluate() const ;

private:

  ClassDef(RooS7CosSwave_Unfold,1) // Your description goes here...
};
 
#endif
