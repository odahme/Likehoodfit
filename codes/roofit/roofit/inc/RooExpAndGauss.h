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

#ifndef ROOEXPANDGAUSS
#define ROOEXPANDGAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

 
class RooExpAndGauss : public RooAbsPdf {
public:
  RooExpAndGauss();
  RooExpAndGauss(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _sh_mean,
	      RooAbsReal& _sh_sigma,
	      RooAbsReal& _sh_trans);
  RooExpAndGauss(const RooExpAndGauss& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpAndGauss(*this,newname); }
  inline virtual ~RooExpAndGauss() { }


protected:

  RooRealProxy x ;
  RooRealProxy sh_mean ;
  RooRealProxy sh_sigma ;
  RooRealProxy sh_trans ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooExpAndGauss,1) // Your description goes here...
};
 
#endif
