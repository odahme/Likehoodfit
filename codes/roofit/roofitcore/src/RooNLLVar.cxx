/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Class RooNLLVar implements a a -log(likelihood) calculation from a dataset
// and a PDF. The NLL is calculated as 
// <pre>
//  Sum[data] -log( pdf(x_data) )
// </pre>
// In extended mode, a (Nexpect - Nobserved*log(NExpected) term is added
// END_HTML
//

#include "RooFit.h"
#include "Riostream.h"

#include "RooNLLVar.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooCmdConfig.h"
#include "RooMsgService.h"
#include "RooAbsDataStore.h"
#include "RooRealMPFE.h"

#include "RooRealVar.h"


ClassImp(RooNLLVar)
;

RooArgSet RooNLLVar::_emptySet ;


//_____________________________________________________________________________
RooNLLVar::RooNLLVar(const char *name, const char* title, RooAbsPdf& pdf, RooAbsData& indata,
		     const RooCmdArg& arg1, const RooCmdArg& arg2,const RooCmdArg& arg3,
		     const RooCmdArg& arg4, const RooCmdArg& arg5,const RooCmdArg& arg6,
		     const RooCmdArg& arg7, const RooCmdArg& arg8,const RooCmdArg& arg9) :
  RooAbsOptTestStatistic(name,title,pdf,indata,
			 *(const RooArgSet*)RooCmdConfig::decodeObjOnTheFly("RooNLLVar::RooNLLVar","ProjectedObservables",0,&_emptySet
									    ,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9),
			 RooCmdConfig::decodeStringOnTheFly("RooNLLVar::RooNLLVar","RangeWithName",0,"",arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9),
			 RooCmdConfig::decodeStringOnTheFly("RooNLLVar::RooNLLVar","AddCoefRange",0,"",arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9),
			 RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","NumCPU",0,1,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9),
			 RooFit::BulkPartition,
			 RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","Verbose",0,1,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9),
			 RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","SplitRange",0,0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9),
			 RooCmdConfig::decodeIntOnTheFly("RooNLLVar::RooNLLVar","CloneData",0,1,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9))             
{
  // Construct likelihood from given p.d.f and (binned or unbinned dataset)
  //
  //  Extended()     -- Include extended term in calculation
  //  NumCPU()       -- Activate parallel processing feature
  //  Range()        -- Fit only selected region
  //  SumCoefRange() -- Set the range in which to interpret the coefficients of RooAddPdf components 
  //  SplitRange()   -- Fit range is split by index catory of simultaneous PDF
  //  ConditionalObservables() -- Define conditional observables 
  //  Verbose()      -- Verbose output of GOF framework classes
  //  CloneData()    -- Clone input dataset for internal use (default is kTRUE)

  RooCmdConfig pc("RooNLLVar::RooNLLVar") ;
  pc.allowUndefined() ;
  pc.defineInt("extended","Extended",0,kFALSE) ;

  pc.process(arg1) ;  pc.process(arg2) ;  pc.process(arg3) ;
  pc.process(arg4) ;  pc.process(arg5) ;  pc.process(arg6) ;
  pc.process(arg7) ;  pc.process(arg8) ;  pc.process(arg9) ;

  _extended = pc.getInt("extended") ;
  _weightSq = kFALSE ;
  _first = kTRUE ;
  _offset = 0 ;

}



//_____________________________________________________________________________
RooNLLVar::RooNLLVar(const char *name, const char *title, RooAbsPdf& pdf, RooAbsData& indata,
		     Bool_t extended, const char* rangeName, const char* addCoefRangeName,
		     Int_t nCPU, RooFit::MPSplit interleave, Bool_t verbose, Bool_t splitRange, Bool_t cloneData) : 
  RooAbsOptTestStatistic(name,title,pdf,indata,RooArgSet(),rangeName,addCoefRangeName,nCPU,interleave,verbose,splitRange,cloneData),
  _extended(extended),
  _weightSq(kFALSE),
  _first(kTRUE)
{
  // Construct likelihood from given p.d.f and (binned or unbinned dataset)
  // For internal use.

}



//_____________________________________________________________________________
RooNLLVar::RooNLLVar(const char *name, const char *title, RooAbsPdf& pdf, RooAbsData& indata,
		     const RooArgSet& projDeps, Bool_t extended, const char* rangeName,const char* addCoefRangeName, 
		     Int_t nCPU,RooFit::MPSplit interleave,Bool_t verbose, Bool_t splitRange, Bool_t cloneData) : 
  RooAbsOptTestStatistic(name,title,pdf,indata,projDeps,rangeName,addCoefRangeName,nCPU,interleave,verbose,splitRange,cloneData),
  _extended(extended),
  _weightSq(kFALSE),
  _first(kTRUE)
{
  // Construct likelihood from given p.d.f and (binned or unbinned dataset)
  // For internal use.  


}



//_____________________________________________________________________________
RooNLLVar::RooNLLVar(const RooNLLVar& other, const char* name) : 
  RooAbsOptTestStatistic(other,name),
  _extended(other._extended),
  _weightSq(other._weightSq),
  _first(kTRUE)
{
  // Copy constructor
}




//_____________________________________________________________________________
RooNLLVar::~RooNLLVar()
{
  // Destructor
}




//_____________________________________________________________________________
void RooNLLVar::applyWeightSquared(Bool_t flag) 
{ 
  if (_gofOpMode==Slave) {
    _weightSq = flag ; 
    setValueDirty() ; 

  } else if ( _gofOpMode==MPMaster) {

    for (Int_t i=0 ; i<_nCPU ; i++) {
      _mpfeArray[i]->applyNLLWeightSquared(flag) ;
    }    

  } else if ( _gofOpMode==SimMaster) {

    for (Int_t i=0 ; i<_nGof ; i++) {
      ((RooNLLVar*)_gofArray[i])->applyWeightSquared(flag) ;
    }

  }
} 



//_____________________________________________________________________________
Double_t RooNLLVar::evaluatePartition(Int_t firstEvent, Int_t lastEvent, Int_t stepSize) const 
{
  // Calculate and return likelihood on subset of data from firstEvent to lastEvent
  // processed with a step size of 'stepSize'. If this an extended likelihood and
  // and the zero event is processed the extended term is added to the return
  // likelihood.

  // Throughout the calculation, we use Kahan's algorithm for summing to
  // prevent loss of precision - this is a factor four more expensive than
  // straight addition, but since evaluating the PDF is usually much more
  // expensive than that, we tolerate the additional cost...
  Int_t i ;
  Double_t result(0), carry(0);
  
  RooAbsPdf* pdfClone = (RooAbsPdf*) _funcClone ;

  // cout << "RooNLLVar::evaluatePartition(" << GetName() << ") projDeps = " << (_projDeps?*_projDeps:RooArgSet()) << endl ;

  _dataClone->store()->recalculateCache( _projDeps, firstEvent, lastEvent, stepSize ) ;

  Double_t sumWeight(0), sumWeightCarry(0);
  for (i=firstEvent ; i<lastEvent ; i+=stepSize) {
    
    // get the data values for this event
    //Double_t wgt = _dataClone->weight(i) ;
    //if (wgt==0) continue ;

    _dataClone->get(i) ;
    //cout << "NLL - now loading event #" << i << endl ;
//     _funcObsSet->Print("v") ;
    

    if (!_dataClone->valid()) continue;

    Double_t eventWeight = _dataClone->weight();
    if (0. == eventWeight * eventWeight) continue ;
    if (_weightSq) eventWeight = _dataClone->weightSquared() ;

    Double_t term = -eventWeight * pdfClone->getLogVal(_normSet);

    Double_t y = eventWeight - sumWeightCarry;
    Double_t t = sumWeight + y;
    sumWeightCarry = (t - sumWeight) - y;
    sumWeight = t;

    y = term - carry;
    t = result + y;
    carry = (t - result) - y;
    result = t;
  }
  
  // include the extended maximum likelihood term, if requested
  if(_extended && _setNum==_extSet) {
    if (_weightSq) {

      // Calculate sum of weights-squared here for extended term
      Double_t sumW2(0), sumW2carry(0);
      for (i=0 ; i<_dataClone->numEntries() ; i++) {
	_dataClone->get(i);
	Double_t y = _dataClone->weightSquared() - sumW2carry;
	Double_t t = sumW2 + y;
	sumW2carry = (t - sumW2) - y;
	sumW2 = t;
      }
      Double_t y = pdfClone->extendedTerm(sumW2 , _dataClone->get()) - carry;
      Double_t t = result + y;
      carry = (t - result) - y;
      result = t;
    } else {
      Double_t y = pdfClone->extendedTerm(_dataClone->sumEntries(), _dataClone->get()) - carry;
      Double_t t = result + y;
      carry = (t - result) - y;
      result = t;
    }
  }    

  // If part of simultaneous PDF normalize probability over 
  // number of simultaneous PDFs: -sum(log(p/n)) = -sum(log(p)) + N*log(n) 
  if (_simCount>1) {
    Double_t y = sumWeight*log(1.0*_simCount) - carry;
    Double_t t = result + y;
    carry = (t - result) - y;
    result = t;
  }

  //timer.Stop() ;
  //cout << "RooNLLVar::evalPart(" << GetName() << ") SET=" << _setNum << " first=" << firstEvent << ", last=" << lastEvent << ", step=" << stepSize << ") result = " << result << " CPU = " << timer.CpuTime() << endl ;

  // At the end of the first full calculation, wire the caches
  if (_first) {
    _first = kFALSE ;
    _funcClone->wireAllCaches() ;
  }
  

  // Check if value offset flag is set.
  if (_doOffset) {
    
    // If no offset is stored enable this feature now
    if (_offset==0 && result !=0 ) {
      coutI(Minimization) << "RooNLLVar::evaluatePartition(" << GetName() << ") first = "<< firstEvent << " last = " << lastEvent << " Likelihood offset now set to " << result << endl ;
      _offset = result ;
      _offsetCarry = carry;
    }

    // Substract offset
    Double_t y = -_offset - (carry + _offsetCarry);
    Double_t t = result + y;
    carry = (t - result) - y;
    result = t;
  }
    
  return result ;
}




