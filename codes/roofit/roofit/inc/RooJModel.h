#ifndef ROOJMODEL_H 
#define ROOJMODEL_H 1

// Include files
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "AngularFunc.h"

#include <vector> 



/** @class RooJModel RooJModel.h
 *  
 *
 *  @author Thomas Blake
 *  @date   2013-10-27
 */


class RooJModel :  public RooAbsPdf  {
public: 
  /// Standard constructor
  RooJModel( const char *name, 
             const char *title, 
             RooAbsReal& costhetal, 
             RooAbsReal& costhetak, 
             RooAbsReal& phi, 
             const RooArgList& J );

  RooJModel( const RooJModel& other,
             const char* name=0 ) ;
  
  virtual ~RooJModel( ); ///< Destructor

  virtual TObject* clone( const char* newname ) const {
    return new RooJModel(*this, newname );
  }

  Int_t getAnalyticalIntegral( RooArgSet& allVars,
                               RooArgSet& analVars,
                               const char* /*rangeName*/) const ;
  
  Double_t analyticalIntegral( Int_t code, const char* rangeName ) const;
  
protected:
  Double_t evaluate() const ;
  
  RooRealProxy costhetal_;
  RooRealProxy costhetak_;
  RooRealProxy phi_;

  RooListProxy J_;
  
  std::vector< AngularFunc::AngularDependence > terms_;

private:
  void BuildTerms() ;
  
  ClassDef(RooJModel, 1)
};
#endif // ROOAMPLITUDEMODEL_H
