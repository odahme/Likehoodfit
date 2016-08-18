#ifndef ROOJMODELSW_H 
#define ROOJMODELW_H 1

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
 *  @author Marcin Chrzaszcz,based on Thomas Blake
 *  @date   2014-06-21
 */


class RooJModel_Swave :  public RooAbsPdf  {
public: 
  /// Standard constructor
  RooJModel_Swave( const char *name, 
             const char *title, 
             RooAbsReal& costhetal, 
             RooAbsReal& costhetak, 
             RooAbsReal& phi, 
             const RooArgList& J );

  RooJModel_Swave( const RooJModel_Swave& other,
             const char* name=0 ) ;
  
  virtual ~RooJModel_Swave( ); ///< Destructor

  virtual TObject* clone( const char* newname ) const {
    return new RooJModel_Swave(*this, newname );
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
  
  ClassDef(RooJModel_Swave, 1)
};
#endif // ROOAMPLITUDEMODEL_H
