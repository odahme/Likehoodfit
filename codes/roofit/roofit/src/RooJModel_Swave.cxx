// Include files 

#include "RooJModel_Swave.h"
#include "TIterator.h"
#include "TMath.h"

#include <iostream> 

using namespace AngularFunc::Functors ;

ClassImp( RooJModel_Swave )

//-----------------------------------------------------------------------------
// Implementation file for class : RooJModel_Swave
//
// 2013-10-27 : Thomas Blake
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================

RooJModel_Swave::RooJModel_Swave( const char *name, 
                      const char *title,
                      RooAbsReal& costhetal, 
                      RooAbsReal& costhetak,
                      RooAbsReal& phi,  
                      const RooArgList& J ) : 
  RooAbsPdf( name , title ),
  costhetal_( "costhetal", "costhetal", this, costhetal ),
  costhetak_( "costhetak", "costhetak", this, costhetak ),
  phi_( "phi", "phi", this, phi ),
  J_("J","J",this)
{
  assert( 17 == J.getSize() );

  TIterator* iter = J.createIterator() ;
  RooAbsArg* arg  = 0;
  
  while ((arg = (RooAbsArg*) iter->Next())) { J_.add( *arg ); }

  BuildTerms() ;
}



RooJModel_Swave::RooJModel_Swave( const RooJModel_Swave& other,
                      const char* name ) : 
  RooAbsPdf(other,name),
  costhetal_("costhetal",this,other.costhetal_),
  costhetak_("costhetak",this,other.costhetak_),
  phi_      ("phi",this,other.phi_),
  J_        ("J",this,other.J_)
{
  BuildTerms();
}

void RooJModel_Swave::BuildTerms()
{
  AngularFunc::Functors::Instantiate();
  ///////////////////////////////////////////      Theta_l,     theta_k,     phi
  terms_.push_back( AngularFunc::AngularDependence( FLAT,       SINSQTHETA, FLAT ) );
  terms_.push_back( AngularFunc::AngularDependence( FLAT,       COSSQTHETA, FLAT ) ); // J1
  terms_.push_back( AngularFunc::AngularDependence( COS2THETA,  SINSQTHETA, FLAT ) ); // J2s
  terms_.push_back( AngularFunc::AngularDependence( COS2THETA,  COSSQTHETA, FLAT ) ); // J2c 
  terms_.push_back( AngularFunc::AngularDependence( SINSQTHETA, SINSQTHETA, COS2PHI ) ); // J3
  terms_.push_back( AngularFunc::AngularDependence( SIN2THETA,  SIN2THETA,  COSPHI ) );  // J4
  terms_.push_back( AngularFunc::AngularDependence( SINTHETA,   SIN2THETA,  COSPHI ) ); // J5
  terms_.push_back( AngularFunc::AngularDependence( COSTHETA,   SINSQTHETA, FLAT ) );   // J6
  terms_.push_back( AngularFunc::AngularDependence( SINTHETA,   SIN2THETA,  SINPHI ) ); // J7
  terms_.push_back( AngularFunc::AngularDependence( SIN2THETA,  SIN2THETA,  SINPHI ) ); // J8
  terms_.push_back( AngularFunc::AngularDependence( SINSQTHETA, SINSQTHETA, SIN2PHI ) ); // J9
  terms_.push_back( AngularFunc::AngularDependence( SINSQTHETA, FLAT, FLAT ) ); //Fs
  terms_.push_back( AngularFunc::AngularDependence( SINSQTHETA, COSTHETA, FLAT ) ); //As
  terms_.push_back( AngularFunc::AngularDependence( SIN2THETA, SINTHETA , COSPHI ) );  //As4 
  terms_.push_back( AngularFunc::AngularDependence( SINTHETA  , SINTHETA , COSPHI ) );  //As5
  terms_.push_back( AngularFunc::AngularDependence( SINTHETA  , SINTHETA , SINPHI ) );  //As7
  terms_.push_back( AngularFunc::AngularDependence( SIN2THETA  , SINTHETA , SINPHI ) );  //As8 





}



//=============================================================================
// Destructor
//=============================================================================

RooJModel_Swave::~RooJModel_Swave() {} 

Double_t RooJModel_Swave::evaluate() const 
{  
  Double_t val = 0.0;
  
  for ( int i = 0; i < J_.getSize(); ++i )
  {
    val += (((RooAbsReal&)J_[i]).getVal())*(terms_[i].evaluate( costhetal_, costhetak_, phi_));
  }

  return (9.*val)/(32.*TMath::Pi());
}



Int_t RooJModel_Swave::getAnalyticalIntegral( RooArgSet& allVars,
                                                RooArgSet& analVars,
                                                const char* /*rangeName*/) const
{
  
  
  // Match all 3 
  if (matchArgs(allVars,analVars,costhetal_,costhetak_,phi_)) return 7;
    
  // Match 2 of 3
  if (matchArgs(allVars,analVars,costhetal_,costhetak_)) return 3;
  if (matchArgs(allVars,analVars,costhetal_,phi_)) return 5;
  if (matchArgs(allVars,analVars,costhetak_,phi_)) return 6;
  
  // Match 1 of 3
  if (matchArgs(allVars,analVars,costhetal_)) return 1 ;
  if (matchArgs(allVars,analVars,costhetak_)) return 2;
  if (matchArgs(allVars,analVars,phi_)) return 4;
  
  return 0 ;
}


Double_t RooJModel_Swave::analyticalIntegral( Int_t code,
                                        const char* rangeName) const
{
  Double_t minCTL = costhetal_.min(rangeName);
  Double_t maxCTL = costhetal_.max(rangeName);
  Double_t minCTK = costhetak_.min(rangeName);
  Double_t maxCTK = costhetak_.max(rangeName);
  Double_t minPHI = phi_.min(rangeName);
  Double_t maxPHI = phi_.max(rangeName);
  
  Double_t val = 0;
  
  for ( int i = 0; i < J_.getSize(); ++i )
  {
    Double_t func = 1.0;
    
    if ( code & 0x01 ) { 
      func *= terms_[i].CTL_->integral( minCTL, maxCTL );
    } else {
      func *= (*terms_[i].CTL_)(costhetal_);
    }
    
    if ( code & 0x02 ) { 
      func *= terms_[i].CTK_->integral( minCTK, maxCTK );
    } else {
      func *= (*terms_[i].CTK_)(costhetak_);
    }
    
    if ( code & 0x04 ) { 
      func *= terms_[i].PHI_->integral( minPHI, maxPHI );
    } else {
      func *= (*terms_[i].PHI_)(phi_);
    }
    
    val += ((RooAbsReal&)J_[i]).getVal()*func;
  }
  
  return (9.*val)/(32.*TMath::Pi());
}

//=============================================================================
