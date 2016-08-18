#ifndef ANGULARFUNC_H 
#define ANGULARFUNC_H 1

#include <math.h> 
#include <cmath>

/* 
   Implementation of integrals associated to the angular distribution 

   Terms inherit from IAngularDependence interface, which defines an 
   evaluation and an integration method.
*/ 

namespace AngularFunc {
  
  class IAngularDependence 
  {
  public :
    virtual double operator()( const double costheta ) const = 0 ;
    
    virtual double integral( const double min, const double max ) const = 0;
  };
  
  class AngularDependence 
  {
  public:
    AngularDependence() : 
      CTL_( 0 ), CTK_( 0 ), PHI_( 0 ){} ;
    
    AngularDependence( IAngularDependence* CTL, 
                       IAngularDependence* CTK,
                       IAngularDependence* PHI ) : 
      CTL_( CTL ), CTK_( CTK ), PHI_( PHI ) {};
    
    double evaluate( const double costhetal,
                     const double costhetak,
                     const double phi ) const 
    {
      return ((*CTL_)( costhetal ))*((*CTK_)( costhetak ))*((*PHI_)( phi ));
    }
    

    ~AngularDependence(){};
    
    const IAngularDependence* CTL_;
    const IAngularDependence* CTK_;
    const IAngularDependence* PHI_;
  };
  
      

  class Flat : virtual public IAngularDependence 
  {
    public:
    double operator()( const double /*phi*/ ) const { return 1.0; }
    
    double integral( const double min, const double max ) const 
    {
      return max - min ;
    }
  };
  
  class CosPhi : virtual public IAngularDependence
  {
  public:
    double operator()( const double phi ) const { return cos(phi); }
    
    double integral( const double min, const double max ) const 
    {
      return sin( max ) - sin( min );
    }
  };
  
  class SinPhi : virtual public IAngularDependence
  {
  public:
    double operator()( const double phi ) const { return sin(phi); }
    
    double integral( const double min, const double max )  const
    {
      return cos( min ) - cos( max );
    }
  };
  
  class Cos2Phi : virtual public IAngularDependence
  {
  public:
    double operator()( const double phi ) const { return cos(2.*phi); }
    
    double integral( const double min, const double max ) const
    {
      return 0.5*( sin( 2.*max ) - sin( 2.*min ) );
    }
  };
  
  class Sin2Phi : virtual public IAngularDependence
  {
  public:
    double operator()( const double phi ) const { return sin(2.*phi); }
    
    double integral( const double min, const double max ) const
    {
      return 0.5*( cos(2.*min) - cos(2.*max) );
    }
  };
  
  
  class CosTheta  : virtual public IAngularDependence
  { 
  public: 
    double operator()( const double costheta ) const { return costheta; }
    
    double integral( const double min, const double max ) const
    {
      return 0.5*(max*max - min*min);
    }
  };
  
  class SinTheta : virtual public IAngularDependence 
  {
  public:
    double operator()( const double costheta ) const { return sqrt( 1. - std::pow( costheta , 2 ) ); }
    
    double integral( const double min, const double max ) const
    {
      return 0.5*( asin( max ) + max*sqrt( 1. - std::pow( max , 2 ) ) - 
                   asin( min ) - min*sqrt( 1. - std::pow( min , 2 ) ) ) ;
    }
  };
  
  class CosSqTheta  : virtual public IAngularDependence
  {
  public: 
    double operator()( const double costheta ) const { return costheta*costheta; } 

    double integral( const double min, const double max )  const
    {
      return (1./3.)*( std::pow(max,3) - std::pow(min,3));
    }
  };

  class SinSqTheta  : virtual public IAngularDependence
  {
  public: 
    double operator()( const double costheta ) const { return 1. - costheta*costheta; }
    
    double integral( const double min, const double max )  const
    {
      return (max - min) - (1./3.)*( std::pow(max,3) - std::pow(min,3));
    }
  };  
  
  class Sin2Theta  : virtual public IAngularDependence
  {
  public:
    
    double operator()( const double costheta ) const { return 2.*sqrt( 1. - std::pow( costheta , 2 ) )*costheta ; }
    
    double integral( const double min, const double max ) const
    {
      return (2./3.)*( std::pow(1. - min*min,3./2.) - std::pow( 1. - max*max , 3./2. ) );
    }
  };
  
  class Cos2Theta : virtual public IAngularDependence
  {
  public:
    double operator()( const double costheta ) const { return 2.*std::pow(costheta,2) - 1. ;}
    
    double integral( const double min, const double max ) const
    {
      return (2./3.)*(std::pow( max, 3 ) - std::pow( min, 3 )) - max + min ;
    }
  };

  namespace Functors 
  {
    extern IAngularDependence* FLAT;
    extern IAngularDependence* COSPHI;
    extern IAngularDependence* COS2PHI;
    extern IAngularDependence* SINPHI;
    extern IAngularDependence* SIN2PHI;
    extern IAngularDependence* COSTHETA;
    extern IAngularDependence* SINTHETA;
    extern IAngularDependence* COSSQTHETA;
    extern IAngularDependence* SINSQTHETA;
    extern IAngularDependence* COS2THETA;
    extern IAngularDependence* SIN2THETA;
    
    void Instantiate();

  }
};



#endif // ANGULARFUNC_H
