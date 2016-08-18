// Include files 



// local
#include "AngularFunc.h"

namespace AngularFunc 
{
  namespace Functors 
  {
    IAngularDependence* FLAT;
    IAngularDependence* COSPHI;
    IAngularDependence* COS2PHI;
    IAngularDependence* SINPHI;
    IAngularDependence* SIN2PHI;
    IAngularDependence* COSTHETA;
    IAngularDependence* SINTHETA;
    IAngularDependence* COSSQTHETA;
    IAngularDependence* SINSQTHETA;
    IAngularDependence* COS2THETA;
    IAngularDependence* SIN2THETA;
    
    void Instantiate()
    {
      if ( !FLAT )       FLAT       = new Flat();
      if ( !COSPHI )     COSPHI     = new CosPhi();
      if ( !SINPHI )     SINPHI     = new SinPhi();
      if ( !COS2PHI )    COS2PHI    = new Cos2Phi();
      if ( !SIN2PHI )    SIN2PHI    = new Sin2Phi();
      if ( !COSTHETA )   COSTHETA   = new CosTheta();
      if ( !SINTHETA )   SINTHETA   = new SinTheta();
      if ( !COS2THETA )  COS2THETA  = new Cos2Theta();
      if ( !SIN2THETA )  SIN2THETA  = new Sin2Theta();
      if ( !COSSQTHETA ) COSSQTHETA = new CosSqTheta();
      if ( !SINSQTHETA ) SINSQTHETA = new SinSqTheta();
    } 
  }
}
