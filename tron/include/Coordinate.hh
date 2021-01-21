#ifndef Tron_Coordinate_hh
#define Tron_Coordinate_hh

#include "TVector3.h"
#include "Types.hh"

namespace Tron {

  class Coordinate {
  protected:
    TVector3 fOrigin;
    TVector3 fXaxis;
    TVector3 fYaxis;
    TVector3 fZaxis;
    
  public:
    Coordinate();
    virtual ~Coordinate();

    virtual       TVector3  GetLocalPosition(const TVector3& globalPosition) const;
    virtual       TVector3  GetGlobalPosition(const TVector3& localPosition) const;

    virtual const TVector3& GetOrigin() const;
    virtual       void      SetOrigin(const TVector3& origin);
    virtual       void      ShiftOrigin(const TVector3 dx);
    virtual       void      ResetOrigin();

    virtual const TVector3& GetXaxis() const;
    virtual const TVector3& GetYaxis() const;
    virtual const TVector3& GetZaxis() const;
    virtual       void      SetXaxis(const TVector3& xaxis);
    virtual       void      SetYaxis(const TVector3& yaxis);
    virtual       void      SetZaxis(const TVector3& zaxis);
    virtual       void      RotateAxes(Double_t theta, const TVector3& axis);
    virtual       void      ResetAxes();
    
  };

}

#endif
