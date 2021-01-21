#ifndef Tron_Ellipse_hh
#define Tron_Ellipse_hh

#include <string>
#include "Types.hh"

class TGraph;
class TF1;
class TF2;
class TMarker;

namespace Tron {

  class Ellipse {
  protected:
    Double_t    fMuX;
    Double_t    fMuY;
    Double_t    fSigmaX;
    Double_t    fSigmaY;
    Double_t    fCorrelation;
    Int_t       fNpx = 120;
  
  public:
    Ellipse();
    virtual ~Ellipse() = default;

    virtual        void     SetMuX(Double_t mu);
    virtual        void     SetMuY(Double_t mu);
    virtual        void     SetSigmaX(Double_t sigma);
    virtual        void     SetSigmaY(Double_t sigma);
    virtual        void     SetCorrelation(Double_t correlation);

    virtual inline Double_t GetMuX() const { return fMuX; }
    virtual inline Double_t GetMuY() const { return fMuY; }
    virtual inline Double_t GetSigmaX() const { return fSigmaX; }
    virtual inline Double_t GetSigmaY() const { return fSigmaY; }
    virtual inline Double_t GetCorrelation() const { return fCorrelation; };

    virtual        TGraph*  CreateGraph(const std::string& name);
    virtual        TF1*     CreateFuncX(const std::string& name);
    virtual        TF1*     CreateFuncY(const std::string& name);
    virtual        TF2*     CreateFunc (const std::string& name);
    virtual        TMarker* CreateCenterMarker();

    static Double_t  GetSemiAxisX(Double_t sigmaX, Double_t sigmaY, Double_t correlation);
    static Double_t  GetSemiAxisY(Double_t sigmaX, Double_t sigmaY, Double_t correlation);
    static Double_t  GetRotationAngle(Double_t sigmaX, Double_t sigmaY, Double_t correlation);
    static Double_t  GetX(Double_t theta, Double_t semiAxisX, Double_t semiAxisY, Double_t rotationAngle);
    static Double_t  GetY(Double_t theta, Double_t semiAxisX, Double_t semiAxisY, Double_t rotationAngle);

  };

}

#endif
