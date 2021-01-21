#ifndef Tron_HoughTransformator_hh
#define Tron_HoughTransformator_hh

#include <vector>
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "Types.hh"

namespace Tron {

  struct HoughLine {
    static Double_t R(Double_t* x, Double_t* p);
    static Double_t Y(Double_t* x, Double_t* p);
    static Double_t Distance(Double_t x, Double_t y, Double_t theta, Double_t r);
  };

  class HoughTransformator {
  private:
    TH2*                  fTransResult;
    TF2*                  fChiSquare;

    std::vector<Double_t> fX;
    std::vector<Double_t> fY;
    std::vector<Double_t> fW;

    Double_t              fTheta;
    Double_t              fR;
  
  public:
    HoughTransformator(const Char_t* name,
                       Int_t nbinsr, Double_t rlow, Double_t rup,
                       Int_t nbinst, Double_t tlow, Double_t tup);
    ~HoughTransformator();

    inline Double_t GetTheta() const { return fTheta; }
    inline Double_t GetR()     const { return fR;     }

    Int_t           GetN() const;
    void            SetPoints(Int_t n, const Double_t* x, const Double_t* y, const Double_t* w);
    void            SetPoints(Int_t n, const Double_t* x, const Double_t* y, Double_t w = 1.0);

    const TH2*      Transform(Int_t smoothTimes = 1);
    Int_t           RemoveFarPoints(Double_t theta, Double_t r, Double_t threshold);
    const TF2*      MinimizeChiSquare();

    TF1*            CreateHoughLine(const Char_t* name, Double_t tmin, Double_t tmax);
    TF1*            CreateTrack(const Char_t* name, Double_t xmin, Double_t xmax);

  private:
    void            Vote(Double_t x, Double_t y, Double_t w);
    Double_t        ChiSquare(Double_t* x, Double_t* p);
    Double_t        EvalChiSquare(Double_t theta, Double_t r);

  };

}

#endif
