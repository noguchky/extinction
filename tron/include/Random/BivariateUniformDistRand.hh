#ifndef Tron_BivariateUniformDistRand_hh
#define Tron_BivariateUniformDistRand_hh

#include "ContinuousUniformDistRand.hh"
#include "Vector2.hh"

namespace Tron {

  /// 二変数一様分布
  class BivariateUniformDistRand {
  private:
    ContinuousUniformDistRand fUniformDistRand;
    double         fMeanX;
    double         fMeanY;
    double         fSigmaX;
    double         fSigmaY;
    double         fCorrelation;

  public:
    using ResultType = Vector2<double>;
  
  public:
    BivariateUniformDistRand();
    BivariateUniformDistRand(unsigned int seed);

    void SetParams(double meanX, double meanY, double sigmaX, double sigmaY, double correlation);
  
    ResultType operator()();
  };

}

#endif
