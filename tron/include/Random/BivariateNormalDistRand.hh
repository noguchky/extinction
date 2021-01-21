#ifndef Tron_BivariateNormalDistRand_hh
#define Tron_BivariateNormalDistRand_hh

#include "NormalDistRand.hh"
#include "Vector2.hh"

namespace Tron {

  /// 二変数正規分布
  class BivariateNormalDistRand {
  private:
    NormalDistRand fNormalDistRand;
    double         fMeanX;
    double         fMeanY;
    double         fSigmaX;
    double         fSigmaY;
    double         fCorrelation;

  public:
    using ResultType = Vector2<double>;
  
  public:
    BivariateNormalDistRand();
    BivariateNormalDistRand(unsigned int seed);

    void SetParams(double meanX, double meanY, double sigmaX, double sigmaY, double correlation);
  
    ResultType operator()();
  };

}

#endif
