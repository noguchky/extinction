#ifndef Tron_NormalDistRand_hh
#define Tron_NormalDistRand_hh

#include <random>
#include "ContinuousUniformDistRand.hh"

namespace Tron {

  /// 正規分布
  class NormalDistRand {
  private:
    std::mt19937               fMersenneTwister;
    std::normal_distribution<> fNormalDistribution;

  public:
    using ResultType = double;
  
  public:
    NormalDistRand();
    NormalDistRand(unsigned int seed);

    void SetParams(double mean, double sigma);
  
    ResultType operator()();
  };

  /// 正規分布 (連続一様分布から生成)
  class NormalDistRand0 {
  private:
    ContinuousUniformDistRand fContinuousUniformDistRand;
    double                    fMean;
    double                    fSigma;

  public:
    using ResultType = double;

  public:
    NormalDistRand0();
    NormalDistRand0(unsigned int seed);

    void SetParams(double mean, double sigma);

    ResultType operator()();
  };
  
}

#endif
