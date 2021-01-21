#ifndef Tron_ContinuousUniformDistRand_hh
#define Tron_ContinuousUniformDistRand_hh

#include <random>

namespace Tron {

  /// 連続一様分布
  class ContinuousUniformDistRand {
  private:
    std::mt19937                     fMersenneTwister;
    std::uniform_real_distribution<> fUniformRealDistribution;

  public:
    using ResultType = double;
  
  public:
    ContinuousUniformDistRand();
    ContinuousUniformDistRand(unsigned int seed);

    void SetParams(double xmin, double xmax);
    
    ResultType operator()();
  };

}

#endif
