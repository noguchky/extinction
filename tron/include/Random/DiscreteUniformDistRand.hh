#ifndef Tron_DiscreteUniformDistRand_hh
#define Tron_DiscreteUniformDistRand_hh

#include <random>

namespace Tron {

  /// 離散一様分布
  class DiscreteUniformDistRand {
  private:
    std::mt19937                    fMersenneTwister;
    std::uniform_int_distribution<> fUniformIntDistribution;

  public:
    using ResultType = double;
  
  public:
    DiscreteUniformDistRand();
    DiscreteUniformDistRand(unsigned int seed);
    
    void SetParams(int xmin, int xmax);
    
    ResultType operator()();
  };

}

#endif
