#ifndef Tron_CosmicDistRand_hh
#define Tron_CosmicDistRand_hh

#include "NormalDistRand.hh"
#include "Vector3.hh"

namespace Tron {

  /// 原点周辺に降る宇宙線の分布
  class CosmicDistRand {
  private:
    NormalDistRand fNormalDistRand;
    double         fSourceRadius;
    double         fSourceDistance;
    
  public:
    struct ResultType {
      Vector3<double> Position;
      Vector3<double> Direction;
    };

  public:
    CosmicDistRand();
    CosmicDistRand(unsigned int seed);

    void SetParams(double radius, double distance);

    ResultType operator()();
  };
  
}

#endif
