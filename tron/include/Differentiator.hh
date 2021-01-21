#ifndef Tron_Differentiator_hh
#define Tron_Differentiator_hh

#include <functional>
#include "Types.hh"
#include "Singleton.hh"

namespace Tron {

  namespace Differentiator {

    using Func_t = std::function<Double_t(Double_t)>;

    class Differentiator {
      TRON_SINGLETONIZE(Differentiator);
    public:
      virtual Double_t Differentiate(const Func_t& func, Double_t x, Double_t dx) const;
    };

    class  NewtonsDifference : public Differentiator {
      TRON_SINGLETONIZE(NewtonsDifference);
    public:
      virtual Double_t Differentiate(const Func_t& func, Double_t x, Double_t dx) const override;
    };

    class SymmetricDifference : public Differentiator {
      TRON_SINGLETONIZE(SymmetricDifference);
    public:
      virtual Double_t Differentiate(const Func_t& func, Double_t x, Double_t dx) const override;
    };

    class FivePointStencil : public Differentiator {
      TRON_SINGLETONIZE(FivePointStencil);
    public:
      virtual Double_t Differentiate(const Func_t& func, Double_t x, Double_t dx) const override;
    };
    
  }

}

#endif
