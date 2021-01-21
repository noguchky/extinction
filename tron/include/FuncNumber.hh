#ifndef Tron_FuncNumber_hh
#define Tron_FuncNumber_hh

#include <iostream>
#include <functional>
#include "Types.hh"

namespace Tron {

  template <typename Number_t> 
  class FuncNumber_t {
  private:
    Number_t                          fValue;
    Number_t                          fDefaultArg;
    std::function<Number_t(Number_t)> fFunc;
    std::function<Number_t(Number_t)> fConst = [this](Number_t) { return fValue; };

  public:
    FuncNumber_t() noexcept : fFunc() {
    }

    FuncNumber_t(std::nullptr_t f) noexcept : fFunc(f) {
    }

    FuncNumber_t(const FuncNumber_t& f) : fFunc(f.fFunc) {
    }

    FuncNumber_t(FuncNumber_t&& f) : fFunc(f.fFunc) {
    }

    FuncNumber_t(const std::function<Number_t(Number_t)>& f) : fFunc(f) {
    }

    FuncNumber_t(std::function<Number_t(Number_t)>&& f) : fFunc(f) {
    }

    template <class F>
    FuncNumber_t(F f) : fFunc(f) {
    }

    FuncNumber_t(Number_t value) : fValue(value), fFunc(fConst) {
    }

    FuncNumber_t& operator=(const FuncNumber_t& f) {
      fFunc = f.fFunc;
      return *this;
    }

    FuncNumber_t& operator=(FuncNumber_t&& f) {
      fFunc = f.fFunc;
      return *this;
    }

    FuncNumber_t& operator=(const std::function<Number_t(Number_t)>& f) {
      fFunc = f;
      return *this;
    }

    FuncNumber_t& operator=(std::function<Number_t(Number_t)>&& f) {
      fFunc = f;
      return *this;
    }

    FuncNumber_t& operator=(std::nullptr_t f) {
      fFunc = f;
      return *this;
    }

    template<class F>
    FuncNumber_t& operator=(F&& f) {
      fFunc = f;
      return *this;
    }

    template<class F>
    FuncNumber_t& operator=(std::reference_wrapper<F> f) noexcept {
      fFunc = f;
      return *this;
    }

    FuncNumber_t& operator=(Number_t value) {
      fValue = value;
      fFunc = fConst;
      return *this;
    }

    inline Number_t operator()(Number_t value) const {
      return fFunc ? fFunc(value) : Number_t();
    }

    inline operator Number_t() {
      return fFunc ? fFunc(fDefaultArg) : Number_t();
    }

    inline void SetDefaultArgument(Number_t value) {
      fDefaultArg = value;
    }

    inline Number_t GetDefaultArgument() const {
      return fDefaultArg;
    }
    
  };

}

#endif
