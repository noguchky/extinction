#ifndef Tron_Param_hh
#define Tron_Param_hh

#include "Types.hh"
#include "Math.hh"

namespace Tron {

  template <typename Number_t>
  struct Param_t {
    Number_t Value;
    Number_t Error;

    Param_t() = default;

    Param_t(const Param_t& param) {
      Value = param.Value;
      Error = param.Error;
    }

    Param_t(const Number_t& value, const Number_t& error) {
      Value = value;
      Error = error;
    }

    inline Bool_t operator==(const Param_t& param) const {
      return Value == param.Value && Error == param.Error;
    }

    inline Bool_t operator!=(const Param_t& param) const {
      return !(*this != param);
    }

    inline Bool_t operator<(const Param_t& param) const {
      return Value < param.Value;
    }

    inline Bool_t operator<=(const Param_t& param) const {
      return Value < param.Value;
    }

    inline Bool_t operator>(const Param_t& param) const {
      return Value > param.Value;
    }

    inline Bool_t operator>=(const Param_t& param) const {
      return Value >= param.Value;
    }

    inline const Param_t operator+() const {
      return Param_t(+Value, +Error);
    }

    inline const Param_t operator-() const {
      return Param_t(-Value, +Error);
    }

    inline const Param_t operator+(const Param_t& param) const {
      return Param_t(Value + param.Value, Math::SqrtOfSumOfSquared(Error, param.Error));
    }

    inline const Param_t operator-(const Param_t& param) const {
      return Param_t(Value - param.Value, Math::SqrtOfSumOfSquared(Error, param.Error));
    }

    inline const Param_t operator*(const Number_t& number) const {
      return Param_t(Value * number, Error * number);
    }

    inline const Param_t operator/(const Number_t& number) const {
      return Param_t(Value / number, Error / number);
    }

    inline const Param_t operator*(const Param_t& param) const {
      const Number_t totalError = Math::SqrtOfSumOfSquared(Error * param.Value, param.Error * Value);
      return Param_t(Value * param.Value, totalError);
    }

    inline const Param_t operator/(const Param_t& param) const {
      const Number_t totalError = Math::SqrtOfSumOfSquared(Error / param.Value, param.Error * Value / param.Value / param.Value);
      return Param_t(Value / param.Value, totalError);
    }

    inline Param_t& operator=(const Param_t& param) {
      Value = param.Value;
      Error = param.Error;
      return *this;
    }
    
    inline Param_t& operator+=(const Param_t& param) {
      Value += param.Value;
      Error  = Math::SqrtOfSumOfSquared(Error, param.Error);
      return *this;
    }

    inline Param_t& operator-=(const Param_t& param) {
      Value -= param.Value;
      Error  = Math::SqrtOfSumOfSquared(Error, param.Error);
      return *this;
    }

    inline Param_t& operator*=(const Number_t& number) {
      Value *= number;
      Error *= number;
      return *this;
    }

    inline Param_t& operator/=(const Number_t& number) {
      Value /= number;
      Error /= number;
      return *this;
    }

    inline Param_t& operator*=(const Param_t& param) {
      Value *= param.Value;
      Error  = Math::SqrtOfSumOfSquared(Error * param.Value, param.Error * Value);
      return *this;
    }

    inline Param_t& operator/=(const Param_t& param) {
      Value /= param.Value;
      Error  = Math::SqrtOfSumOfSquared(Error / param.Value, param.Error * Value / param.Value / param.Value);
      return *this;
    }

    inline operator Number_t() const {
      return Value;
    }
  };

  template <typename Number_t>
  inline const Param_t<Number_t> operator*(const Number_t& number, const Param_t<Number_t>& param) {
    return Param_t<Number_t>(param.Value * number, param.Error * number);
  }

  template <typename Number_t>
  inline const Param_t<Number_t> operator/(const Number_t& number, const Param_t<Number_t>& param) {
    return Param_t<Number_t>(param.Value / number, param.Error / number);
  }

  using ParamD_t = Param_t<Double_t>;
  
}

#endif
