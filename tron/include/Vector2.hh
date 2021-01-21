#ifndef Tron_Vector2_hh
#define Tron_Vector2_hh

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "TVector2.h"

inline TVector2 operator+(const TVector2& vector) {
  return vector;
}

namespace Tron {

  template <typename Number_t>
  struct Vector2 {
    Number_t X;
    Number_t Y;

    Vector2() = default;

    Vector2(Number_t x, Number_t y) {
      X = x;
      Y = y;
    }

    Vector2(const Vector2<Number_t>& vector) {
      X = vector.X;
      Y = vector.Y;
    }

    inline Number_t& operator[](int i) {
      switch (i) {
      case 0: return X;
      case 1: return Y;
      }
      std::cout << "Vector2::operator[]() [error] out of range, " << i << std::endl;
      return Y;
    }

    inline const Number_t& operator[](int i) const {
      switch (i) {
      case 0: return X;
      case 1: return Y;
      }
      std::cout << "Vector2::operator[]() [error] out of range, " << i << std::endl;
      return Y;
    }

    inline bool operator==(const Vector2& vector) const {
      return X == vector.X && Y == vector.Y;
    }

    inline bool operator!=(const Vector2& vector) const {
      return !(*this == vector);
    }

    inline Vector2 operator+() const {
      return Vector2(+X, +Y);
    }

    inline Vector2 operator-() const {
      return Vector2(-X, -Y);
    }

    inline Vector2 operator+(const Vector2& vector) const {
      return Vector2(X + vector.X, Y + vector.Y);
    }

    inline Vector2 operator-(const Vector2& vector) const {
      return Vector2(X - vector.X, Y - vector.Y);
    }

    inline Vector2 operator*(const Number_t& number) const {
      return Vector2(X * number, Y * number);
    }

    inline Vector2 operator/(const Number_t& number) const {
      return Vector2(X / number, Y / number);
    }

    inline Vector2& operator=(const Vector2& vector) {
      X = vector.X;
      Y = vector.Y;
      return *this;
    }
    
    inline Vector2& operator+=(const Vector2& vector) {
      X += vector.X;
      Y += vector.Y;
      return *this;
    }

    inline Vector2& operator-=(const Vector2& vector) {
      X -= vector.X;
      Y -= vector.Y;
      return *this;
    }

    inline Vector2& operator*=(const Number_t& number) {
      X *= number;
      Y *= number;
      return *this;
    }

    inline Vector2& operator/=(const Number_t& number) {
      X /= number;
      Y /= number;
      return *this;
    }

    inline Number_t Length() const {
      // Ignore: Number_t is IMAGINARY NUMBER
      return X * X + Y * Y;
    }

    inline Vector2& Rotate(double theta) {
      const double c = std::cos(theta), s = std::sin(theta);
      const double x = X * c - Y * s;
      const double y = X * s + Y * c;
      X = x;
      Y = y;
      return *this;
    }

    inline static Number_t InnerProduct(const Vector2& v1, const Vector2& v2) {
      return v1.X * v2.X + v1.Y * v2.Y;
    }

    inline Number_t InnerProduct(const Vector2& vector) const {
      return InnerProduct(*this, vector);
    }

  };

  template <typename Number_t>
  inline Vector2<Number_t> operator*(const Number_t& number, const Vector2<Number_t>& vector) {
    return vector * number;
  }

}

#endif
