#ifndef Tron_Vector3_hh
#define Tron_Vector3_hh

#include <iostream>
#include <stdexcept>
#include <cmath>
#include "TVector3.h"

inline TVector3 operator+(const TVector3& vector) {
  return vector;
}

inline TVector3 operator/(const TVector3& vector, Double_t value) {
  return vector * (1.0 / value);
}

inline TVector3& operator/=(TVector3& vector, Double_t value) {
  vector = vector / value;
  return vector;
}

namespace Tron {

  template <typename Number_t>
  struct Vector3 {
    Number_t X;
    Number_t Y;
    Number_t Z;

    Vector3() = default;

    Vector3(Number_t x, Number_t y, Number_t z) {
      X = x;
      Y = y;
      Z = z;
    }

    Vector3(const Vector3& vector) {
      X = vector.X;
      Y = vector.Y;
      Z = vector.Z;
    }

    inline Number_t& operator[](int i) {
      switch (i) {
      case 0: return X;
      case 1: return Y;
      case 2: return Z;
      }
      std::cout << "Vector3::operator[]() [error] out of range, " << i << std::endl;
      return Z;
    }

    inline const Number_t& operator[](int i) const {
      switch (i) {
      case 0: return X;
      case 1: return Y;
      case 2: return Z;
      }
      std::cout << "Vector3::operator[]() [error] out of range, " << i << std::endl;
      return Z;
    }

    inline bool operator==(const Vector3& vector) const {
      return X == vector.X && Y == vector.Y && Z == vector.Z;
    }

    inline bool operator!=(const Vector3& vector) const {
      return !(*this == vector);
    }

    inline Vector3 operator+() const {
      return Vector3(+X, +Y, +Z);
    }

    inline Vector3 operator-() const {
      return Vector3(-X, -Y, -Z);
    }

    inline Vector3 operator+(const Vector3& vector) const {
      return Vector3(X + vector.X, Y + vector.Y, Z + vector.Z);
    }

    inline Vector3 operator-(const Vector3& vector) const {
      return Vector3(X - vector.X, Y - vector.Y, Z - vector.Z);
    }

    inline Vector3 operator*(const Number_t& number) const {
      return Vector3(X * number, Y * number, Z * number);
    }

    inline Vector3 operator/(const Number_t& number) const {
      return Vector3(X / number, Y / number, Z / number);
    }

    inline Vector3& operator=(const Vector3& vector) {
      X = vector.X;
      Y = vector.Y;
      Z = vector.Z;
      return *this;
    }
    
    inline Vector3& operator+=(const Vector3& vector) {
      X += vector.X;
      Y += vector.Y;
      Z += vector.Z;
      return *this;
    }

    inline Vector3& operator-=(const Vector3& vector) {
      X -= vector.X;
      Y -= vector.Y;
      Z -= vector.Z;
      return *this;
    }

    inline Vector3& operator*=(const Number_t& number) {
      X *= number;
      Y *= number;
      Z *= number;
      return *this;
    }

    inline Vector3& operator/=(const Number_t& number) {
      X /= number;
      Y /= number;
      Z /= number;
      return *this;
    }

    inline Number_t Length() const {
      // Ignore: Number_t is IMAGINARY NUMBER
      return X * X + Y * Y + Z * Z;
    }

    inline Vector3& RotateX(double theta) {
      const double c = std::cos(theta), s = std::sin(theta);
      const double y = Y * c - Z * s;
      const double z = Y * s + Z * c;
      Y = y;
      Z = z;
      return *this;
    }

    inline Vector3& RotateY(double theta) {
      const double c = std::cos(theta), s = std::sin(theta);
      const double z = Z * c - X * s;
      const double x = Z * s + X * c;
      Z = z;
      X = x;
      return *this;
    }

    inline Vector3& RotateZ(double theta) {
      const double c = std::cos(theta), s = std::sin(theta);
      const double x = X * c - Y * s;
      const double y = X * s + Y * c;
      X = x;
      Y = y;
      return *this;
    }

    inline Vector3& Rotate(double theta, const Vector3& axis) {
      const double c = std::cos(theta), s = std::sin(theta), cc = 1.0 - c;
      const double nx = axis.X, ny = axis.Y, nz = axis.Z;
      const double x = X * (nx * nx * cc +      c) + Y * (nx * ny * cc - nz * s) + Z * (nz * nx * cc - ny * s);
      const double y = X * (nx * ny * cc - nz * s) + Y * (ny * ny * cc +      c) + Z * (ny * nz * cc - nx * s);
      const double z = X * (nz * nx * cc - ny * s) + Y * (ny * nz * cc - nx * s) + Z * (nz * nz * cc +      c);
      X = x;
      Y = y;
      Z = z;
      return *this;
    }

    inline static Number_t InnerProduct(const Vector3& v1, const Vector3& v2) {
      return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
    }

    inline Number_t InnerProduct(const Vector3& vector) const {
      return InnerProduct(*this, vector);
    }
    
    inline static Vector3 OuterProduct(const Vector3& v1, const Vector3& v2) {
      return Vector3(v1.Y * v2.Z - v1.Z * v2.Y,
                     v1.Z * v2.X - v1.X * v2.Z,
                     v1.Z * v2.X - v1.X * v2.Z);
    }

    inline Vector3 OuterProduct(const Vector3& vector) const {
      return OuterProduct(*this, vector);
    }

  };

  template <typename Number_t>
  inline Vector3<Number_t> operator*(const Number_t& number, const Vector3<Number_t>& vector) {
    return vector * number;
  }

}

#endif
