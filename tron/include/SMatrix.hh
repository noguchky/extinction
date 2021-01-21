#ifndef Tron_SMatrix_hh
#define Tron_SMatrix_hh

#include "Math/SMatrix.h"
#include "Math/Dfact.h"
#include "VMatrix.hh"

namespace Tron {

  template <typename Number_t, UInt_t D1, UInt_t D2 = D1>
  class SMatrix {
  public:
    static const UInt_t kRows = D1;
    static const UInt_t kCols = D2;
    static const UInt_t kSize = D1 * D2;

    using value_type     = Number_t;
    using iterator       = Number_t*;
    using const_iterator = const Number_t*;

  protected:
    Number_t* fValues;
    
  public:
    inline Number_t& operator[](Int_t index) {
      return fValues[index];
    }
    
    inline const Number_t& operator[](Int_t index) const {
      return fValues[index];
    }
    
    inline Number_t& operator()(Int_t row, Int_t col) {
      return fValues[row * kCols + col];
    }

    inline const Number_t& operator()(Int_t row, Int_t col) const {
      return fValues[row * kCols + col];
    }

    inline iterator begin() {
      return fValues;
    }

    inline iterator end() {
      return fValues + kSize;
    }

    inline const_iterator begin() const {
      return fValues;
    }

    inline const_iterator end() const {
      return fValues + kSize;
    }

    SMatrix() {
      fValues = new Number_t[kSize] { };
    }

    SMatrix(const Number_t* values) {
      fValues = new Number_t[kSize] { };
      auto it = values;
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    SMatrix(const std::initializer_list<Number_t>& values) {
      fValues = new Number_t[kSize] { };
      auto it = values.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    SMatrix(const std::vector<Number_t>& values) {
      fValues = new Number_t[kSize] { };
      auto it = values.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    SMatrix(const SMatrix& matrix) {
      fValues = new Number_t[kSize] { };
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    virtual ~SMatrix() {
      delete fValues;
      fValues = nullptr;
    }

    inline iterator Rows(Int_t i) {
      return fValues + i * kCols;
    }

    inline const_iterator Rows(Int_t i) const {
      return fValues + i * kCols;
    }

    inline Bool_t operator==(const SMatrix& matrix) const {
      auto it = matrix.begin();
      for (auto&& value : *this) {
        if (value != *it++) {
          return false;
        }
      }
      return true;
    }

    inline Bool_t operator!=(const SMatrix& matrix) const {
      return !(*this == matrix);
    }

    inline SMatrix operator+() const {
      auto newOne = SMatrix();
      auto it = begin();
      for (auto&& value : newOne) {
        value = *it++;
      }
      return newOne;
    }

    inline SMatrix operator-() const {
      auto newOne = SMatrix();
      auto it = begin();
      for (auto&& value : newOne) {
        value = - *it++;
      }
      return newOne;
    }

    inline SMatrix operator+(const SMatrix& matrix) const {
      auto newOne = SMatrix();
      auto it1 = begin();
      auto it2 = matrix.begin();
      for (auto&& value : newOne) {
        value = *it1++ + *it2++;
      }
      return newOne;
    }

    inline SMatrix operator-(const SMatrix& matrix) const {
      auto newOne = SMatrix();
      auto it1 = begin();
      auto it2 = matrix.begin();
      for (auto&& value : newOne) {
        value = *it1++ - *it2++;
      }
      return newOne;
    }

    inline SMatrix operator*(const Number_t& number) const {
      auto newOne = SMatrix();
      auto it1 = begin();
      for (auto&& value : newOne) {
        value = *it1++ * number;
      }
      return newOne;
    }

    inline SMatrix operator/(const Number_t& number) const {
      auto newOne = SMatrix();
      auto it1 = begin();
      for (auto&& value : newOne) {
        value = *it1++ / number;
      }
      return newOne;
    }

    template <UInt_t D3>
    inline SMatrix<Number_t, D1, D3> operator*(const SMatrix<Number_t, D2, D3>& matrix) const {
      auto newOne = SMatrix<Number_t, D1, D3>();
      for (UInt_t newRow = 0; newRow < D1; ++newRow) {
        for (UInt_t newCol = 0; newCol < D3; ++newCol) {
          for (UInt_t col = 0; col < D2; ++col) {
            newOne(newRow, newCol) += (*this)(newRow, col) * matrix(col, newCol);
          }
        }
      }
      return newOne;
    }

    inline SMatrix& operator=(const SMatrix& matrix) {
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
      return *this;
    }

    inline SMatrix& operator+=(const SMatrix& matrix) {
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value += *it++;
      }
      return *this;
    }

    inline SMatrix& operator-=(const SMatrix& matrix) {
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value -= *it++;
      }
      return *this;
    }

    template <typename ONumber_t>
    inline SMatrix& operator*=(const ONumber_t& number) {
      for (auto&& value : *this) {
        value *= number;
      }
      return *this;
    }

    template <typename ONumber_t>
    inline SMatrix& operator/=(const ONumber_t& number) {
      for (auto&& value : *this) {
        value /= number;
      }
      return *this;
    }

    template <typename ONumber_t>
    inline SMatrix& operator*=(const SMatrix<ONumber_t, D2, D2>& matrix) {
      auto newOne = (*this) * matrix;
      return (*this) = newOne;
    }

    inline operator VMatrix<Number_t>() const {
      return VMatrix<Number_t>(kRows, kCols, fValues);
    }

    inline static SMatrix ZeroMatrix() {
      return SMatrix();
    }

    inline static SMatrix<Number_t, D1> IdentityMatrix() {
      auto newOne = SMatrix<Number_t, D1>();
      for (UInt_t i = 0; i < D1; ++i) {
        newOne(i, i) = 1;
      }
      return newOne;
    }

    inline static SMatrix<Number_t, D2, D1> Transpose(const SMatrix& matrix) {
      auto newOne = SMatrix<Number_t, D2, D1>();
      auto it = matrix.begin();
      for (UInt_t row = 0; row < D1; ++row) {
        for (UInt_t col = 0; col < D2; ++col) {
          newOne(col, row) = *it++;
        }
      }
      return newOne;
    }

    inline SMatrix<Number_t, D2, D1> Transpose() const {
      return Transpose(*this);
    }

    inline SMatrix& TransposeThis() {
      return *this = Transpose(*this);
    }

    inline static SMatrix<Number_t, D1> Inverse(const SMatrix<Number_t, D1>& matrix) {
      auto smatrix = ROOT::Math::SMatrix<Number_t, D1>(matrix.begin(), matrix.end());
      smatrix.Invert();
      return SMatrix<Number_t, D1>(smatrix.Array());
    }

    inline SMatrix<Number_t, D1> Inverse() const {
      return Inverse(*this);
    }

    inline SMatrix& InverseThis() {
      return *this = Inverse(*this);
    }

    inline static Number_t Trace(const SMatrix<Number_t, D1>& matrix) {
      Number_t sum = Number_t();
      for (UInt_t i = 0; i < D1; ++i) {
        sum += matrix(i, i);
      }
      return sum;
    }

    inline Number_t Trace() const {
      return Trace(*this);
    }

    inline static Number_t Determinant(const SMatrix<Number_t, D1>& matrix) {
      auto smatrix = ROOT::Math::SMatrix<Number_t, D1>(matrix.begin(), matrix.end());
      Number_t det = Number_t();
      smatrix.Det(det);
      return det;
    }

    inline Number_t Determinant() const {
      return Determinant(*this);
    }

    inline std::vector<Number_t> MinorRowCol(UInt_t row, UInt_t col) const {
      std::vector<Number_t> newOne;
      if (D1 > 1 && D2 > 1) {
        newOne.assign((D1 - 1) * (D2 - 1), 0.0);
        for (UInt_t i = 0; i < D1; ++i) {
          for (UInt_t j = 0; j < D2; ++j) {
            if      (i < row     && j < col    ) { newOne[(i - 0) * (D2 - 1) + (j - 0)] = (*this)(i, j); }
            else if (i < row     &&     col < j) { newOne[(i - 0) * (D2 - 1) + (j - 1)] = (*this)(i, j); }
            else if (    row < i && j < col    ) { newOne[(i - 1) * (D2 - 1) + (j - 0)] = (*this)(i, j); }
            else if (    row < i &&     col < j) { newOne[(i - 1) * (D2 - 1) + (j - 1)] = (*this)(i, j); }
          }
        }
      }
      return newOne;
    }

    inline std::vector<Number_t> MinorRow(UInt_t row) const {
      std::vector<Number_t> newOne;
      if (D1 > 1) {
        newOne.assign((D1 - 1) * D2, 0.0);
        for (UInt_t i = 0; i < D1; ++i) {
          for (UInt_t j = 0; j < D2; ++j) {
            if      (i < row    ) { newOne[(i - 0) * D2 + j] = (*this)(i, j); }
            else if (    row < i) { newOne[(i - 1) * D2 + j] = (*this)(i, j); }
          }
        }
      }
      return newOne;
    }

    inline std::vector<Number_t> MinorCol(UInt_t col) const {
      std::vector<Number_t> newOne;
      if (D2 > 1) {
        newOne.assign(D1 * (D2 - 1));
        for (UInt_t i = 0; i < D1; ++i) {
          for (UInt_t j = 0; j < D2; ++j) {
            if      (j < col    ) { newOne[i * (D2 - 1) + (j - 0)] = (*this)(i, j); }
            else if (    col < j) { newOne[i * (D2 - 1) + (j - 1)] = (*this)(i, j); }
          }
        }
      }
      return newOne;
    }

  };
  
  template <typename Number_t, UInt_t D1, UInt_t D2>
  inline SMatrix<Number_t, D1, D2> operator*(const Number_t& number, const SMatrix<Number_t, D1, D2>& matrix) {
    return matrix * number;
  }

}

#endif
