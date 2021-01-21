#ifndef Tron_VMatrix_hh
#define Tron_VMatrix_hh

#include <iostream>
#include <vector>
#include "TMatrix.h"

namespace Tron {

  template <typename Number_t>
  class VMatrix {
  public:
    using value_type     = typename std::vector<Number_t>::value_type;
    using iterator       = typename std::vector<Number_t>::iterator;
    using const_iterator = typename std::vector<Number_t>::const_iterator;

  protected:
    UInt_t fRows;
    UInt_t fCols;
    UInt_t fSize;
    std::vector<Number_t> fValues;
    
  public:
    inline Number_t& operator[](Int_t index) {
      return fValues[index];
    }
    
    inline const Number_t& operator[](Int_t index) const {
      return fValues[index];
    }
    
    inline Number_t& operator()(Int_t row, Int_t col) {
      return fValues[row * fCols + col];
    }

    inline const Number_t& operator()(Int_t row, Int_t col) const {
      return fValues[row * fCols + col];
    }

    inline iterator begin() {
      return fValues.begin();
    }

    inline iterator end() {
      return fValues.end();
    }

    inline const_iterator begin() const {
      return fValues.begin();
    }

    inline const_iterator end() const {
      return fValues.end();
    }

    VMatrix()
      : fRows(1),
        fCols(1),
        fSize(fRows * fCols),
        fValues(fSize) {
    }

    VMatrix(int rows, int cols)
     : fRows(rows),
       fCols(cols),
       fSize(fRows * fCols),
       fValues(fSize) {
    }

    VMatrix(int rows)
      : fRows(rows),
        fCols(rows),
        fSize(fRows * fCols),
        fValues(fSize) {
    }

    VMatrix(int rows, int cols, const Number_t* values)
      : fRows(rows),
        fCols(cols),
        fSize(fRows * fCols),
        fValues(fSize) {
      auto it = values;
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    VMatrix(int rows, const Number_t* values)
      : fRows(rows),
        fCols(rows),
        fSize(fRows * fCols),
        fValues(fSize) {
      auto it = values;
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    VMatrix(int rows, const std::initializer_list<Number_t>& values)
      : fRows(rows),
        fCols(rows),
        fSize(fRows * fCols),
        fValues(fSize) {
      auto it = values.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    VMatrix(int rows, int cols, const std::vector<Number_t>& values)
      : fRows(rows),
        fCols(cols),
        fSize(fRows * fCols),
        fValues(fSize) {
      auto it = values.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    VMatrix(int rows, const std::vector<Number_t>& values)
      : fRows(rows),
        fCols(rows),
        fSize(fRows * fCols),
        fValues(fSize) {
      auto it = values.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    VMatrix(const VMatrix& matrix)
      : fRows(matrix.fRows),
        fCols(matrix.fCols),
        fSize(matrix.fSize),
        fValues(fSize) {
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
    }

    virtual ~VMatrix() {
    }
    
    inline iterator Rows(Int_t i) {
      return fValues + i * fCols;
    }

    inline const_iterator Rows(Int_t i) const {
      return fValues + i * fCols;
    }

    inline Bool_t operator==(const VMatrix& matrix) const {
      if (this->fCols != matrix.fCols || this->fRows != matrix.fRows) {
        std::cout << "VMatrix::operator==() [error] matrix dimension mismatch" << std::endl;
        return false;
      }
      auto it = matrix.begin();
      for (auto&& value : *this) {
        if (value != *it++) {
          return false;
        }
      }
      return true;
    }

    inline Bool_t operator!=(const VMatrix& matrix) const {
      return !(*this == matrix);
    }

    inline VMatrix operator+() const {
      auto newOne = VMatrix();
      auto it = begin();
      for (auto&& value : newOne) {
        value = *it++;
      }
      return newOne;
    }

    inline VMatrix operator-() const {
      auto newOne = VMatrix();
      auto it = begin();
      for (auto&& value : newOne) {
        value = - *it++;
      }
      return newOne;
    }

    inline VMatrix operator+(const VMatrix& matrix) const {
      if (this->fCols != matrix.fCols || this->fRows != matrix.fRows) {
        std::cout << "VMatrix::operator+() [error] matrix dimension mismatch" << std::endl;
        return *this;
      }
      auto newOne = VMatrix();
      auto it1 = begin();
      auto it2 = matrix.begin();
      for (auto&& value : newOne) {
        value = *it1++ + *it2++;
      }
      return newOne;
    }

    inline VMatrix operator-(const VMatrix& matrix) const {
      if (this->fCols != matrix.fCols || this->fRows != matrix.fRows) {
        std::cout << "VMatrix::operator-() [error] matrix dimension mismatch" << std::endl;
        return *this;
      }
      auto newOne = VMatrix();
      auto it1 = begin();
      auto it2 = matrix.begin();
      for (auto&& value : newOne) {
        value = *it1++ - *it2++;
      }
      return newOne;
    }

    inline VMatrix operator*(const Number_t& number) const {
      auto newOne = VMatrix();
      auto it1 = begin();
      for (auto&& value : newOne) {
        value = *it1++ * number;
      }
      return newOne;
    }

    inline VMatrix operator/(const Number_t& number) const {
      auto newOne = VMatrix();
      auto it1 = begin();
      for (auto&& value : newOne) {
        value = *it1++ / number;
      }
      return newOne;
    }

    inline VMatrix<Number_t> operator*(const VMatrix<Number_t>& matrix) const {
      if (this->fCols != matrix.fRows) {
        std::cout << "VMatrix::operator*() [error] matrix dimension mismatch" << std::endl;
        return *this;
      }
      auto newOne = VMatrix<Number_t>(this->fRows, matrix.fCols);
      for (UInt_t newRow = 0; newRow < this->fRows; ++newRow) {
        for (UInt_t newCol = 0; newCol < matrix.fCols; ++newCol) {
          for (UInt_t col = 0; col < this->fCols; ++col) {
            newOne(newRow, newCol) += (*this)(newRow, col) * matrix(col, newCol);
          }
        }
      }
      return newOne;
    }

    inline VMatrix& operator=(const VMatrix& matrix) {
      if (this->fCols != matrix.fCols || this->fRows != matrix.fRows) {
        std::cout << "VMatrix::operator=() [error] matrix dimension mismatch" << std::endl;
        return *this;
      }
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value = *it++;
      }
      return *this;
    }

    inline VMatrix& operator+=(const VMatrix& matrix) {
      if (this->fCols != matrix.fCols || this->fRows != matrix.fRows) {
        std::cout << "VMatrix::operator+=() [error] matrix dimension mismatch" << std::endl;
        return *this;
      }
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value += *it++;
      }
      return *this;
    }

    inline VMatrix& operator-=(const VMatrix& matrix) {
      if (this->fCols != matrix.fCols || this->fRows != matrix.fRows) {
        std::cout << "VMatrix::operator-=() [error] matrix dimension mismatch" << std::endl;
        return *this;
      }
      auto it = matrix.begin();
      for (auto&& value : *this) {
        value -= *it++;
      }
      return *this;
    }

    template <typename ONumber_t>
    inline VMatrix& operator*=(const ONumber_t& number) {
      for (auto&& value : *this) {
        value *= number;
      }
      return *this;
    }

    template <typename ONumber_t>
    inline VMatrix& operator/=(const ONumber_t& number) {
      for (auto&& value : *this) {
        value /= number;
      }
      return *this;
    }

    template <typename ONumber_t>
    inline VMatrix& operator*=(const VMatrix<ONumber_t>& matrix) {
      auto newOne = (*this) * matrix;
      return (*this) = newOne;
    }

    inline static VMatrix ZeroVMatrix(UInt_t rows, UInt_t cols) {
      return VMatrix(rows, cols);
    }

    inline static VMatrix<Number_t> IdentityMatrix(UInt_t rows) {
      auto newOne = VMatrix<Number_t>(rows);
      for (UInt_t i = 0; i < rows; ++i) {
        newOne(i, i) = 1;
      }
      return newOne;
    }

    inline static VMatrix<Number_t> Transpose(const VMatrix& matrix) {
      auto newOne = VMatrix<Number_t>(matrix.fCols, matrix.fRows);
      auto it = matrix.begin();
      for (UInt_t row = 0; row < matrix.fRows; ++row) {
        for (UInt_t col = 0; col < matrix.fCols; ++col) {
          newOne(col, row) = *it++;
        }
      }
      return newOne;
    }

    inline VMatrix<Number_t> Transpose() const {
      return Transpose(*this);
    }

    inline VMatrix& TransposeThis() {
      return *this = Transpose(*this);
    }

    inline static VMatrix<Number_t> Inverse(const VMatrix<Number_t>& matrix) {
      if (matrix.fRows != matrix.fCols) {
        std::cout << "VMatrix::Inverse() [error] matrix dimension mismatch" << std::endl;
        return matrix;
      }
      auto tmatrix = TMatrixT<Number_t>(matrix.fRows, matrix.fCols, matrix.fValues.data());
      tmatrix.Invert();
      VMatrix<Number_t> newOne(matrix.fRows);
      for (UInt_t row = 0; row < matrix.fRows; ++row) {
        for (UInt_t col = 0; col < matrix.fCols; ++col) {
          newOne(row, col) = tmatrix(row, col);
        }
      }
      return newOne;
    }

    inline VMatrix<Number_t> Inverse() const {
      return Inverse(*this);
    }

    inline VMatrix& InverseThis() {
      return *this = Inverse(*this);
    }

    inline static Number_t Trace(const VMatrix<Number_t>& matrix) {
      if (matrix.fRows != matrix.fCols) {
        std::cout << "VMatrix::Trace() [error] matrix dimension mismatch" << std::endl;
        return Number_t();
      }
      Number_t sum = Number_t();
      for (UInt_t i = 0; i < matrix.fRows; ++i) {
        sum += matrix(i, i);
      }
      return sum;
    }

    inline Number_t Trace() const {
      return Trace(*this);
    }

    inline static Number_t Determinant(const VMatrix<Number_t>& matrix) {
      auto tmatrix = TMatrixT<Number_t>(matrix.fRows, matrix.fCols, matrix.fValues.data());
      return tmatrix.Determinant();
    }

    inline Number_t Determinant() const {
      return Determinant(*this);
    }

  };
  
  template <typename Number_t>
  inline VMatrix<Number_t> operator*(const Number_t& number, const VMatrix<Number_t>& matrix) {
    return matrix * number;
  }
    
}

#endif
