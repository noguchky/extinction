#ifndef Tron_MultivariateNormalDistRand_hh
#define Tron_MultivariateNormalDistRand_hh

#include <random>
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "NormalDistRand.hh"
#include "SMatrix.hh"

namespace Tron {

  /// 多変量正規分布
  template <unsigned int D>
  class MultivariateNormalDistRand {
  private:
    NormalDistRand fNormalDistRand;
    SMatrix<double, D, 1> fMeans;
    SMatrix<double, D, D> fSigmas;
    TVectorD              fEigenValues;
    TMatrixD              fEigenVectors;
    MultivariateNormalDistRand<D - 1>* fMinor = nullptr;

  public:
    using ResultType = SMatrix<double, D, 1>;

  public:
    MultivariateNormalDistRand();
    MultivariateNormalDistRand(unsigned int seed);
    virtual ~MultivariateNormalDistRand() {
      delete fMinor; fMinor = nullptr;
    }

    inline void SetNormalDistRand(NormalDistRand rand) {
      fNormalDistRand = rand;
    }
    void SetParams(const SMatrix<double, D, 1>& means, const SMatrix<double, D, D>& sigmas);

    ResultType operator()();

  private:
    void UpdateParams();
  };

  template <>
  class MultivariateNormalDistRand<1> {
  private:
    NormalDistRand fNormalDistRand;
    SMatrix<double, 1, 1> fMeans;
    SMatrix<double, 1, 1> fSigmas;
    TVectorD              fEigenValues;
    TMatrixD              fEigenVectors;

  public:
    using ResultType = SMatrix<double, 1, 1>;

  public:
    MultivariateNormalDistRand();
    MultivariateNormalDistRand(unsigned int seed);
    virtual ~MultivariateNormalDistRand() {
    }

    inline void SetNormalDistRand(NormalDistRand rand) {
      fNormalDistRand = rand;
    }
    void SetParams(const SMatrix<double, 1, 1>& means, const SMatrix<double, 1, 1>& sigmas);

    ResultType operator()();

  private:
    void UpdateParams();
  };

  template <unsigned int D>
  inline MultivariateNormalDistRand<D>::MultivariateNormalDistRand() {
    std::random_device seed;
    fNormalDistRand = NormalDistRand(seed());
    fMeans  = SMatrix<double, D, 1>::ZeroMatrix();
    fSigmas = SMatrix<double, D, D>::IdentityMatrix();
    UpdateParams();
  }

  template <unsigned int D>
  inline MultivariateNormalDistRand<D>::MultivariateNormalDistRand(unsigned int seed) {
    fNormalDistRand = NormalDistRand(seed);
    fMeans  = SMatrix<double, D, 1>::ZeroMatrix();
    fSigmas = SMatrix<double, D, D>::IdentityMatrix();
    UpdateParams();
  }

  inline MultivariateNormalDistRand<1>::MultivariateNormalDistRand() {
    std::random_device seed;
    fNormalDistRand = NormalDistRand(seed());
    fMeans  = SMatrix<double, 1, 1>::ZeroMatrix();
    fSigmas = SMatrix<double, 1, 1>::IdentityMatrix();
    UpdateParams();
  }

  inline MultivariateNormalDistRand<1>::MultivariateNormalDistRand(unsigned int seed) {
    fNormalDistRand = NormalDistRand(seed);
    fMeans  = SMatrix<double, 1, 1>::ZeroMatrix();
    fSigmas = SMatrix<double, 1, 1>::IdentityMatrix();
    UpdateParams();
  }

  template <unsigned int D>
  inline void MultivariateNormalDistRand<D>::SetParams(const SMatrix<double, D, 1>& means, const SMatrix<double, D, D>& sigmas) {
    fMeans  = means;
    fSigmas = sigmas;

    // Check triangular or symmetry
    const auto tsigmas = fSigmas.Transpose();
    if (fSigmas != tsigmas) {
      // Convert triangular to symmetric
      fSigmas += tsigmas;
      for (unsigned int i = 0; i < D; ++i) {
        fSigmas(i, i) /= 2.0;
      }
    }

    UpdateParams();

    delete fMinor; fMinor = nullptr;
    for (unsigned int i = 0; i < D; ++i) {
      if (fSigmas(i, i) == 0) {
        fMinor = new MultivariateNormalDistRand<D - 1>();
        std::vector<double> means2  = fMeans.MinorRow(i);
        std::vector<double> sigmas2 = fSigmas.MinorRowCol(i, i);
        fMinor->SetNormalDistRand(fNormalDistRand);
        fMinor->SetParams(means2, sigmas2);
        break;
      }
    }

    if (!fMinor) {
      UpdateParams();
    }
  }

  inline void MultivariateNormalDistRand<1>::SetParams(const SMatrix<double, 1, 1>& means, const SMatrix<double, 1, 1>& sigmas) {
    fMeans  = means;
    fSigmas = sigmas;

    // Check triangular or symmetry
    const auto tsigmas = fSigmas.Transpose();
    if (fSigmas != tsigmas) {
      // Convert triangular to symmetric
      fSigmas += tsigmas;
      for (unsigned int i = 0; i < 1; ++i) {
        fSigmas(i, i) /= 2.0;
      }
    }

    UpdateParams();
  }

  template <unsigned int D>
  inline typename MultivariateNormalDistRand<D>::ResultType MultivariateNormalDistRand<D>::operator()() {
    if (!fMinor) {
      TVectorD rotatedResult(D);
      for (unsigned int i = 0; i < D; ++i) {
        const double variance = fEigenValues[i];
        if (variance > 0) {
          rotatedResult[i] = TMath::Sqrt(variance) * fNormalDistRand();
        }
      }

      TVectorD result = fEigenVectors * rotatedResult;
      result.Add({ D, fMeans.begin() });

      return result.GetMatrixArray();
    } else {
      SMatrix<double, D - 1, 1> minorResult = (*(MultivariateNormalDistRand<D - 1>*)fMinor)();
      SMatrix<double, D, 1> result;
      for (UInt_t i = 0; i < D; ++i) {
        if      (i < D) { result(i, 0) = minorResult(i    , 0); }
        else if (i > D) { result(i, 0) = minorResult(i - 1, 0); }
        else            { result(i, 0) = fMeans(i, 0); }
      }
      return result;
    }
  }

  inline typename MultivariateNormalDistRand<1>::ResultType MultivariateNormalDistRand<1>::operator()() {
    if (fSigmas(0, 0) != 0.0) {
      TVectorD rotatedResult(1);
      for (unsigned int i = 0; i < 1; ++i) {
        const double variance = fEigenValues[i];
        if (variance > 0) {
          rotatedResult[i] = TMath::Sqrt(variance) * fNormalDistRand();
        }
      }

      TVectorD result = fEigenVectors * rotatedResult;
      result.Add({ 1, fMeans.begin() });

      return result.GetMatrixArray();
    } else {
      return { fMeans(0, 0) };
    }
  }

  template <unsigned int D>
  inline void MultivariateNormalDistRand<D>::UpdateParams() {
    auto eigenGenerator = TMatrixDSymEigen(TMatrixDSym(D, fSigmas.begin()));

    fEigenVectors.ResizeTo(D, D);
    fEigenVectors = eigenGenerator.GetEigenVectors();

    fEigenValues.ResizeTo(D);
    fEigenValues = eigenGenerator.GetEigenValues();

    // Check Variance
    for (unsigned int i = 0; i < D; ++i) {
      const double variance = fEigenValues[i];
      if (variance < 0) {
        std::cout << "[warning] variance of sigmas < 0" << std::endl;
        break;
      }
    }
  }
  
  inline void MultivariateNormalDistRand<1>::UpdateParams() {
    auto eigenGenerator = TMatrixDSymEigen(TMatrixDSym(1, fSigmas.begin()));

    fEigenVectors.ResizeTo(1, 1);
    fEigenVectors = eigenGenerator.GetEigenVectors();

    fEigenValues.ResizeTo(1);
    fEigenValues = eigenGenerator.GetEigenValues();

    // Check Variance
    for (unsigned int i = 0; i < 1; ++i) {
      const double variance = fEigenValues[i];
      if (variance < 0) {
        std::cout << "[warning] variance of sigmas < 0" << std::endl;
        break;
      }
    }

  }

}

#endif
