#ifndef Tron_InternalEnumerator_hh
#define Tron_InternalEnumerator_hh

#include <iostream>
#include <algorithm>
#include <functional>
#include <array>
#include <vector>
#include <set>
#include "Types.hh"
#include "Math.hh"
#include "Vector2.hh"

namespace Tron {

  template <UInt_t D>
  class InternalEnumerator {
  public:
    using IPosition = std::array<Int_t, D>;
    using DPosition = std::array<Double_t, D>;
    using Func_t    = std::function<bool(const DPosition&)>;

  private:
    Func_t              fImplicitFunc;
    DPosition           fCenter;
    DPosition           fDifference;
    std::set<IPosition> fInternalGrids;
    std::set<IPosition> fNewFoundGrids;
    UInt_t              fMaxN = 10000;

  public:
    InternalEnumerator(const Func_t& func);

    void SetCenter(const std::initializer_list<Double_t>& list);
    void SetDifference(const std::initializer_list<Double_t>& list);
    void SetMaxN(UInt_t maxN);

    std::vector<DPosition> Enumerate();

  private:
    DPosition GetPosition(const IPosition& grid);
    void      FindNew(IPosition grid);

  };

  namespace CircularEnumeration {
    const typename InternalEnumerator<2>::Func_t Func =
      [](const InternalEnumerator<2>::DPosition& p) -> bool {
        const Double_t& radius = p[0];
        const Double_t& arc    = p[1];
        return
          radius == 0.0 ?
          arc    == 0.0 :
          (0.0 <= radius       && radius       <= 1.0 &&
           0.0 <= arc / radius && arc / radius <= 1.0);
      };
    
    struct Result_t : public Vector2<Double_t> {
      Double_t& R     = X;
      Double_t& Theta = Y;
      Result_t(Double_t r, Double_t theta) : Vector2<Double_t>(r, theta) {
      }
    };
    
    Result_t Convert(const InternalEnumerator<2>::DPosition& p) {
      const Double_t& radius = p[0];
      const Double_t& arc    = p[1];
      return { radius, radius ? arc / radius * twopi : 0.0 };
    }
  }
}

template <UInt_t D>
Tron::InternalEnumerator<D>::InternalEnumerator(const Func_t& func)
  : fImplicitFunc(func) {
  std::fill(fCenter.begin(), fCenter.end(), 0.0);
  std::fill(fDifference.begin(), fDifference.end(), 1.0);
}

template <UInt_t D>
void Tron::InternalEnumerator<D>::SetCenter(const std::initializer_list<Double_t>& list) {
  if (list.size() != D) {
    std::cout << "InternalEnumerator::SetCenter() [error] position dimension mismatch" << std::endl;
    return;
  }
  std::copy(list.begin(), list.end(), fCenter.begin());
}

template <UInt_t D>
void Tron::InternalEnumerator<D>::SetDifference(const std::initializer_list<Double_t>& list) {
  if (list.size() != D) {
    std::cout << "InternalEnumerator::SetDifference() [error] position dimension mismatch" << std::endl;
    return;
  }
  std::copy(list.begin(), list.end(), fDifference.begin());
}

template <UInt_t D>
void Tron::InternalEnumerator<D>::SetMaxN(UInt_t maxN) {
  fMaxN = maxN;
}

template <UInt_t D>
std::vector<typename Tron::InternalEnumerator<D>::DPosition> Tron::InternalEnumerator<D>::Enumerate() {
  fInternalGrids.clear();
  fNewFoundGrids.clear();

  // Create initial grid
  IPosition initialGrid;
  std::fill(initialGrid.begin(), initialGrid.end(), 0);

  // Check initial grid is internal
  DPosition initialPosition = GetPosition(initialGrid);
  if (fImplicitFunc(initialPosition)) {
    fInternalGrids.insert(initialGrid);

    // Find internal grids around initial grid
    FindNew(initialGrid);

    while (fInternalGrids.size() <= fMaxN && fNewFoundGrids.size() != 0) {
      // Clone next grids
      std::set<IPosition> nextGrids(fNewFoundGrids);
      fInternalGrids.insert(fNewFoundGrids.begin(), fNewFoundGrids.end());
      fNewFoundGrids.clear();

      // Find internal grids around each next grids
      for (auto&& nextGrid : nextGrids) {
        FindNew(nextGrid);
      }
    }
    if (fInternalGrids.size() > fMaxN) {
      std::cout << "InternalEnumerator::Enumerate [warning] # of data reaches maxN" << std::endl;
    }
  }

  // Create Results
  std::vector<DPosition> results;
  for (auto&& grid : fInternalGrids) {
    results.push_back(GetPosition(grid));
  }

  return results;
}

template <UInt_t D>
typename Tron::InternalEnumerator<D>::DPosition Tron::InternalEnumerator<D>::GetPosition(const IPosition& grid) {
  DPosition position;
  for (UInt_t i = 0; i < D; ++i) {
    position[i] = fCenter[i] + fDifference[i] * grid[i];
  }
  return position;
}

template <UInt_t D>
void Tron::InternalEnumerator<D>::FindNew(IPosition grid) {
  // Check sorrounding grid
  for (Int_t sign = -1; sign <= 1; sign += 2) {
    for (UInt_t i = 0; i < D; ++i) {
      grid[i] += sign;
  
      // Check douplicate
      if (fInternalGrids.find(grid) == fInternalGrids.end() &&
          fNewFoundGrids.find(grid) == fNewFoundGrids.end()) {
        // Check internal
        DPosition position = GetPosition(grid);
        if (fImplicitFunc(position)) {
          fNewFoundGrids.insert(grid);
        }
      }

      grid[i] -= sign;
    }
  }
}

#endif
