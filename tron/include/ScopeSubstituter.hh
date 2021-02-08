#ifndef Tron_ScopeSubstituter_hh
#define Tron_ScopeSubstituter_hh

namespace Tron {

  template <typename T>
  class ScopeSubstituter {
  private:
    T fOriginalValue;
    T& fValue;
  public:
    ScopeSubstituter(T& value, T subst) : fOriginalValue(value), fValue(value = subst) { }
    ~ScopeSubstituter() { fValue = fOriginalValue; }

    ScopeSubstituter(const ScopeSubstituter&)            = delete;
    ScopeSubstituter(ScopeSubstituter&&)                 = delete;
    ScopeSubstituter& operator=(const ScopeSubstituter&) = delete;
    ScopeSubstituter& operator=(ScopeSubstituter&&)      = delete;
  };
  
}

#endif
