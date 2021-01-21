#ifndef Tron_Nullable_hh
#define Tron_Nullable_hh

#include "Types.hh"

namespace Tron {

  // Value_t = int, char, etc...
  template <typename Value_t>
  class Nullable {
  private:
    Bool_t fIsNull;
    Value_t fValue;

  public:
    Nullable()
      : fIsNull(true), fValue() {
    }
    Nullable(std::nullptr_t)
      : fIsNull(true), fValue() {
    }
    Nullable(const Value_t& value)
      : fIsNull(false), fValue(value) {
    }
    Nullable(const Nullable& nullable)
      : fIsNull(nullable.fIsNull) {
      if (!fIsNull) {
        fValue = nullable.fValue;
      }
    }
    ~Nullable() {
    }

    Nullable& operator=(std::nullptr_t) {
      fIsNull = true;
      return *this;
    }
    Nullable& operator=(const Value_t& value) {
      fIsNull = false;
      fValue = value;
      return *this;
    }
    Nullable& operator=(const Nullable& nullable) {
      fIsNull = nullable.fIsNull;
      if (fIsNull) {
        fValue = nullable.fValue;
      }
    }

    explicit operator Bool_t() const {
      return !fIsNull;
    }

    const Value_t operator|(const Value_t& defaultValue) const { 
      return fIsNull ? defaultValue : fValue;
    }
    Nullable operator|(const Nullable& defaultValue) const {
      return fIsNull ? defaultValue : *this;
    }

    Nullable& operator|=(const Value_t& value) {
      if (fIsNull) {
        fValue = value;
      }
      return *this;
    }
  };

  // Value_t = int*, char*, etc...
  template <typename Value_t>
  class Nullable<Value_t*> {
  private:
    Value_t* fValue;

  public:
    Nullable()
      : fValue(nullptr) {
    }
    Nullable(std::nullptr_t)
      : fValue(nullptr) {
    }
    Nullable(Value_t* value)
      : fValue(value) {
    }
    Nullable(const Nullable& nullable)
      : fValue(nullable.fValue) {
    }
    ~Nullable() {
      fValue = nullptr;
    }

    Nullable& operator=(std::nullptr_t) {
      fValue = nullptr;
      return *this;
    }
    Nullable& operator=(Value_t* value) {
      fValue = value;
      return *this;
    }
    Nullable& operator=(Nullable& nullable) {
      fValue = nullable.fValue;
      return *this;
    }

    explicit operator Bool_t() const {
      return fValue;
    }

    Value_t* operator|(Value_t* defaultValue) { 
      return fValue ? fValue : defaultValue;
    }
    const Value_t* operator|(const Value_t* defaultValue) const { 
      return fValue ? fValue : defaultValue;
    }
    Nullable operator|(const Nullable& defaultValue) const {
      return fValue ? *this : defaultValue;
    }

    Nullable& operator|=(Value_t* value) {
      if (!fValue) {
        fValue = value;
      }
      return *this;
    }
  };

  // Value_t = const int*, const char*, etc...
  template <typename Value_t>
  class Nullable<const Value_t*> {
  private:
    const Value_t* fValue;

  public:
    Nullable()
      : fValue(nullptr) {
    }
    Nullable(std::nullptr_t)
      : fValue(nullptr) {
    }
    Nullable(const Value_t* value)
      : fValue(value) {
    }
    Nullable(const Nullable& nullable)
      : fValue(nullable.fValue) {
    }
    Nullable(const Nullable<Value_t*>& nullable)
      : fValue(nullable.fValue | nullptr) {
    }
    ~Nullable() {
      fValue = nullptr;
    }

    Nullable& operator=(std::nullptr_t) {
      fValue = nullptr;
      return *this;
    }
    Nullable& operator=(const Value_t* value) {
      fValue = value;
      return *this;
    }
    Nullable& operator=(Nullable& nullable) {
      fValue = nullable.fValue;
      return *this;
    }
    Nullable& operator=(const Nullable<Value_t*>& nullable) {
      fValue = nullable | nullptr;
      return *this;
    }

    explicit operator Bool_t() const {
      return fValue;
    }

    const Value_t* operator|(const Value_t* defaultValue) const { 
      return fValue ? fValue : defaultValue;
    }
    Nullable operator|(const Nullable& defaultValue) const {
      return fValue ? *this : defaultValue;
    }
    Nullable operator|(const Nullable<Value_t*>& defaultValue) const {
      return fValue ? *this : defaultValue;
    }

    Nullable& operator|=(const Value_t* value) {
      if (!fValue) {
        fValue = value;
      }
      return *this;
    }
  };

  template <typename Value_t>
  inline Nullable<Value_t> $(Value_t value) {
    return Nullable<Value_t>(value);
  }

}

#endif
