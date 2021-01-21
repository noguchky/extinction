#ifndef Tron_Property_hh
#define Tron_Property_hh

#include <functional>
#include <type_traits>

namespace Tron {
  
  class Setter {
  private:
    struct Base {};
  public:
    struct Default : Base {};
    struct Custom  : Base {};
    struct Remove  : Base {};
    template <typename T>
    using IsSetter = std::is_base_of<Base, T>;
  };

  class Getter {
  private:
    struct Base {};
  public:
    struct Default : Base {};
    struct Custom  : Base {};
    struct Remove  : Base {};
    template <typename T>
    using IsGetter = std::is_base_of<Base, T>;
  };

  template <typename Value_t,
            typename Setter_t = Setter::Default,
            typename Getter_t = Getter::Default>
  struct Property {
    // Types
    template <typename T>
    struct IsValuePointer  { static constexpr bool value = std::is_pointer<Value_t>::value; };

    template <typename T>
    struct IsValidAccessor { static constexpr bool value = Setter::IsSetter<Setter_t>::value &&
                                                           Getter::IsGetter<Getter_t>::value; };
    template <typename F>
    struct IsValidSetter { static constexpr bool value =
        not std::is_same<F,std::nullptr_t>::value &&
        std::is_constructible<std::function<void(Value_t&,Value_t)>,F>::value; };
    template <typename F>
    struct IsValidGetter { static constexpr bool value =
        not std::is_same<F,std::nullptr_t>::value &&
        std::is_constructible<std::function<Value_t(Value_t&)>,F>::value; };

    template <typename T>
    struct IsSetterDefault { static constexpr bool value = std::is_same<Setter_t, Setter::Default>::value; };
    template <typename T>
    struct IsGetterDefault { static constexpr bool value = std::is_same<Getter_t, Getter::Default>::value; };

    template <typename T>
    struct IsSetterCustom  { static constexpr bool value = std::is_same<Setter_t, Setter::Custom >::value; };
    template <typename T>
    struct IsGetterCustom  { static constexpr bool value = std::is_same<Getter_t, Getter::Custom >::value; };

    template <typename T>
    struct IsSetterRemoved { static constexpr bool value = std::is_same<Setter_t, Setter::Remove >::value; };
    template <typename T>
    struct IsGetterRemoved { static constexpr bool value = std::is_same<Getter_t, Getter::Remove >::value; };

    template <bool V>
    using EnableIf = typename std::enable_if<V, std::nullptr_t>::type;

    template <typename T>
    using CanCastToBool = std::is_convertible<Value_t, bool>;

    // Member
    Value_t                                     fValue;
    const std::function<void(Value_t&,Value_t)> fSetter;
    const std::function<Value_t(Value_t&)>      fGetter;

    // Constructor
    template <typename T = void,
              EnableIf<not IsValidAccessor<T>::value> = nullptr>
    Property() = delete;

    template <typename T = void,
              EnableIf<IsSetterDefault<T>::value && IsGetterDefault<T>::value> = nullptr>
    Property() : fValue(), fSetter(nullptr), fGetter(nullptr) { }
    template <typename T = void,
              EnableIf<IsSetterDefault<T>::value && IsGetterDefault<T>::value> = nullptr>
    Property(const Value_t& value) : fValue(value), fSetter(nullptr), fGetter(nullptr) { }

    template <class F, typename T = void,
              EnableIf<IsSetterDefault<T>::value && IsGetterCustom<T>::value &&
                       IsValidGetter<F>::value> = nullptr>
    Property(F getter) : fSetter(nullptr), fGetter(getter) { }
    template <class F, typename T = void,
              EnableIf<IsSetterDefault<T>::value && IsGetterCustom<T>::value &&
                       IsValidGetter<F>::value> = nullptr>
    Property(const Value_t& value, F getter) : fValue(value), fSetter(nullptr), fGetter(getter) { }

    template <typename T = void,
              EnableIf<IsSetterDefault<T>::value && IsGetterRemoved<T>::value> = nullptr>
    Property() : fSetter(nullptr), fGetter(nullptr) { }
    template <typename T = void,
              EnableIf<IsSetterDefault<T>::value && IsGetterRemoved<T>::value> = nullptr>
    Property(const Value_t& value) : fValue(value), fSetter(nullptr), fGetter(nullptr) { }

    template <class F, typename T = void,
              EnableIf<IsSetterCustom<T>::value && IsGetterDefault<T>::value &&
                       IsValidSetter<F>::value> = nullptr>
    Property(F setter) : fSetter(setter), fGetter(nullptr) { }
    template <class F, typename T = void,
              EnableIf<IsSetterCustom<T>::value && IsGetterDefault<T>::value &&
                       IsValidSetter<F>::value> = nullptr>
    Property(const Value_t& value, F setter) : fSetter(setter), fGetter(nullptr) { *this = value; }

    template <class F1, class F2, typename T = void,
              EnableIf<IsSetterCustom<T>::value && IsGetterCustom<T>::value &&
                       IsValidSetter<F1>::value && IsValidGetter<F2>::value> = nullptr>
    Property(F1 setter, F2 getter) : fValue(), fSetter(setter), fGetter(getter) { }
    template <typename F1, typename F2, typename T = void,
              EnableIf<IsSetterCustom<T>::value && IsGetterCustom<T>::value &&
                       IsValidSetter<F1>::value && IsValidGetter<F2>::value> = nullptr>
    Property(const Value_t& value, F1 setter, F2 getter) : fSetter(setter), fGetter(getter) { *this = value; }

    template <class F, typename T = void,
              EnableIf<IsSetterCustom<T>::value && IsGetterRemoved<T>::value &&
                       IsValidSetter<F>::value> = nullptr>
    Property(F setter) : fValue(), fSetter(setter), fGetter(nullptr) { }
    template <class F, typename T = void,
              EnableIf<IsSetterCustom<T>::value && IsGetterRemoved<T>::value &&
                       IsValidSetter<F>::value> = nullptr>
    Property(const Value_t& value, F setter) : fValue(), fSetter(setter), fGetter(nullptr) { *this = value; }

    template <typename T = void,
              EnableIf<IsSetterRemoved<T>::value && IsGetterDefault<T>::value> = nullptr>
    Property(const Value_t& value) : fValue(value), fSetter(nullptr), fGetter(nullptr) { }

    template <class F, typename T = void,
              EnableIf<IsSetterRemoved<T>::value && IsGetterCustom<T>::value &&
                       IsValidGetter<F>::value> = nullptr>
    Property(const Value_t& value, F getter) : fValue(value), fSetter(nullptr), fGetter(getter) { }

    template <class F, typename T = void,
              EnableIf<IsSetterRemoved<T>::value && IsGetterRemoved<T>::value> = nullptr>
    Property() = delete;

    // Setter
    template <typename T = void, EnableIf<IsSetterDefault<T>::value> = nullptr>
    Property& operator=(Value_t value) {
      fValue = value;
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterCustom<T>::value> = nullptr>
    Property& operator=(Value_t value) {
      if (fSetter) { fSetter(fValue, value); } else { fValue = value; }
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterRemoved<T>::value> = nullptr>
    Property& operator=(Value_t value) = delete;

    template <typename T = void, EnableIf<IsSetterDefault<T>::value &&
                                          not IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(const Property& prop) {
      fValue = prop;
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterCustom<T>::value &&
                                          not IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(const Property& prop) {
      if (fSetter) { fSetter(fValue, prop ); } else { fValue = prop ; }
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterRemoved<T>::value ||
                                          IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(const Property& prop) = delete;

    template <typename T = void, EnableIf<IsSetterDefault<T>::value &&
                                          not IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(Property& prop) {
      fValue = prop;
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterCustom<T>::value &&
                                          not IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(Property& prop) {
      if (fSetter) { fSetter(fValue, prop ); } else { fValue = prop ; }
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterRemoved<T>::value ||
                                          IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(Property& prop) = delete;

    template <typename T = void, EnableIf<IsSetterDefault<T>::value &&
                                          not IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(Property&& prop) {
      fValue = prop;
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterCustom<T>::value &&
                                          not IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(Property&& prop) {
      if (fSetter) { fSetter(fValue, prop ); } else { fValue = prop ; }
      return *this;
    }
    template <typename T = void, EnableIf<IsSetterRemoved<T>::value ||
                                          IsGetterRemoved<T>::value> = nullptr>
    Property& operator=(Property&& prop) = delete;

    // Getter
    template <typename T = void, EnableIf<IsGetterDefault<T>::value &&
                                          IsValuePointer <T>::value> = nullptr>
    Value_t operator->() const {
      return fValue;
    }
    template <typename T = void, EnableIf<IsGetterCustom <T>::value &&
                                          IsValuePointer <T>::value> = nullptr>
    Value_t operator->() const {
      return fGetter ? fGetter(fValue) : fValue;
    }

    template <typename T = void, EnableIf<IsGetterDefault<T>::value> = nullptr>
    operator Value_t&() {
      return fValue;
    }
    template <typename T = void, EnableIf<IsGetterDefault<T>::value> = nullptr>
    operator const Value_t&() const {
      return fValue;
    }
    template <typename T = void, EnableIf<IsGetterCustom<T>::value> = nullptr>
    operator Value_t() const {
      return fGetter ? fGetter(fValue) : fValue;
    }

    template <typename T = void, EnableIf<IsGetterDefault<T>::value &&
                                          CanCastToBool  <T>::value> = nullptr>
    explicit operator bool() const {
      return fValue;
    }
    template <typename T = void, EnableIf<IsSetterCustom<T>::value &&
                                          CanCastToBool <T>::value> = nullptr>
    explicit operator bool() const {
      return fGetter ? fGetter(fValue) : fValue;
    }

  };

}

#endif
