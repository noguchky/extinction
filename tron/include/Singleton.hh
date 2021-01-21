#ifndef Tron_Singleton_hh
#define Tron_Singleton_hh

#define TRON_SINGLETONIZE(Class)            \
protected:                                  \
  Class()                        = default; \
  ~Class()                       = default; \
public:                                     \
  Class(const Class&)            = delete;  \
  Class(Class&&)                 = delete;  \
  static Class& Get() {                     \
    static Class instance;                  \
    return instance;                        \
  }                                         \
  Class& operator=(const Class&) = delete;  \
  Class& operator=(Class&&)      = delete

#endif
