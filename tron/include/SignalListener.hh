#ifndef Tron_SignalListener_hh
#define Tron_SignalListener_hh

#include <functional>
#include <map>
#include <iostream>

namespace Tron {

  class SignalListener {
  public:
    using Func_t  = std::function<void(int)>;

  private:
    static bool                  gEnabled;
    static bool                  gInterrupted;
    static std::map<int, Func_t> gInterruptedFuncs;
    SignalListener() = delete;
    
  public:
    static void SetEnabled(bool enabled);
    static bool IsEnabled();
    static bool IsInterrupted();
    static void Add(int id, const Func_t& func);
    static void Remove(int id);
    
  private:
    static void SignalInterrupted(int);
    
  };

}

#endif
