#include <csignal>
#include <iostream>
#include "SignalListener.hh"

bool                                        Tron::SignalListener::gInterrupted = false;
bool                                        Tron::SignalListener::gEnabled     = false;
std::map<int, Tron::SignalListener::Func_t> Tron::SignalListener::gInterruptedFuncs;

void Tron::SignalListener::SetEnabled(bool enabled) {
  if (gEnabled != enabled) {
    gInterrupted = false;
    if (enabled) {
      std::signal(SIGINT, SignalInterrupted);
      gEnabled = true;
    } else {
      std::signal(SIGINT, SIG_DFL);
      gEnabled = false;
    }
  }
}

bool Tron::SignalListener::IsEnabled() {
  return gEnabled;
}

bool Tron::SignalListener::IsInterrupted() {
  return gInterrupted;
}

void Tron::SignalListener::Add(int id, const Func_t& func) {
  gInterruptedFuncs[id] = func;
}

void Tron::SignalListener::Remove(int id) {
  gInterruptedFuncs.erase(id);
}

void Tron::SignalListener::SignalInterrupted(int signal) {
  if (gEnabled) {
    std::cout << "Signal Interrupted" << std::endl;
    gInterrupted = true;
    
    for (auto&& idfunc : gInterruptedFuncs) {
      idfunc.second(signal);
    }
  }
}
