#ifndef Tron_Console_hh
#define Tron_Console_hh

#include <iostream>
#include <sstream>
#include "TFormula.h"
#include "Types.hh"

namespace Tron {
  
  namespace Console {

    std::string Prompt(const std::string& appName, const std::string& funcName);
    std::string Prompt(const std::string& appName, const std::string& hostName, const std::string& funcName);

    inline void Clear() {
      std::cin.clear();
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(),' ');
    }
    
    template <typename Value_t>
    Bool_t ReadTerm(Value_t& value);

    Bool_t ReadLine(std::string& value);
    
    template <typename Value_t>
    Bool_t ReadLine(Value_t& value);

    template <typename Value_t, typename Conv_t/*Int_t(std::string& line, Value_t& value)*/>
    Bool_t ReadLine(Value_t& value, Conv_t converter);

    Bool_t Confirm(const std::string& message = "OK?", Bool_t defaultAnswer = false);
    
    void   Wait(const std::string& message = "");

    Bool_t KeyboardHit();

  }

}

template <typename Value_t>
inline Bool_t Tron::Console::ReadTerm(Value_t& value) {
  std::string buff;
  while (true) {
    std::cout << "> " << std::flush;
    if (!std::getline(std::cin, buff)) {
      Clear();
      continue;
    } else if (buff.empty()) {
      continue;
    }
    std::stringstream line(buff);
    if (!(line >> value)) {
      continue;
    }
    return true;
  };
}

inline Bool_t Tron::Console::ReadLine(std::string& value) {
  while (true) {
    std::cout << "> " << std::flush;
    if (!std::getline(std::cin, value)) {
      Clear();
      continue;
    } else if (value.empty()) {
      continue;
    }
    return true;
  }
}

template <typename Value_t>
inline Bool_t Tron::Console::ReadLine(Value_t& value) {
  static TFormula formula;
  std::string buff;
  while (true) {
    std::cout << "> " << std::flush;
    if (!std::getline(std::cin, buff)) {
      Clear();
      continue;
    } else if (buff.empty()) {
      continue;
    } else if (formula.Compile((buff + " + x").data())) {
      continue;
    }
    value = formula.Eval(0.0);
    return true;
  }
}

template <typename Value_t, typename Conv_t>
inline Bool_t Tron::Console::ReadLine(Value_t& value, Conv_t converter) {
  std::string buff;
  while (true) {
    std::cout << "> " << std::flush;
    if (!std::getline(std::cin, buff)) {
      Clear();
      continue;
    } else if (buff.empty()) {
      continue;
    } else if (!converter(buff, value)) {
      continue;
    }
    return true;
  }
}

#endif
