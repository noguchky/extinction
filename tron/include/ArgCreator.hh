#ifndef Tron_ArgCreator_hh
#define Tron_ArgCreator_hh

#include <string>
#include <vector>
#include "Types.hh"

namespace Tron {

  class ArgCreator {
  private:
    std::vector<std::string> fArgs;

  public:
    ArgCreator(const std::string& appName);

    void AddArg(const std::string& arg);

    Int_t    GetArgc() const;
    Char_t** GetArgv() const;

  };
  
}

#endif
