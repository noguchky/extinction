#include <cstring>
#include "ArgCreator.hh"

Tron::ArgCreator::ArgCreator(const std::string& appName) {
  fArgs.push_back(appName);
}

void Tron::ArgCreator::AddArg(const std::string& arg) {
  fArgs.push_back(arg);
}

int Tron::ArgCreator::GetArgc() const {
  return fArgs.size();
}

char** Tron::ArgCreator::GetArgv() const {
  char** result = nullptr;

  if (fArgs.empty()) {
    return result;
  }

  result = new char*[fArgs.size()];
  for (int i = 0, n = fArgs.size(); i < n; ++i) {
    const std::string::size_type length = fArgs[i].size();
    result[i] = new char[length + 1];
    std::memcpy(result[i], fArgs[i].data(), length);
    result[i][length] = '\0';
  }

  return result;
}
