#ifndef Tron_TreeData_hh
#define Tron_TreeData_hh

#include <vector>
#include <map>
#include "String.hh"
#include "Linq.hh"

namespace Tron {
  
  class TreeData {
  public:
    static const TreeData gEmpty;
    std::map<std::string, TreeData> fSubBranch;
    std::string                     fData;
    
  public:
    TreeData() = default;

    TreeData(const TreeData& tree)
      : fSubBranch(tree.fSubBranch), fData(tree.fData){
    }

    inline TreeData& operator=(const TreeData& tree) {
      fSubBranch = tree.fSubBranch;
      fData      = tree.fData;
      return *this;
    }

    inline TreeData& operator=(const std::string& data) {
      fData = data;
      return *this;
    }

    inline TreeData& operator[](const std::string& name) {
      return fSubBranch[name];
    }

    inline const TreeData& operator[](const std::string& name) const {
      if (fSubBranch.find(name) == fSubBranch.end()) {
        return gEmpty;
      }
      return fSubBranch.at(name);
    }

    inline const TreeData& Get(const std::string& name) {
      return fSubBranch[name];
    }

    inline const TreeData& Get(const std::string& name) const {
      if (fSubBranch.find(name) == fSubBranch.end()) {
        return gEmpty;
      }
      return fSubBranch.at(name);
    }

    inline operator std::string() const {
      return fData;
    }
    
    template <typename Value_t>
    inline Value_t As(const char* option = "") const {
      return String::Convert<Value_t>(fData, option);
    }

    template <typename Value_t>
    inline std::vector<Value_t> AsVector(const char* option = "") const {
      return String::Converts<Value_t>(fData, option);
    }

    inline std::vector<std::string> GetNames() const {
      return Linq::From(fSubBranch)
        .Select([](std::pair<std::string, TreeData> pair) {
                  return pair.first;
                })
        .ToVector();
    }

  };

}

#endif
