#ifndef Tron_ArgReader_hh
#define Tron_ArgReader_hh

#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <getopt.h>
#include "Types.hh"
#include "String.hh"

namespace Tron {

  class ArgReader {
  private:
    struct BaseConfig {
      std::string              Key;
      std::string              Description;
      std::string              DefaultValue;
      std::vector<std::string> Values;
    };
    
    struct ArgConfig : BaseConfig {
      Bool_t IsRequired;
      Bool_t IsMultiple;
    };
    
    struct OptConfig : BaseConfig {
      enum RequirementType {
                            NoArg       = no_argument,
                            RequiredArg = required_argument,
                            OptionalArg = optional_argument,
      };
      Int_t           ShortOpt;
      std::string     LongOpt;
      RequirementType Requirement;
      std::string     GetOptString() const;
      option          GetOption() const;
    };
    
    std::string                        fAppName;
    std::map<std::string, BaseConfig*> fConfigs;
    std::vector<ArgConfig*>            fArgConfigs;
    std::vector<OptConfig*>            fOptConfigs;

  public:
    ArgReader(const std::string& appName);
    virtual ~ArgReader();

    template <typename Value_t>
    void                        AddArg (const std::string& key,
                                        std::string description, const Char_t* defaultValue = nullptr);

    template <typename Value_t>
    void                        AddArgs(const std::string& key,
                                        std::string description, const Char_t* defaultValue = nullptr);
    
    void                        AddOpt (const std::string& key,
                                        Int_t shortOpt, const std::string& longOpt,
                                        const std::string& description);

    template <typename Value_t>
    void                        AddOpt (const std::string& key,
                                        Int_t shortOpt, const std::string& longOpt,
                                        const std::string& description, const Char_t* defaultValue = nullptr);
    
    Bool_t                      Parse(int argc, Char_t** argv, const Char_t* option = "");

    Bool_t                      HasUnsetRequired() const;
    Bool_t                      IsSet(const std::string& key) const;

    void                        Set(const std::string& key, const std::string& value = "");
    void                        Add(const std::string& key, const std::string& value);
    void                        Unset(const std::string& key);

    void                        ShowUsage(std::ostream& stream = std::cout) const;
    
    template <typename Value_t = std::string>
    inline std::vector<Value_t> GetValues(const std::string& key, const Char_t* optoin = "") const;

    template <typename Value_t = std::string>
    inline Value_t              GetValue(const std::string& key, const Char_t* option = "") const;
  };
  
}

template <typename Value_t>
void Tron::ArgReader::AddArg (const std::string& key,
                              std::string description, const Char_t* defaultValue) {
  ArgConfig* config = new ArgConfig();
  config->Key          = key;
  config->Description  = description;
  config->IsMultiple   = false;
  if (defaultValue == nullptr) {
    config->IsRequired   = true;
  } else {
    config->IsRequired   = false;
    config->DefaultValue = defaultValue;
  }

  fConfigs[key] = config;
  fArgConfigs.push_back(config);
}
    
template <typename Value_t>
void Tron::ArgReader::AddArgs(const std::string& key,
                              std::string description, const Char_t* defaultValue) {
  ArgConfig* config = new ArgConfig();
  config->Key          = key;
  config->Description  = description;
  config->IsMultiple   = true;
  if (defaultValue == nullptr) {
    config->IsRequired   = true;
  } else {
    config->IsRequired   = false;
    config->DefaultValue = defaultValue;
  }

  fConfigs[key] = config;
  fArgConfigs.push_back(config);
}

inline void Tron::ArgReader::AddOpt (const std::string& key,
                                     Int_t shortOpt, const std::string& longOpt,
                                     const std::string& description) {
  OptConfig* config = new OptConfig();
  config->Key          = key;
  config->Description  = description;
  config->ShortOpt     = shortOpt;
  config->LongOpt      = longOpt;
  config->Requirement  = OptConfig::NoArg;

  fConfigs[key] = config;
  fOptConfigs.push_back(config);
}

template <typename Value_t>
void Tron::ArgReader::AddOpt (const std::string& key,
                              Int_t shortOpt, const std::string& longOpt,
                              const std::string& description, const Char_t* defaultValue) {
  OptConfig* config = new OptConfig();
  config->Key          = key;
  config->Description  = description;
  config->ShortOpt     = shortOpt;
  config->LongOpt      = longOpt;
  if (defaultValue == nullptr) {
    config->Requirement = OptConfig::RequiredArg;
  } else {
    config->Requirement  = OptConfig::OptionalArg;
    config->DefaultValue = defaultValue;
  }

  fConfigs[key] = config;
  fOptConfigs.push_back(config);
}

template <typename Value_t>
inline std::vector<Value_t> Tron::ArgReader::GetValues(const std::string& key, const Char_t* option) const {
  // std::cout << "Tron::ArgReader::GetValues" << std::endl;  
  std::vector<Value_t> values;
  
  if (fConfigs.find(key) != fConfigs.end()) {
    const auto* config = fConfigs.at(key);
    if (!config->Values.empty()) {
      for (auto&& value : config->Values) { 
        values.push_back(String::Convert<Value_t>(value, "a"));
      }
    } else {
      values.push_back(String::Convert<Value_t>(config->DefaultValue, "a"));
    }
  } else {
    std::cout << "[warning] key was not found, " << String::Wrap(key) << std::endl;
  }
  
  std::string opt = String::ToLower(option);
  if (String::Contains(opt, "v")) {
    std::cout << "ArgReader::GetValues() "
              << "... key = " << String::Wrap(key)
              << ", value = " << String::Wrap(String::Join(values, "\", \"")) << std::endl;
  }
  
  return values;
}

template <typename Value_t>
inline Value_t Tron::ArgReader::GetValue(const std::string& key, const Char_t* option) const {
  // std::cout << "Tron::ArgReader::GetValue" << std::endl;  
  Value_t value {};

  if (fConfigs.find(key) != fConfigs.end()) {
    const auto* config = fConfigs.at(key);
    if (!config->Values.empty()) {
      value = String::Convert<Value_t>(config->Values.front());
    } else {
      value = String::Convert<Value_t>(config->DefaultValue, "a");
    }
  } else {
    std::cout << "[warning] key was not found, " << String::Wrap(key) << std::endl;
  }

  std::string opt = String::ToLower(option);
  if (String::Contains(opt, "v")) {
    std::cout << "ArgReader::GetValue() "
              << "... key = " << String::Wrap(key)
              << ", value = " << String::Wrap(value) << std::endl;
  }

  return value;
}

#endif
