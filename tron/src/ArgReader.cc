#include "ArgReader.hh"

std::string Tron::ArgReader::OptConfig::GetOptString() const {
  std::string str;
  str += (char)ShortOpt;
  if (Requirement == RequiredArg || Requirement == OptionalArg) {
    str += ":";
  }
  return str;
}

option Tron::ArgReader::OptConfig::GetOption() const {
  return { LongOpt.data(), Requirement, nullptr, ShortOpt };
}

Tron::ArgReader::ArgReader(const std::string& appName)
  : fAppName(appName), fConfigs() {
}

Tron::ArgReader::~ArgReader() {
  for (auto&& config : fConfigs) {
    delete config.second;
    config.second = nullptr;
  }
  for (auto&& config : fArgConfigs) {
    config = nullptr;
  }
  for (auto&& config : fOptConfigs) {
    config = nullptr;
  }
  fConfigs.clear();
  fArgConfigs.clear();
  fOptConfigs.clear();
}

bool Tron::ArgReader::Parse(int argc, char** argv) {
  // Load Options With Option
  std::string optstring;
  std::vector<option> longopts;
  for (auto&& config : fOptConfigs) {
    optstring += config->GetOptString();
    longopts.push_back(config->GetOption());
  }
  
  int opt;
  int longindex;
  while ((opt = getopt_long(argc, argv, optstring.data(), longopts.data(), &longindex)) != -1) {
    std::vector<OptConfig*> configs;
    for (auto&& config : fOptConfigs) {
      if (config->ShortOpt == opt) {
        configs.push_back(config);
      }
    }
    
    if (configs.empty()) {
      std::cout << "ArgReader::Parse [error] unknown option " << String::Wrap((char)optopt) << std::endl;
      return false;
    }
      
    BaseConfig* config   = configs.front();
    std::string argValue = (optarg != nullptr ? optarg : config->DefaultValue);
    config->Values.push_back(argValue);
  }

  // Load Options With Required/Optional
  int    argc2 = argc - optind;
  char** argv2 = argv + optind;

  int iconfig = 0;
  for (int iarg = 0, nconfig = fArgConfigs.size();
       iarg < argc2 && iconfig < nconfig; ++iarg) {
    ArgConfig* config = fArgConfigs[iconfig];

    std::string argValue = argv2[iarg];
    config->Values.push_back(argValue);
    
    if (!config->IsMultiple) {
      ++iconfig;
    }
  }
    
  return true;
}

bool Tron::ArgReader::HasUnsetRequired() const {
  std::vector<ArgConfig*> notSatisfied;
  for (auto&& config : fArgConfigs) {
    if (config->IsRequired && config->Values.empty()) {
      notSatisfied.push_back(config);
      break;
    }
  }
  if (notSatisfied.empty()) {
    return false;
  }
  
  std::cout << "ArgReader::Parse [error] do not meet requirement for " << String::Wrap(notSatisfied.front()->Key) << std::endl;
  return true;
}

bool Tron::ArgReader::IsSet(const std::string& name) const {
  return (fConfigs.find(name) != fConfigs.end()) && (!fConfigs.at(name)->Values.empty());
}

void Tron::ArgReader::Set(const std::string& key, const std::string& value) {
  if (fConfigs.find(key) == fConfigs.end()) {
    std::cout << "ArgReader::Set [error] do not exist " << String::Wrap(key) << std::endl;
    return;
  }
  
  auto values = fConfigs[key]->Values;
  values.clear();
  values.push_back(value);
}

void Tron::ArgReader::Add(const std::string& key, const std::string& value) {
  if (fConfigs.find(key) == fConfigs.end()) {
    std::cout << "ArgReader::Set [error] do not exist " << String::Wrap(key) << std::endl;
    return;
  }
  
  auto values = fConfigs[key]->Values;
  values.push_back(value);
}

void Tron::ArgReader::Unset(const std::string& key) {
  if (fConfigs.find(key) == fConfigs.end()) {
    std::cout << "ArgReader::Unset [error] do not exist " << String::Wrap(key) << std::endl;
    return;
  }
  
  auto values = fConfigs[key]->Values;
  values.clear();
}

void Tron::ArgReader::ShowUsage() const {
  std::stringstream message, help1, help2;

  message << "Usage: " << fAppName;
  help2 << "Mandatory arguments to long options are mandatory for short options too." << std::endl;

  if        (fOptConfigs.size() == 1) {
    message << " [option]";
  } else if (fOptConfigs.size() >  1) {
    message << " [option]...";
  }
    
  const int FullLength = 30;
  for (auto&& config : fArgConfigs) {
    message << " ";
    if (config->IsRequired) {
      message << config->Key;
    } else {
      message << "[" << config->Key << "]";
    }
    if (config->IsMultiple) {
      message << "...";
    }
        
    int length = FullLength;
    help1 << "  " << config->Key;
    length -= 2 + config->Key.length();

    if (length > 0) {
      help1 << std::string(length, ' ');
    } else {
      help1 << std::endl
            << std::string(FullLength, ' ');
    }
    help1 << config->Description << std::endl;
  }

  for (auto&& config : fOptConfigs) {
    int length = FullLength;
    help2 << "  -" << (char)config->ShortOpt;
    length -= 4;
        
    if (!config->LongOpt.empty()) {
      help2 << ", --" << config->LongOpt;
      length -= 4 + config->LongOpt.length();
    }
        
    if        (config->Requirement == OptConfig::RequiredArg) {
      help2 << "=" << config->Key;
      length -= 1 + config->Key.length();
    } else if (config->Requirement == OptConfig::OptionalArg) {
      help2 << "=[" << config->Key << "]";
      length -= 3 + config->Key.length();
    }

    if (length > 0) {
      help2 << std::string(length, ' ');
    } else {
      help2 << std::endl
            << std::string(FullLength, ' ');
    }
    help2 << config->Description << std::endl;
  }

  std::cout << message.str() << std::endl;
  if (!fArgConfigs.empty()) {
    std::cout << help1.str();
  }
  if (!fOptConfigs.empty()) {
    std::cout << std::endl
              << help2.str();
  }
}
