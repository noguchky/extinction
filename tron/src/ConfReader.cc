#include <fstream>
#include <regex>
#include "TFormula.h"
#include "ConfReader.hh"
#include "String.hh"
#include "TreeData.hh"

Tron::ConfReader::ConfReader(const std::string& filename, const char* option) {
  // std::cout << "ConfReader::ConfReader" << std::endl;
  Clear();
  if (!filename.empty()) {
    Read(filename.data(), option);
  }
}

Tron::ConfReader::ConfReader(const ConfReader& conf)
  : fFilenames (conf.fFilenames),
    fIndexMap  (conf.fIndexMap),
    fContents  (conf.fContents),
    fParsedLine(conf.fParsedLine) {
}

void Tron::ConfReader::Clear() {
  // std::cout << "ConfReader::Clear" << std::endl;
  fFilenames.clear();
  fIndexMap .clear();
  fContents .clear();
  fParsedLine = 0;
}

void Tron::ConfReader::SetDefaultFormat(const std::string& format) {
  fDefaultFormat = format;
}

int Tron::ConfReader::Read(const std::string& filename, const char* option) {
  // std::cout << "ConfReader::Read" << std::endl;
  const std::string opt = String::ToLower(option);
  if (!String::Contains(opt, "a")) {
    Clear();
  }

  if (!ReadFile(filename)) {
    return 0;
  }

  Parse(fParsedLine);
  fParsedLine = fContents.size();
  
  //--- Verbose Contents
  if (String::Contains(opt, "v")) {
    std::cout << "ConfReader::Read() ... " << std::endl;
    ShowContents();
  }

  return 1;
}

bool Tron::ConfReader::IsOpen() const {
  // std::cout << "ConfReader::IsOpen" << std::endl;
  return !fFilenames.empty();
}

bool Tron::ConfReader::IsOpen(const std::string& filename) const {
  // std::cout << "ConfReader::IsOpen" << std::endl;
  return std::find(fFilenames.begin(), fFilenames.end(), filename) != fFilenames.end();
}

void Tron::ConfReader::ShowContents() const {
  const int kKeyLength = 40;
  
  for (auto&& content : fContents) {
    const std::string keyPart   = GetKeyPart(content);
    const std::string valuePart = GetValuePart(content);
    const int spaceLength = kKeyLength - keyPart.length();
    if (spaceLength > 0) { 
      std::cout << keyPart
                << std::string(spaceLength, ' ')
                << String::Wrap(valuePart) << std::endl;
    } else {
      std::cout << keyPart << std::endl
                << std::string(kKeyLength, ' ')
                << String::Wrap(valuePart) << std::endl;
    }
  }
}

std::vector<std::string> Tron::ConfReader::GetFilenames() const {
  // std::cout << "ConfReader::GetFilename" << std::endl;
  return fFilenames;
}

std::vector<std::string> Tron::ConfReader::GetKeys() const {
  // std::cout << "ConfReader::GetKeys" << std::endl;
  std::vector<std::string> keys;
  for (auto&& keyindex : fIndexMap) {
    keys.push_back(keyindex.first);
  }
  return keys;
}

std::vector<std::string> Tron::ConfReader::GetLines() const {
  // std::cout << "ConfReader::GetLines" << std::endl;
  return fContents;
}

std::string Tron::ConfReader::GetLine(const std::string& key) const {
  // std::cout << "ConfReader::GetLine" << std::endl;
  return Exists(key) ? fContents[fIndexMap.at(key)] : "";
}

bool Tron::ConfReader::Exists(const std::string& key, const char*) const {
  // std::cout << "ConfReader::Exists" << std::endl;
  return fIndexMap.find(key) != fIndexMap.end();
}

void Tron::ConfReader::SetValue(const std::string& key, const std::string& value) {
  const std::string line = key + " " + value;

  if (Exists(key)) {
    fContents[fIndexMap.at(key)] = line;
  } else {
    fIndexMap[key] = fContents.size();
    fContents.push_back(line);
    fParsedLine = fContents.size();
  }
}

void Tron::ConfReader::SetValue(const std::string& key, Double_t value) {
  SetValue(key, Form(fDefaultFormat.data(), value));
}

void Tron::ConfReader::SetValues(const std::string& key, const std::vector<std::string>& values) {
  SetValue(key, Tron::String::Join(values, " "));
}

void Tron::ConfReader::SetValues(const std::string& key, const std::vector<Double_t>& values) {
  SetValue(key, Tron::String::Join(Tron::Linq::From(values).Select([this](Double_t value) { return Form(fDefaultFormat.data(), value); }), " "));
}

int Tron::ConfReader::ReadFile(const std::string& filename) {
  //--- Open File
  std::ifstream file(filename);
  if (!file || !file.is_open()) {
    std::cerr << "ConfReader::ReadFile() [error] file is not open, " << String::Wrap(filename) << std::endl;
    return 0;
  }
  fFilenames.push_back(filename);

  const std::string basepath = GetPath(filename);
  
  //--- Get Contents
  std::string buf;
  bool  isContinuationLine = false;
  bool wasContinuationLine = false;
  while (file && getline(file, buf)) {
    //--- Remove Space
    buf = RemoveSpace(buf);

    //--- Remove Comments
    buf = RemoveComment(buf);

    //--- Remove Escape Sequence of Continuation Line
    buf = RemoveContinuationLine(buf, isContinuationLine);

    //--- Remove Escape Sequence
    buf = RemoveEscapeSequence(buf);

    //--- Replace with Included Contents
    buf = IncludeFile(buf, basepath);

    //--- Append Contents
    if (!String::IsEmptyOrWhiteSpace(buf)) {
      if (wasContinuationLine) {
        fContents.back().append(buf);
      } else {
        fContents.push_back(buf);
      }
    }

    wasContinuationLine = isContinuationLine;
  }
  file.close();

  return 1;
}

void Tron::ConfReader::Parse(int ifrom) {
  //--- Get Key & Index
  std::string key;
  for (unsigned int iline = ifrom; iline < fContents.size(); ++iline) {
    //--- Replace with Alias/Environment/Evaluated
    ReplaceAliasEnvEval(iline);
    
    //--- Replace with Aliases For Each Key
    ReplaceEachKey(iline);
    
    //--- Replace with Aliases For Each Value
    ReplaceEachValue(iline);

    //--- Insert/Update Key & Index
    key = GetKeyPart(fContents[iline]);
    
    if (!String::IsEmptyOrWhiteSpace(key)) {
      fIndexMap[key] = iline;
    }
  }
}

std::string Tron::ConfReader::IncludeFile(const std::string& buf, const std::string& basepath) {
  // std::cout << "ConfReader::ReplaceInclude" << std::endl;
  static const std::regex kIncludeRegex(R"(\$INCLUDE\{([^{}]+)\})");
  std::cmatch match;

  //--- Include File
  std::string str = buf;
  while (std::regex_search(str.data(), match, kIncludeRegex)) {
    const std::string name = String::RemoveSpace(match[1].str());

    //--- Read File
    if (IsAbsolutePath(name)) {
      ReadFile(name);
    } else {
      ReadFile(basepath + name);
    }

    //--- Replace
    str = "";
  }
  return str;
}

std::string Tron::ConfReader::GetKeyPart(const std::string& str) {
  std::stringstream line(str);
  std::string key;
  line >> key;
  return key;
}

std::string Tron::ConfReader::GetValuePart(const std::string& str) { 
  std::stringstream line(str);
  std::string key;
  line >> key;
  return RemoveSpace(str.substr((int)line.tellg() + 1));
}

bool Tron::ConfReader::IsAbsolutePath(const std::string &str) {
  std::string buf = RemoveSpace(str);
  return String::StartWith(buf, "/");
}

std::string Tron::ConfReader::GetPath(const std::string& str) {
  static const std::regex kPathRegex(R"([^/]*$)");
  return std::regex_replace(str, kPathRegex, "");
}

std::string Tron::ConfReader::RemoveSpace(const std::string& str) {
  return String::RemoveSpace(str);
}

std::string Tron::ConfReader::RemoveComment(const std::string& str) {
  static const std::regex kCommentRegex(R"((^|[^\\])#.*)");
  return std::regex_replace(str, kCommentRegex, "");
}

std::string Tron::ConfReader::RemoveEscaped(const std::string& str) {
  static const std::regex kEscapedBackSlash(R"(\\\\)");
  return std::regex_replace(str, kEscapedBackSlash, "\\");
}

std::string Tron::ConfReader::RemoveContinuationLine(const std::string& str, bool& isMatch) {
  static const std::regex kContinuationLine(R"(\\$)");
  isMatch = std::regex_search(str.data(), kContinuationLine);
  return std::regex_replace(str, kContinuationLine, "");
}

std::string Tron::ConfReader::RemoveEscapeSequence(const std::string& str) {
  static const std::regex kEscapeSequence(R"(\\(.))");
  return std::regex_replace(str, kEscapeSequence, "$1");
}

void Tron::ConfReader::ReplaceAliasEnvEval(int iline) {
  // std::cout << "ConfReader::ReplaceAliasEnvEval" << std::endl;
  static const std::regex kAliasOrEnvRegex(R"(\$(ALIAS|ENV|EVAL|)\{([^{}]+)\}(?:\$FMT\{([^{}]+)\})?)");
  std::cmatch match;

  std::string str = fContents[iline];
  while (std::regex_search(str.data(), match, kAliasOrEnvRegex)) {
    const int sample = match.position();
    const int length = match.length();
    const std::string type   = String::RemoveSpace(match[1].str());
    const std::string name   = String::RemoveSpace(match[2].str());
    const std::string format = String::RemoveSpace(match[3].str());

    //--- Replace
    std::string value;
    if        (type == "EVAL") {
      TFormula formula;
      if (formula.Compile((name + " + x").data())) {
        value = 0.0;
      } else {
        const double number = formula.Eval(0.0);
        value = String::ToString(number, (!format.empty() ? format : fDefaultFormat));
      }

    } else if (type == "ALIAS") {
      value = GetValue<std::string>(name.data());

    } else if (type == "ENV") {
      if (const char* env = getenv(name.data())) {
        value = env;
      } else {
        value = "";
      }

    } else {
      value = GetValue<std::string>(name.data());
      if (value.empty()) {
        if (const char* env = getenv(name.data())) {
          value = env;
        } else {
          value = "";
        }
      }

    }
    str.replace(sample, length, value);
  }
  fContents[iline] = str;
}

void Tron::ConfReader::ReplaceEachKey(int iline) {
  // std::cout << "ConfReader::ReplaceEachKey" << std::endl;
  static const std::regex kEachRegex(R"(\$EACH\{([^{}]+)\})");
  std::cmatch match;

  //--- Replace Key
  while (true) {
    const std::string term = GetKeyPart(fContents[iline]);
    if (std::regex_search(term.data(), match, kEachRegex)) {
      const int sample = match.position();
      const int length = match.length();
      const std::string name = String::RemoveSpace(match[1].str());

      //--- Get Replace Values
      auto values = GetValues<std::string>(name.data());
      if (values.empty()) {
        values.push_back("");
      }

      //--- Replace
      const std::string originLine = fContents[iline];
      for (unsigned int ivalue = 0; ivalue < values.size(); ++ivalue) {
        std::string eachLine = originLine;
        eachLine.replace(sample, length, values[ivalue]);
        if (ivalue) {
          fContents.insert(fContents.begin() + iline + ivalue, eachLine);
        } else {
          fContents[iline] = eachLine;
        }
      }
    } else {
      break;
    }
  }
}

void Tron::ConfReader::ReplaceEachValue(int iline) {
  // std::cout << "ConfReader::ReplaceEachValue" << std::endl;
  //--- Replace Value
  std::string key   = GetKeyPart  (fContents[iline]);
  std::string terms = GetValuePart(fContents[iline]);
  ReplaceEachValueWhile(terms);
  fContents[iline] = key + " " + terms;
}

void Tron::ConfReader::ReplaceEachValueWhile(std::string& terms) {
  // std::cout << "ConfReader::ReplaceEachValueWhile" << std::endl;
  static const std::regex kEachRegex(R"(\$EACH\{([^{}]+)\})");
  std::cmatch match;

  std::stringstream line(terms);
  std::stringstream newTerms;
  for (std::string term; line >> term; newTerms << term << " ") {
    if (std::regex_search(term.data(), match, kEachRegex)) {
      const int sample = match.position();
      const int length = match.length();
      const std::string name = String::RemoveSpace(match[1].str());

      //--- Get Replace Values
      auto values = GetValues<std::string>(name.data());
      if (values.empty()) {
        values.push_back("");
      }

      //--- Replace
      std::stringstream buffer;
      for (auto&& value : values) {
        std::string newTerm = term;
        newTerm.replace(sample, length, value);
        buffer << newTerm << " ";
      }
      term = buffer.str();

      //--- Replace Value
      ReplaceEachValueWhile(term);
    }
  }

  terms = newTerms.str();
}

Tron::TreeData Tron::ConfReader::GetAsTree(const char* delim) const {
  // std::cout << "ConfReader::GetAsTree()" << std::endl;
  TreeData tree;

  auto keys = GetKeys();
  for (auto&& key : keys) {
    auto paths = String::Split(key, delim);
    TreeData* branch = &tree;
    for (auto&& path : paths) {
      branch = &(*branch)[path];
    }
    *branch = GetValue(key, "a");
  }

  return tree;
}

Tron::TreeData Tron::ConfReader::GetAsTree(const std::string& basepath, const char* delim) const {
  // std::cout << "ConfReader::GetAsTree()" << std::endl;
  TreeData tree;
  tree = GetValue(basepath, "a");

  auto header = basepath + delim;
  auto keys   = std::vector<std::pair<std::string, std::vector<std::string>>>();
  for (auto&& key : GetKeys()) {
    if (String::StartWith(key, header)) {
      auto paths = String::Split(key.substr(header.length()), delim);

      // Seek End Branch
      TreeData* branch = &tree;
      for (auto&& path : paths) {
        branch = &(*branch)[path];
      }

      // Set Value into End Branch
      *branch = GetValue(key, "a");
    }
  }

  return tree;
}
