#ifndef Tron_ConfReader_hh
#define Tron_ConfReader_hh

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include "Types.hh"
#include "String.hh"

namespace Tron {  

  class TreeData;

  class ConfReader {
  private:
    std::vector<std::string>     fFilenames;
    std::map<std::string, Int_t> fIndexMap;
    std::vector<std::string>     fContents;
    UInt_t                       fParsedLine;
    std::string                  fDefaultFormat = "%23.15e";

  public:
    ConfReader() = default;
    ConfReader(const std::string& filename, const Char_t* option = "");
    ConfReader(const ConfReader& conf);
    ~ConfReader() { }

    void                     SetDefaultFormat(const std::string& format);
    Int_t                    Read(const std::string& filename, const Char_t* option = "");
    Bool_t                   IsOpen() const;
    Bool_t                   IsOpen(const std::string& filename) const;
    void                     ShowContents() const;
    void                     Clear();
    std::string              GetDefaultFormat() const { return fDefaultFormat; }
    std::vector<std::string> GetFilenames() const;
    std::vector<std::string> GetKeys() const;
    std::vector<std::string> GetLines() const;
    std::string              GetLine(const std::string& key) const;
    Bool_t                   Exists(const std::string& key, const Char_t* option = "") const;
    TreeData                 GetAsTree(const Char_t* delim = ".") const;
    TreeData                 GetAsTree(const std::string& name, const Char_t* delim = ".") const;

    template <typename Value_t = std::string>
    Value_t                  GetValue(const std::string& key, const Char_t* option = "") const;

    template <typename Value_t = std::string>
    std::vector<Value_t>     GetValues(const std::string& key, const Char_t* option = "") const;

    void                     SetValue (const std::string& key, const std::string& value);
    void                     SetValue (const std::string& key, Double_t value);
    void                     SetValues(const std::string& key, const std::vector<std::string>& values);
    void                     SetValues(const std::string& key, const std::vector<Double_t>& values);

  private:
    Int_t                    ReadFile(const std::string& filename);
    void                     Parse                (Int_t ifrom);
    std::string              IncludeFile          (const std::string& buf, const std::string& basepath);

    static std::string       GetKeyPart            (const std::string& str);
    static std::string       GetValuePart          (const std::string& str);
    static Bool_t            IsAbsolutePath        (const std::string& str);
    static std::string       GetPath               (const std::string& str);
    static std::string       RemoveSpace           (const std::string& str);
    static std::string       RemoveComment         (const std::string& str);
    static std::string       RemoveEscaped         (const std::string& str);
    static std::string       RemoveContinuationLine(const std::string& str, Bool_t& isMatch);
    static std::string       RemoveEscapeSequence  (const std::string& str);
    void                     ReplaceAliasEnvEval   (Int_t iline);
    void                     ReplaceEachKey        (Int_t iline);
    void                     ReplaceEachValue      (Int_t iline);
    void                     ReplaceEachValueWhile (std::string& term);
  };

}

template <typename Value_t>
Value_t Tron::ConfReader::GetValue(const std::string& key, const Char_t* option) const {
  // std::cout << "ConfReader::GetValue" << std::endl;
  auto linestr   = GetLine(key);
  auto valuepart = GetValuePart(linestr);
  auto value     = String::Convert<Value_t>(valuepart, option);
  if (linestr.empty()) {
    std::cerr << "[warning] key is not found, " << key << std::endl;
  }

  std::string opt = String::ToLower(option);
  if (String::Contains(opt, "v")) {
    std::cout << "ConfReader::GetValue() "
              << "... key = " << String::Wrap(key)
              << ", value = " << String::Wrap(value) << std::endl;
  }
  return value;
}

template <typename Value_t>
std::vector<Value_t> Tron::ConfReader::GetValues(const std::string& key, const Char_t* option) const {
  // std::cout << "ConfReader::GetValues" << std::endl;
  auto linestr   = GetLine(key);
  auto valuepart = GetValuePart(linestr);
  auto values    = String::Converts<Value_t>(valuepart, option);
  if (linestr.empty()) {
    std::cerr << "[warning] key is not found, " << key << std::endl;
  }

  std::string opt = String::ToLower(option);
  if (String::Contains(opt, "v")) {
    std::cout << "ConfReader::GetValues() "
              << "... key = " << String::Wrap(key)
              << ", value = " << String::Wrap(String::Join(values, "\", \"")) << std::endl;
  }
  return values;
}

#endif
