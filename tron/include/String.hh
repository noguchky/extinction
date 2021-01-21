#ifndef Tron_String_hh
#define Tron_String_hh

#include <string>
#include <sstream>
#include <cstdio>
#include <vector>

namespace Tron {

  namespace String {

    bool                     StartWith(const std::string& str, const std::string& pattern);
    bool                     EndWith(const std::string& str, const std::string& pattern);

    bool                     Contains(const std::string& str, const std::string& pattern);

    bool                     IsEmptyOrWhiteSpace(const std::string& str);

    int                      FindCloseParen(const std::string& str, int index, char openParen, char closeParen);

    std::vector<std::string> Split(const std::string& str, const std::string& delim);

    std::string              Replace(const std::string& str, const std::string& pattern, const std::string& replacement);

    std::string              Wrap(const std::string& str, const char* prefix = "\"", const char* suffix = nullptr);

    std::string              Wrap(char chara, const char* prefix = "\'", const char* suffix = nullptr);

    std::string              ToLower(const std::string& str);

    std::string              ToUpper(const std::string& str);

    template <typename ... Args>
    std::string              Format(const std::string& format, Args... args);

    std::string              ToString(double number, const std::string& format);

    template <typename Iterator>
    std::string              Join(Iterator begin, Iterator end, const std::string& delim = "");

    template <typename Iterable>
    std::string              Join(const Iterable& iterable, const char* delim = "");

    template <typename Value_t>
    std::string              Join(const std::initializer_list<Value_t>& list, const char* delim = "");

    template <typename Value_t>
    Value_t                  Convert(const std::string& str, const char* option = "");

    template <>
    std::string              Convert<std::string>(const std::string& str, const char* option);

    template <typename Value_t>
    std::vector<Value_t>     Converts(const std::string& str, const char* option = "");

    double                   EvalFormula(const std::string& str);

    std::string              GetFormatSpecifiers(const std::string& str);

    std::string              RemoveSpace(const std::string& str);

  }

}

template <typename ... Args>
std::string Tron::String::Format(const std::string& format, Args... args) {
  const std::size_t length = std::snprintf( nullptr, 0, format.data(), args...);
  std::vector<char> buf(length + 1);
  std::snprintf(&buf[0], length + 1, format.data(), args...);
  return std::string(&buf[0], &buf[0] + length);
}

template <typename Iterator>
std::string Tron::String::Join(Iterator begin, Iterator end, const std::string& delim) {
  std::stringstream buffer;
  auto it = begin;
  if (it != end) {
    for (buffer << *it++; it != end; buffer << delim << *it++);
  }
  return buffer.str();
}

template <typename Iterable>
std::string Tron::String::Join(const Iterable& iterable, const char* delim) {
  return Join(iterable.begin(), iterable.end(), delim);
}

template <typename Value_t>
std::string Tron::String::Join(const std::initializer_list<Value_t>& list, const char* delim) {
  return Join(list.begin(), list.end(), delim);
}

template <typename Value_t>
Value_t Tron::String::Convert(const std::string& str, const char*) {
  std::stringstream line(str);
  Value_t value;
  line >> value;
  return value;
}

template <typename Value_t>
std::vector<Value_t> Tron::String::Converts(const std::string& str, const char*) {
  std::stringstream line(str);
  std::vector<Value_t> values;
  Value_t value;
  while (line >> value) {
    values.push_back(value);
  }
  return values;
}

#endif
