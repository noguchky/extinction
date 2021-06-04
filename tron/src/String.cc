#include <algorithm>
#include <iostream>
#include <regex>
#include "TFormula.h"
#include "String.hh"

bool Tron::String::StartWith(const std::string& str, const std::string& pattern) {
  // std::cout << "String::StartWith" << std::endl;
  // return str.find(pattern, 0) == 0; // If str was too long, this method would get slow.
  return str.length() > pattern.length() && std::equal(str.begin(), str.begin() + pattern.length(), pattern.begin());
}

bool Tron::String::EndWith(const std::string& str, const std::string& pattern) {
  // std::cout << "String::StartWith" << std::endl;
  return str.length() > pattern.length() && std::equal(str.end() - pattern.length(), str.end(), pattern.begin());
}

bool Tron::String::Contains(const std::string& str, const std::string& pattern) {
  // std::cout << "String::Contains" << std::endl;
  return str.find(pattern, 0) != std::string::npos;
}

bool Tron::String::IsEmptyOrWhiteSpace(const std::string& str) {
  // std::cout << "String::IsEmptyOrWhiteSpace" << std::endl;
  std::string buf = RemoveSpace(str);
  return buf.empty();
}

int Tron::String::FindCloseParen(const std::string& str, int index, char openParen, char closeParen) {
  // std::cout << "String::FindCloseParen" << std::endl;
  int i = index;
  int hierarchy = 1;
  for (int n = str.size(); i < n; ++i) {
    if        (str[i] ==  openParen) {
      ++hierarchy;
    } else if (str[i] == closeParen) {
      --hierarchy;
    }
    if (hierarchy == 0) {
      break;
    }
  }  
  return i;
}

std::vector<std::string> Tron::String::Split(const std::string& str, const std::string& delim) {
  // std::cout << "String::Split" << std::endl;
  std::string::size_type delimLength = delim.length();
  std::vector<std::string> results;

  if (delimLength == 0) {
    results.push_back(str);
  } else {
    std::string::size_type pos = 0, offset = pos;
    for (;;) {
      pos = str.find(delim, offset);
      if (pos == std::string::npos) {
        // Regist all the rest
        results.push_back(str.substr(offset));
        break;
      } else {
        // Regist before delimiter
        results.push_back(str.substr(offset, pos - offset));
        offset = pos + delimLength;
      }
    }
  }
  return results;
}

std::string Tron::String::Wrap(const std::string& str, const char* prefix, const char* suffix) {
  return (prefix ? prefix : "") + str + (suffix ? suffix : prefix ? prefix : "");
}

std::string Tron::String::Wrap(char chara, const char* prefix, const char* suffix) {
  return (prefix ? prefix : "") + std::string(1, chara) + (suffix ? suffix : prefix ? prefix : "");
}

std::string Tron::String::Replace(const std::string& str, const std::string& pattern, const std::string& replacement) {
  // std::cout << "String::Replace" << std::endl;
  std::string buf = str;

  const std::string::size_type patternLength = pattern.length();
  const std::string::size_type replaceLength = replacement.length();
  if (patternLength != 0) {
    std::string::size_type pos = 0;
    while ((pos = buf.find(pattern, pos)) != std::string::npos) {
      buf.replace(pos, patternLength, replacement);
      pos += replaceLength;
    }
  }
  return buf;
}

std::string Tron::String::ToLower(const std::string& str) {
  // std::cout << "String::ToLower" << std::endl;
  std::string buf = str;
  std::transform(buf.begin(), buf.end(), buf.begin(), ::tolower);
  return buf;
}

std::string Tron::String::ToUpper(const std::string& str) {
  // std::cout << "String::ToUpper" << std::endl;
  std::string buf = str;
  std::transform(buf.begin(), buf.end(), buf.begin(), ::toupper);
  return buf;
}

std::string Tron::String::ToString(double number, const std::string& format) {
  // std::cout << "String::ToString" << std::endl; 
  if (format.empty()) {
    // std::cerr << "String::ToString [trace] format is empty" << std::endl;
    return std::to_string(number);
  }

  std::string formatSpecifier = GetFormatSpecifiers(format);
  if (formatSpecifier.empty()) {
    // std::cout << "String::ToString [trace] format specifier is empty, "
    //           << String::Wrap(format) << std::endl;
    return format;
  } else {
    // std::cout << "String::ToString [trace] format specifier is "
    //           << String::Wrap(format) << std::endl;
  }

  if        (formatSpecifier == "c" ||
             formatSpecifier == "C" ||
             formatSpecifier == "d" ||
             formatSpecifier == "i") {
    return Format(format, (int)number);

  } else if (formatSpecifier == "u" ||
             formatSpecifier == "o" ||
             formatSpecifier == "x" ||
             formatSpecifier == "X") {
    return Format(format, (unsigned int)number);

  } else if (formatSpecifier == "hd" ||
             formatSpecifier == "hi") {
    return Format(format, (short)number);

  } else if (formatSpecifier == "hu") {
    return Format(format, (unsigned short)number);
    
  } else if (formatSpecifier == "ld" ||
             formatSpecifier == "li") {
    return Format(format, (long)number);
    
  } else if (formatSpecifier == "lu" ||
             formatSpecifier == "lo" ||
             formatSpecifier == "lx" ||
             formatSpecifier == "lX") {
    return Format(format, (unsigned long)number);
    
  } else if (formatSpecifier == "s" ||
             formatSpecifier == "S" ||
             formatSpecifier == "n" ||
             formatSpecifier == "p") {
    return format;
    
  }
  
  return Format(format, number);
}

template <>
std::string Tron::String::Convert<std::string>(const std::string& str, const char* option) {
  std::string opt = ToLower(option);
  if (String::Contains(opt, "a")) {
    return str;
  }
  std::stringstream line(str);
  std::string value;
  line >> value;
  return value;
}

std::string Tron::String::GetFormatSpecifiers(const std::string& str) {
  // std::cout << "String::GetFormatSpecifiers" << std::endl;
  static const std::regex kFormatSpecifiersRegex(R"((?:^|[^%])(?:%%)*%[-+0-9.]*([cCsSdiuoxXfFaAeEgGnp]|h[diu]|[lL][diuoxXfeE]))");
  std::cmatch match;

  //--- Include File
  std::string buf = str;
  if (std::regex_search(str.data(), match, kFormatSpecifiersRegex)) {
    buf = match[1].str();
  } else {
    buf.clear();
  }
  return buf;
}

std::string Tron::String::RemoveSpace(const std::string& str) {
  // std::cout << "String::RemoveSpace" << std::endl;
  static const std::regex kSpaceRegex(R"(^\s+|\s+$)");
  return std::regex_replace(str, kSpaceRegex, "");
}

double Tron::String::EvalFormula(const std::string& str) {
  double value = 0.0;
  TFormula form;
  if (!form.Compile((str + " + x").data())) {
    value = form.Eval(0.0);
  }
  return value;
}
