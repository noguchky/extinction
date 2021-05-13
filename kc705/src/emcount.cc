#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <ctime>
#include <regex>
#include "TFile.h"
#include "TTree.h"
#include "Kc705.hh"
#include "Units.hh"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"     ,                    "Set rawdata filename");
  args->AddOpt<Int_t>      ("Board"     , 'b', "board"     , "Set board id");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const Int_t       board     = args->GetValue<Int_t>("Board");

  // Open Input File
  std::istream* istr = nullptr;
  std::ifstream ifile;
  if (ifilename == "std::cin") {
    istr = &std::cin;
  } else {
    ifile.open(ifilename, std::ios::binary);
    if (!ifile) {
      std::cerr << "[error] input file is not opened, " << ifilename << std::endl;
      return 1;
    }
    istr = &ifile;
  }

  Extinction::Kc705::Decoder decoder;

  // Read file
  while (decoder.Read(*istr)) {
    if (decoder.Data.IsFooter()) {
      break;
    }
  }
  ifile.close();

  const Int_t emcount = decoder.Data.EMCount;
  if (emcount < 0) {
    // warning messages were sent from DecodeEventMatchNumber
    std::cerr << "[warning] invalid emcount, " << ifilename << std::endl;
    return 1;
  }

  if (decoder.Data.Date) {
    if (auto tm = std::localtime((time_t*)&decoder.Data.Date)) {
      const Int_t year   = tm->tm_year + 1900;
      const Int_t month  = tm->tm_mon + 1;
      const Int_t day    = tm->tm_mday;
      const Int_t hour   = tm->tm_hour;
      const Int_t minute = tm->tm_min;
      const Int_t second = tm->tm_sec;
      const std::string strdate = Form("%04d-%02d-%02d_%02d:%02d:%02d", year, month, day, hour, minute, second);
      std::cout << emcount   << "\t"
                << board     << "\t"
                << strdate   << "\t"
                << ifilename << std::endl;
    } else {
      std::cerr << "[error] time struct is not obtained" << std::endl;
    }
  } else {
    const std::regex pattern(R"((\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2}))");
    std::cmatch match;
    if (std::regex_search(ifilename.data(), match, pattern)) {
      const Int_t year   = Tron::String::Convert<Int_t>(match[1].str());
      const Int_t month  = Tron::String::Convert<Int_t>(match[2].str());
      const Int_t day    = Tron::String::Convert<Int_t>(match[3].str());
      const Int_t hour   = Tron::String::Convert<Int_t>(match[4].str());
      const Int_t minute = Tron::String::Convert<Int_t>(match[5].str());
      const Int_t second = Tron::String::Convert<Int_t>(match[6].str());
      const std::string strdate = Form("%04d-%02d-%02d_%02d:%02d:%02d", year, month, day, hour, minute, second);
      std::cout << emcount   << "\t"
                << board     << "\t"
                << strdate   << "\t"
                << ifilename << std::endl;
    } else {
      std::cerr << "[error] date information was not obtained" << std::endl;
    }
  }

  return 0;
}
