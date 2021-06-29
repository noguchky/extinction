#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include "TFile.h"
#include "TTree.h"
#include "Fct.hh"
#include "Units.hh"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"    ,                   "Set rawdata filename");
  args->AddOpt             ("Help"     , 'h', "help"     , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename  = args->GetValue("Input");

  std::cout << "=== Open Input File" << std::endl;
  std::istream* istr = nullptr;
  std::ifstream ifile;
  if (ifilename == "std::cin") {
    istr = &std::cin;
  } else {
    ifile.open(ifilename, std::ios::binary);
    if (!ifile) {
      std::cout << "[error] input file is not opened, " << ifilename << std::endl;
      return 1;
    }
    istr = &ifile;
  }

  Extinction::Fct::Decoder decoder;

  std::cout << "=== Initialize Variables" << std::endl;
  Extinction::Fct::Packet_t packet = 0;

  std::cout << "=== Decode" << std::endl;
  std::size_t count = 0UL;
  for (; decoder.Read(*istr, &packet); ++count) {
    std::cout << count << "\t"
              << Form("0x%08x", packet) << "\t"
              << (decoder.Data.Type == 0 ? "None" :
                  decoder.Data.Type == 1 ? "Header" :
                  decoder.Data.Type == 2 ? "GateStart" :
                  decoder.Data.Type == 3 ? "GateEnd" :
                  decoder.Data.Type == 4 ? "Carry" :
                  decoder.Data.Type == 5 ? "Data" :
                  decoder.Data.Type == 6 ? "Error" :
                  "Unknown") << "\t"
              << Form("%02d", decoder.Data.Channel) << "\t"
              << Form("%06x", decoder.Data.Tdc) << "\t"
              << Form("%02x", decoder.Data.Carry) << std::endl;
  }
  std::cout << "# of data record = " << count << std::endl;

  std::cout << "=== Close Files" << std::endl;
  ifile.close();

  return 0;
}
