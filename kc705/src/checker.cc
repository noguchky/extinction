#include <iostream>
#include <fstream>
#include <string>
#include "Kc705.hh"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  // if (argc < 2) {
  //   std::cout << "Usage: " << argv[0] << " [Input] [(optional)Output] [(optional)Entries]" << std::endl;
  //   return 1;
  // }

  // const std::string ifilename(argv[1]);
  // const std::string ofilename(argc > 2 ? argv[2] : "checker.pdf");
  // const std::size_t amount          = 1000UL;

  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"     ,                    "A rawdata filename");
  args->AddOpt<std::size_t>("Start"     , 's', "start"     , "Read start position [record]", "0");
  args->AddOpt<std::size_t>("Amount"    , 'a', "amount"    , "Scroll amount [record]", "1000");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename       = args->GetValue("Input");
  const std::size_t start           = args->GetValue<std::size_t>("Start");
  const std::size_t amount          = args->GetValue<std::size_t>("Amount");

  std::cout << "=== Open Input File" << std::endl;
  std::ifstream ifile(ifilename, std::ios::binary);
  if (!ifile) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  std::cout << "=== Read" << std::endl;
  Extinction::Kc705::Packet_t packet;
  ULong64_t entry = 0ULL;
  for (; ifile.read((char*)packet, sizeof(packet)); ++entry) {
    if (entry < start) {
      if ((entry + 1) % 100000UL == 0) {
        std::cout << ">> " << entry << std::endl;
      }
    } else {
      for (std::size_t i = 0, n = sizeof(packet); i < n; ++i) {
        printf("%02X ", packet[i]);
      }
      if ((entry - start + 1) % amount == 0) {
        while (std::getchar() != '\n');
      } else {
        puts("");
      }
    }
  }
  std::cout << "[info] file end, entries = " << entry << std::endl;

  ifile.close();

  return 0;
}
