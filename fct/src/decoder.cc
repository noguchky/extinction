#include <iostream>
#include <string>
#include <getopt.h>
#include "TFile.h"
#include "TParameter.h"
#include "Fct.hh"
#include "Units.hh"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  // if (argc < 2) {
  //   std::cout << "Usage: " << argv[0] << " [Input] [(optional)Output]" << std::endl;
  //   return 1;
  // }

  // const std::string ifilename(argv[1]);
  // const std::string ofilename(argc > 2 ? argv[2] : "");
  // const Bool_t      isVerbose = false;

  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"     ,                    "A rawdata filename");
  args->AddOpt<std::string>("Output"    , 'o', "output"    , "A output filename", "");
  args->AddOpt             ("Verbose"   , 'v', "verbose"   , "Output filled values");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const std::string ofilename = args->GetValue("Output");
  const Bool_t      isVerbose = args->IsSet("Verbose");

  std::string ofilename2;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    if (buff.EndsWith(".dat") ||
        buff.EndsWith(".txt")) {
      buff.Replace(buff.Length() - 4, 4, "");
    } else if (buff.EndsWith(".data")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofilename2 = buff + ".root";
  } else {
    ofilename2 = ofilename;
  }

  if (isVerbose) std::cout << "=== Open Input File" << std::endl;
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

  if (isVerbose) std::cout << "=== Create Output File" << std::endl;
  std::cout << ofilename2 << std::endl;
  TFile* ofile = new TFile(ofilename2.data(), "RECREATE");
  if (!ofile->IsOpen()) {
    std::cout << "[error] output file is not opened, " << ofilename2 << std::endl;
    return 1;
  }

  Extinction::Fct::Decoder decoder;

  if (isVerbose) std::cout << "=== Initialize Tree" << std::endl;
  decoder.InitializeTree();

  if (isVerbose) std::cout << "=== Decode" << std::endl;
  std::size_t count = 0UL;
  for (; decoder.Read(*istr); ++count) {
    if (isVerbose && count % 100000 == 0) {
      std::cout << ">> " << count << std::endl;
    }
    if (decoder.Data.IsData()) {
      decoder.Tree->Fill();
    }
  }
  std::cout << "# of data record = " << count << std::endl;

  if (isVerbose) std::cout << "=== Write Objects" << std::endl;
  std::cout << decoder.Tree->GetName() << std::endl;
  decoder.Tree->Write();

  if (isVerbose) std::cout << "=== Close Files" << std::endl;
  ifile.close();
  ofile->Close();

  return 0;
}
