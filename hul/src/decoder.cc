#include <iostream>
#include <string>
#include <getopt.h>
#include "TFile.h"
#include "TParameter.h"
#include "Hul.hh"
#include "Units.hh"

Int_t main(Int_t argc, Char_t** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " [Input] [(option)Output]" << std::endl;
    return 1;
  }

  const std::string ifilename(argv[1]);
  const std::string ofilename(argc > 2 ? argv[2] : "");
  const Double_t    clock     = 1.04 * Extinction::GHz;
  const Bool_t      isVerbose = true;

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
  std::ifstream ifile(ifilename, std::ios::binary);
  if (!ifile) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  if (isVerbose) std::cout << "=== Create Output File" << std::endl;
  std::cout << ofilename2 << std::endl;
  TFile* ofile = new TFile(ofilename2.data(), "RECREATE");
  if (!ofile->IsOpen()) {
    std::cout << "[error] output file is not opened, " << ofilename2 << std::endl;
    return 1;
  }

  Extinction::Hul::Decoder decoder;

  if (isVerbose) std::cout << "=== Initialize Tree" << std::endl;
  decoder.InitializeTree();

  if (isVerbose) std::cout << "=== Decode" << std::endl;
  std::size_t count = 0UL;
  for (; decoder.Read(ifile); ++count) {
    if (isVerbose && count % 100000 == 0) {
      std::cout << ">> " << count << std::endl;
    }
    if (decoder.Data.Type == Extinction::Hul::DataType::Data) {
      decoder.Tree->Fill();
    }
  }
  std::cout << "# of data record = " << count << std::endl;

  if (isVerbose) std::cout << "=== Write Objects" << std::endl;
  std::cout << decoder.Tree->GetName() << std::endl;
  decoder.Tree->Write();
  auto param = new TParameter<Double_t>("ClockFreq", clock);
  std::cout << param->GetName() << " = " << param->GetVal() << std::endl;
  param->Write();

  if (isVerbose) std::cout << "=== Close Files" << std::endl;
  ifile.close();
  ofile->Close();

  return 0;
}
