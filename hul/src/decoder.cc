#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "Hul.hh"
#include "Units.hh"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"     ,                    "Set rawdata filename");
  args->AddOpt<std::string>("Output"    , 'o', "output"    , "Set output filename", "");
  args->AddOpt<Int_t>      ("EMChannel" , 'e', "cmchannel" , "Set channel of event match", "2");
  args->AddOpt<Int_t>      ("EMCount"   , 'c', "cmcount"   , "Set default count of event match", "-1");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename  = args->GetValue("Input");
  const std::string ofilename  = args->GetValue("Output");
  const Int_t       emChannel  = args->GetValue<Int_t>("EMChannel");
  const Int_t       emDefCount = args->GetValue<Int_t>("EMCount");
  const Double_t    clock      = 1.04 * Extinction::GHz;

  std::string ofilenameRoot;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    if (buff.EndsWith(".dat") ||
        buff.EndsWith(".txt")) {
      buff.Replace(buff.Length() - 4, 4, "");
    } else if (buff.EndsWith(".data")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofilenameRoot = buff + ".root";
  } else {
    ofilenameRoot = ofilename;
  }

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

  std::cout << "=== Create Output File" << std::endl;
  std::cout << ofilenameRoot << std::endl;
  TFile* ofile = new TFile(ofilenameRoot.data(), "RECREATE");
  if (!ofile->IsOpen()) {
    std::cout << "[error] output file is not opened, " << ofilenameRoot << std::endl;
    return 1;
  }

  Extinction::Hul::Decoder decoder;
  decoder.Data.ClockFreq = clock;
  std::vector<Extinction::TdcData> emdata;

  std::cout << "=== Initialize Tree" << std::endl;
  decoder.InitializeTree();

  std::cout << "=== Decode" << std::endl;
  std::vector<std::pair<Long64_t, Int_t>> emCount;
  Int_t nextEmCount = emDefCount;
  emCount.push_back({ std::numeric_limits<Long64_t>::max(), nextEmCount });
  std::size_t count = 0UL;
  for (; decoder.Read(*istr); ++count) {
    if (count % 100000 == 0) {
      std::cout << ">> " << count << std::endl;
    }
    if (decoder.Data.IsData()) {
      decoder.Tree->Fill();
      if (decoder.Data.Channel == emChannel) {
        emdata.push_back(decoder.Data.GetTdcData().front());
      }
    } else if (decoder.Data.IsFooter()) {
      std::cout << "end of spill " << decoder.Data.Spill << std::endl;
      emCount.back() = { decoder.Tree->GetEntries(), decoder.Data.DecodeEventMatchNumber(emdata) };
      if (emCount.back().second < 0) {
        emCount.back().second = nextEmCount < 0 ? nextEmCount : nextEmCount++;
      } else {
        nextEmCount = emCount.back().second + 1;
      }
      emdata.clear();
      emCount.push_back({ std::numeric_limits<Long64_t>::max(), nextEmCount });
    }
  }
  std::cout << "# of data record = " << count << std::endl;

  TBranch* emBranch = decoder.Data.AddEMBranch(decoder.Tree);
  for (Long64_t entry = 0, entries = decoder.Tree->GetEntries(), iem = 0; entry < entries; ++entry) {
    if (entry < emCount[iem].first || entry < emCount[++iem].first) {
      decoder.Data.EMCount = emCount[iem].second;
    } else {
      std::cout << "[error] unexpected scenario" << std::endl;
      return 1;
    }
    emBranch->Fill();
  }

  std::cout << "=== Write Objects" << std::endl;
  std::cout << decoder.Tree->GetName() << std::endl;
  decoder.Tree->Write();
  auto param = new TParameter<Double_t>("ClockFreq", clock);
  std::cout << param->GetName() << " = " << param->GetVal() << std::endl;
  param->Write();

  std::cout << "=== Close Files" << std::endl;
  ifile.close();
  ofile->Close();

  return 0;
}
