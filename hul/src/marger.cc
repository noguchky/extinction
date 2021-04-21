#include <iostream>
#include <fstream>
#include <map>
#include <regex>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "Linq.hh"
#include "ArgReader.hh"
#include "Units.hh"
#include "TreeMarger.hh"
#include "Hul.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                "Set configure filename");
  args->AddArg<std::string>("Boards"      ,                "Set comma separated board numbers");
  args->AddArg<std::string>("Input"       ,                "Set comma separated root filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output", "Set output filename", "");
  args->AddOpt             ("Help"        , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const auto confFilename = args->GetValue("ConfFilename");
  const auto boards       = Tron::Linq::From(Tron::String::Split(args->GetValue("Boards"), ","))
    .Select([](const std::string& board) { return Tron::String::Convert<Int_t>(board); })
    .ToVector();
  const auto ifilenames   = Tron::Linq::From(Tron::String::Split(args->GetValue("Input" ), ","))
    .Join (boards.begin())
    .ToMap([&](std::pair<std::string, Int_t> pair) { return pair.second; },
           [&](std::pair<std::string, Int_t> pair) { return pair.first;  });
  const auto ofilename    = args->GetValue("Output");

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  Extinction::Hul::ChannelMapWithBoard::Load(conf, boards);

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize tree marger" << std::endl;
  Extinction::Hul::HulData defaultProvider;
  auto marger = new Extinction::Analyzer::TreeMarger(&defaultProvider);

  std::map<Int_t, Double_t> timePerTdc;
  do {
    const std::string ifilename = conf->GetValue("TimePerTdc");
    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[warning] time resolution file was not opened, " << ifilename << std::endl;
      break;
    }

    std::string buff;
    Int_t board; Double_t resolution;
    while (std::getline(ifile, buff)) {
      std::stringstream line(buff);
      if (line >> board >> resolution) {
        timePerTdc[board] = resolution * Extinction::nsec;
      }
    }
  } while (false);

  marger->SetReadBufferSize   (conf->GetValue<std::size_t>("ReadBufferSize"   ));
  marger->SetReadBufferMargin (conf->GetValue<std::size_t>("ReadBufferMargin" ));
  marger->SetTimePerTdc       (timePerTdc                                      );

  auto providers = Tron::Linq::From(boards)
    .ToMap([](Int_t board) { return board; },
           [](Int_t) -> Extinction::ITdcDataProvider* { return new Extinction::Hul::HulData(); });

  std::cout << "--- Output tree" << std::endl;

  return 0;
}
