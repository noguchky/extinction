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
#include "TimelineCoincidence.hh"
#include "Fct.hh"
#include "TApplication.h"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                    "Set configure filename");
  args->AddArg<std::string>("Boards"      ,                    "Set comma separated board numbers");
  args->AddArg<std::string>("Input"       ,                    "Set comma separated root filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output"    , "Set prefix of output filename", "");
  args->AddOpt             ("Efficiency"  , 'e', "efficiency", "Execute efficiency analysis");
  args->AddOpt             ("Debug"       , 'd', "debug"     , "Debug mode");
  args->AddOpt             ("Help"        , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  if (args->IsSet("Debug")) {
    new TApplication("app", nullptr, nullptr);
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
  const auto effAna       = args->IsSet("Efficiency");

  std::string ofileprefix;
  if (ofilename.empty()) {
    TString buff = ifilenames.begin()->second.data();
    if (buff.EndsWith(".root")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofileprefix = buff;
  } else {
    TString buff = ofilename.data();
    if (buff.EndsWith(".pdf")) {
      buff.Replace(buff.Length() - 4, 4, "");
    }
    ofileprefix = buff;
  }

  const std::string ofilenameRoot          = ofileprefix + "_coin.root";
  const std::string ofilenamePdf           = ofileprefix + "_coin.pdf";
  // std::cout << "ofilenameRoot          " << ofilenameRoot          << std::endl;
  // std::cout << "ofilenamePdf           " << ofilenamePdf           << std::endl;

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  Extinction::Fct::ChannelMapWithBoard::Load(conf, boards);

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(520, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize histgram generator" << std::endl;
  Extinction::Fct::FctData defaultProvider;
  auto generator = new Extinction::Analyzer::TimelineCoincidence(&defaultProvider);

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

  std::map<Int_t, Double_t> mrSyncInterval;
  do {
    const std::string ifilename = conf->GetValue("MrSyncInterval");
    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[warning] MR Sync interval file was not opened, " << ifilename << std::endl;
      break;
    }

    std::string buff;
    Int_t board; Double_t interval;
    while (std::getline(ifile, buff)) {
      std::stringstream line(buff);
      if (line >> board >> interval) {
        mrSyncInterval[board] = interval;
      }
    }
  } while (false);

  {
    const std::string ifilename = conf->GetValue("Offset");
    if (generator->LoadOffset(ifilename) == 0) {
      std::cerr << "[warning] offset file was not found, " << ifilename << std::endl;
    }
  }

  generator->SetMrSyncOffset     (conf->GetValue <Double_t   >("MrSyncOffset"     ));
  generator->SetCoinTimeWidth    (conf->GetValue <Double_t   >("CoinTimeWidth"    ));
  generator->SetTimePerTdc       (                              timePerTdc         );
  generator->SetMrSyncInterval   (                              mrSyncInterval     );
  generator->SetCoincidenceTarget(conf->GetValues<Int_t      >("CoincidenceTarget"));

  Extinction::Analyzer::TimelineCoincidence::PlotsProfiles profile;
  profile.TimeInSpill   .NbinsX   = conf->GetValue<Double_t>("TimeInSpill.NbinsX" );
  profile.TimeInSpill   .Xmin     = conf->GetValue<Double_t>("TimeInSpill.Xmin"   );
  profile.TimeInSpill   .Xmax     = conf->GetValue<Double_t>("TimeInSpill.Xmax"   );
  profile.TimeInSync    .Xmin     = conf->GetValue<Double_t>("TimeInSync.Xmin"    );
  profile.TimeInSync    .Xmax     = conf->GetValue<Double_t>("TimeInSync.Xmax"    );
  profile.TimeInSync    .BinWidth = conf->GetValue<Double_t>("TimeInSync.BinWidth");
  profile.MrSyncInterval.Mean     = conf->GetValue<Double_t>("MrSyncInterval.Mean");
  profile.MrSyncInterval.Xmin     = conf->GetValue<Double_t>("MrSyncInterval.Xmin");
  profile.MrSyncInterval.Xmax     = conf->GetValue<Double_t>("MrSyncInterval.Xmax");
  profile.TimeDiff      .Xmin     = conf->GetValue<Double_t>("TimeDiff.Xmin"      );
  profile.TimeDiff      .Xmax     = conf->GetValue<Double_t>("TimeDiff.Xmax"      );
  generator->InitializePlots(profile);

  auto providers = Tron::Linq::From(boards)
    .ToMap([](Int_t board) { return board; },
           [](Int_t) -> Extinction::ITdcDataProvider* { return new Extinction::Fct::FctData(); });

  auto parser =
    [](const std::string& filename) {
      static const std::regex pattern(R"((\d+)-(\d+)-(\d+)_(\d+)-(\d+)-(\d+))");

      std::cmatch match;
      if (std::regex_search(filename.data(), match, pattern)) {
        const Int_t year   = Tron::String::Convert<Int_t>(match[1].str());
        const Int_t month  = Tron::String::Convert<Int_t>(match[2].str());
        const Int_t day    = Tron::String::Convert<Int_t>(match[3].str());
        const Int_t hour   = Tron::String::Convert<Int_t>(match[4].str());
        const Int_t minute = Tron::String::Convert<Int_t>(match[5].str());
        const Int_t second = Tron::String::Convert<Int_t>(match[6].str());
        return TDatime(year, month, day, hour, minute, second);
      }

      return TDatime(0U);
    };

  if (effAna) {
    std::cout << "--- Generate efficiency" << std::endl;
    generator->GenerateEfficiency(providers, ifilenames, "tree", parser);
  } else {
    std::cout << "--- Generate hists" << std::endl;
    generator->GeneratePlots(providers, ifilenames, "tree", parser);
  }

  generator->DrawPlots          (ofilenamePdf   );
  generator->WritePlots         (ofilenameRoot  );

  return 0;
}
