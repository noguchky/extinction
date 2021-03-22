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
#include "HistGenerator.hh"
#include "Fct.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                "Set configure filename");
  args->AddArg<std::string>("Boards"      ,                "Set comma separated board numbers");
  args->AddArg<std::string>("Input"       ,                "Set comma separated root filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output", "Set prefix of output filename", "");
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

  const std::string ofilenameRoot   = ofileprefix + "_hists.root";
  const std::string ofilenamePdf    = ofileprefix + "_hists.pdf";
  const std::string ofilenameSpill  = ofileprefix + "_spill.root";
  const std::string ofilenameTPT    = ofileprefix + "_timePerTdc.dat";
  const std::string ofilenameMrSync = ofileprefix + "_mrSync.dat";
  const std::string ofilenameBunch  = ofileprefix + "_bunch.dat";
  const std::string ofilenameOffset = ofileprefix + "_offset.dat";
  // std::cout << "ofilenameRoot   " << ofilenameRoot   << std::endl;
  // std::cout << "ofilenamePdf    " << ofilenamePdf    << std::endl;
  // std::cout << "ofilenameSpill  " << ofilenameSpill  << std::endl;
  // std::cout << "ofilenameTPT    " << ofilenameTPT    << std::endl;
  // std::cout << "ofilenameMrSync " << ofilenameMrSync << std::endl;
  // std::cout << "ofilenameBunch  " << ofilenameBunch  << std::endl;
  // std::cout << "ofilenameOffset " << ofilenameOffset << std::endl;

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
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize histgram generator" << std::endl;
  Extinction::Fct::FctData defaultProvider;
  auto generator = new Extinction::Analyzer::HistGenerator(&defaultProvider);

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

  Double_t bunchCenters[Extinction::SpillData::kNofBunches] = { 0 };
  Double_t bunchWidths [Extinction::SpillData::kNofBunches] = { 0 };
  do {
    const std::string ifilename = conf->GetValue("BunchProfile");
    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[warning] bunch profile file was not opened, " << ifilename << std::endl;
      break;
    }

    std::string buff;
    Int_t bunch; Double_t center, width;
    while (std::getline(ifile, buff)) {
      std::stringstream line(buff);
      if (line >> bunch >> buff >> center >> buff >> width) {
        bunchCenters[bunch] = center * Extinction::nsec;
     // bunchWidths [bunch] = width  * Extinction::nsec;
        bunchWidths [bunch] = 250.0  * Extinction::nsec;
      }
    }
  } while (false);

  {
    const std::string ifilename = conf->GetValue("Offset");
    if (generator->LoadOffset(ifilename) == 0) {
      std::cerr << "[warning] offset file was not found, " << ifilename << std::endl;
    }
  }

  generator->SetCyclicCoincidence(conf->GetValue<Int_t      >("CyclicCoincidence"));
  generator->SetHistoryWidth     (conf->GetValue<Double_t   >("HistoryWidth"     ));
  generator->SetCoinTimeWidth    (conf->GetValue<Double_t   >("CoinTimeWidth"    ));
  generator->SetReadBufferSize   (conf->GetValue<std::size_t>("ReadBufferSize"   ));
  generator->SetReadBufferMargin (conf->GetValue<std::size_t>("ReadBufferMargin" ));
  generator->SetMrSyncRefInterval(conf->GetValue<Double_t   >("MrSyncRefInterval"));
  generator->SetMrSyncRefSize    (conf->GetValue<std::size_t>("MrSyncRefSize"    ));
  generator->SetTimePerTdc       (                             timePerTdc         );
  generator->SetMrSyncInterval   (                             mrSyncInterval     );
  generator->SetBunchCenters     (                             bunchCenters       );
  generator->SetBunchWidths      (                             bunchWidths        );

  Extinction::Analyzer::HistGenerator::PlotsProfiles profile;
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

  generator->InitializeSpillSummary(ofilenameSpill);

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

  std::cout << "--- Generate hists" << std::endl;
  generator->GeneratePlots(providers, ifilenames, "tree", "", parser);

  generator->DrawPlots          (ofilenamePdf   );

  generator->WritePlots         (ofilenameRoot  );
  generator->WriteSpillSummary  (               );
  generator->WriteTimePerTdc    (ofilenameTPT   );
  generator->WriteMrSyncInterval(ofilenameMrSync);
  generator->WriteBunchProfile  (ofilenameBunch );
  generator->WriteOffset        (ofilenameOffset);

  return 0;
}
