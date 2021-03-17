#include <iostream>
#include <map>
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
  args->AddArg<std::string>("ConfFilename",                "A configure filename");
  args->AddArg<std::string>("Boards"      ,                "A comma separated board numbers");
  args->AddArg<std::string>("Input"       ,                "A comma separated root filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output", "A output pdf filename", "");
  args->AddOpt             ("Help"        , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  for (Int_t i = 0; i < argc; ++i) {
    std::cout << argv[i] << "  ";
  }
  std::cout << std::endl;

  std::cout << "ConfFilename  " << args->GetValue("ConfFilename") << std::endl;
  std::cout << "Boards        " << args->GetValue("Boards"      ) << std::endl;
  std::cout << "Input         " << args->GetValue("Input"       ) << std::endl;
  std::cout << "Output        " << args->GetValue("Output"      ) << std::endl;
  
  const auto confFilename = args->GetValue("ConfFilename");
  const auto boards       =
    Tron::Linq::From(Tron::String::Split(args->GetValue("Boards"), ","))
    .Select([](const std::string& board) { return Tron::String::Convert<Int_t>(board); })
    .ToVector();
  const auto ifilenames   =
    Tron::Linq::From(Tron::String::Split(args->GetValue("Input" ), ","))
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
    ofileprefix = buff + "_hists";
  } else {
    TString buff = ofilename.data();
    if (buff.EndsWith(".pdf")) {
      buff.Replace(buff.Length() - 4, 4, "");
    }
    ofileprefix = buff;
  }

  const std::string ofilenameRoot   = ofileprefix + ".root";
  const std::string ofilenamePdf    = ofileprefix + ".pdf";
  const std::string ofilenameMrSync = ofileprefix + "_mrSync.dat";
  const std::string ofilenameOffset = ofileprefix + "_offset.dat";
  std::cout << "ofilenameRoot   " << ofilenameRoot   << std::endl;
  std::cout << "ofilenamePdf    " << ofilenamePdf    << std::endl;
  std::cout << "ofilenameMrSync " << ofilenameMrSync << std::endl;
  std::cout << "ofilenameOffset " << ofilenameOffset << std::endl;

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  const std::map<Int_t, Double_t> timePerTdc = Tron::Linq::From(boards)
    .ToMap([&](int board) { return board; },
           [&](int board) { return conf->GetValue<Double_t>(Form("TimePerTdc.%d", board)) * Extinction::nsec; });
  const std::map<Int_t, Double_t> mrSyncInterval = Tron::Linq::From(boards)
    .ToMap([&](int board) { return board; },
           [&](int board) { return conf->GetValue<Double_t>(Form("MrSyncInterval.%d", board)); });

  Extinction::Fct::ChannelMapWithBoard::Load(conf, boards);

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize monitor window" << std::endl;
  Extinction::Fct::FctData defaultProvider;
  auto generator = new Extinction::Analyzer::HistGenerator(&defaultProvider);

  if (generator->LoadOffset(conf->GetValue("Offset")) == 0) {
    // std::cerr << "[error] offset file was not found" << std::endl;
    // exit(1);
    std::cerr << "[warning] offset file was not found" << std::endl;
  }

  generator->SetHistoryWidth    (conf->GetValue<Double_t   >("HistoryWidth"));
  generator->SetCoinTimeWidth   (conf->GetValue<Double_t   >("CoinTimeWidth"));
  generator->SetReadBufferSize  (conf->GetValue<std::size_t>("ReadBufferSize"));
  generator->SetReadBufferMargin(conf->GetValue<std::size_t>("ReadBufferMargin"));
  generator->SetTimePerTdc      (timePerTdc);
  generator->SetMrSyncInterval  (mrSyncInterval);

  Extinction::Analyzer::HistGenerator::PlotsProfiles profile;
  profile.TimeInSpill   .NbinsX   = conf->GetValue<Double_t>("TimeInSpill.NbinsX");
  profile.TimeInSpill   .Xmin     = conf->GetValue<Double_t>("TimeInSpill.Xmin");
  profile.TimeInSpill   .Xmax     = conf->GetValue<Double_t>("TimeInSpill.Xmax");
  profile.TimeInSpill   .Unit     = conf->GetValue<Double_t>("TimeInSpill.Unit");
  profile.TimeInSync    .Xmin     = conf->GetValue<Double_t>("TimeInSync.Xmin");
  profile.TimeInSync    .Xmax     = conf->GetValue<Double_t>("TimeInSync.Xmax");
  profile.TimeInSync    .BinWidth = conf->GetValue<Double_t>("TimeInSync.BinWidth");
  profile.MrSyncInterval.Mean     = conf->GetValue<Double_t>("MrSyncInterval.Mean");
  profile.MrSyncInterval.Xmin     = conf->GetValue<Double_t>("MrSyncInterval.Xmin");
  profile.MrSyncInterval.Xmax     = conf->GetValue<Double_t>("MrSyncInterval.Xmax");
  profile.TimeDiff      .Xmin     = conf->GetValue<Double_t>("TimeDiff.Xmin");
  profile.TimeDiff      .Xmax     = conf->GetValue<Double_t>("TimeDiff.Xmax");
  generator->InitializePlots(profile);

  auto providers = Tron::Linq::From(boards)
    .ToMap([](Int_t board) { return board; },
           [](Int_t) -> Extinction::ITdcDataProvider* { return new Extinction::Fct::FctData(); });

  generator->UpdatePlots(providers, ifilenames, "tree");

  gStyle->SetOptStat(111111);
  generator->DrawPlots(ofilenamePdf);

  generator->WritePlots(ofilenameRoot);
  generator->WriteMrSyncInterval(ofilenameMrSync);
  generator->WriteCoinDiffs(ofilenameOffset);

  return 0;
}
