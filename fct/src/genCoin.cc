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
#include "MargedReader.hh"
#include "Fct.hh"
#include "TApplication.h"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                      "Set configure filename");
  args->AddArg<std::string>("Boards"      ,                      "Set comma separated board numbers");
  args->AddArg<std::string>("Input"       ,                      "Set comma separated root filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output"      , "Set prefix of output filename", ""); 
  args->AddOpt             ("Timeline"    , 't', "timeline"    , "Draw timeline of coincidence");
  args->AddOpt<Int_t      >("MSCount"     , 'm', "mscount"     , "Draw timeline in mr sync count", "-1");
  args->AddOpt             ("Efficiency"  , 'e', "efficiency"  , "Execute efficiency analysis");
  args->AddOpt             ("Help"        , 'h', "help"        , "Show usage");

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
  const auto ofilename        = args->GetValue("Output");
  const auto drawCoinTimeline = args->IsSet("Timeline");
  const auto mscountSelection = args->GetValue<Int_t>("MSCount");
  const auto efficiency       = args->IsSet("Efficiency");

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

  const std::string ofilenameRoot     = ofileprefix + "_coin.root";
  const std::string ofilenamePdf      = ofileprefix + "_coin.pdf";
  const std::string ofilenameCoinTree = ofileprefix + "_ctree.root";

  if (drawCoinTimeline || mscountSelection > 0) {
    new TApplication("app", nullptr, nullptr);
  }
  
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

  std::cout << "--- Initialize providers" << std::endl;
  Extinction::Fct::FctData defaultProvider;

  auto providers = Tron::Linq::From(boards)
    .ToMap([](Int_t board) { return board; },
           [](Int_t) -> Extinction::ITdcDataProvider* { return new Extinction::Fct::FctData(); });

  std::cout << "--- Initialize marged reader" << std::endl;
  auto reader = new Extinction::Analyzer::MargedReader(&defaultProvider);
  reader->SetMrSyncTimeOffset(conf->GetValue<Double_t>("MrSyncTimeOffset"));

  {
    const std::string ifilename = conf->GetValue("TdcOffsets");
    if (reader->LoadTdcOffsets(ifilename) == 0) {
      std::cerr << "[warning] offset is empty" << std::endl;
    }
  }

  // for (auto&& ifilename : conf->GetValues("TdcOffsets.Add")) {
  //   if (reader->AddOffset(ifilename) == 0) {
  //     std::cerr << "[warning] offset is empty" << std::endl;
  //   }
  // }

  std::cout << "--- Initialize histogram profile" << std::endl;
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

  std::cout << "--- Initialize coincidence generator" << std::endl;
  auto generator = new Extinction::Analyzer::TimelineCoincidence(&defaultProvider);
  generator->SetCoinTimeWidth    (conf->GetValue <Double_t>("CoinTimeWidth"    ));
  generator->SetCoincidenceTarget(conf->GetValues<Int_t   >("CoincidenceTarget"));

  generator->InitializePlots(profile);
  generator->InitializeCoinTree(ofilenameCoinTree);

  if (efficiency) {
    std::cout << "--- Generate efficiency" << std::endl;
    generator->GenerateEfficiency(reader, providers, ifilenames, "tree");
  } else {
    std::cout << "--- Generate hists" << std::endl;

    if (reader->Open(providers, ifilenames, "tree")) {
      exit(1);
    }

    generator->GeneratePlots(reader, drawCoinTimeline, mscountSelection);

    reader->Close();
  }

  generator->DrawPlots(ofilenamePdf);

  generator->WritePlots(ofilenameRoot);
  generator->WriteCoinTree();

  return 0;
}
