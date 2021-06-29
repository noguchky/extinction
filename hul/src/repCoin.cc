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
#include "Hul.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                    "Set configure filename");
  args->AddArg<std::string>("Input"       ,                    "Set input hist filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output"    , "Set prefix of output filename", "");
  args->AddOpt             ("Efficiency"  , 'e', "efficiency", "Execute efficiency analysis");
  args->AddOpt             ("Keyword"     , 'k', "keyword"   , "Set keyword for filename");
  args->AddOpt             ("Help"        , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const auto confFilename = args->GetValue("ConfFilename");
  const auto ifilename    = args->GetValue("Input" );
  const auto ofilename    = args->GetValue("Output");
  const auto efficiency   = args->IsSet("Efficiency");
  const auto keyword      = args->GetValue("Keyword");
  const std::string keywordSuffix = keyword.empty() ? "" : ("_" + keyword);

  std::string ofileprefix;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    const std::string coinSuffix = "_coin"       + keywordSuffix + ".root";
    const std::string effSuffix  = "_efficiency" + keywordSuffix + ".root";
    if       (!efficiency && buff.EndsWith(coinSuffix.data())) {
      buff.Replace(buff.Length() - coinSuffix.size(), coinSuffix.size(), "");
    } else if (efficiency && buff.EndsWith( effSuffix.data())) {
      buff.Replace(buff.Length() -  effSuffix.size(),  effSuffix.size(), "");
    } else if (buff.EndsWith(".root")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofileprefix = buff;
  } else {
    TString buff = ofilename.data();
    const std::string coinSuffix = "_coin"       + keywordSuffix + ".pdf";
    const std::string effSuffix  = "_efficiency" + keywordSuffix + ".pdf";
    if       (!efficiency && buff.EndsWith(coinSuffix.data())) {
      buff.Replace(buff.Length() - coinSuffix.size(), coinSuffix.size(), "");
    } else if (efficiency && buff.EndsWith( effSuffix.data())) {
      buff.Replace(buff.Length() -  effSuffix.size(),  effSuffix.size(), "");
    } else if (buff.EndsWith(".pdf")) {
      buff.Replace(buff.Length() - 4, 4, "");
    }
    ofileprefix = buff;
  }

  const std::string ofilenamePdf             = ofileprefix + "_coin"       + keywordSuffix + ".pdf";
  const std::string ofilenamePdf_Efficiency  = ofileprefix + "_efficiency" + keywordSuffix + ".pdf";
  const std::string ofilenameBunch           = ofileprefix + "_cbunch"     + keywordSuffix + ".dat";
  // std::cout << "ofilenamePdf             " << ofilenamePdf             << std::endl;
  // std::cout << "ofilenamePdf_Efficiency  " << ofilenamePdf_Efficiency  << std::endl;
  // std::cout << "ofilenameBunch           " << ofilenameBunch           << std::endl;

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(520, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize histgram generator" << std::endl;
  Extinction::Hul::HulData defaultProvider;
  auto generator = new Extinction::Analyzer::TimelineCoincidence(&defaultProvider);
  generator->SetBunchEdgeMargin(conf->GetValue<Double_t>("BunchEdgeMargin"));

  if (generator->ReadPlots(ifilename)) {
    return 1;
  }

  if (efficiency) {
    generator->DrawPlots(ofilenamePdf);

  } else {
    generator->CalcBunchProfile();

    generator->DrawPlots(ofilenamePdf);

    generator->WriteBunchProfile(ofilenameBunch);
  }

  return 0;
}
