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
#include "Hul.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"       ,                "Set input hist filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output", "Set prefix of output filename", "");
  args->AddOpt             ("Help"        , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const auto ifilename = args->GetValue("Input" );
  const auto ofilename = args->GetValue("Output");

  std::string ofileprefix;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    if (buff.EndsWith("_hists.root")) {
      buff.Replace(buff.Length() - 11, 11, "");
    } else if (buff.EndsWith(".root")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofileprefix = buff;
  } else {
    TString buff = ofilename.data();
    if (buff.EndsWith("_hists.pdf")) {
      buff.Replace(buff.Length() - 10, 10, "");
    } else if (buff.EndsWith(".pdf")) {
      buff.Replace(buff.Length() - 4, 4, "");
    }
    ofileprefix = buff;
  }

  const std::string ofilenamePdf           = ofileprefix + "_hists.pdf";
  const std::string ofilenamePdf_Offset    = ofileprefix + "_offset.pdf";
  const std::string ofilenameSpill         = ofileprefix + "_spill.root";
  const std::string ofilenameMrSync        = ofileprefix + "_mrSync.dat";
  const std::string ofilenameOffset        = ofileprefix + "_offset.dat";
  // std::cout << "ofilenamePdf           " << ofilenamePdf           << std::endl;
  // std::cout << "ofilenamePdf_Offset    " << ofilenamePdf_Offset    << std::endl;
  // std::cout << "ofilenameSpill         " << ofilenameSpill         << std::endl;
  // std::cout << "ofilenameMrSync        " << ofilenameMrSync        << std::endl;
  // std::cout << "ofilenameOffset        " << ofilenameOffset        << std::endl;

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(520, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize histgram generator" << std::endl;
  Extinction::Hul::HulData defaultProvider;
  auto generator = new Extinction::Analyzer::HistGenerator(&defaultProvider);

  if (generator->ReadPlots(ifilename)) {
    return 1;
  }

  generator->CalcMrSyncInterval();
  generator->CalcTdcOffsets();

  generator->DrawPlots(ofilenamePdf, ofilenamePdf_Offset);

  generator->WriteMrSyncInterval(ofilenameMrSync);
  generator->WriteTdcOffsets    (ofilenameOffset);

  return 0;
}
