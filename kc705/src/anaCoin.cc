#include <iostream>
#include <map>
#include "Units.hh"
#include "Kc705.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "ArgReader.hh"
#include "AnaCoin.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input" ,                "A root filename");
  args->AddArg<std::string>("Offset",                "A offset filename format");
  args->AddOpt<std::string>("Output", 'o', "output", "A output pdf filenames", "");
  args->AddOpt             ("Help"  , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const std::string ffilename = args->GetValue("Offset");
  const std::string ofilename = args->GetValue("Output");

  std::string coinPdfFilename;
  std::string coinRootFilename;
  std::string mtnsPdfFilename;
  {
    TString buff0 = ifilename.data();
    if (buff0.EndsWith(".root")) {
      buff0.Replace(buff0.Length() - 5, 5, "");
    }
    TString buff1 = ofilename.data();
    if (buff1.EndsWith(".pdf")) {
      buff1.Replace(buff1.Length() - 4, 4, "");
    }

    if (ofilename.empty()) {
      coinPdfFilename = buff0 + ".pdf";
      coinRootFilename = buff0 + ".root";
      mtnsPdfFilename = buff0 + "_mtns.pdf";
    } else {
      coinPdfFilename = buff1 + ".pdf";
      coinRootFilename = buff1 + ".root";
      mtnsPdfFilename = buff1 + "_mtns.pdf";
    }
  }

  std::cout << "=== Open Input File" << std::endl;
  TFile* file = new TFile(ifilename.data());
  if (!file->IsOpen()) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  std::cout << "=== Get Tree" << std::endl;
  TTree* tree = dynamic_cast<TTree*>(file->Get("tree"));
  if (!tree) {
    std::cout << "[error] input tree is not found" << std::endl;
    return 1;
  }

  const Long64_t entries = tree->GetEntries();
  std::cout << "entries = " << entries << std::endl;
  if (!entries) {
    std::cout << "[error] input tree is empty" << std::endl;
    return 1;
  }

  std::cout << "=== Set Style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);

  std::cout << "=== Set Branch Address" << std::endl;
  Extinction::Kc705::Kc705Data data;

  std::cout << "=== Get Offset" << std::endl;
  std::map<std::size_t, std::map<std::size_t, Extinction::Analyzer::AnaCoin::CoinInfo>> coinInfo;
  for (std::size_t extCh = 0; extCh < Extinction::ExtinctionDetector::NofChannels; ++extCh) {
    Extinction::Analyzer::AnaCoin::CoinInfo info;
    std::ifstream ffile(Form(ffilename.data(), extCh));
    if (ffile) {
      while (info.Read(ffile)) {
        coinInfo[extCh][info.Index] = info;
      }
    }
  }

  std::cout << "=== Analyze Coincidence" << std::endl;
  gStyle->SetOptStat(111111);
  auto results = Extinction::Analyzer::AnaCoin::Execute(tree, &data, coinInfo);
  {
    std::cout << "=== Draw Hists" << std::endl;
    results.Print(coinPdfFilename);
    results.PrintMtns(mtnsPdfFilename);

    std::cout << "=== Write Hists" << std::endl;
    if (coinRootFilename.size()) {
      std::cerr << "Info in <TFile::TFile>: root file " << coinRootFilename << " has been created" << std::endl;
      TFile* ofile = new TFile(coinRootFilename.data(), "RECREATE");
      if (ofile->IsOpen()) {
        results.Write(ofile);
        ofile->Close();
      }
      delete ofile;
    }
  }
  
  return 0;
}
