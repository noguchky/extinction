#include <iostream>
#include <map>
#include "Units.hh"
#include "Hul.hh"
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
#include "AnaBeam.hh"
#include "AnaTimeOffset.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input" ,                "A root filename");
  args->AddOpt<std::string>("Output", 'o', "output", "A output filenames (split with ':')", "");
  args->AddOpt             ("Help"  , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const std::string ofilename = args->GetValue("Output");

  std::string ofilename1; // For AnaBeam       (pdf)
  std::string ofilename2; // For AnaTimeOffset (pdf)
  std::string ofilename3; // For AnaTimeOffset (dat)
  if (!ofilename.empty()) {
    std::size_t pos = 0, prepos = 0;
    pos = ofilename.find(":", prepos);
    ofilename1 = ofilename.substr(0, pos);
    if (pos != std::string::npos) {
      prepos = pos + 1;
      pos = ofilename.find(":", prepos);
      ofilename2 = ofilename.substr(prepos, pos);
      if (pos != std::string::npos) {
        prepos = pos + 1;
        pos = ofilename.find(":", prepos);
        ofilename3 = ofilename.substr(prepos, pos);
      }
    }
  }

  if (ofilename1.empty() ||
      ofilename2.empty() ||
      ofilename3.empty()) {
    TString buff0 = ifilename.data();
    if (buff0.EndsWith(".root")) {
      buff0.Replace(buff0.Length() - 5, 5, "");
    }
    TString buff1 = ofilename1.data();
    if (buff1.EndsWith(".pdf")) {
      buff1.Replace(buff1.Length() - 4, 4, "");
    }
    TString buff2 = ofilename2.data();
    if (buff2.EndsWith(".pdf")) {
      buff2.Replace(buff2.Length() - 4, 4, "");
    }
    if (ofilename1.empty()) {
      ofilename1 = buff0 + ".pdf";
      if (ofilename2.empty()) {
        ofilename2 = buff0 + "_timeOffset%ld.pdf";
        ofilename3 = buff0 + "_timeOffset%ld.dat";
      } else if (ofilename3.empty()) {
        ofilename3 = buff2 + ".dat";
      }
    } else if (ofilename2.empty()) {
      ofilename2 = buff1 + "_timeOffset%ld.pdf";
      ofilename3 = buff1 + "_timeOffset%ld.dat";
    } else if (ofilename3.empty()) {
      ofilename3 = buff2 + ".dat";
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

  std::cout << "=== Get Clock Frequency" << std::endl;
  TParameter<Double_t>* clockFreq = dynamic_cast<TParameter<Double_t>*>(file->Get("ClockFreq"));
  if (!clockFreq) {
    std::cout << "[error] parameter of clock frequency was not found" << std::endl;
    return 1;
  } else if (!clockFreq->GetVal()) {
    std::cout << "[warning] clock frequency is 0, it set to 1.0 GHz" << std::endl;
    clockFreq->SetVal(1.0 * Extinction::GHz);
  }
  std::cout << "clock freq = " << clockFreq->GetVal() << std::endl;

  std::cout << "=== Set Branch Address" << std::endl;
  Extinction::Hul::HulData data;
  data.ClockFreq = clockFreq->GetVal();

  std::cout << "=== Analyze Beam" << std::endl;
  gStyle->SetOptStat(111111);
  auto results = Extinction::Analyzer::AnaBeam::Execute(tree, &data);

  std::cout << "=== Draw Hists" << std::endl;
  results.Print(ofilename1, data.GetTimePerTdc());

  std::cout << "=== Analyze Time Offset" << std::endl;
  gStyle->SetOptStat(1111);
  for (std::size_t ch = 0; ch < Extinction::ExtinctionDetector::NofChannels; ++ch) {
    if (results.hExtTdcCoin[ch]->GetEntries()) {
      const std::string ofilename2_2 = Form(ofilename2.data(), ch);
      const std::string ofilename3_2 = Form(ofilename3.data(), ch);

      auto results2 = Extinction::Analyzer::AnaTimeOffset::Execute(results.hExtTdcCoin[ch], ofilename3_2);
      if (!results2.hExtTdcCoinY) {
        continue;
      }

      std::cout << "=== Draw Hists" << std::endl;
      results2.Print(ofilename2_2);
    }
  }
  
  return 0;
}
