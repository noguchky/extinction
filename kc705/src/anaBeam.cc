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
#include "AnaBeam.hh"
#include "AnaTimeOffset.hh"
#include "AnaMrSyncInterval.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input" ,                "A root filename");
  args->AddOpt<std::string>("Output", 'o', "output", "A output pdf filenames", "");
  args->AddOpt             ("Help"  , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const std::string ofilename = args->GetValue("Output");

  std::string beamPdfFilename;
  std::string beamRootFilename;
  std::string timeOffsetPdfFilename;
  std::string timeOffsetDatFilename;
  std::string timeOffsetRootFilename;
  std::string mrSyncIntervalPdfFilename;
  std::string mrSyncIntervalDatFilename;
  std::string mrSyncIntervalRootFilename;

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
      beamPdfFilename = buff0 + ".pdf";
      beamRootFilename = buff0 + ".root";
      timeOffsetPdfFilename = buff0 + "_timeOffset%d.pdf";
      timeOffsetDatFilename = buff0 + "_timeOffset%d.dat";
      timeOffsetRootFilename = buff0 + "_timeOfffset%d.root";
      mrSyncIntervalPdfFilename = buff0 + "_mrSyncTdcInterval.pdf";
      mrSyncIntervalDatFilename = buff0 + "_mrSyncTdcInterval.dat";
      mrSyncIntervalRootFilename = buff0 + "_mrSyncTdcInterval.root";
    } else {
      beamPdfFilename = buff1 + ".pdf";
      beamRootFilename = buff1 + ".root";
      timeOffsetPdfFilename = buff1 + "_timeOffset%d.pdf";
      timeOffsetDatFilename = buff1 + "_timeOffset%d.dat";
      timeOffsetRootFilename = buff1 + "_timeOfffset%d.root";
      mrSyncIntervalPdfFilename = buff1 + "_mrSyncTdcInterval.pdf";
      mrSyncIntervalDatFilename = buff1 + "_mrSyncTdcInterval.dat";
      mrSyncIntervalRootFilename = buff1 + "_mrSyncTdcInterval.root";
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

  std::cout << "=== Analyze Beam" << std::endl;
  gStyle->SetOptStat(111111);
  auto results = Extinction::Analyzer::AnaBeam::Execute(tree, &data);
  {
    std::cout << "=== Draw Hists" << std::endl;
    results.Print(beamPdfFilename);

    std::cout << "=== Write Hists" << std::endl;
    if (beamRootFilename.size()) {
      std::cerr << "Info in <TFile::TFile>: root file " << beamRootFilename << " has been created" << std::endl;
      TFile* ofile = new TFile(beamRootFilename.data(), "RECREATE");
      if (ofile->IsOpen()) {
        results.Write(ofile);
        ofile->Close();
      }
      delete ofile;
    }
  }

  std::cout << "=== Analyze Time Offset" << std::endl;
  gStyle->SetOptStat(1111);
  for (std::size_t ch = 0; ch < Extinction::ExtinctionDetector::NofChannels; ++ch) {
    if (results.hExtTdcCoin[ch]->GetEntries()) {
      const std::string pdfFilename  = Form(timeOffsetPdfFilename .data(), ch);
      const std::string datFilename  = Form(timeOffsetDatFilename .data(), ch);
      const std::string rootFilename = Form(timeOffsetRootFilename.data(), ch);

      auto results2 = Extinction::Analyzer::AnaTimeOffset::Execute(results.hExtTdcCoin[ch], datFilename);
      if (results2.hExtTdcCoinY) {
        std::cout << "=== Draw Hists" << std::endl;
        results2.Print(pdfFilename);

        std::cout << "=== Write Hists" << std::endl;
        std::cerr << "Info in <TFile::TFile>: root file " << rootFilename << " has been created" << std::endl;
        TFile* ofile = new TFile(rootFilename.data(), "RECREATE");
        if (ofile->IsOpen()) {
          results2.Write(ofile);
          ofile->Close();
        }
        delete ofile;
      }
    }
  }

  std::cout << "=== Analyzer MR Sync Interval" << std::endl;
  gStyle->SetOptStat(111111);
  auto results3 = Extinction::Analyzer::AnaMrSyncInterval::Execute(results.hMrSyncTdcInterval2, mrSyncIntervalDatFilename);
  {
    std::cout << "=== Draw Hists" << std::endl;
    results3.Print(mrSyncIntervalPdfFilename);

    std::cout << "=== Write Hists" << std::endl;
    std::cerr << "Info in <TFile::TFile>: root file " << mrSyncIntervalRootFilename << " has been created" << std::endl;
    TFile* ofile = new TFile(mrSyncIntervalRootFilename.data(), "RECREATE");
    if (ofile->IsOpen()) {
      results3.Write(ofile);
      ofile->Close();
    }
    delete ofile;
  }

  return 0;
}
