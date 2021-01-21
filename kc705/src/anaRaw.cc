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

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input" ,                "A root filename");
  args->AddOpt<std::string>("Output", 'o', "output", "A output filename", "");
  args->AddOpt             ("Help"  , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const std::string ofilename = args->GetValue("Output");

  std::string ofilename2;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    if (buff.EndsWith(".root")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofilename2 = buff + ".pdf";
  } else {
    ofilename2 = ofilename;
  }

  std::cout << "=== Open Input File" << std::endl;
  TFile* file = new TFile(ifilename.data());
  if (!file->IsOpen()) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  std::cout << "=== Get Tree" << std::endl;
  TTree* tree = dynamic_cast<TTree*>(file->Get("FCTTDC"));
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

  std::cout << "=== Set Branch Address" << std::endl;
  const std::string tdcName = Extinction::Kc705::Name;
  Extinction::Kc705::Kc705Data data;
  data.SetBranchAddress(tree);

  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);

  // Raw channel hit count
  TH1* hKc705MppcEntryByCh = new TH1D("hKc705MppcEntryByCh",
                                      Form("%s, MPPC Entries by Channel;"
                                           "%s Channel;"
                                           "Entries", tdcName.data(), tdcName.data()),
                                      64, 0.0 - 0.5, 64.0 - 0.5);

  TH2* hKc705MppcEntryVsCh = new TH2D("hKc705MppcEntryVsCh",
                                      Form("%s, MPPC Entries vs Channel;"
                                           "%s Channel;"
                                           "Entries", tdcName.data(), tdcName.data()),
                                      64, 0.0 - 0.5, 64.0 - 0.5,
                                      200, 0, entries);

  TH1* hKc705SubEntryByCh = new TH1D("hKc705SubEntryByCh",
                                     Form("%s, Sub Entries by Channel;"
                                          "%s Channel;"
                                          "Entries", tdcName.data(), tdcName.data()),
                                     24, 0.0 - 0.5, 24.0 - 0.5);

  TH2* hKc705SubEntryVsCh = new TH2D("hKc705SubEntryVsCh",
                                     Form("%s, Sub Entries vs Channel;"
                                          "%s Channel;"
                                          "Entries", tdcName.data(), tdcName.data()),
                                     24, 0.0 - 0.5, 24.0 - 0.5,
                                     200, 0, entries);

  std::cout << "=== Get Entry " << std::endl;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (entry % 1000000 == 0) {
      std::cout << ">> " << entry << std::endl;
    }

    tree->GetEntry(entry);

    for (auto&& mppcCh : data.GetMppcHitChannels()) {
      hKc705MppcEntryByCh->Fill(mppcCh);
      hKc705MppcEntryVsCh->Fill(mppcCh, entry);
    }
    for (auto&& subCh : data.GetSubHitChannels()) {
      hKc705SubEntryByCh->Fill(subCh);
      hKc705SubEntryVsCh->Fill(subCh, entry);
    }
  }

  std::cout << "=== Draw Hists" << std::endl;
  TCanvas::MakeDefCanvas();
  gPad->SetGrid();

  gPad->Print((ofilename2 + "[").data());

  // Raw Hit
  gPad->SetLogy(true);
  hKc705MppcEntryByCh->Draw();
  gPad->Print(ofilename2.data());
  gPad->SetLogy(false);

  hKc705MppcEntryVsCh->Draw("colz");
  gPad->Print(ofilename2.data());

  hKc705MppcEntryVsCh->Draw("box");
  gPad->Print(ofilename2.data());

  gPad->SetLogy(true);
  hKc705SubEntryByCh->Draw();
  gPad->Print(ofilename2.data());
  gPad->SetLogy(false);

  hKc705SubEntryVsCh->Draw("colz");
  gPad->Print(ofilename2.data());

  hKc705SubEntryVsCh->Draw("box");
  gPad->Print(ofilename2.data());

  gPad->Print((ofilename2 + "]").data());

  return 0;
}
