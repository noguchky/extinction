#include <iostream>
#include <map>
#include "Units.hh"
#include "Fct.hh"
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

  std::cout << "=== Set Branch Address" << std::endl;
  const std::string tdcName = Extinction::Fct::Name;
  Extinction::Fct::FctData data;
  data.SetBranchAddress(tree);

  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);

  // Raw channel hit count
  TH1* hFctEntryByCh = new TH1D("hFctEntryByCh",
                                Form("%s, Entries by Channel;"
                                     "%s Channel;"
                                     "Entries", tdcName.data(), tdcName.data()),
                                32, 0.0 - 0.5, 32.0 - 0.5);

  TH2* hFctEntryVsCh = new TH2D("hFctEntryVsCh",
                                Form("%s, Entries vs Channel;"
                                     "%s Channel;"
                                     "Entries", tdcName.data(), tdcName.data()),
                                32, 0.0 - 0.5, 32.0 - 0.5,
                                200, 0, entries);

  TH1** hTdcInSpill = new TH1*[32];
  TH1** hTdcInSync  = new TH1*[32];
  TH2** hTdcMount   = new TH2*[32];
  for (std::size_t ch = 0; ch < 32; ++ch) {
    hTdcInSpill[ch] = new TH1D(Form("hTdcInSpill_%03lu", ch),
                               Form("hTdcUnSync_%03lu;ms", ch),
                               170, -100, 1600);
    hTdcInSync [ch] = new TH1D(Form("hTdcInSync_%03lu", ch),
                               Form("hTdcUnSync_%03lu;TDC", ch),
                               200, -100, 900);
    hTdcMount  [ch] = new TH2D(Form("hTdcMount_%03lu", ch),
                               Form("hTdcUnSync_%03lu;TDC;ms", ch),
                               200, -100, 900,
                               170, -100, 1600);
  }

  const Int_t mrSyncCh = 26;

  Long64_t syncTdc = 0;
  
  std::cout << "=== Get Entry " << std::endl;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (entry % 1000000 == 0) {
      std::cout << ">> " << entry << std::endl;
    }

    tree->GetEntry(entry);

    hFctEntryByCh->Fill(data.Channel);
    hFctEntryVsCh->Fill(data.Channel, entry);

    if (data.Channel == mrSyncCh) {
      syncTdc = data.Tdc;
    }
    hTdcInSpill[data.Channel]->Fill(data.GetTime() / Extinction::msec);
    hTdcInSync [data.Channel]->Fill(data.Tdc - syncTdc);
    hTdcMount  [data.Channel]->Fill(data.Tdc - syncTdc,
                                    data.GetTime() / Extinction::msec);
  }

  std::cout << "=== Draw Hists" << std::endl;
  TCanvas::MakeDefCanvas();
  gPad->SetGrid();

  gPad->Print((ofilename2 + "[").data());

  // Raw Hit
  gPad->SetLogy(true);
  hFctEntryByCh->Draw();
  gPad->Print(ofilename2.data());
  gPad->SetLogy(false);

  hFctEntryVsCh->Draw("colz");
  gPad->Print(ofilename2.data());

  hFctEntryVsCh->Draw("box");
  gPad->Print(ofilename2.data());

  for (std::size_t ch = 0; ch < 32; ++ch) {
    hTdcInSpill[ch]->Draw();
    gPad->Print(ofilename2.data());

    hTdcInSync [ch]->Draw();
    gPad->Print(ofilename2.data());

    hTdcMount  [ch]->Draw("colz");
    gPad->Print(ofilename2.data());
  }

  gPad->Print((ofilename2 + "]").data());

  return 0;
}
