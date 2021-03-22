#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TSystem.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "Fct.hh"
#include "Units.hh"
#include "Detector.hh"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                    "Set configure filename");
  args->AddOpt<Int_t      >("Board"       , 'b', "board"     , "Set used board", "0");
  args->AddOpt             ("Help"        , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string confFilename  = args->GetValue("ConfFilename");
  const Int_t       board         = args->GetValue<Int_t>("Board");

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  Extinction::Fct::ChannelMapWithBoard::Load(conf, { board });

  std::cout << "=== Initialize hitmap" << std::endl;
  TGraph* gHitCh = new TGraph();

  TH2*    hExtHitMap     = Extinction::ExtinctionDetector::CreateHitMap("hExtHitMap");
  TList*  lExtBorderLine = Extinction::ExtinctionDetector::CreateBorderLine();

  TCanvas* c1 = new TCanvas();
  gHitCh->Draw("AP");

  TCanvas* c2 = new TCanvas();
  hExtHitMap->Draw("colz");
  lExtBorderLine->Draw();

  std::cout << "=== Decode" << std::endl;
  Extinction::Fct::Decoder decoder;

  std::size_t count = 0UL;
  for (; decoder.Read(std::cin); ++count) {
    if (count % 100000 == 0) {
      std::cout << ">> " << count << std::endl;
    }
    if (decoder.Data.IsData()) {
      gHitCh->SetPoint(count & 0x7FF, count, decoder.Data.Channel);
      for (auto&& tdcdata : decoder.Data.GetTdcData(board)) {
        if (Extinction::ExtinctionDetector::Contains(tdcdata.Channel)) {
          Extinction::ExtinctionDetector::Fill(hExtHitMap, Extinction::ExtinctionDetector::GetChannel(tdcdata.Channel));
        }
      }
    } else if (decoder.Data.IsFooter()) {
      std::cout << "end of spill " << decoder.Data.Spill << std::endl;
      c1->Modified(); c1->Update();
      c2->Modified(); c2->Update();
      gSystem->ProcessEvents();
      hExtHitMap->Clear();
    }
  }

  return 0;
}
