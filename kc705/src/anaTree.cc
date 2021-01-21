#include <iostream>
#include <map>
#include "Units.hh"
#include "Kc705.hh"
#include "TApplication.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "ArgReader.hh"
#include "AnaTree.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input" ,                "A root filename");
  args->AddOpt<Long64_t   >("Start" , 's', "start" , "A entry starts from", "0");
  args->AddOpt             ("Help"  , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename = args->GetValue("Input");
  const Long64_t    start     = args->GetValue<Long64_t>("Start");

  std::cout << "=== Start Application" << std::endl;
  TApplication* app = new TApplication("anaTree", nullptr, nullptr);
  
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

  std::cout << "=== Set Branch Address" << std::endl;
  Extinction::Kc705::Kc705Data data;

  std::cout << "=== Analyze Beam" << std::endl;
  gStyle->SetOptStat(110011);
  Extinction::Analyzer::AnaTree::Execute(tree, &data, start);
  
  // Finalize
  app->Terminate();
  delete app;
  app = nullptr;

  return 0;
}
