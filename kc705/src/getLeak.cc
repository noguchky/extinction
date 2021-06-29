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
#include "ConfReader.hh"
#include "ArgReader.hh"
#include "Units.hh"
#include "Spill.hh"
#include "MargedReader.hh"
#include "HistGenerator.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("ConfFilename",                "Set configure filename");
  args->AddArg<std::string>("Input"       ,                "Set coin root filenames");
  args->AddOpt<std::string>("Output"      , 'o', "output", "Set prefix of output filename", "");
  args->AddOpt             ("Help"        , 'h', "help"  , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const auto confFilename = args->GetValue("ConfFilename");
  const auto ifilename    = args->GetValue("Input");
  const auto ofilename    = args->GetValue("Output");

  std::string ofileprefix;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    if (buff.EndsWith("_ctree.root")) {
      buff.Replace(buff.Length() - 11, 11, "");
    } else if (buff.EndsWith(".root")) {
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

  const std::string ofilenameRoot          = ofileprefix + "_leak.root";
  const std::string ofilenamePdf           = ofileprefix + "_leak.pdf";
  // std::cout << "ofilenameRoot          " << ofilenameRoot          << std::endl;
  // std::cout << "ofilenamePdf           " << ofilenamePdf           << std::endl;

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  std::cout << "--- Get bunch position" << std::endl;
  Double_t BunchCenters [Extinction::kNofBunches];
  // Double_t BunchWidths  [Extinction::kNofBunches];
  Double_t BunchMinEdges[Extinction::kNofBunches];
  Double_t BunchMaxEdges[Extinction::kNofBunches];

  {
    const std::string ifilename = conf->GetValue("BunchProfile");
    std::cout << ifilename << std::endl;

    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[error] bunch profile is not opened, " << ifilename << std::endl;
      return 1;
    }

    Int_t bunch;
    Double_t center, sigma;
    Long64_t minEdge, maxEdge;
    while (ifile >> bunch >> center >> sigma >> minEdge >> maxEdge) {
      BunchCenters [bunch] = center;
      // BunchWidths  [bunch] = sigma;
      BunchMinEdges[bunch] = minEdge;
      BunchMaxEdges[bunch] = maxEdge;
    }
  }

  auto isInBunch =
    [&] (Long64_t dtdc) {
      for (std::size_t bunch = 0; bunch < Extinction::kNofBunches; ++bunch) {
        if (!BunchCenters[bunch]) {
          return false;
        }

        if        (dtdc <  BunchMinEdges[bunch]) {
          return false;
        } else if (dtdc <= BunchMaxEdges[bunch]) {
          return true;
        }
      }
      return false;
    };

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(520, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Open input file" << std::endl;
  TFile* ifile = new TFile(ifilename.data());
  if (!ifile->IsOpen()) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  std::cout << "--- Get coincidence tree" << std::endl;
  TTree* itree = dynamic_cast<TTree*>(ifile->Get("ctree"));
  if (!itree) {
    std::cout << "[error] input tree is not found" << std::endl;
    return 1;
  }

  const Long64_t entries = itree->GetEntries();

  std::cout << "--- Set branch address" << std::endl;
  ULong64_t fDate          = 0;
  Int_t     fEMCount       = 0;
  UInt_t    fCoinTarget    = 0;
  Int_t     fMrSyncCount   = 0;
  Long64_t  fTdcFromMrSync = 0;
  Int_t     fBlurWidth     = 0;
  Int_t     fCoinWidth     = 0;
  Int_t     fNofHits       = 0;
  Int_t     fHitChannels      [Extinction::ExtinctionDetector::NofChannels] = { };
  // Int_t     fHitMrSyncCounts  [Extinction::ExtinctionDetector::NofChannels] = { };
  Long64_t  fHitTdcFromMrSyncs[Extinction::ExtinctionDetector::NofChannels] = { };
  Int_t     fHitTots          [Extinction::ExtinctionDetector::NofChannels] = { };

  itree->SetBranchAddress("date"    , &fDate             );
  itree->SetBranchAddress("emcount" , &fEMCount          );
  itree->SetBranchAddress("ctarget" , &fCoinTarget       );
  itree->SetBranchAddress("mscount" , &fMrSyncCount      );
  itree->SetBranchAddress("dtdc"    , &fTdcFromMrSync    );
  itree->SetBranchAddress("bwidth"  , &fBlurWidth        );
  itree->SetBranchAddress("cwidth"  , &fCoinWidth        );
  itree->SetBranchAddress("nofhits" , &fNofHits          );
  itree->SetBranchAddress("hitchs"  ,  fHitChannels      );
  // itree->SetBranchAddress("mscounts",  fHitMrSyncCounts  );
  itree->SetBranchAddress("dtdcs"   ,  fHitTdcFromMrSyncs);
  itree->SetBranchAddress("tots"    ,  fHitTots          );

  std::cout << "--- Open output file" << std::endl;
  TFile* ofile = new TFile(ofilenameRoot.data(), "RECREATE");

  std::cout << "--- Create output tree" << std::endl;
  TTree* otree = new TTree(itree->GetName(), itree->GetTitle());

  otree->Branch("date"    , &fDate             , "date"           "/l");
  otree->Branch("emcount" , &fEMCount          , "emcount"        "/I");
  otree->Branch("ctarget" , &fCoinTarget       , "ctarget"        "/i");
  otree->Branch("mscount" , &fMrSyncCount      , "mscount"        "/I");
  otree->Branch("dtdc"    , &fTdcFromMrSync    , "dtdc"           "/L");
  otree->Branch("bwidth"  , &fBlurWidth        , "bwidth"         "/I");
  otree->Branch("cwidth"  , &fCoinWidth        , "cwidth"         "/I");
  otree->Branch("nofhits" , &fNofHits          , "nofhits"        "/I");
  otree->Branch("hitchs"  ,  fHitChannels      , "hitchs""[nofhits]/I");
  // otree->Branch("mscounts",  fHitMrSyncCounts  , "mscounts[nofhits]/I");
  otree->Branch("dtdcs"   ,  fHitTdcFromMrSyncs, "dtdcs" "[nofhits]/L");
  otree->Branch("tots"    ,  fHitTots          , "tots"  "[nofhits]/I");

  std::cout << "--- Loop" << std::endl;
  for (Long64_t entry = 0; entry < entries; ++entry) {
    if (entry % 100000 == 0) {
      std::cout << ">> " << entry << "(" << otree->GetEntries() << ")" << std::endl;
    }

    if (itree->GetEntry(entry)) {
      if (!isInBunch(fTdcFromMrSync)) {
        otree->Fill();
      }
    }
  }

  std::cout << "--- Write tree" << std::endl;
  std::cout << otree->GetName() << "\t" << otree->GetEntries() << std::endl;

  ofile->cd();
  otree->Write();

  std::cout << "--- Close" << std::endl;
  ofile->Close();
  ifile->Close();

  return 0;
}
