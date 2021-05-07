#ifndef Extinction_TreeMarger_hh
#define Extinction_TreeMarger_hh

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TDatime.h"
#include "TF1.h"
#include "TParameter.h"

#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "Spill.hh"

#include "Math.hh"
#include "Linq.hh"
#include "String.hh"
#include "ObjectHelper.hh"

namespace Extinction {

  namespace Analyzer {

    class TreeMarger {
    private:
      ITdcDataProvider*            fProvider               = nullptr;
      TdcData                      fTdcData;

      TFile*                       fMargedFile             = nullptr;
      TTree*                       fMargedTree             = nullptr;

      std::map<Int_t, Double_t>    fStdTimePerTdc;

      std::size_t                  fBufferSize             = 5000;
      std::size_t                  fBufferMargin           =  100;

      std::map<ULong64_t, TdcData> fTdcBuffer;

    public:
      TreeMarger(ITdcDataProvider* provider);
      ~TreeMarger();

      void                 SetTimePerTdc(const std::map<Int_t, Double_t>& map);

      inline void          SetReadBufferSize(std::size_t size) {
        std::cout << "SetReadBufferSize ... " << size << std::endl;
        fBufferSize = size;
      }
      inline std::size_t   GetReadBufferSize() const { return fBufferSize; }

      inline void          SetReadBufferMargin(std::size_t size) {
        std::cout << "SetReadBufferMargin ... " << size << std::endl;
        fBufferMargin = size;
      }
      inline std::size_t   GetReadBufferMargin() const { return fBufferMargin; }

      Int_t                ReadPlots(const std::string& ifilename);
      Int_t                InitializeTree(const std::string& filename,
                                          const std::string& treename);
      void                 WriteTree();

      Int_t                MargeTree(std::map<Int_t, ITdcDataProvider*> providers,
                                     const std::map<Int_t, std::string>& ifilenames,
                                     const std::string& treename);

    private:
      void                 ClearLastSpill();

      template <typename T, typename V>
      inline bool          IsAllOfSecondsTrue(const std::map<T, V> map) {
        return std::all_of(map.begin(), map.end(), [](std::pair<T, V> p) { return p.second; });
      }
      template <typename T, typename V>
      inline void          FillToSeconds(std::map<T, V>* map, V value) {
        for (auto&& pair : *map) {
          pair.second = value;
        }
      }
    };

    TreeMarger::TreeMarger(ITdcDataProvider* provider)
      : fProvider(provider) {
    }

    TreeMarger::~TreeMarger() {
      WriteTree();
    }

    void TreeMarger::SetTimePerTdc(const std::map<Int_t, Double_t>& map) {
      std::cout << "SetTimePerTdc" << std::endl;
      for (auto&& pair : map) {
        std::cout << pair.first << "\t" << pair.second << std::endl;
      }
      fStdTimePerTdc = map;
    }

    Int_t TreeMarger::InitializeTree(const std::string& filename,
                                     const std::string& treename) {
      std::cout << "Initialize tree" << std::endl;
      fMargedFile = new TFile(filename.data(), "RECREATE");
      if (!fMargedFile->IsOpen()) {
        std::cout << "[error] output file is not opened, " << filename << std::endl;
        delete fMargedFile; fMargedFile = nullptr;
        return 1;
      }

      fMargedTree = new TTree(treename.data(), "");
      fTdcData.CreateBranch(fMargedTree);
      fMargedTree->AutoSave();
      return 0;
    }

    void TreeMarger::WriteTree() {
      if (fMargedFile && fMargedTree) {
        fMargedFile->cd();
        fMargedTree->Write();
        fMargedFile->Close();
        delete fMargedFile; fMargedFile = nullptr; fMargedTree = nullptr;
      }
    }

    void TreeMarger::ClearLastSpill() {
      fTdcBuffer.clear();
    }

    Int_t TreeMarger::MargeTree(std::map<Int_t, ITdcDataProvider*> providers,
                                const std::map<int, std::string>& ifilenames,
                                const std::string& treename) {
      const clock_t startClock = clock();

      std::cout << "Initialize decoder" << std::endl;
      for (auto&& pair : ifilenames) {
        const Int_t board = pair.first;
        if (fStdTimePerTdc[board] == 0) {
          providers[board]->SetTimePerTdc(fProvider->GetTimePerTdc());
        } else {
          providers[board]->SetTimePerTdc(fStdTimePerTdc[board]);
        }
      }

      std::cout << "Initialzie history" << std::endl;
      ClearLastSpill();

      std::cout << "Open file" << std::endl;
      std::map<Int_t, TFile*>   ifiles;
      std::map<Int_t, TTree*>   itrees;
      std::map<Int_t, Long64_t> entrieses;
      for (auto&& pair : ifilenames) {
        const Int_t       board     = pair.first;
        const std::string ifilename = pair.second;
        std::cout << " - " << ifilename << std::endl;

        ifiles[board] = new TFile(ifilename.data(), "READ");
        if (!ifiles[board]->IsOpen()) {
          std::cerr << "[error] input file is not opened, " << ifilename << std::endl;
          return 1;
        }

        itrees[board] = dynamic_cast<TTree*>(ifiles[board]->Get(treename.data()));
        if (!itrees[board]) {
          std::cerr << "[error] input tree is not found, " << ifilename << " - " << treename << std::endl;
          return 1;
        }

        entrieses[board] = itrees[board]->GetEntries();
        std::cout << " + " << entrieses[board] << std::endl;

        providers[board]->SetBranchAddress(itrees[board]);
      }
      std::cout << " = " << Tron::Linq::From(entrieses).Sum([](std::pair<Int_t, Long64_t> p) { return p.second; }) << std::endl;

      {
        Int_t targetBoard = 0;
        std::map<Int_t, Long64_t>             entries;
        std::map<Int_t, Long64_t>             lastSpills;
        std::map<Int_t, Bool_t>               spillEnded;
        std::map<Int_t, Bool_t>               fileEnded;
        std::map<Int_t, ULong64_t>            lastTdcTags;
        std::map<Int_t, std::vector<TdcData>> firstData;
        for (auto&& pair : ifilenames) {
          targetBoard = pair.first;
          entries         [pair.first] = 0;
          lastSpills      [pair.first] = -1;
          spillEnded      [pair.first] = false;
          fileEnded       [pair.first] = false;
          lastTdcTags     [pair.first] = 0;
          firstData       [pair.first] = { };
        }
        // std::cout << "targetBoard = " << targetBoard << std::endl; 

        Long64_t lastEMCount = -1;
        for (std::size_t count = 0UL;;) {
          {
            ITdcDataProvider* provider   = providers  [targetBoard];
            ULong64_t&        lastTdcTag = lastTdcTags[targetBoard];
            Long64_t&         lastSpill  = lastSpills [targetBoard];
            TTree*            itree      = itrees     [targetBoard];
            Long64_t&         entry      = entries    [targetBoard];
            const Long64_t    ientries   = entrieses  [targetBoard];
            for (std::size_t ibuf = 0; ibuf < fBufferSize; ++ibuf, ++count, ++entry) {
              if (count % 1000000UL == 0) {
                std::cout << ">> " << count << std::endl;
              }

              if (entry >= ientries) {
                std::cout << "[info] detect file end @ " << targetBoard << std::endl;
                fileEnded[targetBoard] = true;
                break;

              } else if (!itree->GetEntry(entry)) {
                std::cout << "[info] detect file end @ " << targetBoard << " (TTree::GetEntry)" << std::endl;
                fileEnded[targetBoard] = true;
                break;

              }

              if (!provider->IsData()) {
                continue;

              } else if (lastSpill != -1 && lastSpill != provider->GetSpill()) {
                std::cout << "[info] detect spill end @ " << targetBoard << std::endl;
                lastSpill               = provider->GetSpill();
                spillEnded[targetBoard] = true;
                firstData [targetBoard] = provider->GetTdcData(targetBoard);
                break;

              } else {
                // std::cout << "[info] detect data @ " << targetBoard << std::endl;
                if (lastSpill != -1 && lastEMCount != provider->GetEMCount()) {
                  std::cout << "[error] conflict EMCount, " << lastEMCount << " <--> " << provider->GetEMCount() << std::endl;
                  exit(1);
                }
                lastSpill   = provider->GetSpill();
                lastEMCount = provider->GetEMCount();
                const std::vector<TdcData> tdcData = provider->GetTdcData(targetBoard);
                for (auto&& data : tdcData) {
                  fTdcBuffer[data.GetTdcTag()] = data;
                }
                if (ibuf < fBufferSize - fBufferMargin) {
                  lastTdcTag = tdcData.back().GetTdcTag();
                }

              }
            }
          }

          {
            Int_t boardNotLoadedYet = -1;
            for (auto&& pair : lastTdcTags) {
              const Int_t board = pair.first;
              if (pair.second == 0 &&
                  !(spillEnded[board] || fileEnded[board])) {
                boardNotLoadedYet = board;
                break;
              }
            }
            if (boardNotLoadedYet >= 0) {
              // std::cout << "[info] there is a board not loaded yet @ " << boardNotLoadedYet << std::endl;
              targetBoard = boardNotLoadedYet;
              // std::cout << "targetBoard = " << targetBoard << std::endl; 
              continue;
            }
          }

          for (auto&& itr = fTdcBuffer.begin(); itr != fTdcBuffer.end(); itr = fTdcBuffer.begin()) {
            fTdcData = itr->second;
            fMargedTree->Fill();
            fTdcBuffer.erase(itr);
          }

          if (IsAllOfSecondsTrue(spillEnded) || IsAllOfSecondsTrue(fileEnded)) {
            std::cout << "[info] end of spill" << std::endl;
            for (auto&& pair : firstData) {
              auto& board      = pair.first;
              auto& tdcData    = pair.second;
              auto& lastTdcTag = lastTdcTags[board];
              for (auto&& data : tdcData) {
                fTdcBuffer[data.GetTdcTag()] = data;
                lastTdcTag = data.GetTdcTag();
              }
              tdcData.clear();
            }

            if (IsAllOfSecondsTrue(fileEnded)) {
              break;
            }

            // Clear buffers
            ClearLastSpill();
            FillToSeconds(&spillEnded , false);
            FillToSeconds(&lastTdcTags, 0ULL );
          }
        }
      }

      std::cout << "Close files" << std::endl;
      for (auto&& pair : ifiles) {
        pair.second->Close();
      }

      const clock_t stopClock = clock();
      std::cout << "time: " << (double)(stopClock - startClock) / CLOCKS_PER_SEC << " sec\n";

      return 0;
    }

  }

}

#endif
