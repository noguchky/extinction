#ifndef Extinction_MargedReader_hh
#define Extinction_MargedReader_hh

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TDatime.h"
#include "TF1.h"
#include "TParameter.h"
#include "TEfficiency.h"

#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"

#include "Linq.hh"
#include "String.hh"

namespace Extinction {

  namespace Analyzer {

    class MargedReader {
    public:
      using TdcOffsets_t = std::map<std::size_t/*globalChannel*/, Long64_t>;
      
    private:
      ITdcDataProvider*                   fProvider               = nullptr;

      BoardMap_t<ITdcDataProvider*>       fProviders;
      BoardMap_t<TFile*>                  fIfiles;
      BoardMap_t<TTree*>                  fItrees;
      BoardMap_t<Long64_t>                fEntrieses;
      BoardMap_t<Long64_t>                fEntries;
      BoardMap_t<Bool_t>                  fSpillEnded;
      BoardMap_t<Bool_t>                  fFileEnded;
      BoardMap_t<SortedTdcData_t>         fTdcBuffers;
      BoardMap_t<TdcData>                 fLastMrSync;
      BoardMap_t<TdcData>                 fNextMrSync;
      BoardMap_t<TdcData>                 fNext2MrSync;
      BoardMap_t<std::deque<TdcData>>     fRsvdMrSync;

      std::size_t                         fReadCount              =  0;
      ULong64_t                           fDate                   =  0;
      Int_t                               fEMCount                = -1;
      Int_t                               fMrSyncCount            =  0;
      clock_t                             fLastClock;
      Bool_t                              fReadFirst              = false;

      Long64_t                            fMrSyncTdcOffset        = 0;
      std::size_t                         fRefExtChannel          = ExtinctionDetector::NofChannels / 2;
      TdcOffsets_t                        fTdcOffsets;

      TFile*                              fMargedFile             = nullptr;
      TTree*                              fMargedTree             = nullptr;

    public:
      MargedReader(ITdcDataProvider* provider);
      ~MargedReader();

      // Timing tuning
      inline void          SetMrSyncTimeOffset(Double_t offset) {
        std::cout << "SetMrSyncTimeOffset ... " << offset / nsec << " nsec" << std::endl;
        fMrSyncTdcOffset = offset / fProvider->GetTimePerTdc();
      }
      inline Double_t      GetMrSyncTdcOffset() const { return fMrSyncTdcOffset; }

      Int_t                LoadTdcOffsets(const std::string& ffilename);
      Int_t                AddTdcOffsets(const std::string& ffilename);

      void                 InitializeMargedTree(const std::string& filename, const std::string& treename = "mtree");
      void                 WriteMargedTree();

      // File operation
      Int_t                Open(BoardMap_t<ITdcDataProvider*>  providers,
                                const BoardMap_t<std::string>& ifilenames,
                                const std::string&             itreename);

      Int_t                Read(SortedTdcData_t& tdcDataInMrSync);

      Int_t                Close();

      void                 ClearLastSpill();

      // Status getters
      inline ULong64_t     GetDate() const {
        return fDate;
      }
      inline ULong64_t     GetEMCount() const {
        return fEMCount;
      }
      inline Int_t         GetMrSyncCount() const {
        return fMrSyncCount;
      }
      inline Bool_t        IsSpillEnded() const {
        return IsAllOfSecondsTrue(fSpillEnded);
      }
      inline Bool_t        IsFileEnded() const {
        return IsAllOfSecondsTrue(fFileEnded);
      }

    private:
      Int_t ReadUntillNextMrSync(Int_t board);


      template <typename T, typename V>
      inline bool          IsAllOfSecondsTrue(const std::map<T, V> map) const {
        return std::all_of(map.begin(), map.end(), [](std::pair<T, V> p) { return p.second; });
      }
      template <typename T, typename V>
      inline void          FillToSeconds(std::map<T, V>* map, V value) {
        for (auto&& pair : *map) {
          pair.second = value;
        }
      }
    };

    MargedReader::MargedReader(ITdcDataProvider* provider)
      : fProvider(provider) {
    }

    MargedReader::~MargedReader() {
    }

    Int_t MargedReader::LoadTdcOffsets(const std::string& ffilename) {
      std::cout << "Load tdc offsets" << std::endl;
      if (auto ffile = std::ifstream(ffilename)) {
        Int_t gch; Long64_t offset;
        while (ffile >> gch >> offset) {
          fTdcOffsets[gch] = offset;
          // std::cout << gch << "\t" << fTdcOffsets[gch] << std::endl;
        }
      } else {
        std::cerr << "[warning] offset file was not opened, " << ffilename << std::endl;
      }
      return fTdcOffsets.size();
    }

    Int_t MargedReader::AddTdcOffsets(const std::string& ffilename) {
      std::cout << "Add offset tdcs" << std::endl;
      if (auto ffile = std::ifstream(ffilename)) {
        Int_t gch; Double_t offset;
        while (ffile >> gch >> offset) {
          fTdcOffsets[gch] += offset;
          // std::cout << gch << "\t" << fTdcOffsets[gch] << std::endl;
        }
      } else {
        std::cerr << "[warning] offset file was not opened, " << ffilename << std::endl;
      }
      return fTdcOffsets.size();
    }

    void MargedReader::InitializeMargedTree(const std::string& filename, const std::string& treename) {
      std::cout << "Initialize marged tree" << std::endl;

      fMargedFile = new TFile(filename.data(), "RECREATE");
      if (!fMargedFile->IsOpen()) {
        std::cout << "[error] marged file is not opened, " << filename << std::endl;
        return;
      }

      fMargedTree = new TTree(treename.data(), "Marged tree");




      fMargedTree->AutoSave();
    }

    void MargedReader::WriteMargedTree() {
      if (!fMargedFile || !fMargedTree) {
        std::cout << "[error] output marged file is not initialized" << std::endl;
        return;
      } else if (!fMargedFile->IsOpen()) {
        std::cout << "[error] output marged file is not opened, " << fMargedFile->GetName() << std::endl;
        return;
      }

      fMargedFile->cd();
      fMargedTree->Write();
      fMargedFile->Close();
    }

    Int_t MargedReader::Open(BoardMap_t<ITdcDataProvider*>  providers,
                             const BoardMap_t<std::string>& ifilenames,
                             const std::string&             itreename) {
      // Check arguments
      {
        std::set<Int_t> boardP;
        for (auto&& pair : providers) {
          boardP.insert(pair.first);
        }
        std::set<Int_t> boardI;
        for (auto&& pair : ifilenames) {
          boardI.insert(pair.first);
        }
        if (boardP != boardI) {
          std::cout << "[error] invalid arguments, unmatch board number between providers and ifilenames" << std::endl;
          return 1;
        }
      }

      std::cout << "Initialize provider" << std::endl;
      fProviders = providers;

      std::cout << "Open file" << std::endl;
      for (auto&& pair : ifilenames) {
        const Int_t       board     = pair.first;
        const std::string ifilename = pair.second;
        std::cout << " - " << ifilename << std::endl;

        fIfiles[board] = new TFile(ifilename.data(), "READ");
        if (!fIfiles[board]->IsOpen()) {
          std::cerr << "[error] input file is not opened, " << ifilename << std::endl;
          return 1;
        }

        fItrees[board] = dynamic_cast<TTree*>(fIfiles[board]->Get(itreename.data()));
        if (!fItrees[board]) {
          std::cerr << "[error] input tree is not found, " << ifilename << " - " << itreename << std::endl;
          return 1;
        }

        fEntrieses[board] = fItrees[board]->GetEntries();
        std::cout << " + " << fEntrieses[board] << std::endl;

        fProviders[board]->SetBranchAddress(fItrees[board]);
      }
      std::cout << " = " << Tron::Linq::From(fEntrieses).Sum([](std::pair<Int_t, Long64_t> p) { return p.second; }) << std::endl;

      std::cout << "Initialize parameters" << std::endl;
      ClearLastSpill();
      FillToSeconds(&fFileEnded , false);

      for (auto&& pair : fProviders) {
        const Int_t board = pair.first;
        fEntries  [board] = 0;
        fFileEnded[board] = false;
      }

      fReadCount = 0;
      fLastClock = clock();

      return 0;
    }

    Int_t MargedReader::ReadUntillNextMrSync(Int_t board) {
      ITdcDataProvider* provider   = fProviders[board];
      TTree*            itree      = fItrees   [board];
      Long64_t&         entry      = fEntries  [board];
      const Long64_t    entries    = fEntrieses[board];

      // Shift MrSync
      // std::cout << "[debug] Shift MrSync " << board << std::endl;
      fLastMrSync [board] = fNextMrSync [board];
      fNextMrSync [board] = fNext2MrSync[board];
      if (fRsvdMrSync[board].size()) {
        fNext2MrSync[board] = fRsvdMrSync[board].front();
        fRsvdMrSync [board].pop_front();
      } else {
        fNext2MrSync[board] = { };
      }

      // std::cout << "[debug] Read board " << board << std::endl;
      if (fSpillEnded[board] || fFileEnded[board]) {
        return 1;
      }

      for (; !fNext2MrSync[board].Tdc; ++fReadCount, ++entry) {
        if (fReadCount % 200000UL == 0) {
          // std::cout << ">> " << fReadCount << std::endl;
          const clock_t currentClock = clock();
          std::cout << ">> " << fReadCount << "  (" << (double)(currentClock - fLastClock) / CLOCKS_PER_SEC << " sec)" << std::endl;
          fLastClock = currentClock;
        }

        // if (entry >= entries - 3) {
        //   std::cout << "[debug] entry = " << entry << " @ " << board << std::endl;
        // }

        // std::cout << "[debug] GetEntry" << std::endl;
        if (entry >= entries) {
          std::cout << "[info] detect file end @ " << board << std::endl;
          fFileEnded[board] = true;
          break;

        } else if (!itree->GetEntry(entry)) {
          std::cout << "[info] detect file end @ " << board << " (TTree::GetEntry)" << std::endl;
          fFileEnded[board] = true;
          break;

        }

        // if (entry >= entries - 3) {
        //   std::cout << "[debug] spill = " << (Int_t)provider->GetSpill() << std::endl;
        //   std::cout << "[debug] IsData = " << (Int_t)provider->IsData() << ", IsFooter = " << (Int_t)provider->IsFooter() << std::endl;
        // }

        // std::cout << "[debug] Data Check" << std::endl;
        if (provider->GetSpill() < 0) {
          continue;

        } else if (provider->IsFooter()) {
          std::cout << "[info] detect spill end @ " << board << std::endl;
          ++entry;
          fSpillEnded[board] = true;
          break;

        } else if (fEMCount != -1 && fEMCount != provider->GetEMCount()) {
          std::cout << "[error] conflict EMCount, " << fEMCount << " <--> " << provider->GetEMCount() << std::endl;
          exit(1);

        } else if (!provider->IsData()) {
          std::cout << "[info] detect non data entry of " << entry << " @ " << board << std::endl;
          auto tdcData = provider->GetTdcData(board);

          // std::cout << "[debug] data process" << std::endl;
          for (auto&& data : tdcData) {
            // std::cout << "[debug] push data into tdc buffers" << std::endl;
            Tag_t tag = { 0, data.Tdc, data.Channel };
            fTdcBuffers[board].emplace_hint(fTdcBuffers[board].end(), tag, data);
          }

        } else {
          // std::cout << "[info] detect data @ " << board << std::endl;
          fDate        = provider->GetDate();
          fEMCount     = provider->GetEMCount();
          auto tdcData = provider->GetTdcData(board);

          // std::cout << "[debug] data process" << std::endl;
          for (auto&& data : tdcData) {
            if (MrSync::Contains(data.Channel)) {
              // std::cout << "[debug] mr sync " << data.Channel << std::endl;
              data.Tdc -= fMrSyncTdcOffset;
              if        (fNextMrSync[board].Tdc) {
                // std::cout << "[debug] get next next mrsync" << std::endl;
                const Double_t dmsNext = data              .Tdc - fNextMrSync[board].Tdc;
                const Double_t dmsLast = fNextMrSync[board].Tdc - fLastMrSync[board].Tdc;
                const Int_t step = TMath::Nint(dmsNext / dmsLast);
                // const Int_t step = fMrSyncCount <= 1 ? 1 : TMath::Nint(dmsNext / dmsLast);
                // std::cout << "[debug] step == " << step << std::endl;
                if (step == 1) {
                  fNext2MrSync[board] = data;
                } else if (step > 1) {
                  TdcData rsvdData = data;
                  const Long64_t dtdc = (rsvdData.Tdc - fNextMrSync[board].Tdc) / step;
                  for (Int_t i = step; i > 1; --i) {
                    fRsvdMrSync[board].push_front(rsvdData);
                    rsvdData.Tdc -= dtdc;
                  }
                  fNext2MrSync[board] = rsvdData;
                } else {
                  std::cout << "[error] detect time inverted mr sync"                                         << std::endl
                            << "   Board               = " << board                                           << std::endl
                            << "   File                = " << fIfiles[board]->GetName()                       << std::endl
                            << "   Entry               = " << entry                                           << std::endl
                            << "   Last (ms0)          = " << fLastMrSync[board].Tdc                          << std::endl
                            << "   Next (ms1)          = " << fNextMrSync[board].Tdc                          << std::endl
                            << "   data (ms2)          = " << data              .Tdc                          << std::endl
                            << "   ms1 - ms0           = " << fNextMrSync[board].Tdc - fLastMrSync[board].Tdc << std::endl
                            << "   ms2 - ms1           = " << data              .Tdc - fNextMrSync[board].Tdc << std::endl
                            << "   (ms2-ms1)/(ms1-ms0) = " << step                                            << std::endl
                    ;
                  exit(1);
                }

              } else if (fLastMrSync[board].Tdc) {
                // std::cout << "[debug] get next mrsync" << std::endl;
                fNextMrSync[board] = data;
              } else {
                // std::cout << "[debug] get first mrsync" << std::endl;
                fLastMrSync[board] = data;
              }

            } else if (BeamlineHodoscope::Contains(data.Channel)) {
              // std::cout << "[debug] beamline hodoscope " << data.Channel << std::endl;
              data.Tdc -= fTdcOffsets[data.Channel];

            } else if (Hodoscope::Contains(data.Channel)) {
              // std::cout << "[debug] hodoscope " << data.Channel << std::endl;
              data.Tdc -= fTdcOffsets[data.Channel];

            } else if (ExtinctionDetector::Contains(data.Channel)) {
              // std::cout << "[debug] extinction detector " << data.Channel << std::endl;
              data.Tdc -= fTdcOffsets[data.Channel];

            } else if (TimingCounter::Contains(data.Channel)) {
              // std::cout << "[debug] timing counter " << data.Channel << std::endl;
              data.Tdc -= fTdcOffsets[data.Channel];

            } else if (Veto::Contains(data.Channel)) {
              // std::cout << "[debug] veto " << data.Channel << std::endl;
              data.Tdc -= fMrSyncTdcOffset;

            } else if (EventMatch::Contains(data.Channel)) {
              // std::cout << "[debug] event match " << data.Channel << std::endl;
              data.Tdc -= fMrSyncTdcOffset;

            } else {
              // std::cout << "[debug] other " << data.Channel << std::endl;
              continue;
            }

            // std::cout << "[debug] push data into tdc buffers" << std::endl;
            Tag_t tag = { 0, data.Tdc, data.Channel };
            fTdcBuffers[board].emplace_hint(fTdcBuffers[board].end(), tag, data);
          }
        }
      }

      return 0;
    }

    Int_t MargedReader::Read(SortedTdcData_t& tdcDataInMrSync) {
      ++fMrSyncCount;

      for (auto&& pair : fProviders) {
        const Int_t board = pair.first;
        ReadUntillNextMrSync(board);
      }

      // Incorrect mr sync scan
      if (fReadFirst) {
        std::size_t retry;
        for (retry = 0; retry < 100; ++retry) {
          Int_t    minBoard     = fLastMrSync.begin()->first;
          Long64_t minMrSyncTdc = fLastMrSync.begin()->second.Tdc;
          Int_t    maxBoard     = minBoard;
          Long64_t maxMrSyncTdc = minMrSyncTdc;
          for (auto&& pair : fLastMrSync) {
            const Long64_t tdc = pair.second.Tdc;
            if        (tdc < minMrSyncTdc) {
              minBoard     = pair.first;
              minMrSyncTdc = tdc;
            } else if (tdc > maxMrSyncTdc) {
              maxBoard     = pair.first;
              maxMrSyncTdc = tdc;
            }
          }

          if ((maxMrSyncTdc - minMrSyncTdc) * fProvider->GetTimePerTdc() > 1.0 * usec) {
            std::cout << "[warning] first mr syncs are so different, "
                      << "min tdc = " << minMrSyncTdc << " @ " << minBoard << ", "
                      << "max tdc = " << maxMrSyncTdc << " @ " << maxBoard << std::endl;
            auto itr = fTdcBuffers[minBoard].find({ 0, fLastMrSync[minBoard].Tdc, fLastMrSync[minBoard].Channel });
            if (itr == fTdcBuffers[minBoard].end()) {
              std::cerr << "[error] mr sync was not found on removing till incorrect mr sync" << std::endl;
              exit(1);
            }
            fTdcBuffers[minBoard].erase(fTdcBuffers[minBoard].begin(), ++itr);
            ReadUntillNextMrSync(minBoard);
          } else {
            std::cout << "[info] first mr syncs are not so different, "
                      << "min tdc = " << minMrSyncTdc << " @ " << minBoard << ", "
                      << "max tdc = " << maxMrSyncTdc << " @ " << maxBoard << std::endl;
            break;
          }
        }
        if (retry == 100) {
          std::cerr << "[error] too many retry, incorrect mr sync scan" << std::endl;
          exit(1);
        }

        fReadFirst = false;
      }

      // std::cout << "[debug] check end" << std::endl;
      if (IsSpillEnded()) {
        std::cout << "[info] spill ended" << std::endl;
        return 0;
      }
      if (IsFileEnded()) {
        std::cout << "[info] file ended" << std::endl;
        return 0;
      }

      // Extract tdc data
      for (auto&& pair : fProviders) {
        const Int_t    board      = pair.first;
        const TdcData& lastMrSync = fLastMrSync[board];
        const TdcData& nextMrSync = fNextMrSync[board];
        SortedTdcData_t::iterator itr, end;
        for (itr = fTdcBuffers[board].begin(), end = fTdcBuffers[board].end(); itr != end; ++itr) {
          TdcData& data = itr->second;
          if        (data.Tdc > nextMrSync.Tdc) {
            break;
          } else if (data.Tdc > lastMrSync.Tdc) {
            data.LastMrSyncCount = fMrSyncCount;
            data.LastMrSyncTdc   = lastMrSync.Tdc;
            data.NextMrSyncTdc   = nextMrSync.Tdc;
            data.TdcFromMrSync   = data.Tdc - data.LastMrSyncTdc;
            const Tag_t tag { data.LastMrSyncCount, data.TdcFromMrSync, data.Channel };
            tdcDataInMrSync.emplace(tag, data);
          } else {
            // Nothing to do
          }
        }
        // std::cout << "[debug] erase buffer" << std::endl;
        fTdcBuffers[board].erase(fTdcBuffers[board].begin(), itr);
      }
      // std::cout << "[debug] end of read" << std::endl;

      return 1;
    }

    Int_t MargedReader::Close() {
      std::cout << "Close files" << std::endl;
      for (auto&& pair : fIfiles) {
        pair.second->Close();
      }

      return 0;
    }

    void MargedReader::ClearLastSpill() {
      fTdcBuffers.clear();
      fDate        =  0;
      fEMCount     = -1;
      fMrSyncCount =  0;
      fReadFirst   = true;

      FillToSeconds(&fSpillEnded, false);

      for (auto&& pair : fProviders) {
        const Int_t board = pair.first;
        fLastMrSync [board].Clear();
        fNextMrSync [board].Clear();
        fNext2MrSync[board].Clear();
        fRsvdMrSync [board].clear();
      }
    }
  }

}

#endif
