#ifndef Extinction_Fct_MonitorWindow_hh
#define Extinction_Fct_MonitorWindow_hh

#include <iostream>
#include <fstream>
#include <time.h>
#include <thread>

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

#include "Linq.hh"
#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "Fct.hh"
#include "AnaCoin.hh"

namespace Extinction {

  namespace Fct {

    class MonitorWindow {
    public:
      enum class ChannelAlign {
        Raw, Projection,
      };

      enum class MonitorMode {
        Channel, Ext, Hod, Offset, Coincidence,
      };

      struct PlotsProfiles {
        struct {
          Double_t NbinsX, Xmin, Xmax, Unit;
        } TimeInSpill;
        struct {
          Double_t Xmin, Xmax, BinWidth;
          Double_t Xwidth() const { return Xmax - Xmin; }
        } TimeInSync;
        struct {
          Double_t Mean, Xmin, Xmax;
          Double_t Xwidth() const { return Xmax - Xmin; }
        } MrSyncInterval;
        struct {
          Double_t Xmin, Xmax;
          Double_t Xwidth() const { return Xmax - Xmin; }
        } TimeDiff;
      };

    private:
      static const std::size_t kHistLimit = 10000;

    private:
      using CoinOffset = Analyzer::AnaTimeOffset::CoinOffset;
      using CoinInfo   = Analyzer::AnaCoin::CoinInfo;
      using CoinDiffs  = std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Double_t>>;
      using Contains   = std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Bool_t>>;

      struct CoinOffsetX {
        enum {
          BH  = 0,
          Hod = 0 + BeamlineHodoscope::NofChannels,
          TC  = 0 + BeamlineHodoscope::NofChannels + Hodoscope::NofChannels,
          N   = 0 + BeamlineHodoscope::NofChannels + Hodoscope::NofChannels + TimingCounter::NofChannels,
        };
      };

      struct ContainsOffset {
        enum : unsigned long long {
          Max = (unsigned long long)-1,
          TC  = (unsigned long long)-1 - TimingCounter::NofChannels,
          Hod = (unsigned long long)-1 - TimingCounter::NofChannels - 1,
          BH  = (unsigned long long)-1 - TimingCounter::NofChannels - 1 - BeamlineHodoscope::NofChannels,
        };
      };

      FctData                      fProvider;
      CoinDiffs                    fCoinDiffs;
      Contains                     fContains;

      TApplication*                fApplication         = nullptr;

      TCanvas*                     fCanvas              = nullptr;
      TPad*                        fPadBh1Tdc           = nullptr;
      TPad*                        fPadBh2Tdc           = nullptr;
      TPad*                        fPadExtTdc           = nullptr;
      TPad*                        fPadHodTdc           = nullptr;
      TPad*                        fPadTc1Tdc           = nullptr;
      TPad*                        fPadTc2Tdc           = nullptr;
      TPad*                        fPadHodHitMap        = nullptr;
      TPad*                        fPadHodHitCnt        = nullptr;
      TPad*                        fPadExtHitMap        = nullptr;
      TPad*                        fPadExtHitCnt        = nullptr;
      TPad*                        fPadMountain         = nullptr;
      TPad*                        fPadTdcSync          = nullptr;
      TPad*                        fPadHit              = nullptr;

      TH1*                         hExtTdcInSpill_Any   = nullptr;
      TH1*                         hHodTdcInSpill_Any   = nullptr;
      TH1**                        hTcTdcInSpill        = nullptr;
      TH1**                        hBhTdcInSpill        = nullptr;
      TH2*                         hExtHitMap           = nullptr;
      TList*                       lExtBorderLine       = nullptr;
      TH1*                         hExtEntryByCh        = nullptr;
      TH1*                         hExtEntryByChBottom  = nullptr;
      TH1*                         hExtEntryByChCenter1 = nullptr;
      TH1*                         hExtEntryByChCenter2 = nullptr;
      TH1*                         hExtEntryByChTop     = nullptr;
      TH2*                         hHodHitMap           = nullptr;
      TList*                       lHodBorderLine       = nullptr;
      TH1*                         hHodEntryByCh        = nullptr;
      TH2*                         hExtMountain_Any     = nullptr;
      TH1*                         hExtTdcInSync_Any    = nullptr;
      TGraphErrors*                gHitInSpill          = nullptr;
      TH2*                         hMountain            = nullptr;
      TH1*                         hTdcInSync           = nullptr;
      TH1*                         hTdcInSpill          = nullptr;
      TH1**                        hMrSyncInterval      = nullptr;
      TH2**                        hExtTdcOffset        = nullptr;

      TDatime                      fDate;
      Long64_t                     fSpillCount          = 0;
      Long64_t                     fCoinCount           = 0;
      Int_t                        fEventMatchNumber    = -1;

      std::map<Int_t, Double_t>    fTimePerTdc;
      std::map<Int_t, Double_t>    fMrSyncInterval;
      Double_t                     fMrSyncIntervalAverage;
      std::map<ULong64_t, TdcData> fTdcBuffer;

      std::vector<TdcData>         fLastExtData;
      std::vector<TdcData>         fLastHodData;
      std::vector<TdcData>         fLastTcData;
      std::vector<TdcData>         fLastBhData;
      std::map<Int_t, TdcData>     fLastMrSyncData;
      std::map<Int_t, std::vector<TdcData>> fEventMatchData;

      Bool_t                       fIsTerminated        = false;
      ChannelAlign                 fChannelAlign        = ChannelAlign::Raw;
      MonitorMode                  fMonitorMode         = MonitorMode::Coincidence;
      Int_t                        fMonitorChannel      = -1; // global channel
      Int_t                        fMonitorBoard        = -1;
      Bool_t                       fCyclicCoincidence   = false; //true;

      Double_t                     fHistoryWidth        = 600.0 * nsec;
      Double_t                     fCoinTimeWidth       =  10.0 * nsec;
      Long64_t                     fSpillLimit          = 1 * 60 * 60 / 5;
      std::size_t                  fBufferSize          = 5000;
      std::size_t                  fBufferMargin        =  100;

      std::string                  fMrSyncIntervalFilename;
      std::string                  fCoinDiffsFilename;

    public:
      MonitorWindow();
      ~MonitorWindow();

      inline void          SetTimePerTdc(const std::map<Int_t, Double_t>& map) { fTimePerTdc = map; }
      inline void          SetMrSyncInterval(const std::map<Int_t, Double_t>& map) {
        fMrSyncInterval        = map;
        fMrSyncIntervalAverage = Tron::Linq::From(map)
          .Select([](const std::pair<Int_t, Double_t>& pair) { return pair.second; })
          .Average();
        std::cout << "SetMrSyncInterval" << std::endl;
        for (auto&& pair : fMrSyncInterval) {
          std::cout << pair.first << "\t" << pair.second << std::endl;
        }
        std::cout << "Ave: " << fMrSyncIntervalAverage << std::endl;
      }
      inline void          SetChannelAlign(ChannelAlign align) { fChannelAlign = align; }
      inline ChannelAlign  GetChannelAlign() const { return fChannelAlign; }
      inline Bool_t        SetMonitorMode(MonitorMode mode, Int_t option = -1) {
        switch (mode) {
        case MonitorMode::Channel:
          if (option >= 0) {
            fMonitorMode    = mode;
            fMonitorChannel = option;
            return true;
          } else {
            std::cout << "[warning] invalid channel number, " << option << std::endl;
          }
          break;
        case MonitorMode::Offset:
          if (ExtinctionDetector::Contains(option)) {
            fMonitorMode    = mode;
            fMonitorChannel = option;
            fMonitorBoard   = ChannelMapWithBoard::Board[option];
            return true;
          } else {
            std::cout << "[warning] invalid channel number, " << option << std::endl;
          }
          break;
        case MonitorMode::Ext:
        case MonitorMode::Hod:
        case MonitorMode::Coincidence:
          fMonitorMode = mode;
          return true;
        }
        return false;
      }
      inline MonitorMode   GetMonitorMode() const { return fMonitorMode; }
      inline Int_t         GetMonitorChannel() const { return fMonitorChannel; }

      void                 SetHistoryWidth(Double_t width) { fHistoryWidth = width; }
      Double_t             GetHistoryWidth() const { return fHistoryWidth; }

      void                 SetCoinTimeWidth(Double_t width) { fCoinTimeWidth = width; }
      Double_t             GetCoinTimeWidth() const { return fCoinTimeWidth; }

      void                 SetSpillLimit(Long64_t limit) { fSpillLimit = limit; }
      Long64_t             GetSpillLimit() const { return fSpillLimit; }

      void                 SetReadBufferSize(std::size_t size)  { fBufferSize = size; }
      std::size_t          GetReadBufferSize() const { return fBufferSize; }

      void                 SetReadBufferMargin(std::size_t size)  { fBufferMargin = size; }
      std::size_t          GetReadBufferMargin() const { return fBufferMargin; }

      void                 SetMrSyncIntervalFilename(const std::string& filename) {
        fMrSyncIntervalFilename = filename;
      }
      void                 SetCoinDiffsFilename(const std::string& filename) {
        fCoinDiffsFilename = filename;
      }
      std::string          GetMrSyncIntervalFilename() const { return fMrSyncIntervalFilename; }
      std::string          GetCoinDiffsFilename() const { return fCoinDiffsFilename; }

      Int_t                LoadOffset(const std::string& ffilename);
      void                 InitializeWindow(Int_t width = 1600, Int_t height = 1200);
      void                 InitializePlots(const PlotsProfiles& profile);

      void                 DrawPlots();

      Int_t                UpdatePlots(const std::map<Int_t, std::string>& ifilenames,
                                       const std::function<TDatime(const std::string&)>& parser = nullptr);

      inline void          Run() { fApplication->Run(true); }

      inline Bool_t        IsClosed() const { return !fCanvas->GetCanvasImp(); }
      inline Bool_t        IsTerminated() const { return fIsTerminated; }

      inline void          Terminate() {
        fIsTerminated = true;
        gSystem->ExitLoop();
      }

    private:
      void                 ClearLastSpill();
      void                 FillCoincidence(const TdcData& tdcData);
      void                 FillCoincidences(const std::vector<TdcData>& tdcData);
      void                 FillCoincidence2(const TdcData& tdcData);
      std::vector<TdcData> CollectCoinExtData(const TdcData& tdcData, std::size_t i);
      std::size_t          RemoveOldTdc(std::vector<TdcData>* lastData, Double_t time);

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

    MonitorWindow::MonitorWindow() {
      fApplication = new TApplication("monitor", nullptr, nullptr);
      fLastExtData.reserve(100);
      fLastHodData.reserve(100);
      fLastTcData .reserve(100);
      fLastBhData .reserve(100);
    }

    MonitorWindow::~MonitorWindow() {
      delete fApplication;
      fApplication = nullptr;
    }

    Int_t MonitorWindow::LoadOffset(const std::string& ffilename) {
      std::cout << "Load offset" << std::endl;
      if (Tron::String::GetFormatSpecifiers(ffilename).empty()) {
        std::ifstream ffile(ffilename);
        if (ffile) {
          Int_t ch, index; Double_t mean;
          while (ffile >> ch >> index >> mean) {
            fCoinDiffs[ch][index] = mean * fProvider.GetTimePerTdc();
          }
        }
      } else {
        for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
          CoinInfo info;
          std::ifstream ffile(Form(ffilename.data(), extCh));
          if (ffile) {
            while (info.Read(ffile)) {
              fCoinDiffs[extCh][info.Index] = info.FitMean * fProvider.GetTimePerTdc();
            }
          }
        }
      }

      for (std::size_t ch = ContainsOffset::BH; ch < (std::size_t)ContainsOffset::Max; ++ch) {
        const std::size_t index = ch - ContainsOffset::BH;
        // std::cout << "++ ch (index) " << ch << " (" << index << ")" << std::endl;
        std::map<std::size_t/*index*/, Double_t> sum, cnt;
        for (auto&& info : fCoinDiffs) {
          const std::size_t extCh = info.first;
          if (ExtinctionDetector::Contains(extCh)) {
            for (auto&& subInfo : info.second) {
              const std::size_t index2 = subInfo.first;
              // std::cout << "     " << extCh << "," << index2 << "," << index << "    " << (fCoinDiffs[extCh][index2] - fCoinDiffs[extCh][index]) / nsec << " nsec" << std::endl;
              sum[index2] += fCoinDiffs[extCh][index2] - fCoinDiffs[extCh][index];
              cnt[index2]++;
            }
          }
        }
        for (auto&& pair : sum) {
          fCoinDiffs[ch][pair.first] = sum[pair.first] / cnt[pair.first];
          // std::cout << " [" << pair.first << "] " << fCoinDiffs[ch][pair.first] / nsec << " nsec" << std::endl;
        }
      }

      for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
        for (std::size_t index = 0; index < CoinOffset::N; ++index) {
          fContains[extCh][index] = false;
        }
      }
      for (std::size_t ch = ContainsOffset::BH; ch < (std::size_t)ContainsOffset::Max; ++ch) {
        for (std::size_t index = 0; index < CoinOffset::N; ++index) {
          fContains[ch][index] = false;
        }
      }
      for (auto&& info : fCoinDiffs) {
        const std::size_t extCh = info.first;
        for (auto&& subInfo : info.second) {
          const std::size_t index = subInfo.first;
          fContains[extCh][index] = true;
        }
      }

      return fCoinDiffs.size();
    }

    void MonitorWindow::InitializeWindow(Int_t width, Int_t height) {
      std::cout << "Initialize window" << std::endl;
      fCanvas = new TCanvas("c1", (fProvider.GetName() + " | Semi-online Monitor").data(), width, height);

      std::vector<TPad*> pads;
      const Double_t seg = 1.0 / 6.0;

      pads.push_back(fPadBh1Tdc    = new TPad("fPadBh1Tdc"   , "", 0.0 * seg, 5.0 * seg, 1.0 * seg, 6.0 * seg));
      pads.push_back(fPadBh2Tdc    = new TPad("fPadBh2Tdc"   , "", 0.0 * seg, 4.0 * seg, 1.0 * seg, 5.0 * seg));
      pads.push_back(fPadExtTdc    = new TPad("fPadExtTdc"   , "", 0.0 * seg, 3.0 * seg, 1.0 * seg, 4.0 * seg));
      pads.push_back(fPadHodTdc    = new TPad("fPadHodTdc"   , "", 0.0 * seg, 2.0 * seg, 1.0 * seg, 3.0 * seg));
      pads.push_back(fPadTc1Tdc    = new TPad("fPadTc1Tdc"   , "", 0.0 * seg, 1.0 * seg, 1.0 * seg, 2.0 * seg));
      pads.push_back(fPadTc2Tdc    = new TPad("fPadTc2Tdc"   , "", 0.0 * seg, 0.0 * seg, 1.0 * seg, 1.0 * seg));
      fPadBh1Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadBh1Tdc   ->SetGrid(); fPadBh1Tdc   ->SetLogy();
      fPadBh2Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadBh2Tdc   ->SetGrid(); fPadBh2Tdc   ->SetLogy();
      fPadExtTdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadExtTdc   ->SetGrid(); fPadExtTdc   ->SetLogy();
      fPadHodTdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHodTdc   ->SetGrid(); fPadHodTdc   ->SetLogy();
      fPadTc1Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTc1Tdc   ->SetGrid(); fPadTc1Tdc   ->SetLogy();
      fPadTc2Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTc2Tdc   ->SetGrid(); fPadTc2Tdc   ->SetLogy();

      pads.push_back(fPadHodHitMap = new TPad("fPadHodHitMap", "", 1.0 * seg, 4.5 * seg, 3.0 * seg, 6.0 * seg));
      pads.push_back(fPadHodHitCnt = new TPad("fPadHodHitCnt", "", 1.0 * seg, 3.0 * seg, 3.0 * seg, 4.5 * seg));
      pads.push_back(fPadExtHitMap = new TPad("fPadExtHitMap", "", 1.0 * seg, 1.5 * seg, 3.0 * seg, 3.0 * seg));
      pads.push_back(fPadExtHitCnt = new TPad("fPadExtHitCnt", "", 1.0 * seg, 0.0 * seg, 3.0 * seg, 1.5 * seg));
      fPadHodHitMap->SetMargin(0.06, 0.06, 0.10, 0.10);
      fPadHodHitCnt->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHodHitCnt->SetGrid(); fPadHodHitCnt->SetLogy();
      fPadExtHitMap->SetMargin(0.06, 0.06, 0.10, 0.10);
      fPadExtHitCnt->SetMargin(0.06, 0.06, 0.10, 0.10); fPadExtHitCnt->SetGrid(); fPadExtHitCnt->SetLogy();

      if (fMonitorMode == MonitorMode::Offset) {
        pads.push_back(fPadMountain  = new TPad("fPadMountain" , "", 3.0 * seg, 4.0 * seg, 6.0 * seg, 6.0 * seg));
        pads.push_back(fPadTdcSync   = new TPad("fPadTdcSync"  , "", 3.0 * seg, 0.0 * seg, 6.0 * seg, 4.0 * seg));
        fPadMountain ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadMountain ->SetGrid(); fPadMountain ->SetLogy();
        fPadTdcSync  ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTdcSync  ->SetGrid();
      } else if (fMonitorMode == MonitorMode::Channel ||
                 fMonitorMode == MonitorMode::Ext     ||
                 fMonitorMode == MonitorMode::Hod     ) {
        pads.push_back(fPadMountain  = new TPad("fPadMountain" , "", 3.0 * seg, 4.0 * seg, 6.0 * seg, 6.0 * seg));
        pads.push_back(fPadTdcSync   = new TPad("fPadTdcSync"  , "", 3.0 * seg, 2.0 * seg, 6.0 * seg, 4.0 * seg));
        pads.push_back(fPadHit       = new TPad("fPadHit"      , "", 3.0 * seg, 0.0 * seg, 6.0 * seg, 2.0 * seg));
        fPadMountain ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadMountain ->SetGrid();
        fPadTdcSync  ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTdcSync  ->SetGrid(); fPadTdcSync  ->SetLogy();
        fPadHit      ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHit      ->SetGrid(); fPadHit      ->SetLogy();
      } else if (fMonitorMode == MonitorMode::Coincidence) {
        pads.push_back(fPadMountain  = new TPad("fPadMountain" , "", 3.0 * seg, 4.0 * seg, 6.0 * seg, 6.0 * seg));
        pads.push_back(fPadTdcSync   = new TPad("fPadTdcSync"  , "", 3.0 * seg, 2.0 * seg, 6.0 * seg, 4.0 * seg));
        pads.push_back(fPadHit       = new TPad("fPadHit"      , "", 3.0 * seg, 0.0 * seg, 6.0 * seg, 2.0 * seg));
        fPadMountain ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadMountain ->SetGrid();
        fPadTdcSync  ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTdcSync  ->SetGrid(); fPadTdcSync  ->SetLogy();
        fPadHit      ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHit      ->SetGrid();
      }

      for (Int_t ipad = 0, npad = pads.size(); ipad < npad; ++ipad) {
        fCanvas->cd();
        pads[ipad]->Draw();
        pads[ipad]->SetNumber(ipad + 1);
      }
    }

    void MonitorWindow::InitializePlots(const PlotsProfiles& profile) {
      std::cout << "Initialize plots" << std::endl;
      const std::string tdcName    = fProvider.GetName();
      const Double_t    timePerTdc = fProvider.GetTimePerTdc();

      const Double_t  xminInSpill = profile.TimeInSpill.Xmin / profile.TimeInSpill.Unit;
      const Double_t  xmaxInSpill = profile.TimeInSpill.Xmax / profile.TimeInSpill.Unit;
      const Int_t    xbinsInSpill = profile.TimeInSpill.NbinsX;

      const Double_t  xminInSync  = (Int_t)(profile.TimeInSync.Xmin     / timePerTdc) - 0.5;
      const Int_t    xbinsInSync  = (Int_t)(profile.TimeInSync.Xwidth() / timePerTdc) / profile.TimeInSync.BinWidth;
      const Double_t  xmaxInSync  = xminInSync + xbinsInSync * profile.TimeInSync.BinWidth;

      const Double_t  xminInt     = (Int_t)(profile.MrSyncInterval.Xmin     / timePerTdc) - 0.5;
      const Int_t    xbinsInt     = (Int_t)(profile.MrSyncInterval.Xwidth() / timePerTdc);
      const Double_t  xmaxInt     = xminInt + xbinsInt;

      const Double_t  xminInDiff  = (Int_t)(profile.TimeDiff.Xmin     / timePerTdc) - 0.5;
      const Int_t    xbinsInDiff  = (Int_t)(profile.TimeDiff.Xwidth() / timePerTdc);
      const Double_t  xmaxInDiff  = xminInDiff + xbinsInDiff;

      // Extinction Detector TDC in spill
      hExtTdcInSpill_Any = new TH1D("hExtTdcInSpill_Any",
                                    Form("%s, Extinction Detector TDC in Spill;"
                                         "Time [ms]", tdcName.data()),
                                    xbinsInSpill, xminInSpill, xmaxInSpill);
      hExtTdcInSpill_Any->SetStats(false);

      // Hodoscope TDC in spill
      hHodTdcInSpill_Any = new TH1D("hHodTdcInSpill_Any",
                                    Form("%s, Hodoscope TDC in Spill;"
                                         "Time [ms]", tdcName.data()),
                                    xbinsInSpill, xminInSpill, xmaxInSpill);
      hHodTdcInSpill_Any->SetStats(false);

      // Timing Counter TDC in spill
      hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSpill[ch] = new TH1D(Form("hTcTdcInSpill_%03ld", ch),
                                     Form("%s, Timing Counter %ld TDC in Spill;"
                                          "Time [ms]", tdcName.data(), ch + 1),
                                     xbinsInSpill, xminInSpill, xmaxInSpill);
        hTcTdcInSpill[ch]->SetStats(false);
      }

      // Beamline Hodoscope TDC in spill
      hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSpill[ch] = new TH1D(Form("hBhTdcInSpill_%03ld", ch),
                                     Form("%s, BH%ld TDC in Spill;"
                                          "Time [ms]", tdcName.data(), ch + 1),
                                     xbinsInSpill, xminInSpill, xmaxInSpill);
        hBhTdcInSpill[ch]->SetStats(false);
      }

      // Extinction detector hit map
      hExtHitMap     = ExtinctionDetector::CreateHitMap("hExtHitMap");
      lExtBorderLine = ExtinctionDetector::CreateBorderLine("lExtBorderLine", kBlack, kSolid, 1);

      // Extinction detector hit count
      hExtEntryByCh = new TH1D("hExtEntryByCh",
                               Form("%s, Extinction Detector Entries by Channel;"
                                    "Channel;"
                                    "", tdcName.data()),
                               ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5);
      hExtEntryByCh->SetStats(false);

      hExtEntryByChBottom  = hExtHitMap->ProjectionX("hExtEntryByChBottom" , 1, 1);
      hExtEntryByChCenter1 = hExtHitMap->ProjectionX("hExtEntryByChCenter1", 2, 2);
      hExtEntryByChCenter2 = hExtHitMap->ProjectionX("hExtEntryByChCenter2", 3, 3);
      hExtEntryByChTop     = hExtHitMap->ProjectionX("hExtEntryByChTop"    , 4, 4);
      hExtEntryByChBottom ->SetLineColor(kBlue   + 1);
      hExtEntryByChCenter1->SetLineColor(kRed    + 1);
      hExtEntryByChCenter2->SetLineColor(kOrange + 1);
      hExtEntryByChTop    ->SetLineColor(kGreen  + 1);
      hExtEntryByChBottom ->SetStats(false);
      hExtEntryByChCenter1->SetStats(false);
      hExtEntryByChCenter2->SetStats(false);
      hExtEntryByChTop    ->SetStats(false);

      // Hodoscope hit map
      hHodHitMap     = Hodoscope::CreateHitMap("hHodHitMap");
      lHodBorderLine = Hodoscope::CreateBorderLine("lHodBorderLine", kBlack, kSolid, 1);

      // Hodoscope hit count
      hHodEntryByCh = new TH1D("hHodEntryByCh",
                               Form("%s, Hodoscope Entries by Channel;"
                                    "Channel;"
                                    "", tdcName.data()),
                               Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5);
      hHodEntryByCh->SetStats(false);

      // Extinction Detector Mountain Plot (coincidence)
      hExtMountain_Any = new TH2D("hExtMountain_Any",
                                  Form("%s, Extinction Detector Mountain Plot;"
                                       "TDC [count];"
                                       "Time [ms]", tdcName.data()),
                                  xbinsInSync / 2, xminInSync, xmaxInSync,
                                  xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hExtMountain_Any->SetStats(false);

      // Extinction Detector TDC in sync (coincidence)
      hExtTdcInSync_Any = new TH1D("hExtTdcInSync_Any",
                                   Form("%s, Extinction Detector TDC in MR Sync;"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInSync, xminInSync, xmaxInSync);

      // Hit in spill (coincidence)
      gHitInSpill = new TGraphErrors();
      gHitInSpill->SetNameTitle("gHitInSpill",
                                Form("%s, # of Hits in Spill;"
                                     "Spill", tdcName.data()));
      gHitInSpill->SetMarkerStyle(kPlus);
      gHitInSpill->SetMarkerColor(kBlue + 1);
      gHitInSpill->SetPoint(0, 0, 0);

      // Channel monitor hists
      if (fMonitorMode == MonitorMode::Channel) {
        hMountain = new TH2D("hMountain",
                             Form("%s, Mountain Plot @ ch%d;"
                                  "TDC [count];"
                                  "Time [ms]", tdcName.data(), fMonitorChannel),
                             xbinsInSync / 2, xminInSync, xmaxInSync,
                             xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hMountain->SetStats(false);

        hTdcInSync = new TH1D("hTdcInSync",
                              Form("%s, TDC in MR Sync @ ch%d;"
                                   "TDC [count]", tdcName.data(), fMonitorChannel),
                              xbinsInSync, xminInSync, xmaxInSync);

        hTdcInSpill = new TH1D("hTdcInSpill",
                               Form("%s, TDC in Spill @ ch%d;"
                                    "Time [ms]", tdcName.data(), fMonitorChannel),
                               xbinsInSpill, xminInSpill, xmaxInSpill);
        hTdcInSpill->SetStats(false);
      } else if (fMonitorMode == MonitorMode::Ext) {
        hMountain = new TH2D("hMountain",
                             Form("%s, Mountain Plot @ Extinction Detector;"
                                  "TDC [count];"
                                  "Time [ms]", tdcName.data()),
                             xbinsInSync / 2, xminInSync, xmaxInSync,
                             xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hMountain->SetStats(false);

        hTdcInSync = new TH1D("hTdcInSync",
                              Form("%s, TDC in MR Sync @ Extinction Detector;"
                                   "TDC [count]", tdcName.data()),
                              xbinsInSync, xminInSync, xmaxInSync);

        hTdcInSpill = new TH1D("hTdcInSpill",
                               Form("%s, TDC in Spill @ Extinction Detector;"
                                    "Time [ms]", tdcName.data()),
                               xbinsInSpill, xminInSpill, xmaxInSpill);
        hTdcInSpill->SetStats(false);
      } else if (fMonitorMode == MonitorMode::Hod) {
        hMountain = new TH2D("hMountain",
                             Form("%s, Mountain Plot @ Hodoscope;"
                                  "TDC [count];"
                                  "Time [ms]", tdcName.data()),
                             xbinsInSync / 2, xminInSync, xmaxInSync,
                             xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hMountain->SetStats(false);

        hTdcInSync = new TH1D("hTdcInSync",
                              Form("%s, TDC in MR Sync @ Hodoscope;"
                                   "TDC [count]", tdcName.data()),
                              xbinsInSync, xminInSync, xmaxInSync);

        hTdcInSpill = new TH1D("hTdcInSpill",
                               Form("%s, TDC in Spill @ Hodoscope;"
                                    "Time [ms]", tdcName.data()),
                               xbinsInSpill, xminInSpill, xmaxInSpill);
        hTdcInSpill->SetStats(false);
      }

      // Offset monitor hists
      if (fMonitorMode == MonitorMode::Offset) {
        hMrSyncInterval = new TH1*[MrSync::NofChannels];
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncInterval[ch] = new TH1D(Form("hMrSyncInterval_%03ld", ch),
                                         Form("%s, MR Sync TDC Interval @ ch%ld;"
                                              "TDC [count];"
                                              "", tdcName.data(), ch),
                                         xbinsInt, xminInt, xmaxInt);
        }

        hExtTdcOffset = new TH2*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcOffset[ch] = new TH2D(Form("hExtTdcOffset_%03ld", ch),
                                       Form("%s, Extinction Detector TDC Offset @ ch%ld;"
                                            "Channel;"
                                            "TDC [count]", tdcName.data(), ch),
                                       CoinOffsetX::N, 0 - 0.5, CoinOffsetX::N - 0.5,
                                       xbinsInDiff, xminInDiff, xmaxInDiff);
          hExtTdcOffset[ch]->SetStats(false);
        }
      }
    }

    void MonitorWindow::DrawPlots() {
      std::cout << "Draw plots" << std::endl;
      Int_t padnumber = 0;

      if (fEventMatchNumber >= 0) {
        fCanvas->SetTitle((fProvider.GetName() + " | Semi-online Monitor" +
                           Form(" (%04d/%02d/%02d %02d:%02d:%02d - %d)",
                                fDate.GetYear(),
                                fDate.GetMonth(),
                                fDate.GetDay(),
                                fDate.GetHour(),
                                fDate.GetMinute(),
                                fDate.GetSecond(),
                                fEventMatchNumber)).data());
      } else {
        fCanvas->SetTitle((fProvider.GetName() + " | Semi-online Monitor" +
                           Form(" (%04d/%02d/%02d %02d:%02d:%02d)",
                                fDate.GetYear(),
                                fDate.GetMonth(),
                                fDate.GetDay(),
                                fDate.GetHour(),
                                fDate.GetMinute(),
                                fDate.GetSecond()
                                )).data());
      }

      if (fCanvas->cd(++padnumber)) {
        hBhTdcInSpill[0]->Draw();
        hBhTdcInSpill[0]->SetMinimum(0.2);
        hBhTdcInSpill[0]->GetXaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[0]->GetYaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[0]->GetXaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[0]->GetYaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[0]->GetXaxis()->SetTitleOffset(0.6);
        hBhTdcInSpill[0]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hBhTdcInSpill[1]->Draw();
        hBhTdcInSpill[1]->SetMinimum(0.2);
        hBhTdcInSpill[1]->GetXaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[1]->GetYaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[1]->GetXaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[1]->GetYaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[1]->GetXaxis()->SetTitleOffset(0.6);
        hBhTdcInSpill[1]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hHodTdcInSpill_Any->Draw();
        hHodTdcInSpill_Any->SetMinimum(0.2);
        hHodTdcInSpill_Any->GetXaxis()->SetLabelSize(0.07);
        hHodTdcInSpill_Any->GetYaxis()->SetLabelSize(0.07);
        hHodTdcInSpill_Any->GetXaxis()->SetTitleSize(0.07);
        hHodTdcInSpill_Any->GetYaxis()->SetTitleSize(0.07);
        hHodTdcInSpill_Any->GetXaxis()->SetTitleOffset(0.6);
        hHodTdcInSpill_Any->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hExtTdcInSpill_Any->Draw();
        hExtTdcInSpill_Any->SetMinimum(0.2);
        hExtTdcInSpill_Any->GetXaxis()->SetLabelSize(0.07);
        hExtTdcInSpill_Any->GetYaxis()->SetLabelSize(0.07);
        hExtTdcInSpill_Any->GetXaxis()->SetTitleSize(0.07);
        hExtTdcInSpill_Any->GetYaxis()->SetTitleSize(0.07);
        hExtTdcInSpill_Any->GetXaxis()->SetTitleOffset(0.6);
        hExtTdcInSpill_Any->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hTcTdcInSpill[0]->Draw();
        hTcTdcInSpill[0]->SetMinimum(0.2);
        hTcTdcInSpill[0]->GetXaxis()->SetLabelSize(0.07);
        hTcTdcInSpill[0]->GetYaxis()->SetLabelSize(0.07);
        hTcTdcInSpill[0]->GetXaxis()->SetTitleSize(0.07);
        hTcTdcInSpill[0]->GetYaxis()->SetTitleSize(0.07);
        hTcTdcInSpill[0]->GetXaxis()->SetTitleOffset(0.6);
        hTcTdcInSpill[0]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hTcTdcInSpill[1]->Draw();
        hTcTdcInSpill[1]->SetMinimum(0.2);
        hTcTdcInSpill[1]->GetXaxis()->SetLabelSize(0.07);
        hTcTdcInSpill[1]->GetYaxis()->SetLabelSize(0.07);
        hTcTdcInSpill[1]->GetXaxis()->SetTitleSize(0.07);
        hTcTdcInSpill[1]->GetYaxis()->SetTitleSize(0.07);
        hTcTdcInSpill[1]->GetXaxis()->SetTitleOffset(0.6);
        hTcTdcInSpill[1]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hHodHitMap->Draw("col");
        lHodBorderLine->Draw();
        hHodHitMap->GetXaxis()->SetLabelSize(0.06);
        hHodHitMap->GetYaxis()->SetLabelSize(0.06);
        hHodHitMap->GetXaxis()->SetTitleSize(0.06);
        hHodHitMap->GetYaxis()->SetTitleSize(0.06);
        hHodHitMap->GetXaxis()->SetTitleOffset(0.6);
        hHodHitMap->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hHodEntryByCh->Draw();
        hHodEntryByCh->GetXaxis()->SetLabelSize(0.06);
        hHodEntryByCh->GetYaxis()->SetLabelSize(0.06);
        hHodEntryByCh->GetXaxis()->SetTitleSize(0.06);
        hHodEntryByCh->GetYaxis()->SetTitleSize(0.06);
        hHodEntryByCh->GetXaxis()->SetTitleOffset(0.6);
        hHodEntryByCh->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hExtHitMap->Draw("col");
        lExtBorderLine->Draw();
        hExtHitMap->GetXaxis()->SetLabelSize(0.06);
        hExtHitMap->GetYaxis()->SetLabelSize(0.06);
        hExtHitMap->GetXaxis()->SetTitleSize(0.06);
        hExtHitMap->GetYaxis()->SetTitleSize(0.06);
        hExtHitMap->GetXaxis()->SetTitleOffset(0.6);
        hExtHitMap->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        TH1* hist = nullptr;
        switch (fChannelAlign) {
        case ChannelAlign::Raw:
          hist = hExtEntryByCh;
          hExtEntryByCh->Draw();
          break;
        case ChannelAlign::Projection:
          hist = hExtEntryByChBottom;
          hExtEntryByChBottom ->Draw();
          hExtEntryByChCenter1->Draw("same");
          hExtEntryByChCenter2->Draw("same");
          hExtEntryByChTop    ->Draw("same");
          hExtEntryByChBottom->SetMaximum(2.0 * hExtHitMap->GetBinContent(hExtHitMap->GetMaximumBin()));
          break;
        }
        if (hist) {
          hist->GetXaxis()->SetLabelSize(0.05);
          hist->GetYaxis()->SetLabelSize(0.05);
          hist->GetXaxis()->SetTitleSize(0.05);
          hist->GetYaxis()->SetTitleSize(0.05);
          hist->GetXaxis()->SetTitleOffset(0.6);
          hist->GetYaxis()->SetTitleOffset(0.5);
        }
      }

      if (fMonitorMode == MonitorMode::Coincidence) {
        if (fCanvas->cd(++padnumber)) {
          hExtMountain_Any->Draw("col");
          hExtMountain_Any->GetXaxis()->SetLabelSize(0.05);
          hExtMountain_Any->GetYaxis()->SetLabelSize(0.05);
          hExtMountain_Any->GetXaxis()->SetTitleSize(0.05);
          hExtMountain_Any->GetYaxis()->SetTitleSize(0.05);
          hExtMountain_Any->GetXaxis()->SetTitleOffset(0.6);
          hExtMountain_Any->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          hExtTdcInSync_Any->Draw();
          hExtTdcInSync_Any->SetMinimum(0.2);
          hExtTdcInSync_Any->GetXaxis()->SetLabelSize(0.05);
          hExtTdcInSync_Any->GetYaxis()->SetLabelSize(0.05);
          hExtTdcInSync_Any->GetXaxis()->SetTitleSize(0.05);
          hExtTdcInSync_Any->GetYaxis()->SetTitleSize(0.05);
          hExtTdcInSync_Any->GetXaxis()->SetTitleOffset(0.6);
          hExtTdcInSync_Any->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          gHitInSpill->Draw("AP");
          gHitInSpill->GetXaxis()->SetLabelSize(0.05);
          gHitInSpill->GetYaxis()->SetLabelSize(0.05);
          gHitInSpill->GetXaxis()->SetTitleSize(0.05);
          gHitInSpill->GetYaxis()->SetTitleSize(0.05);
          gHitInSpill->GetXaxis()->SetTitleOffset(0.6);
          gHitInSpill->GetYaxis()->SetTitleOffset(0.5);
        }

      } else if (fMonitorMode == MonitorMode::Offset) {
        if (fCanvas->cd(++padnumber)) {
          hMrSyncInterval[fMonitorBoard]->Draw("col");
          hMrSyncInterval[fMonitorBoard]->GetXaxis()->SetLabelSize(0.05);
          hMrSyncInterval[fMonitorBoard]->GetYaxis()->SetLabelSize(0.05);
          hMrSyncInterval[fMonitorBoard]->GetXaxis()->SetTitleSize(0.05);
          hMrSyncInterval[fMonitorBoard]->GetYaxis()->SetTitleSize(0.05);
          hMrSyncInterval[fMonitorBoard]->GetXaxis()->SetTitleOffset(0.6);
          hMrSyncInterval[fMonitorBoard]->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          hExtTdcOffset[fMonitorChannel]->Draw("col");
          hExtTdcOffset[fMonitorChannel]->GetXaxis()->SetLabelSize(0.03);
          hExtTdcOffset[fMonitorChannel]->GetYaxis()->SetLabelSize(0.03);
          hExtTdcOffset[fMonitorChannel]->GetXaxis()->SetTitleSize(0.03);
          hExtTdcOffset[fMonitorChannel]->GetYaxis()->SetTitleSize(0.03);
          hExtTdcOffset[fMonitorChannel]->GetXaxis()->SetTitleOffset(1.0);
          hExtTdcOffset[fMonitorChannel]->GetYaxis()->SetTitleOffset(0.9);
        }

      } else if (fMonitorMode == MonitorMode::Channel ||
                 fMonitorMode == MonitorMode::Ext ||
                 fMonitorMode == MonitorMode::Hod) {
        if (fCanvas->cd(++padnumber)) {
          hMountain->Draw("col");
          hMountain->GetXaxis()->SetLabelSize(0.05);
          hMountain->GetYaxis()->SetLabelSize(0.05);
          hMountain->GetXaxis()->SetTitleSize(0.05);
          hMountain->GetYaxis()->SetTitleSize(0.05);
          hMountain->GetXaxis()->SetTitleOffset(0.6);
          hMountain->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          hTdcInSync->Draw();
          hTdcInSync->SetMinimum(0.2);
          hTdcInSync->GetXaxis()->SetLabelSize(0.05);
          hTdcInSync->GetYaxis()->SetLabelSize(0.05);
          hTdcInSync->GetXaxis()->SetTitleSize(0.05);
          hTdcInSync->GetYaxis()->SetTitleSize(0.05);
          hTdcInSync->GetXaxis()->SetTitleOffset(0.6);
          hTdcInSync->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          hTdcInSpill->Draw();
          hTdcInSpill->GetXaxis()->SetLabelSize(0.05);
          hTdcInSpill->GetYaxis()->SetLabelSize(0.05);
          hTdcInSpill->GetXaxis()->SetTitleSize(0.05);
          hTdcInSpill->GetYaxis()->SetTitleSize(0.05);
          hTdcInSpill->GetXaxis()->SetTitleOffset(0.6);
          hTdcInSpill->GetYaxis()->SetTitleOffset(0.5);
        }
      }

      fCanvas->Modified();
      fCanvas->Update();
      gSystem->ProcessEvents();
    }

    void MonitorWindow::ClearLastSpill() {
      fCoinCount        = 0;
      fEventMatchNumber = -1;

      fTdcBuffer.clear();

      fLastExtData.clear();
      fLastHodData.clear();
      fLastTcData .clear();
      fLastBhData .clear();
      for (auto&& pair: fLastMrSyncData) {
        pair.second.Clear();
      }
      for (auto&& pair: fEventMatchData) {
        pair.second.clear();
      }

      hBhTdcInSpill[0]  ->Reset();
      hBhTdcInSpill[1]  ->Reset();
      hHodTdcInSpill_Any->Reset();
      hExtTdcInSpill_Any->Reset();
      hTcTdcInSpill[0]  ->Reset();
      hTcTdcInSpill[1]  ->Reset();
      hHodHitMap        ->Reset();
      hHodEntryByCh     ->Reset();
      hExtHitMap        ->Reset();
      hExtEntryByCh     ->Reset();
      hExtMountain_Any  ->Reset();
      hExtTdcInSync_Any ->Reset();
      if (fMonitorMode == MonitorMode::Channel ||
          fMonitorMode == MonitorMode::Ext     ||
          fMonitorMode == MonitorMode::Hod     ) {
        hMountain       ->Reset();
        hTdcInSync      ->Reset();
        hTdcInSpill     ->Reset();
      } else if (fMonitorMode == MonitorMode::Offset) {
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncInterval[ch]->Reset();
        }
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcOffset[ch]->Reset();
        }
      }
    }

    void MonitorWindow::FillCoincidence(const TdcData& extData) {
      const std::size_t extCh = ExtinctionDetector::GetChannel(extData.Channel);
      const Double_t    time  = extData.Time;
      const Int_t       board = extData.Board;

      Bool_t coincidence[CoinOffset::N];
      if (fCyclicCoincidence) {
        for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
          const std::size_t i = bhCh + CoinOffset::BH;
          coincidence[i] = !fContains[extCh][i];
        }
        for (auto&& lastData : fLastBhData) {
          const std::size_t bhCh = BeamlineHodoscope::GetChannel(lastData.Channel);
          const std::size_t i = bhCh + CoinOffset::BH;
          if (!coincidence[i]) {
            const Double_t dt0 = lastData.Time - time;
            Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[         board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fMrSyncInterval[         board] *  extData.TimePerTdc,
            //        thre = -0.6 * fMrSyncIntervalAverage *  extData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::BH, coincidence + CoinOffset::BH + BeamlineHodoscope::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          coincidence[i] = !fContains[extCh][i];
        }
        for (auto&& lastData : fLastHodData) {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          if (!coincidence[i]) {
            const Double_t dt0 = lastData.Time - time;
            Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[         board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fMrSyncInterval[         board] *  extData.TimePerTdc,
            //        thre = -0.6 * fMrSyncIntervalAverage *  extData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              break;
            }
          }
        }

        for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
          const std::size_t i = tcCh + CoinOffset::TC;
          coincidence[i] = !fContains[extCh][i];
        }
        for (auto&& lastData : fLastTcData) {
          const std::size_t tcCh = TimingCounter::GetChannel(lastData.Channel);
          const std::size_t i = tcCh + CoinOffset::TC;
          if (!coincidence[i]) {
            const Double_t dt0 = lastData.Time - time;
            Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[         board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fMrSyncInterval[         board] *  extData.TimePerTdc,
            //        thre = -0.6 * fMrSyncIntervalAverage *  extData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean  = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::TC, coincidence + CoinOffset::TC + TimingCounter::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
          ++fCoinCount;
          auto          lastData = fLastMrSyncData[board];
          const Long64_t     tdc = extData.Tdc;
          const Long64_t syncTdc = lastData.Tdc;
          const Long64_t    dtdc = tdc - syncTdc;
          if        (dtdc > +1.05 * fMrSyncInterval[board]) {
            const Long64_t dtdc2 = std::fmod(dtdc, fMrSyncInterval[board]);
            hExtTdcInSync_Any->Fill(dtdc2);
            hExtMountain_Any ->Fill(dtdc2, time / msec);
          } else if (dtdc < -0.05 * fMrSyncInterval[board]) {
            const Long64_t dtdc2 = std::fmod(dtdc, fMrSyncInterval[board]) + fMrSyncInterval[board];
            hExtTdcInSync_Any->Fill(dtdc2);
            hExtMountain_Any ->Fill(dtdc2, time / msec);
          } else {
            hExtTdcInSync_Any->Fill(dtdc);
            hExtMountain_Any ->Fill(dtdc, time / msec);
          }
        }

      } else {
        for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
          const std::size_t i = bhCh + CoinOffset::BH;
          coincidence[i] = !fContains[extCh][i];
        }
        for (auto&& lastData : fLastBhData) {
          const std::size_t bhCh = BeamlineHodoscope::GetChannel(lastData.Channel);
          const std::size_t i = bhCh + CoinOffset::BH;
          if (!coincidence[i]) {
            const Double_t dt   = lastData.Time - time;
            const Double_t mean = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::BH, coincidence + CoinOffset::BH + BeamlineHodoscope::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          coincidence[i] = !fContains[extCh][i];
        }
        for (auto&& lastData : fLastHodData) {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          if (!coincidence[i]) {
            const Double_t dt   = lastData.Time - time;
            const Double_t mean = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              break;
            }
          }
        }

        for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
          const std::size_t i = tcCh + CoinOffset::TC;
          coincidence[i] = !fContains[extCh][i];
        }
        for (auto&& lastData : fLastTcData) {
          const std::size_t tcCh = TimingCounter::GetChannel(lastData.Channel);
          const std::size_t i = tcCh + CoinOffset::TC;
          if (!coincidence[i]) {
            const Double_t dt   = lastData.Time - time;
            const Double_t mean = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::TC, coincidence + CoinOffset::TC + TimingCounter::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
          ++fCoinCount;
          auto          lastData = fLastMrSyncData[board];
          const Long64_t     tdc = extData.Tdc;
          const Long64_t syncTdc = lastData.Tdc;
          hExtTdcInSync_Any->Fill(tdc - syncTdc);
          hExtMountain_Any ->Fill(tdc - syncTdc, time / msec);
        }
      }
    }

    inline void MonitorWindow::FillCoincidences(const std::vector<TdcData>& extData) {
      for (auto&& data : extData) {
        FillCoincidence(data);
      }
    }

    void MonitorWindow::FillCoincidence2(const TdcData& tdcData) {
      const Double_t    time  = tdcData.Time;
      const Int_t       board = tdcData.Board;
      std::size_t tdcCh = 0;
      std::size_t tdcI  = 0;
      if        (BeamlineHodoscope::Contains  (tdcData.Channel)) {
        const std::size_t bhCh = BeamlineHodoscope::GetChannel(tdcData.Channel);
        tdcI  =  bhCh + CoinOffset    ::BH;
        tdcCh =  bhCh + ContainsOffset::BH;
      } else if (Hodoscope        ::Contains  (tdcData.Channel)) {
        const std::size_t hodCh = 0;
        tdcI  = hodCh + CoinOffset    ::Hod;
        tdcCh = hodCh + ContainsOffset::Hod;
      } else if (TimingCounter    ::Contains  (tdcData.Channel)) {
        const std::size_t tcCh = TimingCounter::GetChannel(tdcData.Channel);
        tdcI  =  tcCh + CoinOffset    ::TC;
        tdcCh =  tcCh + ContainsOffset::TC;
      } else {
        return;
      }

      Bool_t coincidence[CoinOffset::N];
      if (fCyclicCoincidence) {
        for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
          const std::size_t i = bhCh + CoinOffset::BH;
          coincidence[i] = !fContains[tdcCh][i];
        }
        for (auto&& lastData : fLastBhData) {
          const std::size_t bhCh = BeamlineHodoscope::GetChannel(lastData.Channel);
          const std::size_t i = bhCh + CoinOffset::BH;
          if (!coincidence[i]) {
            const Double_t dt0 = lastData.Time - time;
            Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[         board] *  tdcData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fMrSyncInterval[         board] *  tdcData.TimePerTdc,
            //        thre = -0.6 * fMrSyncIntervalAverage *  tdcData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fCoinDiffs[tdcCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::BH, coincidence + CoinOffset::BH + BeamlineHodoscope::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          coincidence[i] = !fContains[tdcCh][i];
        }
        for (auto&& lastData : fLastHodData) {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          if (!coincidence[i]) {
            const Double_t dt0 = lastData.Time - time;
            Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[         board] *  tdcData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fMrSyncInterval[         board] *  tdcData.TimePerTdc,
            //        thre = -0.6 * fMrSyncIntervalAverage *  tdcData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fCoinDiffs[tdcCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              break;
            }
          }
        }

        for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
          const std::size_t i = tcCh + CoinOffset::TC;
          coincidence[i] = !fContains[tdcCh][i];
        }
        for (auto&& lastData : fLastTcData) {
          const std::size_t tcCh = TimingCounter::GetChannel(lastData.Channel);
          const std::size_t i = tcCh + CoinOffset::TC;
          if (!coincidence[i]) {
            const Double_t dt0 = lastData.Time - time;
            Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[         board] *  tdcData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fMrSyncInterval[         board] *  tdcData.TimePerTdc,
            //        thre = -0.6 * fMrSyncIntervalAverage *  tdcData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean  = fCoinDiffs[tdcCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::TC, coincidence + CoinOffset::TC + TimingCounter::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {

          std::vector<TdcData> hits;
          for (auto&& extData : fLastExtData) {
            const std::size_t extCh = ExtinctionDetector::GetChannel(extData.Channel);
            const Double_t dt0 = time - extData.Time;
            Double_t dsync = fLastMrSyncData[board].Time - fLastMrSyncData[extData.Board].Time;
            for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                 dsync += fMrSyncInterval[        board] *  tdcData.TimePerTdc);
            for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                 dsync -= fMrSyncInterval[extData.Board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            const Double_t mean  = fCoinDiffs[extCh][tdcI];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              hits.push_back(extData);
            }
          }

          if (hits.size()) {
            ++fCoinCount;
            Double_t tdcSum = 0.0;
            for (auto& hitData : hits) {
              auto          lastData = fLastMrSyncData[board];
              const Long64_t     tdc = hitData.Tdc;
              const Long64_t syncTdc = lastData.Tdc;
              const Long64_t    dtdc = tdc - syncTdc;
              if        (dtdc > +1.05 * fMrSyncInterval[board]) {
                const Long64_t dtdc2 = std::fmod(dtdc, fMrSyncInterval[board]);
                tdcSum += dtdc2;
              } else if (dtdc < -0.05 * fMrSyncInterval[board]) {
                const Long64_t dtdc2 = std::fmod(dtdc, fMrSyncInterval[board]) + fMrSyncInterval[board];
                tdcSum += dtdc2;
              } else {
                tdcSum += dtdc;
              }
            }
            hExtTdcInSync_Any->Fill(tdcSum / hits.size());
            hExtMountain_Any ->Fill(tdcSum / hits.size(), time / msec);
          }
        }

      } else {
        for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
          const std::size_t i = bhCh + CoinOffset::BH;
          coincidence[i] = !fContains[tdcCh][i];
        }
        for (auto&& lastData : fLastBhData) {
          const std::size_t bhCh = BeamlineHodoscope::GetChannel(lastData.Channel);
          const std::size_t i = bhCh + CoinOffset::BH;
          if (!coincidence[i]) {
            const Double_t dt0   = lastData.Time - time;
            const Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            const Double_t dt    = dt0 - dsync;
            const Double_t mean  = fCoinDiffs[tdcCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::BH, coincidence + CoinOffset::BH + BeamlineHodoscope::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          coincidence[i] = !fContains[tdcCh][i];
        }
        for (auto&& lastData : fLastHodData) {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          if (!coincidence[i]) {
            const Double_t dt0   = lastData.Time - time;
            const Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            const Double_t dt    = dt0 - dsync;
            const Double_t mean  = fCoinDiffs[tdcCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              break;
            }
          }
        }

        for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
          const std::size_t i = tcCh + CoinOffset::TC;
          coincidence[i] = !fContains[tdcCh][i];
        }
        for (auto&& lastData : fLastTcData) {
          const std::size_t tcCh = TimingCounter::GetChannel(lastData.Channel);
          const std::size_t i = tcCh + CoinOffset::TC;
          if (!coincidence[i]) {
            const Double_t dt0   = lastData.Time - time;
            const Double_t dsync = fLastMrSyncData[lastData.Board].Time - fLastMrSyncData[board].Time;
            const Double_t dt    = dt0 - dsync;
            const Double_t mean  = fCoinDiffs[tdcCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coincidence[i] = true;
              if (std::all_of(coincidence + CoinOffset::TC, coincidence + CoinOffset::TC + TimingCounter::NofChannels, [](Bool_t b) { return b; })) {
                break;
              }
            }
          }
        }

        if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
          std::vector<TdcData> hits;
          for (auto&& extData : fLastExtData) {
            const std::size_t extCh = ExtinctionDetector::GetChannel(extData.Channel);
            const Double_t dt0   = time - extData.Time;
            const Double_t dsync = fLastMrSyncData[board].Time - fLastMrSyncData[extData.Board].Time;
            const Double_t dt    = dt0 - dsync;
            const Double_t mean  = fCoinDiffs[extCh][tdcI];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              hits.push_back(extData);
            }
          }

          if (hits.size()) {
            ++fCoinCount;
            Double_t tdcSum = 0.0;
            for (auto& hitData : hits) {
              auto          lastData = fLastMrSyncData[board];
              const Long64_t     tdc = hitData.Tdc;
              const Long64_t syncTdc = lastData.Tdc;
              const Long64_t    dtdc = tdc - syncTdc;
              tdcSum += dtdc;
            }
            hExtTdcInSync_Any->Fill(tdcSum / hits.size());
            hExtMountain_Any ->Fill(tdcSum / hits.size(), time / msec);
          }
        }

      }
    }
    
    std::vector<TdcData> MonitorWindow::CollectCoinExtData(const TdcData& tdcData, std::size_t i) {
      const Double_t time = tdcData.Time;

      std::vector<TdcData> coinExtData;
      if (fCyclicCoincidence) {
        for (auto&& lastData : fLastExtData) {
          const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
          if (fContains[extCh][i]) {
              const Double_t dt0 = time - lastData.Time;
              Double_t dsync = fLastMrSyncData[tdcData.Board].Time - fLastMrSyncData[lastData.Board].Time;
              for (Double_t thre = -0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync < thre;
                   dsync += fMrSyncInterval[ tdcData.Board] *  tdcData.TimePerTdc);
              for (Double_t thre = +0.5 * fMrSyncIntervalAverage * fProvider.GetTimePerTdc(); dsync > thre;
                   dsync -= fMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
              const Double_t dt = dt0 - dsync;
            // Double_t dt = (time - fLastMrSyncData[tdcData.Board].Time) - (lastData.Time - fLastMrSyncData[lastData.Board].Time);
            // for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc, thre = -0.6 * fMrSyncIntervalAverage; dt < thre; dt += dint);
            // for (Double_t dint = fMrSyncInterval[ tdcData.Board] *  tdcData.TimePerTdc, thre = +0.6 * fMrSyncIntervalAverage; dt > thre; dt -= dint);
            const Double_t mean  = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coinExtData.push_back(lastData);
            }
          }
        }
      } else {
        for (auto&& lastData : fLastExtData) {
          const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
          if (fContains[extCh][i]) {
            const Double_t dt   = time - lastData.Time;
            const Double_t mean = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coinExtData.push_back(lastData);
            }
          }
        }
      }
      return coinExtData;
    }

    std::size_t MonitorWindow::RemoveOldTdc(std::vector<TdcData>* lastData, Double_t time) {
      for (std::size_t i = 0, n = lastData->size(); i < n; ++i) {
        if (std::abs(time - lastData->at(0).Time) > fHistoryWidth) {
          lastData->erase(lastData->begin());
        } else {
          break;
        }
      }
      return lastData->size();
    }

    Int_t MonitorWindow::UpdatePlots(const std::map<int, std::string>& ifilenames,
                                     const std::function<TDatime(const std::string&)>& parser) {
      const clock_t startClock = clock();

      std::cout << "Initialize decoder" << std::endl;
      Packet_t packet = 0;
      std::map<Int_t, Decoder> decoders;
      for (auto&& pair : ifilenames) {
        const Int_t board = pair.first;
        if (fTimePerTdc[board] == 0) {
          decoders[board].Data.SetTimePerTdc(fProvider.GetTimePerTdc());
        } else {
          decoders[board].Data.SetTimePerTdc(fTimePerTdc[board]);
        }
      }

      std::cout << "Initialzie history" << std::endl;
      ClearLastSpill();

      std::cout << "Open file" << std::endl;
      std::map<Int_t, std::ifstream> ifiles;
      std::vector<TDatime>           datimes;
      for (auto&& pair : ifilenames) {
        const Int_t       board     = pair.first;
        const std::string ifilename = pair.second;
        std::cout << " - " << ifilename << std::endl;
        if (parser) {
          datimes.push_back(parser(ifilename));
        }

        ifiles[board].open(ifilename, std::ios::binary);
        if (!ifiles[board]) {
          std::cerr << "[error] input file is not opened, " << ifilename << std::endl;
          return 1;
        }
      }
      const UInt_t averageTime = Tron::Linq::From(datimes)
        .Select([](TDatime datime) -> ULong64_t { return datime.Convert(); })
        .Average();
      fDate.Set(averageTime);

      {
        Int_t targetBoard = 0;
        std::map<Int_t, std::size_t> mrcount;
        std::map<Int_t, Long64_t>    mrtdc;
        std::map<Int_t, Bool_t>      gateStarted;
        std::map<Int_t, Bool_t>      gateEnded;
        std::map<Int_t, Bool_t>      fileEnded;
        std::map<Int_t, ULong64_t>   lastTdcTags;
        for (auto&& pair : ifilenames) {
          targetBoard = pair.first;
          gateStarted[pair.first] = false;
          gateEnded  [pair.first] = false;
          fileEnded  [pair.first] = false;
          lastTdcTags[pair.first] = 0;
          fLastMrSyncData[pair.first];
          fEventMatchData[pair.first];
        }

        for (std::size_t count = 0UL; !fIsTerminated;) {
          Decoder&   decoder    = decoders   [targetBoard];
          ULong64_t& lastTdcTag = lastTdcTags[targetBoard];
          for (std::size_t ibuf = 0; ibuf < fBufferSize; ++ibuf, ++count) {
            if (count % 1000000UL == 0) {
              std::cout << ">> " << count << std::endl;
            }
            if (count % 3000000UL == 0) {
              DrawPlots();
            }

            if (!decoder.Read(ifiles[targetBoard], &packet)) {
              // std::cout << "[info] detect file end @ " << targetBoard << std::endl;
              fileEnded[targetBoard] = true;
              break;

            } else {
              const Int_t dtype = decoder.Data.Type;
              if (dtype == DataType::Header) {
                // std::cout << "[info] detect header @ " << targetBoard << std::endl;
                break;

              } else if (dtype == DataType::GateStart) {
                // std::cout << "[info] detect gate start @ " << targetBoard << std::endl;
                gateStarted[targetBoard] = true;
                break;

              } else if (dtype == DataType::GateEnd) {
                // std::cout << "[info] detect gate end @ " << targetBoard << std::endl;
                gateEnded[targetBoard] = true;
                break;

              } else {
                // std::cout << "[info] detect data @ " << targetBoard << std::endl;
                auto tdcData = decoder.Data.GetTdcData(targetBoard);

                for (auto&& data : tdcData) {
                  if (MrSync::Contains(data.Channel)) {
                    mrcount[data.Board]++;
                    mrtdc  [data.Board] = data.Tdc;
                  }
                }

                ULong64_t tdcTag = 0;
                for (auto&& data : tdcData) {
                  tdcTag             = data.GetTdcTag(mrcount[data.Board], mrtdc[data.Board]);
                  fTdcBuffer[tdcTag] = data;
                }
                if (ibuf < fBufferSize - fBufferMargin) {
                  lastTdcTag = tdcTag;
                }
              }
            }
          }

          {
            Int_t boardNotLoadedYet = -1;
            for (auto&& pair : lastTdcTags) {
              const Int_t board = pair.first;
              if (pair.second == 0 &&
                  !gateStarted[board] && !gateEnded[board] && !fileEnded[board]) {
                boardNotLoadedYet = board;
                break;
              }
            }
            if (boardNotLoadedYet >= 0) {
              // std::cout << "[info] there is a board not loaded yet @ " << boardNotLoadedYet << std::endl;
              targetBoard = boardNotLoadedYet;
              continue;
            }
          }

          for (auto&& itr = fTdcBuffer.begin(); itr != fTdcBuffer.end(); itr = fTdcBuffer.begin()) {
            auto sortedTdcPair = *itr;
            const ULong64_t tdcTag = sortedTdcPair.first;
            const TdcData&  data   = sortedTdcPair.second;
            fTdcBuffer.erase(itr);

            const Int_t    board         = data.Board;
            const Int_t    globalChannel = data.Channel;
            const Double_t time          = data.Time;

            if ((fMonitorMode == MonitorMode::Channel && fMonitorChannel           == globalChannel ) ||
                (fMonitorMode == MonitorMode::Ext     && ExtinctionDetector::Contains(globalChannel)) ||
                (fMonitorMode == MonitorMode::Hod     && Hodoscope         ::Contains(globalChannel))) {
              const Long64_t     tdc = data.Tdc;
              const Long64_t syncTdc = fLastMrSyncData[board].Tdc;
              hTdcInSpill->Fill(time / msec);
              hTdcInSync ->Fill(tdc - syncTdc);
              hMountain  ->Fill(tdc - syncTdc, time / msec);
            }

            if (ExtinctionDetector::Contains(globalChannel)) {
              const Int_t ch = ExtinctionDetector::GetChannel(globalChannel);

              hExtEntryByCh     ->Fill(ch);
              ExtinctionDetector::Fill(hExtHitMap, ch);
              hExtTdcInSpill_Any->Fill(time / msec);

              if (fMonitorMode == MonitorMode::Offset) {
                const Long64_t tdc = data.Tdc;
                for (auto&& lastData : fLastHodData) {
                  auto lastCh = Hodoscope::GetChannel(lastData.Channel);
                  const Double_t dtdc0  = lastData.Tdc - tdc;
                  const Double_t dsync  = fLastMrSyncData[lastData.Board].Tdc - fLastMrSyncData[board].Tdc;
                  const Double_t dtdc   = dtdc0 - dsync;
                  hExtTdcOffset[ch]->Fill(lastCh + CoinOffsetX::Hod, dtdc);
                }
                for (auto&& lastData : fLastTcData) {
                  auto lastCh = TimingCounter::GetChannel(lastData.Channel);
                  const Double_t dtdc0 = lastData.Tdc - tdc;
                  const Double_t dsync = fLastMrSyncData[lastData.Board].Tdc - fLastMrSyncData[board].Tdc;
                  const Double_t dtdc  = dtdc0 - dsync;
                  hExtTdcOffset[ch]->Fill(lastCh + CoinOffsetX::TC, dtdc);
                }
                for (auto&& lastData : fLastBhData) {
                  auto lastCh = BeamlineHodoscope::GetChannel(lastData.Channel);
                  const Double_t dtdc0 = lastData.Tdc - tdc;
                  const Double_t dsync = fLastMrSyncData[lastData.Board].Tdc - fLastMrSyncData[board].Tdc;
                  const Double_t dtdc  = dtdc0 - dsync;
                  hExtTdcOffset[ch]->Fill(lastCh + CoinOffsetX::BH, dtdc);
                }
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                // FillCoincidence(data);
              }

              if (RemoveOldTdc(&fLastExtData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastExtData reaches " << kHistLimit << std::endl;
                return 1;
              }
              fLastExtData.push_back(data);

            } else if (Hodoscope::Contains(globalChannel)) {
              const Int_t ch = Hodoscope::GetChannel(globalChannel);

              hHodEntryByCh     ->Fill(ch);
              Hodoscope         ::Fill(hHodHitMap, ch);
              hHodTdcInSpill_Any->Fill(time / msec);

              if (fMonitorMode == MonitorMode::Offset) {
                const Long64_t tdc = data.Tdc;
                for (auto&& lastData : fLastExtData) {
                  auto lastCh = ExtinctionDetector::GetChannel(lastData.Channel);
                  const Double_t ddc0  = tdc - lastData.Tdc;
                  const Double_t dsync = fLastMrSyncData[board].Tdc - fLastMrSyncData[lastData.Board].Tdc;
                  const Double_t dtdc  = ddc0 - dsync;
                  hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::Hod, dtdc);
                }
                fLastHodData.push_back(data);
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                fLastHodData.push_back(data);
                // FillCoincidences(CollectCoinExtData(data, CoinOffset::Hod));
                FillCoincidence2(data);
              }

              if (RemoveOldTdc(&fLastHodData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastHodData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (TimingCounter::Contains(globalChannel)) {
              const Int_t ch = TimingCounter::GetChannel(globalChannel);

              hTcTdcInSpill[ch]->Fill(time / msec);

              if (fMonitorMode == MonitorMode::Offset) {
                const Long64_t tdc = data.Tdc;
                for (auto&& lastData : fLastExtData) {
                  auto lastCh = ExtinctionDetector::GetChannel(lastData.Channel);
                  const Double_t dtdc0 = tdc - lastData.Tdc;
                  const Double_t dsync = fLastMrSyncData[board].Tdc - fLastMrSyncData[lastData.Board].Tdc;
                  const Double_t dtdc  = dtdc0 - dsync;
                  hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::TC, dtdc);
                }
                fLastTcData.push_back(data);
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                fLastTcData.push_back(data);
                // FillCoincidences(CollectCoinExtData(data, ch + CoinOffset::TC));
                FillCoincidence2(data);
              }

              if (RemoveOldTdc(&fLastTcData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastTcData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (BeamlineHodoscope::Contains(globalChannel)) {
              const Int_t ch = BeamlineHodoscope::GetChannel(globalChannel);

              hBhTdcInSpill[ch]->Fill(time / msec);

              if (fMonitorMode == MonitorMode::Offset) {
                const Long64_t tdc = data.Tdc;
                for (auto&& lastData : fLastExtData) {
                  auto lastCh   = ExtinctionDetector::GetChannel(lastData.Channel);
                  const Double_t dtdc0 = tdc - lastData.Tdc;
                  const Double_t dsync = fLastMrSyncData[board].Tdc - fLastMrSyncData[lastData.Board].Tdc;
                  const Double_t dtdc  = dtdc0 - dsync;
                  hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::BH, dtdc);
                }
                fLastBhData.push_back(data);
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                fLastBhData.push_back(data);
                // FillCoincidences(CollectCoinExtData(data, ch + CoinOffset::BH));
                FillCoincidence2(data);
              }

              if (RemoveOldTdc(&fLastBhData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastBhData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (MrSync::Contains(globalChannel)) {
              if (fMonitorMode == MonitorMode::Offset) {
                decltype(fLastMrSyncData)::const_iterator itr;
                if ((itr = fLastMrSyncData.find(data.Board)) != fLastMrSyncData.end()) {
                  hMrSyncInterval[data.Board]->Fill(data.Tdc - itr->second.Tdc);
                }
              }
              fLastMrSyncData[board] = data;

            } else if (EventMatch::Contains(globalChannel)) {
              fEventMatchData[board].push_back(data);
            }

            if (tdcTag == lastTdcTags[board] &&
                !gateEnded[board] && !fileEnded[board]) {
              // std::cout << "[info] detect read last buffer @ " << targetBoard << std::endl;
              targetBoard = board;
              break;
            }

          }

          if        (IsAllOfSecondsTrue(gateStarted)) {
              std::cout << "[info] begin of spill " << fSpillCount << std::endl;
              ClearLastSpill();

              FillToSeconds(&gateStarted, false);
              FillToSeconds(&lastTdcTags, 0ULL );

          } else if (IsAllOfSecondsTrue(gateEnded)) {
            std::cout << "[info] end of spill " << fSpillCount << std::endl;

            fEventMatchNumber = fProvider.DecodeEventMatchNumber(fEventMatchData.begin()->second);
            for (auto&& pair : ifilenames) {
              const Int_t board = pair.first;
              const Int_t emcount = fProvider.DecodeEventMatchNumber(fEventMatchData[board]);
              if (emcount != fEventMatchNumber) {
                std::cout << "[warning] conflict EMCount, " << fEventMatchNumber << " <--> " << emcount << std::endl;
              }
            }

            const Int_t np = fSpillCount % fSpillLimit;
            gHitInSpill->SetPoint     (np, fSpillCount,     fCoinCount      );
         // gHitInSpill->SetPointError(np, 0.0, TMath::Sqrt(fCoinCount     ));
            ++fSpillCount;

            for (Int_t xbin = 1, nbinsx = hExtHitMap->GetNbinsX(); xbin <= nbinsx; ++xbin) {
              hExtEntryByChBottom ->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 1));
              hExtEntryByChCenter1->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 2));
              hExtEntryByChCenter2->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 3));
              hExtEntryByChTop    ->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 4));
            }

            if (fMonitorMode == MonitorMode::Offset) {

              std::ofstream mrSyncFile;
              if (fMrSyncIntervalFilename.size()) {
                mrSyncFile.open(fMrSyncIntervalFilename);
                if (!mrSyncFile) {
                  std::cout << "[warning] mrSyncInterval file is not opened, " << fMrSyncIntervalFilename << std::endl;
                } else {
                  std::cout << "[info] write mrSyncInterval to " << fMrSyncIntervalFilename << std::endl;
                }
              }

              if (!mrSyncFile) {
                std::cout << "_____ MrSync _____" << std::endl;
              }
              for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
                if (!hMrSyncInterval[ch]->GetEntries()) {
                  continue;
                }

                const Int_t xbins = hMrSyncInterval[ch]->GetNbinsX();

                Double_t yxsum = 0.0, ysum = 0.0;
                for (Int_t xbin = 1; xbin <= xbins; ++xbin) {
                  const Double_t y = hMrSyncInterval[ch]->GetBinContent(xbin);
                  yxsum += y * hMrSyncInterval[ch]->GetBinCenter(xbin);
                  ysum  += y;
                }
                if (ysum) {
                  const Double_t xmean = yxsum / ysum;
                  if (!mrSyncFile) {
                    std::cout << ch << "\t" << xmean << std::endl;
                  } else {
                    mrSyncFile << ch << "\t" << xmean << std::endl;
                  }
                }
              }
              if (!mrSyncFile) {
                std::cout << "^^^^^^^^^^^^^^^^^^" << std::endl;
              }

              std::ofstream offsetFile;
              if (fCoinDiffsFilename.size()) {
                offsetFile.open(fCoinDiffsFilename);
                if (!offsetFile) {
                  std::cout << "[warning] coinDiffs file is not opened, " << fCoinDiffsFilename << std::endl;
                } else {
                  std::cout << "[info] write coinDiffs to " << fCoinDiffsFilename << std::endl;
                }
              }

              if (!offsetFile) {
                std::cout << "_____ Offset _____" << std::endl;
              }
              for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
                if (!hExtTdcOffset[ch]->GetEntries()) {
                  continue;
                }

                const Int_t ybins = hExtTdcOffset[ch]->GetNbinsY();

                // Beamline Hodoscope 1, 2
                for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
                  const std::size_t i    = bhCh + CoinOffset::BH;
                  const Int_t       xbin = bhCh + CoinOffsetX::BH + 1;
                  Double_t maxsum = 0.0/*, maxtdc = 0.0*/; Int_t maxbin = 0;
                  for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                    const Double_t sum = hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    if (sum > maxsum) {
                      maxbin = ybin;
                      // maxtdc = hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      maxsum = sum;
                    }
                  }
                  if (maxsum) {
                    Double_t ysum = 0.0, yxsum = 0.0;
                    for (Int_t ybin = maxbin - 1; ybin <= maxbin + 1; ++ybin) {
                      yxsum += hExtTdcOffset[ch]->GetBinContent(xbin, ybin) * hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      ysum  += hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    }
                    const Double_t maxtdc = yxsum / ysum;
                    if (!offsetFile) {
                      std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                    } else {
                      offsetFile << ch << "\t" << i << "\t" << maxtdc << std::endl;
                    }
                  }
                }

                // Hodoscope
                {
                  const std::size_t i     = CoinOffset::Hod;
                  const Int_t       xbin1 = CoinOffsetX::Hod + 1;
                  const Int_t       xbin2 = CoinOffsetX::Hod + Hodoscope::NofChannels;
                  Double_t maxsum = 0.0/*, maxtdc = 0.0*/; Int_t maxbin = 0;
                  for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                    Double_t sum = 0.0;
                    for (Int_t xbin = xbin1; xbin <= xbin2; ++xbin) {
                      sum += hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    }
                    if (sum > maxsum) {
                      maxbin = ybin;
                      // maxtdc = hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      maxsum = sum;
                    }
                  }
                  if (maxsum) {
                    Double_t ysum = 0.0, yxsum = 0.0;
                    for (Int_t xbin = xbin1; xbin <= xbin2; ++xbin) {
                      for (Int_t ybin = maxbin - 1; ybin <= maxbin + 1; ++ybin) {
                        yxsum += hExtTdcOffset[ch]->GetBinContent(xbin, ybin) * hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                        ysum  += hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                      }
                    }
                    const Double_t maxtdc = yxsum / ysum;
                    if (!offsetFile) {
                      std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                    } else {
                      offsetFile << ch << "\t" << i << "\t" << maxtdc << std::endl;
                    }
                  }
                }

                // Timing Counter 1, 2
                for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
                  const std::size_t i    = tcCh + CoinOffset::TC;
                  const Int_t       xbin = tcCh + CoinOffsetX::TC + 1;
                  Double_t maxsum = 0.0/*, maxtdc = 0.0*/; Int_t maxbin = 0;
                  for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                    const Double_t sum = hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    if (sum > maxsum) {
                      maxbin = ybin;
                      // maxtdc = hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      maxsum = sum;
                    }
                  }
                  if (maxsum) {
                    Double_t ysum = 0.0, yxsum = 0.0;
                    for (Int_t ybin = maxbin - 1; ybin <= maxbin + 1; ++ybin) {
                      yxsum += hExtTdcOffset[ch]->GetBinContent(xbin, ybin) * hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      ysum  += hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    }
                    const Double_t maxtdc = yxsum / ysum;
                    if (!offsetFile) {
                      std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                    } else {
                      offsetFile << ch << "\t" << i << "\t" << maxtdc << std::endl;
                    }
                  }
                }

              }
              // std::cout << "^^^^^^^^^^^^^^^^^^" << std::endl;
            }

            DrawPlots();

            FillToSeconds(&gateEnded  , false);
            FillToSeconds(&lastTdcTags, 0ULL );

          } else if (IsAllOfSecondsTrue(fileEnded)) {
            // std::cout << "[info] detect file end" << std::endl;
            break;

          }

        }
      }

      std::cout << "Close files" << std::endl;
      for (auto&& pair : ifiles) {
        pair.second.close();
      }

      const clock_t stopClock = clock();
      std::cout << "time: " << (double)(stopClock - startClock) / CLOCKS_PER_SEC << " sec\n";

      return 0;
    }

  }

}

#endif
