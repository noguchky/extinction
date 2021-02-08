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

    enum class ChannelAlign {
      Raw, Projection,
    };

    enum class MonitorMode {
      Channel, Offset, Coincidence,
    };

    class MonitorWindow {
    private:
      const std::size_t kHistLimit     = 10000;
      const Double_t    kLastThreshold = 100.0 * nsec;
      const Double_t    kCoinTimeWidth =  25.0 * nsec;
      const Long64_t    kSpillLimit    = 1 * 60 * 60 / 5;
      const std::size_t kBufferSize    = 1000;

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

      Long64_t                     fSpillCount          = 0;
      Long64_t                     fCoinCount           = 0;

      std::map<Int_t, Double_t>    fTimePerTdc;
      std::map<Int_t, Double_t>    fMrSyncInterval;
      Double_t                     fMrSyncIntervalAverage;
      std::map<ULong64_t, TdcData> fTdcBuffer;

      std::vector<TdcData>         fLastExtData;
      std::vector<TdcData>         fLastHodData;
      std::vector<TdcData>         fLastTcData;
      std::vector<TdcData>         fLastBhData;

      std::map<Int_t, TdcData>     fLastMrSyncData;

      Bool_t                       fIsTerminated        = false;
      ChannelAlign                 fChannelAlign        = ChannelAlign::Raw;
      MonitorMode                  fMonitorMode         = MonitorMode::Coincidence;
      Int_t                        fMonitorChannel      = -1; // global channel
      Int_t                        fMonitorBoard        = -1;
      Bool_t                       fCyclicCoincidence   = true;

    public:
      MonitorWindow();
      ~MonitorWindow();

      inline void          SetTimePerTdc(const std::map<Int_t, Double_t>& map) { fTimePerTdc = map; }
      inline void          SetMrSyncInterval(const std::map<Int_t, Double_t>& map) {
        fMrSyncInterval        = map;
        fMrSyncIntervalAverage = Tron::Linq::From(map)
          .Select([](const std::pair<Int_t, Double_t>& pair) { return pair.second; })
          .Average();
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
          }
          break;
        case MonitorMode::Offset:
          if (ExtinctionDetector::Contains(option)) {
            fMonitorMode    = mode;
            fMonitorChannel = option;
            fMonitorBoard   = ChannelMapWithBoard::Board[option];
            return true;
          }
          break;
        case MonitorMode::Coincidence:
          fMonitorMode = mode;
          return true;
        }
        return false;
      }
      inline MonitorMode   GetMonitorMode() const { return fMonitorMode; }
      inline Int_t         GetMonitorChannel() const { return fMonitorChannel; }

      Int_t                LoadOffset(const std::string& ffilename);
      void                 InitializeWindow(Int_t width = 1600, Int_t height = 1200);
      void                 InitializePlots();

      void                 DrawPlots();

      Int_t                UpdatePlots(const std::map<Int_t, std::string>& ifilenames);

      inline void          Run() { fApplication->Run(true); }

      inline Bool_t        IsClosed() const { return fCanvas->GetCanvasImp(); }
      inline Bool_t        IsTerminated() const { return fIsTerminated; }
      
      inline void          Terminate() {
        fIsTerminated = true;
        gSystem->ExitLoop();
      }

    private:
      void                 ClearLastSpill();
      void                 FillCoincidence(const TdcData& tdcData);
      void                 FillCoincidences(const std::vector<TdcData>& tdcData);
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

      for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
        for (std::size_t index = 0; index < CoinOffset::N; ++index) {
          fContains[extCh][index] = false;
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
      } else {
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

    void MonitorWindow::InitializePlots() {
      std::cout << "Initialize plots" << std::endl;
      const std::string tdcName    = fProvider.GetName();
      const Double_t    timePerTdc = fProvider.GetTimePerTdc();

      const Double_t  xminInSpill =    0; // [msec]
      const Double_t  xmaxInSpill = 3000; // [msec]
      const Int_t    xbinsInSpill =  100;

      const Double_t  xminInSync  = (Int_t)(- 5400 * nsec * 0.2 / timePerTdc) - 0.5; // [count]
      const Int_t    xbinsInSync  = (Int_t)(+ 5400 * nsec * 1.4 / timePerTdc) / 10; // [count/10]
      const Double_t  xmaxInSync  = xminInSync + xbinsInSync * 10; // [count]

      const Double_t meanMrSync   = (5220 * nsec / timePerTdc); // [count]
      const Int_t    xbinsInt     = (Int_t)(375 * nsec / timePerTdc); // [count]
      const Double_t  xminInt     = (Int_t)meanMrSync - xbinsInt / 2 - 0.5; // [count]
      const Double_t  xmaxInt     = xminInt + xbinsInt; // [count]

      const Double_t  xminInDiff  = (Int_t)(-100 * nsec / timePerTdc) - 0.5; // [count]
      const Int_t    xbinsInDiff  = (Int_t)( 300 * nsec / timePerTdc); // [count]
      const Double_t  xmaxInDiff  = xminInDiff + xbinsInDiff; // [count]

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
      lExtBorderLine = ExtinctionDetector::CreateBorderLine(kBlack, kSolid, 1);

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
      lHodBorderLine = Hodoscope::CreateBorderLine(kBlack, kSolid, 1);

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
          hExtTdcOffset[ch] = new TH2D(Form("hhExtTdcOffset_%03ld", ch),
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

      if (fCanvas->cd(++padnumber)) {
        hBhTdcInSpill[0]->Draw();
        hBhTdcInSpill[0]->GetXaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[0]->GetYaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[0]->GetXaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[0]->GetYaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[0]->GetXaxis()->SetTitleOffset(0.6);
        hBhTdcInSpill[0]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hBhTdcInSpill[1]->Draw();
        hBhTdcInSpill[1]->GetXaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[1]->GetYaxis()->SetLabelSize(0.07);
        hBhTdcInSpill[1]->GetXaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[1]->GetYaxis()->SetTitleSize(0.07);
        hBhTdcInSpill[1]->GetXaxis()->SetTitleOffset(0.6);
        hBhTdcInSpill[1]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hHodTdcInSpill_Any->Draw();
        hHodTdcInSpill_Any->GetXaxis()->SetLabelSize(0.07);
        hHodTdcInSpill_Any->GetYaxis()->SetLabelSize(0.07);
        hHodTdcInSpill_Any->GetXaxis()->SetTitleSize(0.07);
        hHodTdcInSpill_Any->GetYaxis()->SetTitleSize(0.07);
        hHodTdcInSpill_Any->GetXaxis()->SetTitleOffset(0.6);
        hHodTdcInSpill_Any->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hExtTdcInSpill_Any->Draw();
        hExtTdcInSpill_Any->GetXaxis()->SetLabelSize(0.07);
        hExtTdcInSpill_Any->GetYaxis()->SetLabelSize(0.07);
        hExtTdcInSpill_Any->GetXaxis()->SetTitleSize(0.07);
        hExtTdcInSpill_Any->GetYaxis()->SetTitleSize(0.07);
        hExtTdcInSpill_Any->GetXaxis()->SetTitleOffset(0.6);
        hExtTdcInSpill_Any->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hTcTdcInSpill[0]->Draw();
        hTcTdcInSpill[0]->GetXaxis()->SetLabelSize(0.07);
        hTcTdcInSpill[0]->GetYaxis()->SetLabelSize(0.07);
        hTcTdcInSpill[0]->GetXaxis()->SetTitleSize(0.07);
        hTcTdcInSpill[0]->GetYaxis()->SetTitleSize(0.07);
        hTcTdcInSpill[0]->GetXaxis()->SetTitleOffset(0.6);
        hTcTdcInSpill[0]->GetYaxis()->SetTitleOffset(0.5);
      }

      if (fCanvas->cd(++padnumber)) {
        hTcTdcInSpill[1]->Draw();
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

      } else if (fMonitorMode == MonitorMode::Channel) {
        if (fCanvas->cd(++padnumber)) {
          hMountain->SetTitle(Form("%s, Mountain Plot @ ch%d", fProvider.GetName().data(), fMonitorChannel));
          hMountain->Draw("col");
          hMountain->GetXaxis()->SetLabelSize(0.05);
          hMountain->GetYaxis()->SetLabelSize(0.05);
          hMountain->GetXaxis()->SetTitleSize(0.05);
          hMountain->GetYaxis()->SetTitleSize(0.05);
          hMountain->GetXaxis()->SetTitleOffset(0.6);
          hMountain->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          hTdcInSync->SetTitle(Form("%s, TDC in Sync @ ch%d", fProvider.GetName().data(), fMonitorChannel));
          hTdcInSync->Draw();
          hTdcInSync->GetXaxis()->SetLabelSize(0.05);
          hTdcInSync->GetYaxis()->SetLabelSize(0.05);
          hTdcInSync->GetXaxis()->SetTitleSize(0.05);
          hTdcInSync->GetYaxis()->SetTitleSize(0.05);
          hTdcInSync->GetXaxis()->SetTitleOffset(0.6);
          hTdcInSync->GetYaxis()->SetTitleOffset(0.5);
        }

        if (fCanvas->cd(++padnumber)) {
          hTdcInSpill->SetTitle(Form("%s, TDC in Spill @ ch%d", fProvider.GetName().data(), fMonitorChannel));
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
      fCoinCount = 0;

      fTdcBuffer.clear();

      fLastExtData.clear();
      fLastHodData.clear();
      fLastTcData .clear();
      fLastBhData .clear();
      for (auto&& pair: fLastMrSyncData) {
        pair.second.Clear();
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
      if (fMonitorMode == MonitorMode::Channel) {
        hMountain       ->Reset();
        hTdcInSync      ->Reset();
        hTdcInSpill     ->Reset();
      } else if (fMonitorMode == MonitorMode::Offset) {
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
          if (fContains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : fLastBhData) {
              Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
              for (Double_t dint = fMrSyncInterval[         board] *  extData.TimePerTdc, thre = -0.6 * fMrSyncIntervalAverage; dt < thre; dt += dint);
              for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc, thre = +0.6 * fMrSyncIntervalAverage; dt > thre; dt -= dint);
              const Double_t mean = fCoinDiffs[extCh][i];
              if (std::abs(dt - mean) < kCoinTimeWidth) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
          }
        }

        {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          if (fContains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : fLastHodData) {
              Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
              for (Double_t dint = fMrSyncInterval[         board] *  extData.TimePerTdc, thre = -0.6 * fMrSyncIntervalAverage; dt < thre; dt += dint);
              for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc, thre = +0.6 * fMrSyncIntervalAverage; dt > thre; dt -= dint);
              const Double_t mean = fCoinDiffs[extCh][i];
              if (std::abs(dt - mean) < kCoinTimeWidth) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
          }
        }

        for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
          const std::size_t i = tcCh + CoinOffset::TC;
          if (fContains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : fLastTcData) {
              // const Double_t dt    = lastData.Time - time;
              Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
              for (Double_t dint = fMrSyncInterval[         board] *  extData.TimePerTdc, thre = -0.6 * fMrSyncIntervalAverage; dt < thre; dt += dint);
              for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc, thre = +0.6 * fMrSyncIntervalAverage; dt > thre; dt -= dint);
              const Double_t mean  = fCoinDiffs[extCh][i];
              if (std::abs(dt - mean) < kCoinTimeWidth) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
          }
        }

        if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
          ++fCoinCount;
          auto          lastData = fLastMrSyncData[board];
          const Long64_t     tdc = extData.Tdc;
          const Long64_t syncTdc = lastData.Tdc;
          const Long64_t    dtdc = tdc - syncTdc;
          if        (dtdc > +1.05 * fMrSyncInterval[board]) {
            const Long64_t dtdc2 = std::fmod(dtdc2, fMrSyncInterval[board]);
            hExtTdcInSync_Any->Fill(dtdc2);
            hExtMountain_Any ->Fill(dtdc2, time / msec);
          } else if (dtdc < -0.05 * fMrSyncInterval[board]) {
            const Long64_t dtdc2 = std::fmod(dtdc2, fMrSyncInterval[board]) + fMrSyncInterval[board];
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
          if (fContains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : fLastBhData) {
              const Double_t dt   = lastData.Time - time;
              const Double_t mean = fCoinDiffs[extCh][i];
              if (std::abs(dt - mean) < kCoinTimeWidth) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
          }
        }

        {
          const std::size_t hodCh = 0;
          const std::size_t i = hodCh + CoinOffset::Hod;
          if (fContains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : fLastHodData) {
              const Double_t dt   = lastData.Time - time;
              const Double_t mean = fCoinDiffs[extCh][i];
              if (std::abs(dt - mean) < kCoinTimeWidth) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
          }
        }

        for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
          const std::size_t i = tcCh + CoinOffset::TC;
          if (fContains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : fLastTcData) {
              const Double_t dt   = lastData.Time - time;
              const Double_t mean = fCoinDiffs[extCh][i];
              if (std::abs(dt - mean) < kCoinTimeWidth) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
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

    std::vector<TdcData> MonitorWindow::CollectCoinExtData(const TdcData& tdcData, std::size_t i) {
      const Double_t time = tdcData.Time;

      std::vector<TdcData> coinExtData;
      if (fCyclicCoincidence) {
        for (auto&& lastData : fLastExtData) {
          const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
          if (fContains[extCh][i]) {
            Double_t dt = (time - fLastMrSyncData[tdcData.Board].Time) - (lastData.Time - fLastMrSyncData[lastData.Board].Time);
            for (Double_t dint = fMrSyncInterval[lastData.Board] * lastData.TimePerTdc, thre = -0.6 * fMrSyncIntervalAverage; dt < thre; dt += dint);
            for (Double_t dint = fMrSyncInterval[ tdcData.Board] *  tdcData.TimePerTdc, thre = +0.6 * fMrSyncIntervalAverage; dt > thre; dt -= dint);
            const Double_t mean  = fCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < kCoinTimeWidth) {
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
            if (std::abs(dt - mean) < kCoinTimeWidth) {
              coinExtData.push_back(lastData);
            }
          }
        }
      }
      return coinExtData;
    }

    std::size_t MonitorWindow::RemoveOldTdc(std::vector<TdcData>* lastData, Double_t time) {
      for (std::size_t i = 0, n = lastData->size(); i < n; ++i) {
        if (std::abs(time - lastData->at(0).Time) > kLastThreshold) {
          lastData->erase(lastData->begin());
        } else {
          break;
        }
      }
      return lastData->size();
    }

    Int_t MonitorWindow::UpdatePlots(const std::map<int, std::string>& ifilenames) {
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
      for (auto&& pair : ifilenames) {
        const Int_t       board     = pair.first;
        const std::string ifilename = pair.second;
        std::cout << " - " << ifilename << std::endl;
        ifiles[board].open(ifilename, std::ios::binary);
        if (!ifiles[board]) {
          std::cerr << "[error] input file is not opened, " << ifilename << std::endl;
          return 1;
        }
      }

      {
        Int_t targetBoard = 0;
        std::map<Int_t, Bool_t>    gateStarted;
        std::map<Int_t, Bool_t>    gateEnded;
        std::map<Int_t, Bool_t>    fileEnded;
        std::map<Int_t, ULong64_t> lastTdcTags;
        for (auto&& pair : ifilenames) {
          targetBoard = pair.first;
          gateStarted[pair.first] = false;
          gateEnded  [pair.first] = false;
          fileEnded  [pair.first] = false;
          lastTdcTags[pair.first] = 0;
        }

        for (std::size_t count = 0UL; !fIsTerminated;) {
          Decoder&   decoder    = decoders   [targetBoard];
          ULong64_t& lastTdcTag = lastTdcTags[targetBoard];
          for (std::size_t ibuf = 0; ibuf < kBufferSize; ++ibuf, ++count) {
            if (count % 1000000UL == 0) {
              std::cout << ">> " << count << std::endl;
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
                const std::vector<TdcData> tdcData = decoder.Data.GetTdcData(targetBoard);
                for (auto&& data : tdcData) {
                  fTdcBuffer[data.GetTdcTag()] = data;
                }
                lastTdcTag = tdcData.back().GetTdcTag();

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

            if (fMonitorMode == MonitorMode::Channel && fMonitorChannel == globalChannel) {
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
                  hExtTdcOffset[ch]->Fill(lastCh + CoinOffsetX::Hod, lastData.Tdc - tdc);
                }
                for (auto&& lastData : fLastTcData) {
                  auto lastCh = TimingCounter::GetChannel(lastData.Channel);
                  hExtTdcOffset[ch]->Fill(lastCh + CoinOffsetX::TC, lastData.Tdc - tdc);
                }
                for (auto&& lastData : fLastBhData) {
                  auto lastCh = BeamlineHodoscope::GetChannel(lastData.Channel);
                  hExtTdcOffset[ch]->Fill(lastCh + CoinOffsetX::BH, lastData.Tdc - tdc);
                }
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                FillCoincidence(data);
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
                  hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::Hod, tdc - lastData.Tdc);
                }
                fLastHodData.push_back(data);
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                fLastHodData.push_back(data);
                FillCoincidences(CollectCoinExtData(data, CoinOffset::Hod));
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
                  hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::TC, tdc - lastData.Tdc);
                }
                fLastTcData.push_back(data);
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                fLastTcData.push_back(data);
                FillCoincidences(CollectCoinExtData(data, ch + CoinOffset::TC));
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
                  hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::BH, tdc - lastData.Tdc);
                }
                fLastBhData.push_back(data);
              } else if (fMonitorMode == MonitorMode::Coincidence) {
                fLastBhData.push_back(data);
                FillCoincidences(CollectCoinExtData(data, ch + CoinOffset::BH));
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
            const Int_t np = fSpillCount % kSpillLimit;
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

              std::cout << "_____ MrSync _____" << std::endl;
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
                  std::cout << ch << "\t" << xmean << std::endl;
                }
              }
              std::cout << "^^^^^^^^^^^^^^^^^^" << std::endl;

              std::cout << "_____ Offset _____" << std::endl;
              for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
                if (!hExtTdcOffset[ch]->GetEntries()) {
                  continue;
                }

                const Int_t ybins = hExtTdcOffset[ch]->GetNbinsY();

                // Beamline Hodoscope 1, 2
                for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
                  const std::size_t i    = bhCh + CoinOffsetX::BH;
                  const Int_t       xbin = bhCh + CoinOffsetX::BH + 1;
                  Double_t maxsum = 0.0, maxtdc = 0.0;
                  for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                    const Double_t sum = hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    if (sum > maxsum) {
                      maxtdc = hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      maxsum = sum;
                    }
                  }
                  if (maxsum) {
                    std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                  }
                }

                // Hodoscope
                {
                  const std::size_t i     = CoinOffsetX::Hod;
                  const Int_t       xbin1 = CoinOffsetX::Hod + 1;
                  const Int_t       xbin2 = CoinOffsetX::Hod + Hodoscope::NofChannels;
                  Double_t maxsum = 0.0, maxtdc = 0.0;
                  for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                    Double_t sum = 0.0;
                    for (Int_t xbin = xbin1; xbin <= xbin2; ++xbin) {
                      sum += hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    }
                    if (sum > maxsum) {
                      maxtdc = hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      maxsum = sum;
                    }
                  }
                  if (maxsum) {
                    std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                  }
                }

                // Timing Counter 1, 2
                for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
                  const std::size_t i    = tcCh + CoinOffsetX::Hod + 1;
                  const Int_t       xbin = tcCh + CoinOffsetX::TC + 1;
                  Double_t maxsum = 0.0, maxtdc = 0.0;
                  for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                    const Double_t sum = hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                    if (sum > maxsum) {
                      maxtdc = hExtTdcOffset[ch]->GetYaxis()->GetBinCenter(ybin);
                      maxsum = sum;
                    }
                  }
                  if (maxsum) {
                    std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                  }
                }

              }
              std::cout << "^^^^^^^^^^^^^^^^^^" << std::endl;
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
