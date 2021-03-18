#ifndef Extinction_HistGenerator_hh
#define Extinction_HistGenerator_hh

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

#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"

#include "Math.hh"
#include "Linq.hh"
#include "String.hh"
#include "ObjectHelper.hh"

namespace Extinction {

  namespace Analyzer {

    class HistGenerator {
    public:
      static const std::size_t kNofBunches = 4;

      struct PlotsProfiles {
        struct {
          Double_t NbinsX, Xmin, Xmax;
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
      using CoinDiffs  = std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Double_t>>;
      using Contains   = std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Bool_t>>;

      struct CoinOffset {
        enum {
          BH  = 0,
          Hod = 0 + BeamlineHodoscope::NofChannels,
          TC  = 0 + BeamlineHodoscope::NofChannels + 1,
          N   = 0 + BeamlineHodoscope::NofChannels + 1 + TimingCounter::NofChannels,
        };
      };

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

      struct CoinInfo {
        std::size_t Index;
        Double_t    Mean;     // gotten using GetMaximumBin
        Double_t    FitMean;  // gotten using Fit
        Double_t    FitSigma; // gotten using Fit

        void Write(std::ofstream& file) {
          file << Index    << "    "
               << Mean     << "    "
               << FitMean  << "    "
               << FitSigma << std::endl;
        }
        std::basic_istream<char>& Read(std::ifstream& file) {
          return file >> Index
                      >> Mean
                      >> FitMean
                      >> FitSigma;
        }
      };

    private:
      ITdcDataProvider*            fProvider            = nullptr;

      TH2*                         hHodHitMap           = nullptr;
      TList*                       lHodBorderLine       = nullptr;
      TH1*                         hHodEntryByCh        = nullptr;
      TH2*                         hExtHitMap           = nullptr;
      TList*                       lExtBorderLine       = nullptr;
      TH1*                         hExtEntryByCh        = nullptr;
      TH1*                         hExtEntryByChBottom  = nullptr;
      TH1*                         hExtEntryByChCenter1 = nullptr;
      TH1*                         hExtEntryByChCenter2 = nullptr;
      TH1*                         hExtEntryByChTop     = nullptr;
      TH1**                        hBhTdcInSpill        = nullptr;
      TH1**                        hHodTdcInSpill       = nullptr;
      TH1*                         hHodTdcInSpill_Any   = nullptr;
      TH1**                        hExtTdcInSpill       = nullptr;
      TH1*                         hExtTdcInSpill_Any   = nullptr;
      TH1**                        hTcTdcInSpill        = nullptr;
      TH1**                        hMrSyncTdcInSpill    = nullptr;
      TH1**                        hEvmTdcInSpill       = nullptr;
      TH1**                        hBhTdcInSync         = nullptr;
      TH1**                        hHodTdcInSync        = nullptr;
      TH1*                         hHodTdcInSync_Any    = nullptr;
      TH1**                        hExtTdcInSync        = nullptr;
      TH1*                         hExtTdcInSync_Any    = nullptr;
      TH1**                        hTcTdcInSync         = nullptr;
      TH2**                        hExtMountain         = nullptr;
      TH2*                         hExtMountain_Any     = nullptr;
      TH1*                         hCoinTdcInSync       = nullptr;
      TH2*                         hCoinMountain        = nullptr;
      TGraphErrors*                gHitInSpill          = nullptr;
      TH1**                        hMrSyncInterval      = nullptr;
      TH2**                        hExtTdcOffset        = nullptr;
      TFile*                       fSpillFile           = nullptr;
      TTree*                       fSpillTree           = nullptr;

      Long64_t                     fSpillCount          = 0;
      Long64_t                     fCoinCount           = 0;
      Int_t                        fYear                = 0;
      Int_t                        fMonth               = 0;
      Int_t                        fDay                 = 0;
      Int_t                        fHour                = 0;
      Int_t                        fMinute              = 0;
      Int_t                        fSecond              = 0;
      Int_t                        fEMCount             = -1;
      Long64_t                     fEntries[GlobalChannel::NofChannels] = { 0 };
      Long64_t                     fLeakage             = 0;
      Long64_t                     fInBunch             = 0;
      Double_t                     fExtinction          = 1;
      Double_t                     fExtinctionErr       = 1;
      Double_t                     fBunchCenters[kNofBunches] = { 0 };
      Double_t                     fBunchWidths [kNofBunches] = { 0 };
      Double_t                     fMrSyncInterval[MrSync::NofChannels] = { 0 };
      CoinDiffs                    fCoinDiffs;

      CoinDiffs                    fStdCoinDiffs;
      Contains                     fContains;
      std::map<Int_t, Double_t>    fTimePerTdc;
      std::map<Int_t, Double_t>    fStdMrSyncInterval;
      Double_t                     fStdMrSyncIntervalAverage;
      Double_t                     fStdBunchCenters[kNofBunches] = { 0 };
      Double_t                     fStdBunchWidths [kNofBunches] = { 0 };

      Bool_t                       fCyclicCoincidence   = true;
      Double_t                     fHistoryWidth        = 600.0 * nsec;
      Double_t                     fCoinTimeWidth       =  10.0 * nsec;
      std::size_t                  fBufferSize          = 5000;
      std::size_t                  fBufferMargin        =  100;

      std::map<ULong64_t, TdcData> fTdcBuffer;
      std::vector<TdcData>         fLastBhData;
      std::vector<TdcData>         fLastHodData;
      std::vector<TdcData>         fLastExtData;
      std::vector<TdcData>         fLastTcData;
      std::map<Int_t, TdcData>     fLastMrSyncData;
      std::vector<TdcData>         fEventMatchData;

      TF1*                         fGauss               = nullptr;

    public:
      HistGenerator(ITdcDataProvider* provider);
      ~HistGenerator();

      Int_t                LoadOffset(const std::string& ffilename);
      inline void          SetTimePerTdc(const std::map<Int_t, Double_t>& map) {
        fTimePerTdc = map;
      }
      inline void          SetMrSyncInterval(const std::map<Int_t, Double_t>& map) {
        fStdMrSyncInterval        = map;
        fStdMrSyncIntervalAverage = Tron::Linq::From(map)
          .Select([](const std::pair<Int_t, Double_t>& pair) { return pair.second; })
          .Average();
      }
      inline void          SetBunchCenters(const Double_t bunchCenters[kNofBunches]) {
        std::memcpy(fStdBunchCenters, bunchCenters, sizeof(fStdBunchCenters));
      }
      inline void          SetBunchWidths(const Double_t bunchWidths[kNofBunches]) {
        std::memcpy(fStdBunchWidths, bunchWidths, sizeof(fStdBunchWidths));
      }

      inline void          SetCyclicCoincidence(Bool_t flag) { fCyclicCoincidence = flag; }
      inline Bool_t        IsCyclicCoincidence() const { return fCyclicCoincidence; }

      inline void          SetHistoryWidth(Double_t width) { fHistoryWidth = width; }
      inline Double_t      GetHistoryWidth() const { return fHistoryWidth; }

      inline void          SetCoinTimeWidth(Double_t width) { fCoinTimeWidth = width; }
      inline Double_t      GetCoinTimeWidth() const { return fCoinTimeWidth; }

      inline void          SetReadBufferSize(std::size_t size)  { fBufferSize = size; }
      inline std::size_t   GetReadBufferSize() const { return fBufferSize; }

      inline void          SetReadBufferMargin(std::size_t size)  { fBufferMargin = size; }
      inline std::size_t   GetReadBufferMargin() const { return fBufferMargin; }

      void                 ReadPlots(const std::string& ifilename);
      void                 InitializePlots(const PlotsProfiles& profile);
      void                 InitializeSpillSummary(const std::string& filename, const std::string& treename = "spilltree");

      void                 DrawPlots(const std::string& ofilename);
      void                 WritePlots(const std::string& ofilename);
      void                 WriteSpillSummary();
      void                 WriteBunchProfile(const std::string& ofilename);
      void                 WriteMrSyncInterval(const std::string& ofilename);
      void                 WriteCoinDiffs(const std::string& ofilename);

      Int_t                GeneratePlots(std::map<Int_t, ITdcDataProvider*> providers,
                                         const std::map<Int_t, std::string>& ifilenames,
                                         const std::string& treename,
                                         const std::string& ofilename = "",
                                         const std::function<TDatime(const std::string&)>& parser = nullptr);

    private:
      void                 ClearLastSpill(Bool_t clearHists);
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

    HistGenerator::HistGenerator(ITdcDataProvider* provider)
      : fProvider(provider) {
      fLastBhData .reserve(100);
      fLastHodData.reserve(100);
      fLastExtData.reserve(100);
      fLastTcData .reserve(100);

      fGauss = new TF1("fGauss", "gaus(0)");
    }

    HistGenerator::~HistGenerator() {
    }

    Int_t HistGenerator::LoadOffset(const std::string& ffilename) {
      std::cout << "Load offset" << std::endl;
      if (Tron::String::GetFormatSpecifiers(ffilename).empty()) {
        std::ifstream ffile(ffilename);
        if (ffile) {
          Int_t ch, index; Double_t mean;
          while (ffile >> ch >> index >> mean) {
            fStdCoinDiffs[ch][index] = mean * fProvider->GetTimePerTdc();
          }
        }
      } else {
        for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
          CoinInfo info;
          std::ifstream ffile(Form(ffilename.data(), extCh));
          if (ffile) {
            while (info.Read(ffile)) {
              fStdCoinDiffs[extCh][info.Index] = info.FitMean * fProvider->GetTimePerTdc();
            }
          }
        }
      }

      for (std::size_t ch = ContainsOffset::BH; ch < (std::size_t)ContainsOffset::Max; ++ch) {
        const std::size_t index = ch - ContainsOffset::BH;
        // std::cout << "++ ch (index) " << ch << " (" << index << ")" << std::endl;
        std::map<std::size_t/*index*/, Double_t> sum, cnt;
        for (auto&& info : fStdCoinDiffs) {
          const std::size_t extCh = info.first;
          if (ExtinctionDetector::Contains(extCh)) {
            for (auto&& subInfo : info.second) {
              const std::size_t index2 = subInfo.first;
              // std::cout << "     " << extCh << "," << index2 << "," << index << "    " << (fCoinDiffs[extCh][index2] - fCoinDiffs[extCh][index]) / nsec << " nsec" << std::endl;
              sum[index2] += fStdCoinDiffs[extCh][index2] - fStdCoinDiffs[extCh][index];
              cnt[index2]++;
            }
          }
        }
        for (auto&& pair : sum) {
          fStdCoinDiffs[ch][pair.first] = sum[pair.first] / cnt[pair.first];
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
      for (auto&& info : fStdCoinDiffs) {
        const std::size_t extCh = info.first;
        for (auto&& subInfo : info.second) {
          const std::size_t index = subInfo.first;
          fContains[extCh][index] = true;
        }
      }

      return fStdCoinDiffs.size();
    }

    void HistGenerator::ReadPlots(const std::string& ifilename) {
      TFile* file = new TFile(ifilename.data(), "READ");
      if (!file->IsOpen()) {
        std::cout << "[error] input file is not opened, " << ifilename << std::endl;
        return;
      }

      {
        hHodHitMap     = dynamic_cast<TH2*>  (file->Get("hExtHitMap"    ));      
        lHodBorderLine = dynamic_cast<TList*>(file->Get("lExtBorderLine"));
      }

      {
        hHodEntryByCh = dynamic_cast<TH1*>(file->Get("hHodEntryByCh"));
      }

      {
        hExtHitMap     = dynamic_cast<TH2*>  (file->Get("hExtHitMap"    ));
        lExtBorderLine = dynamic_cast<TList*>(file->Get("lExtBorderLine"));
      }

      {
        hExtEntryByCh = dynamic_cast<TH1*>(file->Get("hExtEntryByCh"));
      } {
        hExtEntryByChBottom  = dynamic_cast<TH1*>(file->Get("hExtEntryByChBottom" ));
        hExtEntryByChCenter1 = dynamic_cast<TH1*>(file->Get("hExtEntryByChCenter1"));
        hExtEntryByChCenter2 = dynamic_cast<TH1*>(file->Get("hExtEntryByChCenter2"));
        hExtEntryByChTop     = dynamic_cast<TH1*>(file->Get("hExtEntryByChTop"    ));
      }

      hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSpill[ch] = dynamic_cast<TH1*>(file->Get(Form("hBhTdcInSpill_%03lu", ch)));
      }

      hHodTdcInSpill = new TH1*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSpill[ch] = dynamic_cast<TH1*>(file->Get(Form("hHodTdcInSpill_%03lu", ch)));
      } {
        hHodTdcInSpill_Any = dynamic_cast<TH1*>(file->Get("hHodTdcInSpill_Any"));
      }

      hExtTdcInSpill = new TH1*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSpill[ch] = dynamic_cast<TH1*>(file->Get(Form("hExtTdcInSpill_%03lu", ch)));
      } {
        hExtTdcInSpill_Any = dynamic_cast<TH1*>(file->Get("hExtTdcInSpill_Any"));
      }

      hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSpill[ch] = dynamic_cast<TH1*>(file->Get(Form("hTcTdcInSpill_%03lu", ch)));
      }

      hMrSyncTdcInSpill = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncTdcInSpill[ch] = dynamic_cast<TH1*>(file->Get(Form("hMrSyncInSpill_%03lu", ch)));
      }

      hEvmTdcInSpill = new TH1*[EventMatch::NofChannels];
      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        hEvmTdcInSpill[ch] = dynamic_cast<TH1*>(file->Get(Form("hEvmTdcInSpill_%03lu", ch)));
      }

      hBhTdcInSync = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSync[ch] = dynamic_cast<TH1*>(file->Get(Form("hBhTdcInSync_%03lu", ch)));
      }

      hHodTdcInSync = new TH1*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSync[ch] = dynamic_cast<TH1*>(file->Get(Form("hHodTdcInSync_%03lu", ch)));
      } {
        hHodTdcInSync_Any = dynamic_cast<TH1*>(file->Get("hHodTdcInSync_Any"));
      }

      hExtTdcInSync = new TH1*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSync[ch] = dynamic_cast<TH1*>(file->Get(Form("hExtTdcInSync_%03lu", ch)));
      } {
        hExtTdcInSync_Any = dynamic_cast<TH1*>(file->Get("hExtTdcInSync_Any"));
      }

      hTcTdcInSync = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSync[ch] = dynamic_cast<TH1*>(file->Get(Form("hTcTdcInSync_%03lu", ch)));
      }

      hExtMountain = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtMountain[ch] = dynamic_cast<TH2*>(file->Get(Form("hExtMountain_%03lu", ch)));
      } {
        hExtMountain_Any = dynamic_cast<TH2*>(file->Get("hExtMountain_Any"));
      }

      {
        hCoinTdcInSync = dynamic_cast<TH1*>(file->Get("hCoinTdcInSync"));
      }

      {
        hCoinMountain = dynamic_cast<TH2*>(file->Get("hCoinMountain"));
      }

      {
        gHitInSpill = dynamic_cast<TGraphErrors*>(file->Get("gHitInSpill"));
      }

      hMrSyncInterval = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval[ch] = dynamic_cast<TH1*>(file->Get(Form("hMrSyncInterval_%03lu", ch)));
      }

      hExtTdcOffset = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcOffset[ch] = dynamic_cast<TH2*>(file->Get(Form("hExtTdcOffset_%03lu", ch)));
      }

      file->Close();
    }

    void HistGenerator::InitializePlots(const PlotsProfiles& profile) {
      std::cout << "Initialize plots" << std::endl;
      const std::string tdcName    = fProvider->GetName();
      const Double_t    timePerTdc = fProvider->GetTimePerTdc();

      const Double_t  xminInSpill = profile.TimeInSpill.Xmin / msec;
      const Double_t  xmaxInSpill = profile.TimeInSpill.Xmax / msec;
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

      // Beamline Hodoscope TDC in spill
      hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSpill[ch] = new TH1D(Form("hBhTdcInSpill_%03lu", ch),
                                     Form("%s, BH%ld TDC in Spill;"
                                          "Time [ms]", tdcName.data(), ch + 1),
                                     xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      // Hodoscope TDC in spill
      hHodTdcInSpill = new TH1*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSpill[ch] = new TH1D(Form("hHodTdcInSpill_%03lu", ch),
                                      Form("%s, Hodoscope Tdc in Spill @ %lu;"
                                           "Time [ms]", tdcName.data(), ch),
                                      xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      hHodTdcInSpill_Any = new TH1D("hHodTdcInSpill_Any",
                                    Form("%s, Hodoscope TDC in Spill;"
                                         "Time [ms]", tdcName.data()),
                                    xbinsInSpill, xminInSpill, xmaxInSpill);

      // Extinction Detector TDC in spill
      hExtTdcInSpill = new TH1*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSpill[ch] = new TH1D(Form("hExtTdcInSpill_%03lu", ch),
                                      Form("%s, Extinction Detector Tdc in Spill @ %lu;"
                                           "Time [ms]", tdcName.data(), ch),
                                      xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      hExtTdcInSpill_Any = new TH1D("hExtTdcInSpill_Any",
                                    Form("%s, Extinction Detector TDC in Spill;"
                                         "Time [ms]", tdcName.data()),
                                    xbinsInSpill, xminInSpill, xmaxInSpill);

      // Timing Counter TDC in spill
      hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSpill[ch] = new TH1D(Form("hTcTdcInSpill_%03lu", ch),
                                     Form("%s, Timing Counter %ld TDC in Spill;"
                                          "Time [ms]", tdcName.data(), ch + 1),
                                     xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      // MR Sync TDC in spill
      hMrSyncTdcInSpill = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncTdcInSpill[ch] = new TH1D(Form("hMrSyncTdcInSpill_%03lu", ch),
                                         Form("%s, MR Sync %ld TDC in Spill;"
                                              "Time [ms]", tdcName.data(), ch + 1),
                                         xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      // Event Match TDC in spill
      hEvmTdcInSpill = new TH1*[EventMatch::NofChannels];
      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        hEvmTdcInSpill[ch] = new TH1D(Form("hEvmTdcInSpill_%03lu", ch),
                                      Form("%s, Event Match %ld TDC in Spill;"
                                           "Time [ms]", tdcName.data(), ch + 1),
                                      xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      // Beamline Hodoscope TDC in sync
      hBhTdcInSync = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSync[ch] = new TH1D(Form("hBhTdcInSync_%03lu", ch),
                                     Form("%s, BH%ld TDC in MR Sync;"
                                          "TDC [count]", tdcName.data(), ch + 1),
                                     xbinsInSync, xminInSync, xmaxInSync);
      }

      // Hodoscope TDC in sync
      hHodTdcInSync = new TH1*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSync[ch] = new TH1D(Form("hHodTdcInSync_%03lu", ch),
                                      Form("%s, Hodoscope Tdc in MR Sync @ %lu;"
                                           "TDC [count]", tdcName.data(), ch),
                                      xbinsInSync, xminInSync, xmaxInSync);
      }

      hHodTdcInSync_Any = new TH1D("hHodTdcInSync_Any",
                                    Form("%s, Hodoscope TDC in MR Sync;"
                                         "TDC [count]", tdcName.data()),
                                    xbinsInSync, xminInSync, xmaxInSync);

      // Extinction Detector TDC in sync
      hExtTdcInSync = new TH1*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSync[ch] = new TH1D(Form("hExtTdcInSync_%03lu", ch),
                                 Form("%s, Extinction Detector TDC in MR Sync @ %lu;"
                                      "TDC [count]", tdcName.data(), ch),
                                 xbinsInSync, xminInSync, xmaxInSync);
      }

      hExtTdcInSync_Any = new TH1D("hExtTdcInSync_Any",
                                   Form("%s, Extinction Detector TDC in MR Sync;"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInSync, xminInSync, xmaxInSync);

      // Timing Counter TDC in sync
      hTcTdcInSync = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSync[ch] = new TH1D(Form("hTcTdcInSync_%03lu", ch),
                                     Form("%s, Timing Counter %ld TDC in MR Sync;"
                                          "TDC [count]", tdcName.data(), ch + 1),
                                     xbinsInSync, xminInSync, xmaxInSync);
      }

      // Extinction Detector Mountain Plot
      hExtMountain = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtMountain[ch] = new TH2D(Form("hExtMountain_%03lu", ch),
                                    Form("%s, Extinction Detector Mountain Plot @ %lu;"
                                         "TDC [count];"
                                         "Time [ms]", tdcName.data(), ch),
                                    xbinsInSync / 2, xminInSync, xmaxInSync,
                                    xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hExtMountain[ch]->SetStats(false);
      }

      hExtMountain_Any = new TH2D("hExtMountain_Any",
                                  Form("%s, Extinction Detector Mountain Plot;"
                                       "TDC [count];"
                                       "Time [ms]", tdcName.data()),
                                  xbinsInSync / 2, xminInSync, xmaxInSync,
                                  xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hExtMountain_Any->SetStats(false);

      // Extinction Detector TDC in sync (coincidence)
      hCoinTdcInSync = new TH1D("hCoinTdcInSync",
                                Form("%s, TDC in MR Sync;"
                                     "TDC [count]", tdcName.data()),
                                xbinsInSync, xminInSync, xmaxInSync);

      // Extinction Detector Mountain Plot (coincidence)
      hCoinMountain = new TH2D("hCoinMountain",
                               Form("%s, Mountain Plot;"
                                    "TDC [count];"
                                    "Time [ms]", tdcName.data()),
                               xbinsInSync / 2, xminInSync, xmaxInSync,
                               xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hCoinMountain->SetStats(false);

      // Hit in spill (coincidence)
      gHitInSpill = new TGraphErrors();
      gHitInSpill->SetNameTitle("gHitInSpill",
                                Form("%s, # of Hits in Spill;"
                                     "Spill", tdcName.data()));
      gHitInSpill->SetMarkerStyle(kPlus);
      gHitInSpill->SetMarkerColor(kBlue + 1);

      // Offset monitor hists
      hMrSyncInterval = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval[ch] = new TH1D(Form("hMrSyncInterval_%03lu", ch),
                                       Form("%s, MR Sync TDC Interval @ ch%ld;"
                                            "TDC [count];"
                                            "", tdcName.data(), ch),
                                       xbinsInt, xminInt, xmaxInt);
      }

      hExtTdcOffset = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcOffset[ch] = new TH2D(Form("hhExtTdcOffset_%03lu", ch),
                                     Form("%s, Extinction Detector TDC Offset @ ch%ld;"
                                          "Channel;"
                                          "TDC [count]", tdcName.data(), ch),
                                     CoinOffsetX::N, 0 - 0.5, CoinOffsetX::N - 0.5,
                                     xbinsInDiff, xminInDiff, xmaxInDiff);
        hExtTdcOffset[ch]->SetStats(false);
      }
    }

    void HistGenerator::InitializeSpillSummary(const std::string& filename, const std::string& treename) {
      std::cout << "Initialize spill summary" << std::endl;

      fSpillFile = new TFile(filename.data(), "RECREATE");
      if (!fSpillFile->IsOpen()) {
        std::cout << "[error] spill summary file is not opened, " << filename << std::endl;
        return;
      }

      fSpillTree = new TTree(treename.data(), "Spill summary");
      fSpillTree->Branch("year"          , &fYear          , Form("year"                "/I"                            ));
      fSpillTree->Branch("month"         , &fMonth         , Form("month"               "/I"                            ));
      fSpillTree->Branch("day"           , &fDay           , Form("day"                 "/I"                            ));
      fSpillTree->Branch("hour"          , &fHour          , Form("hour"                "/I"                            ));
      fSpillTree->Branch("minute"        , &fMinute        , Form("minute"              "/I"                            ));
      fSpillTree->Branch("second"        , &fSecond        , Form("second"              "/I"                            ));
      fSpillTree->Branch("emcount"       , &fEMCount       , Form("emcount"             "/I"                            ));      
      fSpillTree->Branch("entries"       ,  fEntries       , Form("entries"        "[%lu]/L", GlobalChannel::NofChannels));
      fSpillTree->Branch("coincount"     , &fCoinCount     , Form("coincount"           "/L"                            ));
      fSpillTree->Branch("leakage"       , &fLeakage       , Form("leakage"             "/L"                            ));
      fSpillTree->Branch("inbunch"       , &fInBunch       , Form("inbunch"             "/L"                            ));
      fSpillTree->Branch("extinction"    , &fExtinction    , Form("extinction"          "/D"                            ));
      fSpillTree->Branch("extinctionerr" , &fExtinctionErr , Form("extinctionerr"       "/D"                            ));
      fSpillTree->Branch("bunchcenter"   ,  fBunchCenters  , Form("bunchcenter"    "[%lu]/D",                kNofBunches));
      fSpillTree->Branch("bunchwidth"    ,  fBunchWidths   , Form("bunchwidth"     "[%lu]/D",                kNofBunches));
      fSpillTree->Branch("mrsyncinterval",  fMrSyncInterval, Form("mrsyncinterval" "[%lu]/D", MrSync       ::NofChannels));
    }

    void HistGenerator::DrawPlots(const std::string& ofilename) {
      std::cout << "Draw plots" << std::endl;
      if (!gPad) {
        TCanvas::MakeDefCanvas();
      }
      gPad->SetLogy(false);
      gPad->Print((ofilename + "[").data());

      {
        hHodHitMap->Draw("col");
        lHodBorderLine->Draw();
        hHodHitMap->SetMinimum(-0.001);
        gPad->Print(ofilename.data());
      }

      gPad->SetLogy(true);
      {
        hHodEntryByCh->Draw();
        hHodEntryByCh->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      {
        hExtHitMap->Draw("col");
        lExtBorderLine->Draw();
        hExtHitMap->SetMinimum(-0.001);
        gPad->Print(ofilename.data());
      }

      gPad->SetLogy(true);
      {
        hExtEntryByCh->Draw();
        hExtEntryByCh->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      } {
        hExtEntryByChBottom ->Draw();
        hExtEntryByChCenter1->Draw("same");
        hExtEntryByChCenter2->Draw("same");
        hExtEntryByChTop    ->Draw("same");
        hExtEntryByChBottom->SetMinimum(0.2);
        hExtEntryByChBottom->SetMaximum(2.0 * hExtHitMap->GetBinContent(hExtHitMap->GetMaximumBin()));
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhTdcInSpill[ch]->GetEntries()) {
          hBhTdcInSpill[ch]->Draw();
          hBhTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        if (hHodTdcInSpill[ch]->GetEntries()) {
          hHodTdcInSpill[ch]->Draw();
          hHodTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      } {
        hHodTdcInSpill_Any->Draw();
        hHodTdcInSpill_Any->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcInSpill[ch]->GetEntries()) {
          hExtTdcInSpill[ch]->Draw();
          hExtTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      } {
        hExtTdcInSpill_Any->Draw();
        hExtTdcInSpill_Any->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcInSpill[ch]->GetEntries()) {
          hTcTdcInSpill[ch]->Draw();
          hTcTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncTdcInSpill[ch]->GetEntries()) {
          hMrSyncTdcInSpill[ch]->Draw();
          hMrSyncTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        if (hEvmTdcInSpill[ch]->GetEntries()) {
          hEvmTdcInSpill[ch]->Draw();
          hEvmTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhTdcInSync[ch]->GetEntries()) {
          hBhTdcInSync[ch]->Draw();
          hBhTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        if (hHodTdcInSync[ch]->GetEntries()) {
          hHodTdcInSync[ch]->Draw();
          hHodTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      } {
        hHodTdcInSync_Any->Draw();
        hHodTdcInSync_Any->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcInSync[ch]->GetEntries()) {
          hExtTdcInSync[ch]->Draw();
          hExtTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      } {
        hExtTdcInSync_Any->Draw();
        hExtTdcInSync_Any->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcInSync[ch]->GetEntries()) {
          hTcTdcInSync[ch]->Draw();
          hTcTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtMountain[ch]->GetEntries()) {
          hExtMountain[ch]->Draw("colz");
          hExtMountain[ch]->SetMinimum(0);
          gPad->Print(ofilename.data());
        }
      } {
        hExtMountain_Any->Draw("colz");
        hExtMountain_Any->SetMinimum(0);
        gPad->Print(ofilename.data());
      }

      gPad->SetLogy(true);
      {
        hCoinTdcInSync->Draw();
        hCoinTdcInSync->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      {
        hCoinMountain->Draw("colz");
        hCoinMountain->SetMinimum(0);
        gPad->Print(ofilename.data());
      }

      gPad->SetLogy(true);
      {
        gHitInSpill->Draw("AP");
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncInterval[ch]->GetEntries()) {
          hMrSyncInterval[ch]->Draw();
          hMrSyncInterval[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcOffset[ch]->GetEntries()) {
          hExtTdcOffset[ch]->Draw("col");
          hExtTdcOffset[ch]->SetMinimum(0);
          // hExtTdcOffset[ch]->SetMinimum(-0.001);
          gPad->Print(ofilename.data());
        }
      }

      gPad->Print((ofilename + "]").data());
    }
    void HistGenerator::WritePlots(const std::string& ofilename) {
      std::cout << "Write plots" << std::endl;
      TFile* file = new TFile(ofilename.data(), "RECREATE");
      if (!file->IsOpen()) {
        std::cout << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      {
        hHodHitMap->Write();
        lHodBorderLine->Write(lHodBorderLine->GetName(), TObject::kSingleKey);
      }

      {
        hHodEntryByCh->Write();
      }

      {
        hExtHitMap->Write();
        lExtBorderLine->Write(lExtBorderLine->GetName(), TObject::kSingleKey);
      }

      {
        hExtEntryByCh->Write();
      } {
        hExtEntryByChBottom ->Write();
        hExtEntryByChCenter1->Write();
        hExtEntryByChCenter2->Write();
        hExtEntryByChTop    ->Write();
      }

      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhTdcInSpill[ch]->GetEntries()) {
          hBhTdcInSpill[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        if (hHodTdcInSpill[ch]->GetEntries()) {
          hHodTdcInSpill[ch]->Write();
        }
      } {
        hHodTdcInSpill_Any->Write();
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcInSpill[ch]->GetEntries()) {
          hExtTdcInSpill[ch]->Write();
        }
      } {
        hExtTdcInSpill_Any->Write();
      }

      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcInSpill[ch]->GetEntries()) {
          hTcTdcInSpill[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncTdcInSpill[ch]->GetEntries()) {
          hMrSyncTdcInSpill[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        if (hEvmTdcInSpill[ch]->GetEntries()) {
          hEvmTdcInSpill[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhTdcInSync[ch]->GetEntries()) {
          hBhTdcInSync[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        if (hHodTdcInSync[ch]->GetEntries()) {
          hHodTdcInSync[ch]->Write();
        }
      } {
        hHodTdcInSync_Any->Write();
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcInSync[ch]->GetEntries()) {
          hExtTdcInSync[ch]->Write();
        }
      } {
        hExtTdcInSync_Any->Write();
      }

      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcInSync[ch]->GetEntries()) {
          hTcTdcInSync[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtMountain[ch]->GetEntries()) {
          hExtMountain[ch]->Write();
        }
      } {
        hExtMountain_Any->Write();
      }

      {
        hCoinTdcInSync->Write();
      }

      {
        hCoinMountain->Write();
      }

      {
        gHitInSpill->Write();
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncInterval[ch]->GetEntries()) {
          hMrSyncInterval[ch]->Write();
        }
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcOffset[ch]->GetEntries()) {
          hExtTdcOffset[ch]->Write();
        }
      }

      file->Close();
    }

    void HistGenerator::WriteSpillSummary() {
      std::cout << "Write spill summary" << std::endl;
      if (!fSpillFile) {
        std::cout << "[warning] spill file is not initialized" << std::endl;
        return;
      }

      if (!fSpillTree) {
        std::cout << "[warning] spill tree is not initialized" << std::endl;
        return;
      }
      
      fSpillFile->cd();
      fSpillTree->Write();
    }

    void HistGenerator::WriteBunchProfile(const std::string &ofilename) {
      std::cout << "Write bunch profile" << std::endl;
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cerr << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
        ofile << bunch << "\t"
              << fBunchCenters[bunch] / nsec << "\t"
              << fBunchWidths [bunch] / nsec << std::endl;
      }

      ofile.close();
    }

    void HistGenerator::WriteMrSyncInterval(const std::string& ofilename) {
      std::cout << "Write MR Sync interval" << std::endl;
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cerr << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        ofile << ch << "\t" <<  fMrSyncInterval[ch] << std::endl;
      }

      ofile.close();
    }

    void HistGenerator::WriteCoinDiffs(const std::string& ofilename) {
      std::cout << "Write offset" << std::endl;
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cerr << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      for (auto&& pair : fCoinDiffs) {
        auto& board = pair.first;
        for (auto&& pair2 : pair.second) {
          auto& index = pair2.first;
          auto& value = pair2.second;
          ofile << board << "\t" << index << "\t" << value << std::endl;
        }
      }

      ofile.close();
    }

    void HistGenerator::ClearLastSpill(Bool_t clearHists) {
      fCoinCount     = 0;
      fEMCount       = -1;
      fLeakage       = 0;
      fExtinction    = 1;
      fExtinctionErr = 1;
      std::memset(fEntries    , 0, sizeof(fEntries    ));
      std::memset(fBunchWidths, 0, sizeof(fBunchWidths));

      fTdcBuffer.clear();

      fLastExtData.clear();
      fLastHodData.clear();
      fLastTcData .clear();
      fLastBhData .clear();
      for (auto&& pair: fLastMrSyncData) {
        pair.second.Clear();
      }
      fEventMatchData.clear();

      if (clearHists) {
        hHodHitMap   ->Reset();
        hHodEntryByCh->Reset();

        hExtHitMap   ->Reset();
        hExtEntryByCh->Reset();

        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhTdcInSpill[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodTdcInSpill[ch]->Reset();
        } {
          hHodTdcInSpill_Any->Reset();
        }

        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcInSpill[ch]->Reset();
        } {
          hExtTdcInSpill_Any->Reset();
        }

        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcTdcInSpill[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncTdcInSpill[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
          hEvmTdcInSpill[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhTdcInSync[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodTdcInSync[ch]->Reset();
        } {
          hHodTdcInSync_Any->Reset();
        }

        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcInSync[ch]->Reset();
        } {
          hExtTdcInSync_Any->Reset();
        }

        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcTdcInSync[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtMountain[ch]->Reset();
        } {
          hExtMountain_Any->Reset();
        }

        {
          hCoinTdcInSync->Reset();
        }

        {
          hCoinMountain->Reset();
        }

        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncInterval[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcOffset[ch]->Reset();
        }
      }
    }

    void HistGenerator::FillCoincidence(const TdcData& extData) {
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[         board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fStdMrSyncInterval[         board] *  extData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage *  extData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fStdCoinDiffs[extCh][i];
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[         board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fStdMrSyncInterval[         board] *  extData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage *  extData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fStdCoinDiffs[extCh][i];
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[         board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fStdMrSyncInterval[         board] *  extData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage *  extData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean  = fStdCoinDiffs[extCh][i];
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
          if        (dtdc > +1.05 * fStdMrSyncInterval[board]) {
            const Long64_t dtdc2 = std::fmod(dtdc, fStdMrSyncInterval[board]);
            hCoinTdcInSync->Fill(dtdc2);
            hCoinMountain ->Fill(dtdc2, time / msec);
          } else if (dtdc < -0.05 * fStdMrSyncInterval[board]) {
            const Long64_t dtdc2 = std::fmod(dtdc, fStdMrSyncInterval[board]) + fStdMrSyncInterval[board];
            hCoinTdcInSync->Fill(dtdc2);
            hCoinMountain ->Fill(dtdc2, time / msec);
          } else {
            hCoinTdcInSync->Fill(dtdc);
            hCoinMountain ->Fill(dtdc, time / msec);
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
            const Double_t mean = fStdCoinDiffs[extCh][i];
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
            const Double_t mean = fStdCoinDiffs[extCh][i];
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
            const Double_t mean = fStdCoinDiffs[extCh][i];
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
          hCoinTdcInSync->Fill(tdc - syncTdc);
          hCoinMountain ->Fill(tdc - syncTdc, time / msec);
        }
      }
    }

    inline void HistGenerator::FillCoincidences(const std::vector<TdcData>& extData) {
      for (auto&& data : extData) {
        FillCoincidence(data);
      }
    }

    void HistGenerator::FillCoincidence2(const TdcData& tdcData) {
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[         board] *  tdcData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fStdMrSyncInterval[         board] *  tdcData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage *  tdcData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fStdCoinDiffs[tdcCh][i];
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[         board] *  tdcData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fStdMrSyncInterval[         board] *  tdcData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage *  tdcData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean = fStdCoinDiffs[tdcCh][i];
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[         board] *  tdcData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            // Double_t dt = (lastData.Time - fLastMrSyncData[lastData.Board].Time) - (time - fLastMrSyncData[board].Time);
            // for (Double_t dint = fStdMrSyncInterval[         board] *  tdcData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage *  tdcData.TimePerTdc; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage * lastData.TimePerTdc; dt > thre; dt -= dint);
            const Double_t mean  = fStdCoinDiffs[tdcCh][i];
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
            for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                 dsync += fStdMrSyncInterval[        board] *  tdcData.TimePerTdc);
            for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                 dsync -= fStdMrSyncInterval[extData.Board] *  extData.TimePerTdc);
            const Double_t dt = dt0 - dsync;
            const Double_t mean  = fStdCoinDiffs[extCh][tdcI];
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
              if        (dtdc > +1.05 * fStdMrSyncInterval[board]) {
                const Long64_t dtdc2 = std::fmod(dtdc, fStdMrSyncInterval[board]);
                tdcSum += dtdc2;
              } else if (dtdc < -0.05 * fStdMrSyncInterval[board]) {
                const Long64_t dtdc2 = std::fmod(dtdc, fStdMrSyncInterval[board]) + fStdMrSyncInterval[board];
                tdcSum += dtdc2;
              } else {
                tdcSum += dtdc;
              }
            }
            hCoinTdcInSync->Fill(tdcSum / hits.size());
            hCoinMountain ->Fill(tdcSum / hits.size(), time / msec);
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
            const Double_t mean  = fStdCoinDiffs[tdcCh][i];
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
            const Double_t mean  = fStdCoinDiffs[tdcCh][i];
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
            const Double_t mean  = fStdCoinDiffs[tdcCh][i];
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
            const Double_t mean  = fStdCoinDiffs[extCh][tdcI];
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
            hCoinTdcInSync->Fill(tdcSum / hits.size());
            hCoinMountain ->Fill(tdcSum / hits.size(), time / msec);
          }
        }

      }
    }

    std::vector<TdcData> HistGenerator::CollectCoinExtData(const TdcData& tdcData, std::size_t i) {
      const Double_t time = tdcData.Time;

      std::vector<TdcData> coinExtData;
      if (fCyclicCoincidence) {
        for (auto&& lastData : fLastExtData) {
          const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
          if (fContains[extCh][i]) {
              const Double_t dt0 = time - lastData.Time;
              Double_t dsync = fLastMrSyncData[tdcData.Board].Time - fLastMrSyncData[lastData.Board].Time;
              for (Double_t thre = -0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync < thre;
                   dsync += fStdMrSyncInterval[ tdcData.Board] *  tdcData.TimePerTdc);
              for (Double_t thre = +0.5 * fStdMrSyncIntervalAverage * fProvider->GetTimePerTdc(); dsync > thre;
                   dsync -= fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc);
              const Double_t dt = dt0 - dsync;
            // Double_t dt = (time - fLastMrSyncData[tdcData.Board].Time) - (lastData.Time - fLastMrSyncData[lastData.Board].Time);
            // for (Double_t dint = fStdMrSyncInterval[lastData.Board] * lastData.TimePerTdc,
            //        thre = -0.6 * fStdMrSyncIntervalAverage; dt < thre; dt += dint);
            // for (Double_t dint = fStdMrSyncInterval[ tdcData.Board] *  tdcData.TimePerTdc,
            //        thre = +0.6 * fStdMrSyncIntervalAverage; dt > thre; dt -= dint);
            const Double_t mean  = fStdCoinDiffs[extCh][i];
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
            const Double_t mean = fStdCoinDiffs[extCh][i];
            if (std::abs(dt - mean) < fCoinTimeWidth) {
              coinExtData.push_back(lastData);
            }
          }
        }
      }
      return coinExtData;
    }

    std::size_t HistGenerator::RemoveOldTdc(std::vector<TdcData>* lastData, Double_t time) {
      for (std::size_t i = 0, n = lastData->size(); i < n; ++i) {
        if (std::abs(time - lastData->at(0).Time) > fHistoryWidth) {
          lastData->erase(lastData->begin());
        } else {
          break;
        }
      }
      return lastData->size();
    }

    Int_t HistGenerator::GeneratePlots(std::map<Int_t, ITdcDataProvider*> providers,
                                       const std::map<int, std::string>& ifilenames,
                                       const std::string& treename,
                                       const std::string& ofilename,
                                       const std::function<TDatime(const std::string&)>& parser) {
      const clock_t startClock = clock();

      std::cout << "Initialize decoder" << std::endl;
      for (auto&& pair : ifilenames) {
        const Int_t board = pair.first;
        if (fTimePerTdc[board] == 0) {
          providers[board]->SetTimePerTdc(fProvider->GetTimePerTdc());
        } else {
          providers[board]->SetTimePerTdc(fTimePerTdc[board]);
        }
      }

      std::cout << "Initialzie history" << std::endl;
      ClearLastSpill(true);

      std::cout << "Open file" << std::endl;
      std::map<Int_t, TFile*> ifiles;
      std::map<Int_t, TTree*> itrees;
      std::map<Int_t, Long64_t> entrieses;
      std::map<Int_t, Long64_t> entries;
      std::map<Int_t, Long64_t> lastSpills;
      std::vector<TDatime>      datimes;
      for (auto&& pair : ifilenames) {
        const Int_t       board     = pair.first;
        const std::string ifilename = pair.second;
        std::cout << " - " << ifilename << std::endl;
        if (parser) {
          datimes.push_back(parser(ifilename));
        }

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

        entrieses [board] = itrees[board]->GetEntries();
        entries   [board] = 0;
        lastSpills[board] = -1;

        providers[board]->SetBranchAddress(itrees[board]);
        std::cout << " + " << entrieses[board] << std::endl;
      }
      std::cout << " = " << Tron::Linq::From(entrieses).Sum([](std::pair<Int_t, Long64_t> p) { return p.second; }) << std::endl;
      const UInt_t averageTime = Tron::Linq::From(datimes)
        .Select([](TDatime datime) -> ULong64_t { return datime.Convert(); })
        .Average();
      TDatime averageDatime;
      averageDatime.Set(averageTime);
      fYear   = averageDatime.GetYear();
      fMonth  = averageDatime.GetMonth();
      fDay    = averageDatime.GetDay();
      fHour   = averageDatime.GetHour();
      fMinute = averageDatime.GetMinute();
      fSecond = averageDatime.GetSecond();

      {
        Int_t targetBoard = 0;
        std::map<Int_t, Bool_t>    spillEnded;
        std::map<Int_t, Bool_t>    fileEnded;
        std::map<Int_t, ULong64_t> lastTdcTags;
        std::map<Int_t, std::vector<TdcData>> firstData;
        for (auto&& pair : ifilenames) {
          targetBoard = pair.first;
          spillEnded [pair.first] = false;
          fileEnded  [pair.first] = false;
          lastTdcTags[pair.first] = 0;
        }

        for (std::size_t count = 0UL;;) {
          {
            ITdcDataProvider* provider   = providers  [targetBoard];
            ULong64_t&        lastTdcTag = lastTdcTags[targetBoard];
            Long64_t&         lastSpill  = lastSpills [targetBoard];
            TTree*            itree      = itrees     [targetBoard];
            Long64_t&         entry      = entries    [targetBoard];
            const Long64_t    ientries   = entrieses  [targetBoard];
            for (std::size_t ibuf = 0; ibuf < fBufferSize; ++ibuf, ++count) {
              if (count % 1000000UL == 0) {
                std::cout << ">> " << count << std::endl;
              }

              if (entry >= ientries) {
                std::cout << "[info] detect file end @ " << targetBoard << std::endl;
                fileEnded[targetBoard] = true;
                break;

              } else if (!itree->GetEntry(entry++)) {
                std::cout << "[info] detect file end @ " << targetBoard << " (TTree::GetEntry)" << std::endl;
                fileEnded[targetBoard] = true;
                break;

              } else if (!provider->IsData()) {
                continue;

              } else if (lastSpill != -1 && lastSpill != provider->GetSpill()) {
                std::cout << "[info] detect spill end @ " << targetBoard << std::endl;
                lastSpill               = provider->GetSpill();
                spillEnded[targetBoard] = true;
                firstData [targetBoard] = provider->GetTdcData();
                break;

              } else {
                // std::cout << "[info] detect data @ " << targetBoard << std::endl;
                if (lastSpill != -1 && fEMCount != provider->GetEMCount()) {
                  std::cout << "[error] conflict EMCount, " << fEMCount << " <--> " << provider->GetEMCount() << std::endl;
                  exit(1);
                }
                lastSpill = provider->GetSpill();
                fEMCount  = provider->GetEMCount();
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
            const Long64_t syncTdc       = fLastMrSyncData[board].Tdc;

            if (BeamlineHodoscope::Contains(globalChannel)) {
              const Int_t ch = BeamlineHodoscope::GetChannel(globalChannel);
              const Long64_t tdc = data.Tdc;

              hBhTdcInSpill[ch]->Fill(time / msec);

              hBhTdcInSync[ch]->Fill(tdc - syncTdc);

              for (auto&& lastData : fLastExtData) {
                auto lastCh   = ExtinctionDetector::GetChannel(lastData.Channel);
                hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::BH, tdc - lastData.Tdc);
              }
              fLastBhData.push_back(data);

              fLastBhData.push_back(data);
              // FillCoincidences(CollectCoinExtData(data, ch + CoinOffset::BH));
              FillCoincidence2(data);

              if (RemoveOldTdc(&fLastBhData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastBhData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (Hodoscope::Contains(globalChannel)) {
              const Int_t ch = Hodoscope::GetChannel(globalChannel);
              const Long64_t tdc = data.Tdc;

              hHodEntryByCh     ->Fill(ch);
              Hodoscope         ::Fill(hHodHitMap, ch);

              hHodTdcInSpill[ch]->Fill(time / msec);
              hHodTdcInSpill_Any->Fill(time / msec);

              hHodTdcInSync[ch]->Fill(tdc - syncTdc);
              hHodTdcInSync_Any->Fill(tdc - syncTdc);

              for (auto&& lastData : fLastExtData) {
                auto lastCh = ExtinctionDetector::GetChannel(lastData.Channel);
                hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::Hod, tdc - lastData.Tdc);
              }
              fLastHodData.push_back(data);

              fLastHodData.push_back(data);
              // FillCoincidences(CollectCoinExtData(data, CoinOffset::Hod));
              FillCoincidence2(data);

              if (RemoveOldTdc(&fLastHodData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastHodData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (ExtinctionDetector::Contains(globalChannel)) {
              const Int_t ch = ExtinctionDetector::GetChannel(globalChannel);
              const Long64_t tdc     = data.Tdc;

              hExtEntryByCh     ->Fill(ch);
              ExtinctionDetector::Fill(hExtHitMap, ch);

              hExtTdcInSpill[ch]->Fill(time / msec);
              hExtTdcInSpill_Any->Fill(time / msec);

              hExtTdcInSync[ch]->Fill(tdc - syncTdc);
              hExtTdcInSync_Any->Fill(tdc - syncTdc);

              hExtMountain[ch]->Fill(tdc - syncTdc, time / msec);
              hExtMountain_Any->Fill(tdc - syncTdc, time / msec);

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

              // FillCoincidence(data);

              if (RemoveOldTdc(&fLastExtData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastExtData reaches " << kHistLimit << std::endl;
                return 1;
              }
              fLastExtData.push_back(data);

            } else if (TimingCounter::Contains(globalChannel)) {
              const Int_t ch = TimingCounter::GetChannel(globalChannel);
              const Long64_t tdc = data.Tdc;

              hTcTdcInSpill[ch]->Fill(time / msec);

              hTcTdcInSync[ch]->Fill(tdc - syncTdc);

              for (auto&& lastData : fLastExtData) {
                auto lastCh = ExtinctionDetector::GetChannel(lastData.Channel);
                hExtTdcOffset[lastCh]->Fill(ch + CoinOffsetX::TC, tdc - lastData.Tdc);
              }
              fLastTcData.push_back(data);

              fLastTcData.push_back(data);
              // FillCoincidences(CollectCoinExtData(data, ch + CoinOffset::TC));
              FillCoincidence2(data);

              if (RemoveOldTdc(&fLastTcData, time) > kHistLimit) {
                std::cerr << "[error] size of fLastTcData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (MrSync::Contains(globalChannel)) {
              const Int_t ch = MrSync::GetChannel(globalChannel);

              hMrSyncTdcInSpill[ch]->Fill(time / msec);

              decltype(fLastMrSyncData)::const_iterator itr;
              if ((itr = fLastMrSyncData.find(data.Board)) != fLastMrSyncData.end()) {
                hMrSyncInterval[data.Board]->Fill(data.Tdc - itr->second.Tdc);
              }

              fLastMrSyncData[board] = data;

            } else if (EventMatch::Contains(globalChannel)) {
              const Int_t ch = EventMatch::GetChannel(globalChannel);

              hEvmTdcInSpill[ch]->Fill(time / msec);

              fEventMatchData.push_back(data);
            }

            if (tdcTag == lastTdcTags[board] &&
                !(spillEnded[board] || fileEnded[board])) {
              // std::cout << "[info] detect read last buffer @ " << targetBoard << std::endl;
              targetBoard = board;
              break;
            }

          }

          if (IsAllOfSecondsTrue(spillEnded) || IsAllOfSecondsTrue(fileEnded)) {
            std::cout << "[info] end of spill " << fSpillCount << " (" << fTdcBuffer.size() << ")" << std::endl;
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

            // Set point of hit count
            const Int_t np = fSpillCount;
            gHitInSpill->SetPoint     (np, fSpillCount,     fCoinCount      );
         // gHitInSpill->SetPointError(np, 0.0, TMath::Sqrt(fCoinCount     ));
            ++fSpillCount;

            // Get projections
            for (Int_t xbin = 1, nbinsx = hExtHitMap->GetNbinsX(); xbin <= nbinsx; ++xbin) {
              hExtEntryByChBottom ->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 1));
              hExtEntryByChCenter1->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 2));
              hExtEntryByChCenter2->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 3));
              hExtEntryByChTop    ->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 4));
            }

            // Get entries
            for (std::size_t ch = 0; ch < BeamlineHodoscope ::NofChannels; ++ch) {
              fEntries[ch + BeamlineHodoscope ::GlobalChannelOffset] = hBhTdcInSpill    [ch]->GetEntries();
            }
            for (std::size_t ch = 0; ch < Hodoscope         ::NofChannels; ++ch) {
              fEntries[ch + Hodoscope         ::GlobalChannelOffset] = hHodTdcInSpill   [ch]->GetEntries();
            }
            for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
              fEntries[ch + ExtinctionDetector::GlobalChannelOffset] = hExtTdcInSpill   [ch]->GetEntries();
            }
            for (std::size_t ch = 0; ch < TimingCounter     ::NofChannels; ++ch) {
              fEntries[ch + TimingCounter     ::GlobalChannelOffset] = hTcTdcInSpill    [ch]->GetEntries();
            }
            for (std::size_t ch = 0; ch < MrSync            ::NofChannels; ++ch) {
              fEntries[ch + MrSync            ::GlobalChannelOffset] = hMrSyncTdcInSpill[ch]->GetEntries();
            }
            for (std::size_t ch = 0; ch < EventMatch        ::NofChannels; ++ch) {
              fEntries[ch + EventMatch        ::GlobalChannelOffset] = hEvmTdcInSpill   [ch]->GetEntries();
            }

            // Calc extinction
            fLeakage = 0;
            fInBunch = 0;
            for (Int_t xbin = 1, nbinsx = hCoinTdcInSync->GetNbinsX(); xbin < nbinsx; ++xbin) {
              const Double_t t = hCoinTdcInSync->GetBinCenter(xbin) * fProvider->GetTimePerTdc();
              Bool_t isInBunch = false;
              for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
                if (Tron::Math::Between(t,
                                        fStdBunchCenters[bunch] - fStdBunchWidths[bunch] * 0.5,
                                        fStdBunchCenters[bunch] + fStdBunchWidths[bunch] * 0.5)) {
                  isInBunch = true;
                  break;
                }
              }
              if (isInBunch) {
                fInBunch += hCoinTdcInSync->GetBinContent(xbin);
              } else {
                fLeakage += hCoinTdcInSync->GetBinContent(xbin);
              }
            }
            if (fInBunch) {
              fExtinction = (Double_t)fLeakage / (Double_t)fInBunch;
              if (fLeakage) {
                fExtinctionErr = Tron::Math::SqrtOfSumOfSquared(TMath::Sqrt(fLeakage) / fInBunch, fLeakage / TMath::Power((Double_t)fInBunch, 1.5));
              } else {
                fExtinctionErr = 1.0 / fInBunch;
              }
            } else {
              fExtinction = 1.0;
            }

            // Calc bunch profile
            std::cout << "_____ Bunch Profile _____" << std::endl;
            {
              Int_t xbin = 1, nbinsx = hExtTdcInSync_Any->GetNbinsX();
              const Double_t ymax = hExtTdcInSync_Any->GetBinContent(hExtTdcInSync_Any->GetMaximumBin());
              for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
                Int_t xbin1 = xbin, xbin2 = xbin;
                // std::cout << "... search start edge" << std::endl;
                for (; xbin <= nbinsx && hExtTdcInSync_Any->GetBinContent(xbin) < ymax * 0.10; ++xbin) {
                  // std::cout << xbin << "\t" << hExtTdcInSync_Any->GetBinContent(xbin) << "\t" << ymax * 0.10 << std::endl;
                  xbin1 = xbin;
                }
                // std::cout << "... search end edge" << std::endl;
                for (; xbin <= nbinsx && hExtTdcInSync_Any->GetBinContent(xbin) > ymax * 0.01; ++xbin) {
                  // std::cout << xbin << "\t" << hExtTdcInSync_Any->GetBinContent(xbin) << "\t" << ymax * 0.01 << std::endl;
                  xbin2 = xbin;
                }
                // std::cout << "xbin1 = " << xbin1 << "\txbin2 = " << xbin2 << std::endl;

                if (xbin1 < xbin2) {
                  const Double_t width  = hExtTdcInSync_Any->GetBinCenter(xbin2) - hExtTdcInSync_Any->GetBinCenter(xbin1);
                  const Double_t center = hExtTdcInSync_Any->GetBinCenter((xbin1 + xbin2) / 2);
                  fGauss->SetParameters(ymax, center, width / 2.0);
                  hExtTdcInSync_Any->Fit(fGauss, "+", "", center - 3.0 * width, center + 3.0 * width);

                  fBunchCenters[bunch] = fGauss->GetParameter(1) * fProvider->GetTimePerTdc();
                  fBunchWidths [bunch] = fGauss->GetParameter(2) * fProvider->GetTimePerTdc();
                  std::cout << bunch << "\t"
                            << fBunchCenters[bunch] / nsec << "\t"
                            << fBunchWidths [bunch] / nsec << std::endl;
                } else {
                  std::cout << "[warning] " << (bunch == 0 ? "1st" :
                                                bunch == 1 ? "2nd" :
                                                bunch == 2 ? "3rd" :
                                                Form("%ldth", bunch + 1)) << " bunch was not found" << std::endl;
                  break;
                }
              }
            }

            // Calc MR Sync Interval 
            std::cout << "_____ MR Sync Interval _____" << std::endl;
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
                fMrSyncInterval[ch] = xmean;
              }
            }

            // Calc PMT Offset 
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
                Double_t maxsum = 0.0; Int_t maxbin = 0;
                for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                  const Double_t sum = hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                  if (sum > maxsum) {
                    maxbin = ybin;
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
                  std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                  fCoinDiffs[ch][i] = maxtdc;
                }
              }

              // Hodoscope
              {
                const std::size_t i     = CoinOffsetX::Hod;
                const Int_t       xbin1 = CoinOffsetX::Hod + 1;
                const Int_t       xbin2 = CoinOffsetX::Hod + Hodoscope::NofChannels;
                Double_t maxsum = 0.0; Int_t maxbin = 0;
                for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                  Double_t sum = 0.0;
                  for (Int_t xbin = xbin1; xbin <= xbin2; ++xbin) {
                    sum += hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                  }
                  if (sum > maxsum) {
                    maxbin = ybin;
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
                  std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                  fCoinDiffs[ch][i] = maxtdc;
                }
              }

              // Timing Counter 1, 2
              for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
                const std::size_t i    = tcCh + CoinOffsetX::Hod + 1;
                const Int_t       xbin = tcCh + CoinOffsetX::TC + 1;
                Double_t maxsum = 0.0; Int_t maxbin = 0;
                for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
                  const Double_t sum = hExtTdcOffset[ch]->GetBinContent(xbin, ybin);
                  if (sum > maxsum) {
                    maxbin = ybin;
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
                  std::cout << ch << "\t" << i << "\t" << maxtdc << std::endl;
                  fCoinDiffs[ch][i] = maxtdc;
                }
              }
            }

            // Fill spill summary
            fSpillTree->Fill();

            // Draw plots
            if (ofilename.size()) {
              DrawPlots(ofilename);
            }

            if (IsAllOfSecondsTrue(fileEnded)) {
              break;
            }

            // Clear buffers
            if (ofilename.size()) {
              ClearLastSpill(true);
            } else {
              ClearLastSpill(false);
            }
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
