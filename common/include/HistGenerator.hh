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
#include "TParameter.h"
#include "TEfficiency.h"

#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "Spill.hh"
#include "MargedReader.hh"

#include "Math.hh"
#include "Linq.hh"
#include "String.hh"
#include "ObjectHelper.hh"
#include "ScopeSubstituter.hh"

namespace Extinction {

  namespace Analyzer {

    class HistGenerator {
    public:
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

    public:
      using TdcOffsets_t = MargedReader::TdcOffsets_t;

    private:
      ITdcDataProvider*            fProvider                  = nullptr;

      TH2*                         hHodHitMap                 = nullptr;
      TList*                       lHodBorderLine             = nullptr;
      TH1*                         hHodEntriesByCh            = nullptr;
      TH2*                         hExtHitEvent               = nullptr;
      TH2*                         hExtHitMap                 = nullptr;
      TList*                       lExtBorderLine             = nullptr;
      TH1*                         hExtEntriesByCh            = nullptr;
      TH1*                         hExtEntriesByChBottom      = nullptr;
      TH1*                         hExtEntriesByChCenter1     = nullptr;
      TH1*                         hExtEntriesByChCenter2     = nullptr;
      TH1*                         hExtEntriesByChTop         = nullptr;
      TH2*                         hEntriesInMrSyncByCh       = nullptr;
      TH2*                         hEntriesInMrSyncByDetector = nullptr;
      TGraph*                      gEntriesInMrSyncByDetector = nullptr;
      TH2*                         hExtEntriesInMrSyncInSpill = nullptr;
      TGraph*                      gExtEntriesInMrSyncInSpill = nullptr;
      TH1**                        hBhTdcInSpill              = nullptr;
      TH1**                        hHodTdcInSpill             = nullptr;
      TH1*                         hHodTdcInSpill_Any         = nullptr;
      TH1**                        hExtTdcInSpill             = nullptr;
      TH1*                         hExtTdcInSpill_Any         = nullptr;
      TH1**                        hTcTdcInSpill              = nullptr;
      TH1**                        hMrSyncTdcInSpill          = nullptr;
      TH1**                        hEvmTdcInSpill             = nullptr;
      TH1**                        hVetoTdcInSpill            = nullptr;
      TH1**                        hErrTdcInSpill             = nullptr;
      TH1**                        hBhTdcInSync               = nullptr;
      TH1**                        hHodTdcInSync              = nullptr;
      TH1*                         hHodTdcInSync_Any          = nullptr;
      TH1**                        hExtTdcInSync              = nullptr;
      TH1*                         hExtTdcInSync_Any          = nullptr;
      TH1**                        hTcTdcInSync               = nullptr;
      TH1**                        hVetoTdcInSync             = nullptr;
      TH2**                        hBhMountain                = nullptr;
      TH2**                        hHodMountain               = nullptr;
      TH2*                         hHodMountain_Any           = nullptr;
      TH2**                        hExtMountain               = nullptr;
      TH2*                         hExtMountain_Any           = nullptr;
      TH2**                        hTcMountain                = nullptr;
      TH2**                        hVetoMountain              = nullptr;
      TH2**                        hErrMountain               = nullptr;
      TH1**                        hMrSyncInterval            = nullptr;
      TH2**                        hMrSyncInterval2           = nullptr;
      TH2**                        hBhTdcOffset               = nullptr;
      TH2**                        hTcTdcOffset               = nullptr;
      TH2**                        hExtTdcOffset              = nullptr;

      TFile*                       fSpillFile                 = nullptr;
      TTree*                       fSpillTree                 = nullptr;

      RawSpillData                 fSpillData;
      std::size_t                  fRefExtChannel             = ExtinctionDetector::NofChannels / 2;
      TdcOffsets_t                 fTdcOffsets;

      Double_t                     fHistoryWidth              = 600.0 * nsec;

      std::vector<TdcData>         fLastBhData;
      std::vector<TdcData>         fLastHodData;
      std::vector<TdcData>         fLastExtData;
      std::vector<TdcData>         fLastTcData;
      std::map<Int_t, TdcData>     fLastMrSyncData;

      Double_t                     fBunchCenters [Extinction::kNofBunches] = { };
   // Double_t                     fBunchWidths  [Extinction::kNofBunches] = { };
      Double_t                     fBunchMinEdges[Extinction::kNofBunches] = { };
      Double_t                     fBunchMaxEdges[Extinction::kNofBunches] = { };
      Bool_t                       fOffsetFromBunch           = false;

    public:
      HistGenerator(ITdcDataProvider* provider);
      ~HistGenerator();

      inline void          SetHistoryWidth(Double_t width) {
        std::cout << "SetHistoryWidth ... " << width / nsec << " nsec" << std::endl;
        fHistoryWidth = width;
      }
      inline Double_t      GetHistoryWidth() const { return fHistoryWidth; }

      Int_t                LoadBunchProfile(const std::string& ifilename) {
        std::cout << "LoadBunchProfile ... " << ifilename << std::endl;

        std::ifstream ifile(ifilename);
        if (!ifile) {
          std::cout << "[error] bunch profile is not opened, " << ifilename << std::endl;
          return 1;
        }

        Int_t bunch;
        Double_t center, sigma;
        Long64_t minEdge, maxEdge;
        while (ifile >> bunch >> center >> sigma >> minEdge >> maxEdge) {
          fBunchCenters [bunch] = center;
       // fBunchWidths  [bunch] = sigma;
          fBunchMinEdges[bunch] = minEdge;
          fBunchMaxEdges[bunch] = maxEdge;
        }

        return 0;
      }
      void                 SetOffsetFromBunch(Bool_t flag) { fOffsetFromBunch = flag; }
      Bool_t               IsOffsetFromBunch() const { return fOffsetFromBunch; }

      Int_t                ReadPlots(const std::string& ifilename);
      void                 InitializePlots(const PlotsProfiles& profile);
      void                 InitializeSpillSummary(const std::string& filename, const std::string& treename = "spilltree");

      void                 DrawPlots(const std::string& ofilename, const std::string& ofilename_offset/*, const std::string& ofilename_time*/);
      void                 WritePlots(const std::string& ofilename);
      void                 WriteSpillSummary();
      void                 WriteMrSyncInterval(const std::string& ofilename);
      void                 WriteTdcOffsets(const std::string& ofilename);

      Int_t                GeneratePlots(MargedReader* reader);

      void                 CalcEntries();
      void                 CalcMrSyncInterval();
      void                 CalcTdcOffsets();

    private:
      void                 ClearLastSpill(Bool_t clearHists);
      std::size_t          RemoveOldTdc(std::vector<TdcData>* lastData, const TdcData& tdc);

      Bool_t               IsInBunch(Long64_t dtdc) const {
        for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
          if (!fBunchCenters[bunch]) {
            return false;
          }

          if        (dtdc <  fBunchMinEdges[bunch]) {
            return false;
          } else if (dtdc <= fBunchMaxEdges[bunch]) {
            return true;
          }
        }
        return false;
      }

    };

    HistGenerator::HistGenerator(ITdcDataProvider* provider)
      : fProvider(provider) {
      fLastBhData .reserve(100);
      fLastHodData.reserve(100);
      fLastExtData.reserve(100);
      fLastTcData .reserve(100);
    }

    HistGenerator::~HistGenerator() {
    }

    Int_t HistGenerator::ReadPlots(const std::string& ifilename) {
      TFile* file = new TFile(ifilename.data(), "READ");
      if (!file->IsOpen()) {
        std::cout << "[error] input file is not opened, " << ifilename << std::endl;
        return 1;
      }
      file->ls();

      auto getFromFile =
        [&](const std::string& filename) {
          std::cout << "  - " << filename << std::endl;
          return file->Get(filename.data());
        };

      {
        fSpillData.SetDate((UInt_t)Tron::ObjectHelper::ReadValue<Long64_t>("Time"));
      }

      {
        hHodHitMap     = dynamic_cast<TH2*>  (getFromFile("hHodHitMap"    ));
        lHodBorderLine = dynamic_cast<TList*>(getFromFile("lHodBorderLine"));
        hHodHitMap->SetDirectory(nullptr);
      }

      {
        hHodEntriesByCh = dynamic_cast<TH1*>(getFromFile("hHodEntriesByCh"));
        hHodEntriesByCh->SetDirectory(nullptr);
      }

      {
        hExtHitEvent   = dynamic_cast<TH2*>  (getFromFile("hExtHitEvent"  ));
        hExtHitMap     = dynamic_cast<TH2*>  (getFromFile("hExtHitMap"    ));
        lExtBorderLine = dynamic_cast<TList*>(getFromFile("lExtBorderLine"));
        hExtHitEvent->SetDirectory(nullptr);
        hExtHitMap  ->SetDirectory(nullptr);
      }

      {
        hExtEntriesByCh = dynamic_cast<TH1*>(getFromFile("hExtEntriesByCh"));
        hExtEntriesByCh->SetDirectory(nullptr);
      } {
        hExtEntriesByChBottom  = dynamic_cast<TH1*>(getFromFile("hExtEntriesByChBottom" ));
        hExtEntriesByChCenter1 = dynamic_cast<TH1*>(getFromFile("hExtEntriesByChCenter1"));
        hExtEntriesByChCenter2 = dynamic_cast<TH1*>(getFromFile("hExtEntriesByChCenter2"));
        hExtEntriesByChTop     = dynamic_cast<TH1*>(getFromFile("hExtEntriesByChTop"    ));
        hExtEntriesByChBottom ->SetDirectory(nullptr);
        hExtEntriesByChCenter1->SetDirectory(nullptr);
        hExtEntriesByChCenter2->SetDirectory(nullptr);
        hExtEntriesByChTop    ->SetDirectory(nullptr);
      }

      {
        hEntriesInMrSyncByCh = dynamic_cast<TH2*>(getFromFile("hEntriesInMrSyncByCh"));
        hEntriesInMrSyncByCh->SetDirectory(nullptr);
      }

      {
        hEntriesInMrSyncByDetector = dynamic_cast<TH2*>(getFromFile("hEntriesInMrSyncByDetector"));
        hEntriesInMrSyncByDetector->SetDirectory(nullptr);

        gEntriesInMrSyncByDetector = new TGraph();
        gEntriesInMrSyncByDetector->SetName("gEntriesInMrSyncByDetector");
        gEntriesInMrSyncByDetector->SetMarkerStyle(kPlus);
        gEntriesInMrSyncByDetector->SetMarkerSize(1.0);
      }

      {
        hExtEntriesInMrSyncInSpill = dynamic_cast<TH2*>(getFromFile("hExtEntriesInMrSyncInSpill"));
        hExtEntriesInMrSyncInSpill->SetDirectory(nullptr);

        gExtEntriesInMrSyncInSpill = new TGraph();
        gExtEntriesInMrSyncInSpill->SetName("gExtEntriesInMrSyncInSpill");
        gExtEntriesInMrSyncInSpill->SetMarkerStyle(kPlus);
        gExtEntriesInMrSyncInSpill->SetMarkerSize(0.6);
      }
      
      hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hBhTdcInSpill_%03lu", ch)));
        hBhTdcInSpill[ch]->SetDirectory(nullptr);
      }

      hHodTdcInSpill = new TH1*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hHodTdcInSpill_%03lu", ch)));
        hHodTdcInSpill[ch]->SetDirectory(nullptr);
      } {
        hHodTdcInSpill_Any = dynamic_cast<TH1*>(getFromFile("hHodTdcInSpill_Any"));
        hHodTdcInSpill_Any->SetDirectory(nullptr);
      }

      hExtTdcInSpill = new TH1*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hExtTdcInSpill_%03lu", ch)));
        hExtTdcInSpill[ch]->SetDirectory(nullptr);
      } {
        hExtTdcInSpill_Any = dynamic_cast<TH1*>(getFromFile("hExtTdcInSpill_Any"));
        hExtTdcInSpill_Any->SetDirectory(nullptr);
      }

      hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hTcTdcInSpill_%03lu", ch)));
        hTcTdcInSpill[ch]->SetDirectory(nullptr);
      }

      hVetoTdcInSpill = new TH1*[Veto::NofChannels];
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hVetoTdcInSpill_%03lu", ch)));
        hVetoTdcInSpill[ch]->SetDirectory(nullptr);
      }

      hErrTdcInSpill = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hErrTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hErrTdcInSpill_%03lu", ch)));
        hErrTdcInSpill[ch]->SetDirectory(nullptr);
      }

      hMrSyncTdcInSpill = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hMrSyncTdcInSpill_%03lu", ch)));
        hMrSyncTdcInSpill[ch]->SetDirectory(nullptr);
      }

      hEvmTdcInSpill = new TH1*[EventMatch::NofChannels];
      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        hEvmTdcInSpill[ch] = dynamic_cast<TH1*>(getFromFile(Form("hEvmTdcInSpill_%03lu", ch)));
        hEvmTdcInSpill[ch]->SetDirectory(nullptr);
      }

      hBhTdcInSync = new TH1*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSync[ch] = dynamic_cast<TH1*>(getFromFile(Form("hBhTdcInSync_%03lu", ch)));
        hBhTdcInSync[ch]->SetDirectory(nullptr);
      }

      hHodTdcInSync = new TH1*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSync[ch] = dynamic_cast<TH1*>(getFromFile(Form("hHodTdcInSync_%03lu", ch)));
        hHodTdcInSync[ch]->SetDirectory(nullptr);
      } {
        hHodTdcInSync_Any = dynamic_cast<TH1*>(getFromFile("hHodTdcInSync_Any"));
        hHodTdcInSync_Any->SetDirectory(nullptr);
      }

      hExtTdcInSync = new TH1*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSync[ch] = dynamic_cast<TH1*>(getFromFile(Form("hExtTdcInSync_%03lu", ch)));
        hExtTdcInSync[ch]->SetDirectory(nullptr);
      } {
        hExtTdcInSync_Any = dynamic_cast<TH1*>(getFromFile("hExtTdcInSync_Any"));
        hExtTdcInSync_Any->SetDirectory(nullptr);
      }

      hTcTdcInSync = new TH1*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSync[ch] = dynamic_cast<TH1*>(getFromFile(Form("hTcTdcInSync_%03lu", ch)));
        hTcTdcInSync[ch]->SetDirectory(nullptr);
      }

      hVetoTdcInSync = new TH1*[Veto::NofChannels];
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoTdcInSync[ch] = dynamic_cast<TH1*>(getFromFile(Form("hVetoTdcInSync_%03lu", ch)));
        hVetoTdcInSync[ch]->SetDirectory(nullptr);
      }

      hBhMountain = new TH2*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhMountain[ch] = dynamic_cast<TH2*>(getFromFile(Form("hBhMountain_%03lu", ch)));
        hBhMountain[ch]->SetDirectory(nullptr);
      }

      hHodMountain = new TH2*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodMountain[ch] = dynamic_cast<TH2*>(getFromFile(Form("hHodMountain_%03lu", ch)));
        hHodMountain[ch]->SetDirectory(nullptr);
      } {
        hHodMountain_Any = dynamic_cast<TH2*>(getFromFile("hHodMountain_Any"));
        hHodMountain_Any->SetDirectory(nullptr);
      }

      hExtMountain = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtMountain[ch] = dynamic_cast<TH2*>(getFromFile(Form("hExtMountain_%03lu", ch)));
        hExtMountain[ch]->SetDirectory(nullptr);
      } {
        hExtMountain_Any = dynamic_cast<TH2*>(getFromFile("hExtMountain_Any"));
        hExtMountain_Any->SetDirectory(nullptr);
      }

      hTcMountain = new TH2*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcMountain[ch] = dynamic_cast<TH2*>(getFromFile(Form("hTcMountain_%03lu", ch)));
        hTcMountain[ch]->SetDirectory(nullptr);
      }

      hVetoMountain = new TH2*[Veto::NofChannels];
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoMountain[ch] = dynamic_cast<TH2*>(getFromFile(Form("hVetoMountain_%03lu", ch)));
        hVetoMountain[ch]->SetDirectory(nullptr);
      }

      hErrMountain = new TH2*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hErrMountain[ch] = dynamic_cast<TH2*>(getFromFile(Form("hErrMountain_%03lu", ch)));
        hErrMountain[ch]->SetDirectory(nullptr);
      }

      hMrSyncInterval = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval[ch] = dynamic_cast<TH1*>(getFromFile(Form("hMrSyncInterval_%03lu", ch)));
        hMrSyncInterval[ch]->SetDirectory(nullptr);
      }

      hMrSyncInterval2 = new TH2*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval2[ch] = dynamic_cast<TH2*>(getFromFile(Form("hMrSyncInterval2_%03lu", ch)));
        hMrSyncInterval2[ch]->SetDirectory(nullptr);
      }

      hBhTdcOffset = new TH2*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcOffset[ch] = dynamic_cast<TH2*>(getFromFile(Form("hBhTdcOffset_%03lu", ch)));
        hBhTdcOffset[ch]->SetDirectory(nullptr);
      }

      hTcTdcOffset = new TH2*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcOffset[ch] = dynamic_cast<TH2*>(getFromFile(Form("hTcTdcOffset_%03lu", ch)));
        hTcTdcOffset[ch]->SetDirectory(nullptr);
      }

      hExtTdcOffset = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcOffset[ch] = dynamic_cast<TH2*>(getFromFile(Form("hExtTdcOffset_%03lu", ch)));
        hExtTdcOffset[ch]->SetDirectory(nullptr);
      }

      file->Close();

      return 0;
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
      hHodEntriesByCh = new TH1D("hHodEntriesByCh",
                                 Form("%s, Hodoscope Entries by Channel;"
                                      "Channel;"
                                      "", tdcName.data()),
                                 Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5);
      hHodEntriesByCh->SetStats(false);

      // Extinction detector hit map
      hExtHitEvent   = ExtinctionDetector::CreateHitMap("hExtHitEvent");
      hExtHitMap     = ExtinctionDetector::CreateHitMap("hExtHitMap");
      lExtBorderLine = ExtinctionDetector::CreateBorderLine("lExtBorderLine", kBlack, kSolid, 1);

      // Extinction detector hit count
      hExtEntriesByCh = new TH1D("hExtEntriesByCh",
                                 Form("%s, Extinction Detector Entries by Channel;"
                                      "Channel;"
                                      "", tdcName.data()),
                                 ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5);
      hExtEntriesByCh->SetStats(false);

      hExtEntriesByChBottom  = hExtHitMap->ProjectionX("hExtEntriesByChBottom" , 1, 1);
      hExtEntriesByChCenter1 = hExtHitMap->ProjectionX("hExtEntriesByChCenter1", 2, 2);
      hExtEntriesByChCenter2 = hExtHitMap->ProjectionX("hExtEntriesByChCenter2", 3, 3);
      hExtEntriesByChTop     = hExtHitMap->ProjectionX("hExtEntriesByChTop"    , 4, 4);
      hExtEntriesByChBottom ->SetLineColor(kBlue   + 1);
      hExtEntriesByChCenter1->SetLineColor(kRed    + 1);
      hExtEntriesByChCenter2->SetLineColor(kOrange + 1);
      hExtEntriesByChTop    ->SetLineColor(kGreen  + 1);
      hExtEntriesByChBottom ->SetStats(false);
      hExtEntriesByChCenter1->SetStats(false);
      hExtEntriesByChCenter2->SetStats(false);
      hExtEntriesByChTop    ->SetStats(false);

      // Detectors hit count
      hEntriesInMrSyncByCh = new TH2D("hEntriesInMrSyncByCh",
                                      Form("%s, Detector Entries In MR Sync;"
                                           "Channel;"
                                           "Entries in MR Sync", tdcName.data()),
                                      Detectors::NofChannels, 0.0 - 0.5, Detectors::NofChannels - 0.5,
                                      20, 0 - 0.5, 20 - 0.5);
      hEntriesInMrSyncByCh->SetStats(false);

      hEntriesInMrSyncByDetector = new TH2D("hEntriesInMrSyncByDetector",
                                            Form("%s, Detector Entries In MR Sync;"
                                                 ";"
                                                 "Entries in MR Sync", tdcName.data()),
                                            Detectors::NofTypes, 0.0 - 0.5, Detectors::NofTypes - 0.5,
                                            100, 0 - 0.5, 100 - 0.5);
      hEntriesInMrSyncByDetector->SetStats(false);
      hEntriesInMrSyncByDetector->GetXaxis()->SetBinLabel(Detectors::Bh1 + 1, "BH1");
      hEntriesInMrSyncByDetector->GetXaxis()->SetBinLabel(Detectors::Bh2 + 1, "BH2");
      hEntriesInMrSyncByDetector->GetXaxis()->SetBinLabel(Detectors::Hod + 1, "Hod");
      hEntriesInMrSyncByDetector->GetXaxis()->SetBinLabel(Detectors::Ext + 1, "Ext");
      hEntriesInMrSyncByDetector->GetXaxis()->SetBinLabel(Detectors::Tc1 + 1, "TC1");
      hEntriesInMrSyncByDetector->GetXaxis()->SetBinLabel(Detectors::Tc2 + 1, "TC2");

      gEntriesInMrSyncByDetector = new TGraph();
      gEntriesInMrSyncByDetector->SetName("gEntriesInMrSyncByDetector");
      gEntriesInMrSyncByDetector->SetMarkerStyle(kPlus);
      gEntriesInMrSyncByDetector->SetMarkerSize(1.0);

      // Extinction detector hit count
      hExtEntriesInMrSyncInSpill = new TH2D("hExtEntriesInMrSyncInSpill",
                                            Form("%s, Extinction Detector Entries In MR Sync In Spill;"
                                                 "Time [ms];"
                                                 "Entries in MR Sync", tdcName.data()),
                                            xbinsInSpill, xminInSpill, xmaxInSpill,
                                            100, 0 - 0.5, 100 - 0.5);
      hExtEntriesInMrSyncInSpill->SetStats(false);

      gExtEntriesInMrSyncInSpill = new TGraph();
      gExtEntriesInMrSyncInSpill->SetName("gExtEntriesInMrSyncInSpill");
      gExtEntriesInMrSyncInSpill->SetMarkerStyle(kPlus);
      gExtEntriesInMrSyncInSpill->SetMarkerSize(0.6);

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

      // Veto TDC in spill
      hVetoTdcInSpill = new TH1*[Veto::NofChannels];
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoTdcInSpill[ch] = new TH1D(Form("hVetoTdcInSpill_%03lu", ch),
                                       Form("%s, Veto %ld TDC in Spill;"
                                            "Time [ms]", tdcName.data(), ch + 1),
                                       xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      // Error TDC in spill
      hErrTdcInSpill = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hErrTdcInSpill[ch] = new TH1D(Form("hErrTdcInSpill_%03lu", ch),
                                      Form("%s, Error %ld TDC in Spill;"
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

      // Veto TDC in sync
      hVetoTdcInSync = new TH1*[Veto::NofChannels];
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoTdcInSync[ch] = new TH1D(Form("hVetoTdcInSync_%03lu", ch),
                                      Form("%s, Veto %ld TDC in MR Sync;"
                                           "TDC [count]", tdcName.data(), ch + 1),
                                      xbinsInSync, xminInSync, xmaxInSync);
      }

      // Beamline Hodoscope Mountain Plot
      hBhMountain = new TH2*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhMountain[ch] = new TH2D(Form("hBhMountain_%03lu", ch),
                                   Form("%s, Beamline Hodoscope Mountain Plot @ %lu;"
                                        "TDC [count];"
                                        "Time [ms]", tdcName.data(), ch),
                                   xbinsInSync / 2, xminInSync, xmaxInSync,
                                   xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hBhMountain[ch]->SetStats(false);
      }

      // Hodoscope Mountain Plot
      hHodMountain = new TH2*[Hodoscope::NofChannels];
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodMountain[ch] = new TH2D(Form("hHodMountain_%03lu", ch),
                                    Form("%s, Hodoscope Mountain Plot @ %lu;"
                                         "TDC [count];"
                                         "Time [ms]", tdcName.data(), ch),
                                    xbinsInSync / 2, xminInSync, xmaxInSync,
                                    xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hHodMountain[ch]->SetStats(false);
      }

      hHodMountain_Any = new TH2D("hHodMountain_Any",
                                  Form("%s, Hodoscope Mountain Plot;"
                                       "TDC [count];"
                                       "Time [ms]", tdcName.data()),
                                  xbinsInSync / 2, xminInSync, xmaxInSync,
                                  xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hHodMountain_Any->SetStats(false);

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

      // Timing Counter Mountain Plot
      hTcMountain = new TH2*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hTcMountain[ch] = new TH2D(Form("hTcMountain_%03lu", ch),
                                   Form("%s, Timing Counter Mountain Plot @ %lu;"
                                        "TDC [count];"
                                        "Time [ms]", tdcName.data(), ch),
                                   xbinsInSync / 2, xminInSync, xmaxInSync,
                                   xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hTcMountain[ch]->SetStats(false);
      }

      // Veto Mountain Plot
      hVetoMountain = new TH2*[Veto::NofChannels];
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoMountain[ch] = new TH2D(Form("hVetoMountain_%03lu", ch),
                                     Form("%s, Veto Mountain Plot @ %lu;"
                                          "TDC [count];"
                                          "Time [ms]", tdcName.data(), ch),
                                     xbinsInSync / 2, xminInSync, xmaxInSync,
                                     xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hVetoMountain[ch]->SetStats(false);
      }

      // Error Mountain Plot
      hErrMountain = new TH2*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hErrMountain[ch] = new TH2D(Form("hErrMountain_%03lu", ch),
                                    Form("%s, Error Mountain Plot @ %lu;"
                                         "TDC [count];"
                                         "Time [ms]", tdcName.data(), ch),
                                    xbinsInSync / 2, xminInSync, xmaxInSync,
                                    xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hErrMountain[ch]->SetStats(false);
      }

      // Offset monitor hists
      hMrSyncInterval = new TH1*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval[ch] = new TH1D(Form("hMrSyncInterval_%03lu", ch),
                                       Form("%s, MR Sync TDC Interval @ ch%ld;"
                                            "TDC [count];"
                                            "", tdcName.data(), ch),
                                       xbinsInt, xminInt, xmaxInt);
      }

      hMrSyncInterval2 = new TH2*[MrSync::NofChannels];
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval2[ch] = new TH2D(Form("hMrSyncInterval2_%03lu", ch),
                                        Form("%s, MR Sync TDC Interval @ ch%ld;"
                                             "TDC [count];"
                                             "Time [ms]", tdcName.data(), ch),
                                        xbinsInt, xminInt, xmaxInt,
                                        xbinsInSpill, xminInSpill, xmaxInSpill);
      }

      hBhTdcOffset = new TH2*[BeamlineHodoscope::NofChannels];
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcOffset[ch] = new TH2D(Form("hBhTdcOffset_%03lu", ch),
                                    Form("%s, Beamline Hodoscope TDC Offset @ ch%ld;"
                                         "Channel;"
                                         "TDC [count]", tdcName.data(), ch),
                                    Detectors::NofChannels, 0 - 0.5, Detectors::NofChannels - 0.5,
                                    xbinsInDiff, -xmaxInDiff, -xminInDiff);
        hBhTdcOffset[ch]->SetStats(false);
      }

      hTcTdcOffset = new TH2*[TimingCounter::NofChannels];
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcOffset[ch] = new TH2D(Form("hTcTdcOffset_%03lu", ch),
                                    Form("%s, Timing Counter TDC Offset @ ch%ld;"
                                         "Channel;"
                                         "TDC [count]", tdcName.data(), ch),
                                    Detectors::NofChannels, 0 - 0.5, Detectors::NofChannels - 0.5,
                                    xbinsInDiff, -xmaxInDiff, -xminInDiff);
        hTcTdcOffset[ch]->SetStats(false);
      }

      hExtTdcOffset = new TH2*[ExtinctionDetector::NofChannels];
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcOffset[ch] = new TH2D(Form("hExtTdcOffset_%03lu", ch),
                                     Form("%s, Extinction Detector TDC Offset @ ch%ld;"
                                          "Channel;"
                                          "TDC [count]", tdcName.data(), ch),
                                     Detectors::NofChannels, 0 - 0.5, Detectors::NofChannels - 0.5,
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
      fSpillData.CreateBranch(fSpillTree);
    }

    void HistGenerator::DrawPlots(const std::string& ofilename, const std::string& ofilename_offset) {
      std::cout << "Draw plots" << std::endl;
      if (!gPad) {
        TCanvas::MakeDefCanvas();
      }
      gPad->SetGrid(true, true);
      gPad->SetLogy(false);
      gPad->Print((ofilename        + "[").data());
      gPad->Print((ofilename_offset + "[").data());

      Tron::ScopeSubstituter<Int_t> ss { gErrorIgnoreLevel, kWarning };

      gPad->SetGrid(false, false);
      {
        hHodHitMap->Draw("col");
        lHodBorderLine->Draw();
        hHodHitMap->SetMinimum(-0.001);
        gPad->Print(ofilename.data());
      }
      gPad->SetGrid(true, true);

      gPad->SetLogy(true);
      {
        hHodEntriesByCh->Draw();
        hHodEntriesByCh->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetGrid(false, false);
      {
        hExtHitMap->Draw("col");
        lExtBorderLine->Draw();
        hExtHitMap->SetMinimum(-0.001);
        gPad->Print(ofilename.data());
      }
      gPad->SetGrid(true, true);

      gPad->SetLogy(true);
      {
        hExtEntriesByCh->Draw();
        hExtEntriesByCh->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      } {
        hExtEntriesByChBottom ->Draw("hist");
        hExtEntriesByChCenter1->Draw("histsame");
        hExtEntriesByChCenter2->Draw("histsame");
        hExtEntriesByChTop    ->Draw("histsame");
        hExtEntriesByChBottom->SetMinimum(0.2);
        hExtEntriesByChBottom->SetMaximum(2.0 * hExtHitMap->GetBinContent(hExtHitMap->GetMaximumBin()));
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      gPad->SetLogz(true);
      {
        hEntriesInMrSyncByCh->Draw("colz");
        gPad->Print(ofilename.data());
      } {
        gEntriesInMrSyncByDetector->Set(0);
        for (Int_t xbin = 1, nbinsx = hEntriesInMrSyncByDetector->GetNbinsX(); xbin <= nbinsx; ++xbin) {
          Double_t sumcy = 0.0, sumc = 0.0;
          for (Int_t ybin = 2, nbinsy = hEntriesInMrSyncByDetector->GetNbinsY(); ybin <= nbinsy; ++ybin) {
            const Double_t y       = hEntriesInMrSyncByDetector->GetYaxis()->GetBinCenter(ybin);
            const Double_t content = hEntriesInMrSyncByDetector->GetBinContent(xbin, ybin);
            sumcy += content * y;
            sumc  += content;
          }
          if (sumc) {
            const Double_t x = hEntriesInMrSyncByDetector->GetXaxis()->GetBinCenter(xbin);
            gEntriesInMrSyncByDetector->SetPoint(gEntriesInMrSyncByDetector->GetN(), x, sumcy / sumc);
          }
        }

        hEntriesInMrSyncByDetector->Draw("colz");
        gEntriesInMrSyncByDetector->Draw("P");
        gPad->Print(ofilename.data());
      }
      gPad->SetLogz(true);

      gPad->SetLogz(true);
      {
        gExtEntriesInMrSyncInSpill->Set(0);
        for (Int_t xbin = 1, nbinsx = hExtEntriesInMrSyncInSpill->GetNbinsX(); xbin <= nbinsx; ++xbin) {
          Double_t sumcy = 0.0, sumc = 0.0;
          for (Int_t ybin = 2, nbinsy = hExtEntriesInMrSyncInSpill->GetNbinsY(); ybin <= nbinsy; ++ybin) {
            const Double_t y       = hExtEntriesInMrSyncInSpill->GetYaxis()->GetBinCenter(ybin);
            const Double_t content = hExtEntriesInMrSyncInSpill->GetBinContent(xbin, ybin);
            sumcy += content * y;
            sumc  += content;
          }
          if (sumc) {
            const Double_t x = hExtEntriesInMrSyncInSpill->GetXaxis()->GetBinCenter(xbin);
            gExtEntriesInMrSyncInSpill->SetPoint(gExtEntriesInMrSyncInSpill->GetN(), x, sumcy / sumc);
          }
        }

        hExtEntriesInMrSyncInSpill->Draw("colz");
        gExtEntriesInMrSyncInSpill->Draw("P");
        gPad->Print(ofilename.data());
      }
      gPad->SetLogz(true);

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

      std::cout << "hTcTdcInSpill" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcInSpill[ch]->GetEntries()) {
          hTcTdcInSpill[ch]->Draw();
          hTcTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hMrSyncTdcInSpill" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncTdcInSpill[ch]->GetEntries()) {
          hMrSyncTdcInSpill[ch]->Draw();
          hMrSyncTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hEvmTdcInSpill" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        if (hEvmTdcInSpill[ch]->GetEntries()) {
          hEvmTdcInSpill[ch]->Draw();
          hEvmTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hVetoTdcInSpill" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        if (hVetoTdcInSpill[ch]->GetEntries()) {
          hVetoTdcInSpill[ch]->Draw();
          hVetoTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hErrTdcInSpill" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hErrTdcInSpill[ch]->GetEntries()) {
          hErrTdcInSpill[ch]->Draw();
          hErrTdcInSpill[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hBhTdcInSync" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhTdcInSync[ch]->GetEntries()) {
          hBhTdcInSync[ch]->Draw();
          hBhTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hHodTdcInSync" << std::endl;
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

      std::cout << "hExtTdcInSync" << std::endl;
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

      std::cout << "hTcTdcInSync" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcInSync[ch]->GetEntries()) {
          hTcTdcInSync[ch]->Draw();
          hTcTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hVetoTdcInSync" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        if (hVetoTdcInSync[ch]->GetEntries()) {
          hVetoTdcInSync[ch]->Draw();
          hVetoTdcInSync[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hBhMountain" << std::endl;
      gPad->SetGrid(false, true);
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhMountain[ch]->GetEntries()) {
          hBhMountain[ch]->Draw("colz");
          hBhMountain[ch]->SetMinimum(0);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetGrid(true, true);

      std::cout << "hHodMountain" << std::endl;
      gPad->SetGrid(false, true);
      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        if (hHodMountain[ch]->GetEntries()) {
          hHodMountain[ch]->Draw("colz");
          hHodMountain[ch]->SetMinimum(0);
          gPad->Print(ofilename.data());
        }
      } {
        hHodMountain_Any->Draw("colz");
        hHodMountain_Any->SetMinimum(0);
        gPad->Print(ofilename.data());
      }
      gPad->SetGrid(true, true);

      std::cout << "hExtMountain" << std::endl;
      gPad->SetGrid(false, true);
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
      gPad->SetGrid(true, true);

      std::cout << "hTcMountain" << std::endl;
      gPad->SetGrid(false, true);
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcMountain[ch]->GetEntries()) {
          hTcMountain[ch]->Draw("colz");
          hTcMountain[ch]->SetMinimum(0);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetGrid(true, true);

      std::cout << "hVetoMountain" << std::endl;
      gPad->SetGrid(false, true);
      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        if (hVetoMountain[ch]->GetEntries()) {
          hVetoMountain[ch]->Draw("colz");
          hVetoMountain[ch]->SetMinimum(0);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetGrid(true, true);

      std::cout << "hErrMountain" << std::endl;
      gPad->SetGrid(false, true);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hErrMountain[ch]->GetEntries()) {
          hErrMountain[ch]->Draw("colz");
          hErrMountain[ch]->SetMinimum(0);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetGrid(true, true);

      std::cout << "hMrSyncInterval" << std::endl;
      gPad->SetLogy(true);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncInterval[ch]->GetEntries()) {
          hMrSyncInterval[ch]->Draw();
          hMrSyncInterval[ch]->SetMinimum(0.2);
          gPad->Print(ofilename.data());
        }
      }
      gPad->SetLogy(false);

      std::cout << "hMrSyncInterval2" << std::endl;
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (hMrSyncInterval2[ch]->GetEntries()) {
          hMrSyncInterval2[ch]->Draw("colz");
          gPad->Print(ofilename.data());
        }
      }

      std::cout << "hBhTdcOffset" << std::endl;
      gPad->SetLogz(true);
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        if (hBhTdcOffset[ch]->GetEntries()) {
          hBhTdcOffset[ch]->Draw("colz");
          hBhTdcOffset[ch]->SetMinimum(0);
       // hBhTdcOffset[ch]->SetMinimum(-0.001);
          gPad->Print(ofilename_offset.data());
        }
      }
      gPad->SetLogz(false);

      std::cout << "hTcTdcOffset" << std::endl;
      gPad->SetLogz(true);
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        if (hTcTdcOffset[ch]->GetEntries()) {
          hTcTdcOffset[ch]->Draw("colz");
          hTcTdcOffset[ch]->SetMinimum(0);
       // hTcTdcOffset[ch]->SetMinimum(-0.001);
          gPad->Print(ofilename_offset.data());
        }
      }
      gPad->SetLogz(false);

      std::cout << "hExtTdcOffset" << std::endl;
      gPad->SetLogz(true);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        if (hExtTdcOffset[ch]->GetEntries()) {
          hExtTdcOffset[ch]->Draw("colz");
          hExtTdcOffset[ch]->SetMinimum(0);
       // hExtTdcOffset[ch]->SetMinimum(-0.001);
          gPad->Print(ofilename_offset.data());
        }
      }
      gPad->SetLogz(false);

      gPad->Print((ofilename        + "]").data());
      gPad->Print((ofilename_offset + "]").data());
    }

    void HistGenerator::WritePlots(const std::string& ofilename) {
      std::cout << "Write plots" << std::endl;
      TFile* file = new TFile(ofilename.data(), "RECREATE");
      if (!file->IsOpen()) {
        std::cout << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      {
        Tron::ObjectHelper::WriteValue<Long64_t>(fSpillData.Date.Convert(), "Time");
      }

      {
        hHodHitMap->Write();
        lHodBorderLine->Write(lHodBorderLine->GetName(), TObject::kSingleKey);
      }

      {
        hHodEntriesByCh->Write();
      }

      {
        hExtHitEvent  ->Write();
        hExtHitMap    ->Write();
        lExtBorderLine->Write(lExtBorderLine->GetName(), TObject::kSingleKey);
      }

      {
        hExtEntriesByCh->Write();
      } {
        hExtEntriesByChBottom ->Write();
        hExtEntriesByChCenter1->Write();
        hExtEntriesByChCenter2->Write();
        hExtEntriesByChTop    ->Write();
      }

      {
        hEntriesInMrSyncByCh      ->Write();
        hEntriesInMrSyncByDetector->Write();
      }

      {
        hExtEntriesInMrSyncInSpill->Write();
      }

      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSpill[ch]->Write();
      }

      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSpill[ch]->Write();
      } {
        hHodTdcInSpill_Any->Write();
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSpill[ch]->Write();
      } {
        hExtTdcInSpill_Any->Write();
      }

      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSpill[ch]->Write();
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncTdcInSpill[ch]->Write();
      }

      for (std::size_t ch = 0; ch < EventMatch::NofChannels; ++ch) {
        hEvmTdcInSpill[ch]->Write();
      }

      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoTdcInSpill[ch]->Write();
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hErrTdcInSpill[ch]->Write();
      }

      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcInSync[ch]->Write();
      }

      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodTdcInSync[ch]->Write();
      } {
        hHodTdcInSync_Any->Write();
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcInSync[ch]->Write();
      } {
        hExtTdcInSync_Any->Write();
      }

      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcInSync[ch]->Write();
      }

      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoTdcInSync[ch]->Write();
      }

      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhMountain[ch]->Write();
      }

      for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
        hHodMountain[ch]->Write();
      } {
        hHodMountain_Any->Write();
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtMountain[ch]->Write();
      } {
        hExtMountain_Any->Write();
      }

      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcMountain[ch]->Write();
      }

      for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
        hVetoMountain[ch]->Write();
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hErrMountain[ch]->Write();
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval[ch]->Write();
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hMrSyncInterval2[ch]->Write();
      }

      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        hBhTdcOffset[ch]->Write();
      }

      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        hTcTdcOffset[ch]->Write();
      }

      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hExtTdcOffset[ch]->Write();
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

    void HistGenerator::WriteMrSyncInterval(const std::string& ofilename) {
      std::cout << "Write MR Sync interval" << std::endl;
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cerr << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        ofile << ch << "\t" <<  Form("%23.15e", fSpillData.MrSyncInterval[ch]) << std::endl;
      }

      ofile.close();

      std::cerr << "Info in <HistGenerator::WriteMrSyncInterval>: dat file " << ofilename << " has been created" << std::endl;
    }

    void HistGenerator::WriteTdcOffsets(const std::string& ofilename) {
      std::cout << "Write tdc offset" << std::endl;
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cerr << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      for (auto&& pair : fTdcOffsets) {
        auto& gch   = pair.first;
        auto& value = pair.second;
        ofile << gch << "\t" << value << std::endl;
      }

      ofile.close();

      std::cerr << "Info in <HistGenerator::WriteTdcOffsets>: dat file " << ofilename << " has been created" << std::endl;
    }

    void HistGenerator::ClearLastSpill(Bool_t clearHists) {
      fLastExtData.clear();
      fLastHodData.clear();
      fLastTcData .clear();
      fLastBhData .clear();
      fLastMrSyncData.clear();

      if (clearHists) {
        hHodHitMap     ->Reset();
        hHodEntriesByCh->Reset();

        hExtHitMap            ->Reset();
        hExtEntriesByCh       ->Reset();
        hExtEntriesByChBottom ->Reset();
        hExtEntriesByChCenter1->Reset();
        hExtEntriesByChCenter2->Reset();
        hExtEntriesByChTop    ->Reset();

        hEntriesInMrSyncByCh      ->Reset();
        hEntriesInMrSyncByDetector->Reset();

        hExtEntriesInMrSyncInSpill->Reset();

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

        for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
          hVetoTdcInSpill[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hErrTdcInSpill[ch]->Reset();
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

        for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
          hVetoTdcInSync[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhMountain[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodMountain[ch]->Reset();
        } {
          hHodMountain_Any->Reset();
        }

        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtMountain[ch]->Reset();
        } {
          hExtMountain_Any->Reset();
        }

        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcMountain[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < Veto::NofChannels; ++ch) {
          hVetoMountain[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hErrMountain[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncInterval[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncInterval2[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhTdcOffset[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcTdcOffset[ch]->Reset();
        }

        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcOffset[ch]->Reset();
        }
      }
    }

    std::size_t HistGenerator::RemoveOldTdc(std::vector<TdcData>* lastData, const TdcData& tdc) {
      for (std::size_t i = 0, n = lastData->size(); i < n; ++i) {
        TdcData& lastTdc = lastData->at(0);
        const Double_t dt = TdcData::GetTimeDifference(tdc, lastTdc);
        if (std::abs(dt) > fHistoryWidth) {
          lastData->erase(lastData->begin());
        } else {
          break;
        }
      }
      return lastData->size();
    }

    Int_t HistGenerator::GeneratePlots(MargedReader* reader) {
      const clock_t startClock = clock();

      std::cout << "Initialize history" << std::endl;
      ClearLastSpill(true);

      // std::cout << "[debug] check end of spill" << std::endl;
      // Throw away first mr sync
      std::map<Tag_t, TdcData> tdcDataInMrSync;
      reader->Read(tdcDataInMrSync);
      tdcDataInMrSync.clear();

      fSpillData.SetDate(reader->GetDate());
      fSpillData.EMCount = reader->GetEMCount();

      std::map<std::size_t, Long64_t> entriesInMrSyncByCh;
      std::map<Int_t, Long64_t> entriesInMrSyncByDetector;
      
      while (true) {
        for (; reader->Read(tdcDataInMrSync); tdcDataInMrSync.clear()) {
          // std::cout << "[debug] data process" << std::endl;
          for (auto&& pair : tdcDataInMrSync) {
            // auto& tag  = pair.first;
            auto& data = pair.second;
            const Int_t    board   = data.Board;
            const Int_t    gch     = data.Channel;
            const Double_t time    = data.Time;
            const Long64_t syncTdc = fLastMrSyncData[board].Tdc;

            if (data.Channel < 0) {
              hErrTdcInSpill[board]->Fill(time / msec);
              for (Int_t xbin = 0, nbinsx = hErrMountain[board]->GetNbinsX(); xbin < nbinsx; ++xbin) {
                const Double_t dtdc = hErrMountain[board]->GetXaxis()->GetBinCenter(xbin);
                hErrMountain[board]->Fill(dtdc, time / msec);
              }
              continue;
            }

            if (Detectors::Contains(data.Channel)) {
              ++entriesInMrSyncByCh[data.Channel];
            }

            if (BeamlineHodoscope::Contains(data.Channel)) {
              // std::cout << "[debug] beamline hodoscope" << std::endl;
              const Int_t    ch  = BeamlineHodoscope::GetChannel(gch);
              const Long64_t tdc = data.Tdc;

              ++entriesInMrSyncByDetector[Detectors::Bh1 + ch];

              hBhTdcInSpill[ch]->Fill(time / msec);

              if (syncTdc) {
                hBhTdcInSync[ch]->Fill(tdc - syncTdc);
                hBhMountain [ch]->Fill(tdc - syncTdc, time / msec);

                if (fOffsetFromBunch) {
                  const Bool_t thisIsInBunch = IsInBunch(tdc - syncTdc);

                  for (auto&& lastData : fLastBhData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hBhTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hBhTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastHodData) {
                    auto lastGch     = lastData.Channel;
                    // auto lastCh      = Hodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hBhTdcOffset[ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (tdc - syncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = ExtinctionDetector::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hBhTdcOffset [    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hExtTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = TimingCounter::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hBhTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hTcTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  
                } else {

                  for (auto&& lastData : fLastBhData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hBhTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hBhTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastHodData) {
                    auto lastGch     = lastData.Channel;
                    // auto lastCh      = Hodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hBhTdcOffset[ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (tdc - syncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = ExtinctionDetector::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hBhTdcOffset [    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hExtTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = TimingCounter::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hBhTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hTcTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                }

              }

              fLastBhData.push_back(data);

              if (RemoveOldTdc(&fLastBhData, data) > kHistLimit) {
                std::cerr << "[error] size of fLastBhData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (Hodoscope::Contains(data.Channel)) {
              // std::cout << "[debug] hodoscope" << std::endl;
              const Int_t    ch  = Hodoscope::GetChannel(gch);
              const Long64_t tdc = data.Tdc;

              ++entriesInMrSyncByDetector[Detectors::Hod];

              hHodEntriesByCh   ->Fill(ch);
              Hodoscope         ::Fill(hHodHitMap, ch);

              hHodTdcInSpill[ch]->Fill(time / msec);
              hHodTdcInSpill_Any->Fill(time / msec);

              if (syncTdc) {
                hHodTdcInSync[ch]->Fill(tdc - syncTdc);
                hHodTdcInSync_Any->Fill(tdc - syncTdc);
                hHodMountain [ch]->Fill(tdc - syncTdc, time / msec);
                hHodMountain_Any ->Fill(tdc - syncTdc, time / msec);

                if (fOffsetFromBunch) {
                  // const Bool_t thisIsInBunch = IsInBunch(tdc - syncTdc);

                  for (auto&& lastData : fLastBhData) {
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastData.Channel);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hBhTdcOffset[lastCh]->Fill(gch, (tdc - syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastCh      = ExtinctionDetector::GetChannel(lastData.Channel);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hExtTdcOffset[lastCh]->Fill(gch, (tdc - syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastCh      = TimingCounter::GetChannel(lastData.Channel);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hTcTdcOffset[lastCh]->Fill(gch, (tdc - syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }

                } else {

                  for (auto&& lastData : fLastBhData) {
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastData.Channel);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hBhTdcOffset[lastCh]->Fill(gch, (tdc - syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastCh      = ExtinctionDetector::GetChannel(lastData.Channel);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hExtTdcOffset[lastCh]->Fill(gch, (tdc - syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastCh      = TimingCounter::GetChannel(lastData.Channel);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hTcTdcOffset[lastCh]->Fill(gch, (tdc - syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }

                }
              }

              fLastHodData.push_back(data);

              if (RemoveOldTdc(&fLastHodData, data) > kHistLimit) {
                std::cerr << "[error] size of fLastHodData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (ExtinctionDetector::Contains(data.Channel)) {
              // std::cout << "[debug] extinction detector" << std::endl;
              const Int_t    ch  = ExtinctionDetector::GetChannel(gch);
              const Long64_t tdc = data.Tdc;

              ++entriesInMrSyncByDetector[Detectors::Ext];

              hExtEntriesByCh   ->Fill(ch);
              ExtinctionDetector::Fill(hExtHitMap, ch);

              hExtTdcInSpill[ch]->Fill(time / msec);
              hExtTdcInSpill_Any->Fill(time / msec);

              if (syncTdc) {
                hExtTdcInSync[ch]->Fill(tdc - syncTdc);
                hExtTdcInSync_Any->Fill(tdc - syncTdc);
                hExtMountain [ch]->Fill(tdc - syncTdc, time / msec);
                hExtMountain_Any ->Fill(tdc - syncTdc, time / msec);

                if (fOffsetFromBunch) {
                  const Bool_t thisIsInBunch = IsInBunch(tdc - syncTdc);

                  for (auto&& lastData : fLastBhData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hExtTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hBhTdcOffset [lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastHodData) {
                    auto lastGch     = lastData.Channel;
                    // auto lastCh      = Hodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hExtTdcOffset[ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (tdc - syncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = ExtinctionDetector::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hExtTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hExtTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = TimingCounter::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hExtTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hTcTdcOffset [lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }

                } else {

                  for (auto&& lastData : fLastBhData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hExtTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hBhTdcOffset [lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastHodData) {
                    auto lastGch     = lastData.Channel;
                    // auto lastCh      = Hodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hExtTdcOffset[ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (tdc - syncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = ExtinctionDetector::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hExtTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hExtTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = TimingCounter::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hExtTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hTcTdcOffset [lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }

                }
              }

              if (RemoveOldTdc(&fLastExtData, data) > kHistLimit) {
                std::cerr << "[error] size of fLastExtData reaches " << kHistLimit << std::endl;
                return 1;
              }
              fLastExtData.push_back(data);

            } else if (TimingCounter::Contains(data.Channel)) {
              // std::cout << "[debug] timing counter" << std::endl;
              const Int_t ch = TimingCounter::GetChannel(gch);
              const Long64_t tdc = data.Tdc;

              ++entriesInMrSyncByDetector[Detectors::Tc1 + ch];

              hTcTdcInSpill[ch]->Fill(time / msec);

              if (syncTdc) {
                hTcTdcInSync[ch]->Fill(tdc - syncTdc);
                hTcMountain [ch]->Fill(tdc - syncTdc, time / msec);

                if (fOffsetFromBunch) {
                  const Bool_t thisIsInBunch = IsInBunch(tdc - syncTdc);

                  for (auto&& lastData : fLastBhData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hTcTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hBhTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastHodData) {
                    auto lastGch     = lastData.Channel;
                    // auto lastCh      = Hodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hTcTdcOffset[ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (tdc - syncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = ExtinctionDetector::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hTcTdcOffset [    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hExtTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = TimingCounter::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      if (thisIsInBunch)
                        hTcTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      if (IsInBunch(lastData.Tdc - lastSyncTdc))
                        hTcTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }

                } else {
                
                  for (auto&& lastData : fLastBhData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = BeamlineHodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hTcTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hBhTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastHodData) {
                    auto lastGch     = lastData.Channel;
                    // auto lastCh      = Hodoscope::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hTcTdcOffset[ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (tdc - syncTdc));
                    }
                  }
                  for (auto&& lastData : fLastExtData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = ExtinctionDetector::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hTcTdcOffset [    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hExtTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }
                  for (auto&& lastData : fLastTcData) {
                    auto lastGch     = lastData.Channel;
                    auto lastCh      = TimingCounter::GetChannel(lastGch);
                    auto lastSyncTdc = fLastMrSyncData[lastData.Board].Tdc;
                    if (lastSyncTdc) {
                      hTcTdcOffset[    ch]->Fill(lastGch, (lastData.Tdc - lastSyncTdc) - (         tdc -     syncTdc));
                      hTcTdcOffset[lastCh]->Fill(    gch, (         tdc -     syncTdc) - (lastData.Tdc - lastSyncTdc));
                    }
                  }

                }
              }

              fLastTcData.push_back(data);

              if (RemoveOldTdc(&fLastTcData, data) > kHistLimit) {
                std::cerr << "[error] size of fLastTcData reaches " << kHistLimit << std::endl;
                return 1;
              }

            } else if (Veto::Contains(data.Channel)) {
              // std::cout << "[debug] veto" << std::endl;
              const Int_t ch = Veto::GetChannel(gch);
              const Long64_t tdc = data.Tdc;

              hVetoTdcInSpill[ch]->Fill(time / msec);

              if (syncTdc) {
                hVetoTdcInSync[ch]->Fill(tdc - syncTdc);
                hVetoMountain [ch]->Fill(tdc - syncTdc, time / msec);
              }

            } else if (MrSync::Contains(data.Channel)) {
              // std::cout << "[debug] mrsync" << std::endl;
              const Int_t ch = MrSync::GetChannel(gch);

              hMrSyncTdcInSpill[ch]->Fill(time / msec);

              decltype(fLastMrSyncData)::const_iterator itr;
              if ((itr = fLastMrSyncData.find(data.Board)) != fLastMrSyncData.end()) {
                if (itr->second.Tdc) {
                  hMrSyncInterval [data.Board]->Fill(data.Tdc - itr->second.Tdc);
                  hMrSyncInterval2[data.Board]->Fill(data.Tdc - itr->second.Tdc, data.Time / msec);
                }
              }

              fLastMrSyncData[board] = data;

            } else if (EventMatch::Contains(gch)) {
              const Int_t ch = EventMatch::GetChannel(gch);

              hEvmTdcInSpill[ch]->Fill(time / msec);

            } else {
              // std::cout << "[debug] skip others" << std::endl;
              continue;

            }
          }

          for (std::size_t gch = 0; gch < Detectors::NofChannels; ++gch) {
            hEntriesInMrSyncByCh->Fill(gch, entriesInMrSyncByCh[gch]);
            entriesInMrSyncByCh[gch] = 0;
          }

          hExtEntriesInMrSyncInSpill->Fill(fLastMrSyncData.begin()->second.Time / msec, entriesInMrSyncByDetector[Detectors::Ext]);
          for (std::size_t detector = 0; detector < Detectors::NofTypes; ++detector) {
            hEntriesInMrSyncByDetector->Fill(detector, entriesInMrSyncByDetector[detector]);
            entriesInMrSyncByDetector[detector] = 0;
          }

        }

        // std::cout << "[debug] check end of spill" << std::endl;
        // Check end of spill
        if (reader->IsSpillEnded()) {
          std::cout << "[info] end of spill" << std::endl;

          // Calc spill summary
          CalcEntries();
          CalcMrSyncInterval();
          CalcTdcOffsets();

          // Fill spill summary
          if (fSpillTree) {
            fSpillTree->Fill();
          }

          reader->ClearLastSpill();
          // std::cout << "[debug] Throw away first mr sync" << std::endl;
          // Throw away first mr sync
          reader->Read(tdcDataInMrSync);
          tdcDataInMrSync.clear();

          if (!reader->IsFileEnded()) {
            fSpillData.SetDate(reader->GetDate());
            fSpillData.EMCount = reader->GetEMCount();
          }
        }

        // std::cout << "[debug] check end of file" << std::endl;
        // Check end of file
        if (reader->IsFileEnded()) {
          std::cout << "[info] end of file" << std::endl;

          // Get projections
          for (Int_t xbin = 1, nbinsx = hExtHitMap->GetNbinsX(); xbin <= nbinsx; ++xbin) {
            hExtEntriesByChBottom ->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 1));
            hExtEntriesByChCenter1->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 2));
            hExtEntriesByChCenter2->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 3));
            hExtEntriesByChTop    ->SetBinContent(xbin, hExtHitMap->GetBinContent(xbin, 4));
            hExtEntriesByChBottom ->SetBinError  (xbin, hExtHitMap->GetBinError  (xbin, 1));
            hExtEntriesByChCenter1->SetBinError  (xbin, hExtHitMap->GetBinError  (xbin, 2));
            hExtEntriesByChCenter2->SetBinError  (xbin, hExtHitMap->GetBinError  (xbin, 3));
            hExtEntriesByChTop    ->SetBinError  (xbin, hExtHitMap->GetBinError  (xbin, 4));
          }

          break;
        }
      }

      const clock_t stopClock = clock();
      std::cout << "time: " << (double)(stopClock - startClock) / CLOCKS_PER_SEC << " sec\n";

      return 0;
    }

    void HistGenerator::CalcEntries() {
      for (std::size_t ch = 0; ch < BeamlineHodoscope ::NofChannels; ++ch) {
        fSpillData.Entries[ch + BeamlineHodoscope ::GlobalChannelOffset] = hBhTdcInSpill    [ch]->GetEntries();
      }
      for (std::size_t ch = 0; ch < Hodoscope         ::NofChannels; ++ch) {
        fSpillData.Entries[ch + Hodoscope         ::GlobalChannelOffset] = hHodTdcInSpill   [ch]->GetEntries();
      }
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        fSpillData.Entries[ch + ExtinctionDetector::GlobalChannelOffset] = hExtTdcInSpill   [ch]->GetEntries();
      }
      for (std::size_t ch = 0; ch < TimingCounter     ::NofChannels; ++ch) {
        fSpillData.Entries[ch + TimingCounter     ::GlobalChannelOffset] = hTcTdcInSpill    [ch]->GetEntries();
      }
      for (std::size_t ch = 0; ch < MrSync            ::NofChannels; ++ch) {
        fSpillData.Entries[ch + MrSync            ::GlobalChannelOffset] = hMrSyncTdcInSpill[ch]->GetEntries();
      }
      for (std::size_t ch = 0; ch < EventMatch        ::NofChannels; ++ch) {
        fSpillData.Entries[ch + EventMatch        ::GlobalChannelOffset] = hEvmTdcInSpill   [ch]->GetEntries();
      }
    }

    void HistGenerator::CalcMrSyncInterval() {
      std::cout << "_____ MR Sync Interval _____" << std::endl;
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        if (!hMrSyncInterval[ch]->GetEntries()) {
          fSpillData.MrSyncInterval[ch] = 0.0;
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
          fSpillData.MrSyncInterval[ch] = xmean;
        }
      }
    }

    void HistGenerator::CalcTdcOffsets() {
      std::cout << "_____ TDC Offsets _____" << std::endl;

      const std::size_t targetGch  = TimingCounter::GlobalChannelOffset;
      const TH2*        hTdcOffset = hTcTdcOffset[0];
      const Int_t       ybins      = hTdcOffset->GetNbinsY();

      auto getOffsetBin =
        [&] (Int_t xbin1, Int_t xbin2) {
          Double_t maxsum = 0.0; Int_t maxbin = 0;
          for (Int_t ybin = 1; ybin <= ybins; ++ybin) {
            Double_t sum = 0.0;
            for (Int_t xbin = xbin1; xbin <= xbin2; ++xbin) {
              sum += hTdcOffset->GetBinContent(xbin, ybin);
            }
            if (sum > maxsum) {
              maxbin = ybin;
              maxsum = sum;
            }
          }
          return maxbin;
        };

      // Initialize offset tdc
      fTdcOffsets.clear();
      for (std::size_t gch = 0; gch < Detectors::NofChannels; ++gch) {
        fTdcOffsets[gch] = 0;
      }

      // Extinction Detector
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        const std::size_t gch = ch + ExtinctionDetector::GlobalChannelOffset;
        if (gch != targetGch) {
          const Int_t xbin1 = gch + 1;
          const Int_t xbin2 = gch + 1;
          if (Int_t offsetBin = getOffsetBin(xbin1, xbin2)) {
            const Double_t offset = hTdcOffset->GetYaxis()->GetBinCenter(offsetBin);
            fTdcOffsets[gch] = offset;
          }
        }
      }

      // Hodoscope
      {
        const Int_t xbin1 = Hodoscope::GlobalChannelOffset + 1;
        const Int_t xbin2 = Hodoscope::GlobalChannelOffset + Hodoscope::NofChannels;
        if (Int_t offsetBin = getOffsetBin(xbin1, xbin2)) {
          const Double_t offset = hTdcOffset->GetYaxis()->GetBinCenter(offsetBin);
          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            const std::size_t gch = ch + Hodoscope::GlobalChannelOffset;
            fTdcOffsets[gch] = offset;
          }
        }
      }

      // Timing Counter
      for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
        const std::size_t gch = ch + TimingCounter::GlobalChannelOffset;
        if (gch != targetGch) {
          const Int_t xbin1 = gch + 1;
          const Int_t xbin2 = gch + 1;
          if (Int_t offsetBin = getOffsetBin(xbin1, xbin2)) {
            const Double_t offset = hTdcOffset->GetYaxis()->GetBinCenter(offsetBin);
            fTdcOffsets[gch] = offset;
          }
        }
      }

      // Beamline Hodoscope
      for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
        const std::size_t gch = ch + BeamlineHodoscope::GlobalChannelOffset;
        if (gch != targetGch) {
          const Int_t xbin1 = gch + 1;
          const Int_t xbin2 = gch + 1;
          if (Int_t offsetBin = getOffsetBin(xbin1, xbin2)) {
            const Double_t offset = hTdcOffset->GetYaxis()->GetBinCenter(offsetBin);
            fTdcOffsets[gch] = offset;
          }
        }
      }

      // Adjust offset to reference
      const Long64_t refOffset = fTdcOffsets[fRefExtChannel + ExtinctionDetector::GlobalChannelOffset];
      for (auto&& pair : fTdcOffsets) {
        pair.second -= refOffset;
        std::cout << pair.first << "\t" << pair.second << std::endl;
      }
      
    }

  }

}

#endif
