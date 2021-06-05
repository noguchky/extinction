#ifndef Extinction_TimelineCoincidence_hh
#define Extinction_TimelineCoincidence_hh

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
#include "Spill.hh"

#include "Linq.hh"
#include "String.hh"
#include "ScopeSubstituter.hh"
#include "MargedReader.hh"

namespace Extinction {

  namespace Analyzer {

    class Timeline_t {
    private:
      std::size_t           fSize;
      Long64_t*             pTdcs = nullptr;

    public:
      std::vector<Long64_t> fTdcs;
      std::vector<Int_t>    fMrSyncCounts;
      std::vector<Long64_t> fTdcFromMrSyncs;
      std::vector<Int_t   > fTots;

    public:
      void Resize(std::size_t size) {
        fSize = size;
        fTdcs          .resize(fSize);
        fMrSyncCounts  .resize(fSize);
        fTdcFromMrSyncs.resize(fSize);
        fTots          .resize(fSize);
        pTdcs = fTdcs.data();
        Clear();
      }

      inline std::size_t Size() const {
        return fSize;
      }

      inline void Clear() {
        std::fill(fTdcs          .begin(), fTdcs          .end(), 0);
        std::fill(fMrSyncCounts  .begin(), fMrSyncCounts  .end(), 0);
        std::fill(fTdcFromMrSyncs.begin(), fTdcFromMrSyncs.end(), 0);
        std::fill(fTots          .begin(), fTots          .end(), 0);
      }

      inline Bool_t HasHit(std::size_t tdc) const {
        return pTdcs[tdc];
      }

      inline void FillHit(TdcData& data, Int_t width) {
        const Int_t x = data.TdcFromMrSync;
        const Int_t last = std::min((std::size_t)(x + width + 1), fSize);
        if (x < 0 || last <= x) { return; }
        std::fill(fTdcs          .begin() + x, fTdcs          .begin() + last, data.Tdc            );
        std::fill(fMrSyncCounts  .begin() + x, fMrSyncCounts  .begin() + last, data.LastMrSyncCount);
        std::fill(fTdcFromMrSyncs.begin() + x, fTdcFromMrSyncs.begin() + last, data.TdcFromMrSync  );
        std::fill(fTots          .begin() + x, fTots          .begin() + last, data.Tot            );
      }

      inline void Shift(std::size_t length) {
        const std::size_t shift = fSize - length;
        std::memmove(fTdcs          .data(), fTdcs          .data() + length, sizeof(Long64_t) * shift);
        std::memmove(fMrSyncCounts  .data(), fMrSyncCounts  .data() + length, sizeof(Int_t   ) * shift);
        std::memmove(fTdcFromMrSyncs.data(), fTdcFromMrSyncs.data() + length, sizeof(Long64_t) * shift);
        std::memmove(fTots          .data(), fTots          .data() + length, sizeof(Int_t   ) * shift);
        std::fill(fTdcs          .begin() + shift, fTdcs          .end(), 0);
        std::fill(fMrSyncCounts  .begin() + shift, fMrSyncCounts  .end(), 0);
        std::fill(fTdcFromMrSyncs.begin() + shift, fTdcFromMrSyncs.end(), 0);
        std::fill(fTots          .begin() + shift, fTots          .end(), 0);
      }
    };
    
    class TimelineCoincidence {
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
      ITdcDataProvider*            fProvider               = nullptr;

      Timeline_t                   fBh1Timeline;
      Timeline_t                   fBh2Timeline;
      Timeline_t                   fTc1Timeline;
      Timeline_t                   fTc2Timeline;
      std::vector<Timeline_t>      fExtTimeline;

      TH1I*                        hTmpBh1Timeline         = nullptr;
      TH1I*                        hTmpBh2Timeline         = nullptr;
      TH1I*                        hTmpTc1Timeline         = nullptr;
      TH1I*                        hTmpTc2Timeline         = nullptr;
      std::vector<TH1I*>           hTmpExtTimeline;
      TH1I*                        hTmpExtTimeline_Any     = nullptr;
      TH1I*                        hTmpSumTimeline         = nullptr;
      TLegend*                     hTmpTimelineLegend      = nullptr;

      TH1D*                        hCoinTlTdcInSync        = nullptr;
      TH2D*                        hCoinTlMountain         = nullptr;
      TH1D*                        hCoinExtTdcInSync       = nullptr;
      TH2D*                        hCoinExtMountain        = nullptr;

      TH1*                         hEfficiencyFrame        = nullptr;
      TEfficiency*                 hEfficiency             = nullptr;

      TFile*                       fCoinFile               = nullptr;
      TTree*                       fCoinTree               = nullptr;
      TCanvas*                     fTimelineCanvas         = nullptr;
      TCanvas*                     fHitMapCanvas           = nullptr;

      TH2*                         fExtHitMap              = nullptr;
      TList*                       fExtBorderLine          = nullptr;

      // TFile*                       fSpillFile              = nullptr;
      // TTree*                       fSpillTree              = nullptr;

      // Long64_t                     fSpillCount             = 0;

      Bool_t                       fCoincidenceTargetBh1;
      Bool_t                       fCoincidenceTargetBh2;
      Bool_t                       fCoincidenceTargetHod;
      Bool_t                       fCoincidenceTargetExt;
      Bool_t                       fCoincidenceTargetTc1;
      Bool_t                       fCoincidenceTargetTc2;

      Bool_t                       fEfficiencyTargetBh1;
      Bool_t                       fEfficiencyTargetBh2;
      Bool_t                       fEfficiencyTargetHod;
      Bool_t                       fEfficiencyTargetExt;
      Bool_t                       fEfficiencyTargetTc1;
      Bool_t                       fEfficiencyTargetTc2;
      
      Double_t                     fCoinTimeWidth          =  50.0 * nsec;

      // Coincidence data
      Long64_t                   fDate;
      Long64_t                   fEMCount;
      UInt_t                     fCoinTarget;
      Int_t                      fMrSyncCount;
      Long64_t                   fTdc;
      Long64_t                   fTdcFromMrSync;
      Int_t                      fBlurWidth;
      Int_t                      fCoinWidth;
      Int_t                      fNofHits;
      Int_t                      fHitChannels      [ExtinctionDetector::NofChannels];
      Int_t                      fHitMrSyncCounts  [ExtinctionDetector::NofChannels];
      Long64_t                   fHitTdcFromMrSyncs[ExtinctionDetector::NofChannels];
      Int_t                      fHitTots          [ExtinctionDetector::NofChannels];

      void InitializeCoinData() {
        fDate        = 0;
        fEMCount     = 0;
        fCoinTarget  = 0;
        fMrSyncCount = 0;
        ClearCoinEvent();
      }

      void ClearCoinEvent() {
        fTdc           = 0;
        fTdcFromMrSync = 0;
        fBlurWidth     = 0;
        fCoinWidth     = 0;
        fNofHits       = 0;
        std::memset(fHitChannels      , 0, sizeof(fHitChannels      ));
        std::memset(fHitMrSyncCounts  , 0, sizeof(fHitMrSyncCounts  ));
        std::memset(fHitTdcFromMrSyncs, 0, sizeof(fHitTdcFromMrSyncs));
        std::memset(fHitTots          , 0, sizeof(fHitTots          ));
      }
      
    public:
      TimelineCoincidence(ITdcDataProvider* provider);
      ~TimelineCoincidence();

      void                 SetCoincidenceTarget(const std::vector<Int_t>& flags) {
        std::cout << "SetCoincidenceTarget" << std::endl;
        if (flags.size() < 6) {
          std::cerr << "[error] invalid coincidence target" << std::endl;
          return;
        }
        fCoincidenceTargetBh1 = flags[0];
        fCoincidenceTargetBh2 = flags[1];
        fCoincidenceTargetHod = flags[2];
        fCoincidenceTargetExt = flags[3];
        fCoincidenceTargetTc1 = flags[4];
        fCoincidenceTargetTc2 = flags[5];
        std::cout << "Bh1 " << fCoincidenceTargetBh1 << std::endl;
        std::cout << "Bh2 " << fCoincidenceTargetBh2 << std::endl;
        std::cout << "Hod " << fCoincidenceTargetHod << std::endl;
        std::cout << "Ext " << fCoincidenceTargetExt << std::endl;
        std::cout << "Tc1 " << fCoincidenceTargetTc1 << std::endl;
        std::cout << "Tc2 " << fCoincidenceTargetTc2 << std::endl;
      }

      inline void          SetCoinTimeWidth(Double_t width) {
        std::cout << "SetCoinTimeWidth ... " << width / nsec << " nsec" << std::endl;
        fCoinTimeWidth = width;
      }
      inline Double_t      GetCoinTimeWidth() const { return fCoinTimeWidth; }

      void                 InitializePlots(const PlotsProfiles& profile);
      void                 InitializeCoinTree(const std::string& filename, const std::string& treename = "ctree");

      void                 DrawPlots(const std::string& ofilename);
      void                 WritePlots(const std::string& ofilename);
      void                 WriteCoinTree();

      Int_t                GenerateEfficiency(MargedReader*                  reader,
                                              BoardMap_t<ITdcDataProvider*>  providers,
                                              const BoardMap_t<std::string>& ifilenames,
                                              const std::string&             itreename);

      Int_t                GeneratePlots(MargedReader* reader,
                                         Bool_t        drawCoinTimeline = false,
                                         Int_t         mscountSelection = -1);

    private:
      void                 ClearLastSpill(Bool_t clearHists);

      void                 DrawTmpTimeline(Int_t bin, Int_t range);
    };

    TimelineCoincidence::TimelineCoincidence(ITdcDataProvider* provider)
      : fProvider(provider) {
      fEfficiencyTargetBh1 = false;
      fEfficiencyTargetBh2 = false;
      fEfficiencyTargetHod = false;
      fEfficiencyTargetExt = false;
      fEfficiencyTargetTc1 = false;
      fEfficiencyTargetTc2 = false;
      InitializeCoinData();
    }

    TimelineCoincidence::~TimelineCoincidence() {
    }

    void TimelineCoincidence::InitializePlots(const PlotsProfiles& profile) {
      std::cout << "Initialize plots" << std::endl;
      const std::string tdcName    = fProvider->GetName();
      const Double_t    timePerTdc = fProvider->GetTimePerTdc();

      const Double_t  xminInSpill = profile.TimeInSpill.Xmin / msec;
      const Double_t  xmaxInSpill = profile.TimeInSpill.Xmax / msec;
      const Int_t    xbinsInSpill = profile.TimeInSpill.NbinsX;

      const Double_t  xminInSync  = (Int_t)(profile.TimeInSync.Xmin     / timePerTdc) - 0.5;
      const Int_t    xbinsInSync  = (Int_t)(profile.TimeInSync.Xwidth() / timePerTdc) / profile.TimeInSync.BinWidth;
      const Double_t  xmaxInSync  = xminInSync + xbinsInSync * profile.TimeInSync.BinWidth;

      const Int_t    xbinsTmpTimeline = 1.2 * profile.MrSyncInterval.Mean / fProvider->GetTimePerTdc();
      const Double_t  xminTmpTimeline = 0;
      const Double_t  xmaxTmpTimeline = xbinsTmpTimeline;

      // Canvas
      fTimelineCanvas = new TCanvas();
      fHitMapCanvas   = new TCanvas();

      // Hit map
      fExtHitMap     = ExtinctionDetector::CreateHitMap("fExtHitMap");
      fExtBorderLine = ExtinctionDetector::CreateBorderLine();

      // Timeline
      fBh1Timeline.Resize(xbinsTmpTimeline);
      fBh2Timeline.Resize(xbinsTmpTimeline);
      fTc1Timeline.Resize(xbinsTmpTimeline);
      fTc2Timeline.Resize(xbinsTmpTimeline);
      fExtTimeline.resize(ExtinctionDetector::NofChannels);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        fExtTimeline[ch].Resize(xbinsTmpTimeline);
      }

      // Timeline
      hTmpBh1Timeline  = new TH1I("hTmpBh1Timeline",
                                  Form("%s, Coincidence Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpBh2Timeline  = new TH1I("hTmpBh2Timeline",
                                  Form("%s, Coincidence Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpTc1Timeline  = new TH1I("hTmpTc1Timeline",
                                  Form("%s, Coincidence Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpTc2Timeline  = new TH1I("hTmpTc2Timeline",
                                  Form("%s, Coincidence Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpExtTimeline = std::vector<TH1I*>(ExtinctionDetector::NofChannels);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        hTmpExtTimeline[ch] = new TH1I(Form("hTmpExtTimeline_%03lu", ch),
                                       Form("%s, Coincidence Timeline;"
                                            "TDC [count]", tdcName.data()),
                                       xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      }
      hTmpExtTimeline_Any = new TH1I("hTmpExtORTimeline_Any",
                                     Form("%s, Coincidence Timeline;"
                                          "TDC [count]", tdcName.data()),
                                     xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpSumTimeline  = new TH1I("hTmpSumTimeline",
                                  Form("%s, Coincidence Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);

      auto setTimelineStyle =
        [&] (TH1I* timeline, Int_t width, Int_t color) {
          timeline->SetStats(false);
          timeline->SetLineColor(color);
          timeline->SetLineWidth(width);
          timeline->SetMinimum(0);
          timeline->SetMaximum(5);
        };

      setTimelineStyle(hTmpBh1Timeline    , 5, 51);
   // setTimelineStyle(hTmpBh2Timeline    , 5, 60);
   // setTimelineStyle(hTmpHodTimeline    , 5,  1);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        setTimelineStyle(hTmpExtTimeline[ch], 4, 70);
      }
      setTimelineStyle(hTmpExtTimeline_Any, 4, 70);
      setTimelineStyle(hTmpTc1Timeline    , 3, 91);
      setTimelineStyle(hTmpTc2Timeline    , 2, 99);
      setTimelineStyle(hTmpSumTimeline    , 1,  1);

      hTmpTimelineLegend = new TLegend(0.80, 0.90 - 0.06 * 5, 0.90, 0.90);
      hTmpTimelineLegend->AddEntry(hTmpBh1Timeline    , "Bh1&2", "LF");
   // hTmpTimelineLegend->AddEntry(hTmpBh2Timeline    , "-"    , "LF");
      hTmpTimelineLegend->AddEntry(hTmpExtTimeline_Any, "Ext"  , "LF");
      hTmpTimelineLegend->AddEntry(hTmpTc1Timeline    , "Tc1&2", "LF");
      hTmpTimelineLegend->AddEntry(hTmpTc2Timeline    , "Tc3"  , "LF");
      hTmpTimelineLegend->AddEntry(hTmpSumTimeline    , "Sum"  , "LF");

      // Coincidence TDC in sync (timeline)
      hCoinTlTdcInSync = new TH1D("hCoinTlTdcInSync",
                                  Form("%s, Timeline TDC in MR Sync (Coincidence);"
                                       "TDC [count]", tdcName.data()),
                                  xbinsInSync, xminInSync, xmaxInSync);

      // Coincidence Mountain Plot (timeline)
      hCoinTlMountain = new TH2D("hCoinTlMountain",
                                 Form("%s, Timeline Mountain Plot (Coincidence);"
                                      "TDC [count];"
                                      "Time [ms]", tdcName.data()),
                                 xbinsInSync / 2, xminInSync, xmaxInSync,
                                 xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hCoinTlMountain->SetStats(false);

      // Coincidence TDC in sync (extinction detector)
      hCoinExtTdcInSync = new TH1D("hCoinExtTdcInSync",
                                   Form("%s, Extinction Detector TDC in MR Sync (Coincidence);"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInSync, xminInSync, xmaxInSync);

      // Coincidence Mountain Plot (extinction detector)
      hCoinExtMountain = new TH2D("hCoinExtMountain",
                                  Form("%s, Extinction Detector Mountain Plot (Coincidence);"
                                       "TDC [count];"
                                       "Time [ms]", tdcName.data()),
                                  xbinsInSync / 2, xminInSync, xmaxInSync,
                                  xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hCoinExtMountain->SetStats(false);

      hEfficiencyFrame = new TH1C("hEfficiencyFrame", "Efficiency", 6, 0, 6);
      hEfficiencyFrame->SetStats(false);
      hEfficiencyFrame->SetMinimum(0.00);
      hEfficiencyFrame->SetMaximum(1.05);
      hEfficiencyFrame->GetXaxis()->SetNdivisions(6, 1, 1);
      hEfficiencyFrame->GetYaxis()->SetNdivisions(510);
      hEfficiencyFrame->GetXaxis()->SetBinLabel(1, "BH1&2");
      hEfficiencyFrame->GetXaxis()->SetBinLabel(2, "");
      hEfficiencyFrame->GetXaxis()->SetBinLabel(3, "");
      hEfficiencyFrame->GetXaxis()->SetBinLabel(4, "Ext");
      hEfficiencyFrame->GetXaxis()->SetBinLabel(5, "TC1&2");
      hEfficiencyFrame->GetXaxis()->SetBinLabel(6, "TC3");

      hEfficiency = new TEfficiency("hEfficiency",
                                    "Efficiency",
                                    6, 0, 6);
    }

    void TimelineCoincidence::InitializeCoinTree(const std::string& filename, const std::string& treename) {
      std::cout << "Initialize coincidence tree" << std::endl;

      fCoinFile = new TFile(filename.data(), "RECREATE");
      if (!fCoinFile->IsOpen()) {
        std::cout << "[error] coincidence file is not opened, " << filename << std::endl;
        return;
      }

      fCoinTree = new TTree(treename.data(), "Coincidence resuld of residuals");

      fCoinTree->Branch("date"    , &fDate             , "date"           "/I");
      fCoinTree->Branch("emcount" , &fEMCount          , "emcount"        "/I");
      fCoinTree->Branch("ctarget" , &fCoinTarget       , "ctarget"        "/I");
      fCoinTree->Branch("mscount" , &fMrSyncCount      , "mscount"        "/I");
      fCoinTree->Branch("tdc"     , &fTdc              , "tdc"            "/L");
      fCoinTree->Branch("dtdc"    , &fTdcFromMrSync    , "dtdc"           "/L");
      fCoinTree->Branch("bwidth"  , &fBlurWidth        , "bwidth"         "/I");
      fCoinTree->Branch("cwidth"  , &fCoinWidth        , "cwidth"         "/I");
      fCoinTree->Branch("nofhits" , &fNofHits          , "nofhits"        "/I");
      fCoinTree->Branch("hitchs"  ,  fHitChannels      , "hitchs""[nofhits]/I");
   // fCoinTree->Branch("mscounts",  fHitMrSyncCounts  , "mscounts[nofhits]/I");
      fCoinTree->Branch("dtdcs"   ,  fHitTdcFromMrSyncs, "dtdcs" "[nofhits]/L");
      fCoinTree->Branch("tots"    ,  fHitTots          , "tots"  "[nofhits]/I");

      fCoinTree->AutoSave();
    }

    void TimelineCoincidence::DrawPlots(const std::string& ofilename) {
      std::cout << "Draw plots" << std::endl;
      if (!gPad) {
        TCanvas::MakeDefCanvas();
      }
      gPad->SetGrid(true, true);
      gPad->SetLogy(false);
      gPad->Print((ofilename + "[").data());

      Tron::ScopeSubstituter<Int_t> ss { gErrorIgnoreLevel, kWarning };

      std::cout << "hCoinTlTdcInSync" << std::endl;
      gPad->SetLogy(true);
      {
        hCoinTlTdcInSync->Draw("hist");
        hCoinTlTdcInSync->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      std::cout << "hCoinTlMountain" << std::endl;
      gPad->SetGrid(false, true);
      {
        hCoinTlMountain->Draw("colz");
        hCoinTlMountain->SetMinimum(0);
        gPad->Print(ofilename.data());
      }
      gPad->SetGrid(true, true);

      std::cout << "hCoinExtTdcInSync" << std::endl;
      gPad->SetLogy(true);
      {
        hCoinExtTdcInSync->Draw("hist");
        hCoinExtTdcInSync->SetMinimum(0.2);
        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      std::cout << "hCoinExtMountain" << std::endl;
      gPad->SetGrid(false, true);
      {
        hCoinExtMountain->Draw("colz");
        hCoinExtMountain->SetMinimum(0);
        gPad->Print(ofilename.data());
      }
      gPad->SetGrid(true, true);

      std::cout << "hEfficiency" << std::endl;
      if (hEfficiency->GetTotalHistogram()->GetEntries()) {
        hEfficiencyFrame->Draw();
        hEfficiency->Draw("Psame");
        gPad->Print(ofilename.data());
      }
      
      gPad->Print((ofilename + "]").data());
    }

    void TimelineCoincidence::WritePlots(const std::string& ofilename) {
      std::cout << "Write plots" << std::endl;
      TFile* file = new TFile(ofilename.data(), "RECREATE");
      if (!file->IsOpen()) {
        std::cout << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      hCoinTlTdcInSync ->Write();
      hCoinTlMountain  ->Write();
      hCoinExtTdcInSync->Write();
      hCoinExtMountain ->Write();

      hEfficiencyFrame ->Write();
      hEfficiency      ->Write();

      file->Close();
    }

    void TimelineCoincidence::WriteCoinTree() {
      if (!fCoinFile || !fCoinTree) {
        std::cout << "[error] output coincidence file is not initialized" << std::endl;
        return;
      } else if (!fCoinFile->IsOpen()) {
        std::cout << "[error] output coincidence file is not opened, " << fCoinFile->GetName() << std::endl;
        return;
      }

      fCoinFile->cd();
      fCoinTree->Write();
      fCoinFile->Close();
    }

    void TimelineCoincidence::ClearLastSpill(Bool_t clearHists) {
      if (clearHists) {
        hCoinTlTdcInSync ->Reset();
        hCoinTlMountain  ->Reset();
        hCoinExtTdcInSync->Reset();
        hCoinExtMountain ->Reset();
      }
    }

    Int_t TimelineCoincidence::GenerateEfficiency(MargedReader*                  reader,
                                                  BoardMap_t<ITdcDataProvider*>  providers,
                                                  const BoardMap_t<std::string>& ifilenames,
                                                  const std::string&             itreename) {
      for (std::size_t i = 0; i < 7; ++i) {
        fEfficiencyTargetBh1 = (i == 0);
        fEfficiencyTargetBh2 = (i == 1);
        fEfficiencyTargetHod = (i == 2);
        fEfficiencyTargetExt = (i == 3);
        fEfficiencyTargetTc1 = (i == 4);
        fEfficiencyTargetTc2 = (i == 5);

        if ((fCoincidenceTargetBh1 && fEfficiencyTargetBh1) ||
            (fCoincidenceTargetBh2 && fEfficiencyTargetBh2) ||
            (fCoincidenceTargetHod && fEfficiencyTargetHod) ||
            (fCoincidenceTargetExt && fEfficiencyTargetExt) ||
            (fCoincidenceTargetTc1 && fEfficiencyTargetTc1) ||
            (fCoincidenceTargetTc2 && fEfficiencyTargetTc2) ||
            (i == 6)) {
          reader->Open(providers, ifilenames, itreename);
          GeneratePlots(reader);
          reader->Close();
        }
      }

      return 0;
    }

    Int_t TimelineCoincidence::GeneratePlots(MargedReader* reader,
                                             Bool_t        drawCoinTimeline,
                                             Int_t         mscountSelection) {
      const clock_t startClock = clock();

      if (mscountSelection < 0) {
        mscountSelection = std::numeric_limits<Int_t>::max();
      }

      std::cout << "Initialize history" << std::endl;
      ClearLastSpill(true);

      // std::cout << "[debug] check end of spill" << std::endl;
      // Throw away first mr sync
      SortedTdcData_t tdcDataInMrSync;
      reader->Read(tdcDataInMrSync);
      tdcDataInMrSync.clear();

      fDate    = reader->GetDate();
      fEMCount = reader->GetEMCount();

      {
        Double_t mrSyncTime = 0.0;

        fMrSyncCount = 0;
        fCoinTarget  = 
          fCoincidenceTargetBh1 * 100000 +
          fCoincidenceTargetBh2 *  10000 +
          fCoincidenceTargetHod *   1000 +
          fCoincidenceTargetExt *    100 +
          fCoincidenceTargetTc1 *     10 +
          fCoincidenceTargetTc2 *      1;

        Int_t dtdc = 0;
        const Int_t coinTdcWidth = fCoinTimeWidth / fProvider->GetTimePerTdc();

        while (true) {
          for (; reader->Read(tdcDataInMrSync); tdcDataInMrSync.clear()) {
            fMrSyncCount = reader->GetMrSyncCount();

            // std::cout << "[debug] data process" << std::endl;
            for (auto&& pair : tdcDataInMrSync) {
              // auto& tag  = pair.first;
              auto& data = pair.second;

              if (BeamlineHodoscope::Contains(data.Channel)) {
                // std::cout << "[debug] fill bh timeline" << std::endl;
                const std::size_t bhCh = BeamlineHodoscope::GetChannel(data.Channel);
                if (bhCh == 0) {
                  fBh1Timeline.FillHit(data, coinTdcWidth);
                } else {
                  fBh2Timeline.FillHit(data, coinTdcWidth);
                }

              } else if (ExtinctionDetector::Contains(data.Channel)) {
                // std::cout << "[debug] fill ext timeline" << std::endl;
                const std::size_t extCh = ExtinctionDetector::GetChannel(data.Channel);
                fExtTimeline[extCh].FillHit(data, coinTdcWidth);

              } else if (TimingCounter::Contains(data.Channel)) {
                // std::cout << "[debug] fill tc timeline" << std::endl;
                const std::size_t tcCh = TimingCounter::GetChannel(data.Channel);
                if (tcCh == 0) {
                  fTc1Timeline.FillHit(data, coinTdcWidth);
                } else {
                  fTc2Timeline.FillHit(data, coinTdcWidth);
                }

              } else if (MrSync::Contains(data.Channel)) {
                // std::cout << "[debug] get mrsync" << std::endl;
                mrSyncTime = data.Time;

              } else {
                // std::cout << "[debug] skip others" << std::endl;
                continue;

              }
            }

            {
              static const Int_t xmax = fBh1Timeline.Size();

              static Bool_t hitBh1 = false, hitBh2 = false, hitTc1 = false, hitTc2 = false, hitExt = false;
              auto isCoincident =
                [&] (Int_t i) {
                  {
                    hitBh1 = fBh1Timeline.HasHit(i);
                    hitBh2 = fBh2Timeline.HasHit(i);
                    hitExt = std::any_of(fExtTimeline.begin(), fExtTimeline.end(), [&](Timeline_t& tl) { return tl.HasHit(i); });
                    hitTc1 = fTc1Timeline.HasHit(i);
                    hitTc2 = fTc2Timeline.HasHit(i);
                  }

                  return
                    (hitBh1 || !fCoincidenceTargetBh1 || fEfficiencyTargetBh1) &&
                    (hitBh2 || !fCoincidenceTargetBh2 || fEfficiencyTargetBh2) &&
                 // (hitHod || !fCoincidenceTargetHod || fEfficiencyTargetHod) &&
                    (hitExt || !fCoincidenceTargetExt || fEfficiencyTargetExt) &&
                    (hitTc1 || !fCoincidenceTargetTc1 || fEfficiencyTargetTc1) &&
                    (hitTc2 || !fCoincidenceTargetTc2 || fEfficiencyTargetTc2);
                };

              // std::cout << "[debug] seek timeline" << std::endl;
              for (dtdc = 0; dtdc < xmax; ++dtdc) {

                if (isCoincident(dtdc)) {
                  ClearCoinEvent();

                  fTdcFromMrSync = dtdc;
                  hCoinTlTdcInSync->Fill(dtdc);
                  hCoinTlMountain ->Fill(dtdc, mrSyncTime / msec);

                  if (drawCoinTimeline) {
                    DrawTmpTimeline(dtdc, 2 * coinTdcWidth);
                    gPad->WaitPrimitive();
                  }

                  // std::cout << "[debug] skip dead time" << std::endl;
                  // Skip dead time
                  fCoinWidth = 0;
                  Bool_t filled = false;
                  std::map<std::pair<Int_t, Long64_t>, Long64_t> hitmscounts;
                  std::map<std::pair<Int_t, Long64_t>, Long64_t> hitdtdcs;
                  std::map<std::pair<Int_t, Long64_t>, Long64_t> hittots;
                  for (; dtdc < xmax; ++dtdc) {
                    if (isCoincident(dtdc)) {
                      ++fCoinWidth;

                      // std::cout << "[debug] fill coincidence data" << std::endl;
                      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
                        auto& tl = fExtTimeline[ch];
                        if (tl.HasHit(dtdc)) {
                          std::pair<Int_t, Long64_t> tag { ch, tl.fTdcs[dtdc] };
                          hitmscounts.emplace(tag, tl.fMrSyncCounts  [dtdc]);
                          hitdtdcs   .emplace(tag, tl.fTdcFromMrSyncs[dtdc]);
                          hittots    .emplace(tag, tl.fTots          [dtdc]);
                        }
                      }

                      if (fEfficiencyTargetBh1 && hitBh1 && !filled) { hEfficiency->Fill(true, 0); filled = true; }
                      if (fEfficiencyTargetBh2 && hitBh2 && !filled) { hEfficiency->Fill(true, 1); filled = true; }
                   // if (fEfficiencyTargetHod && hitHod && !filled) { hEfficiency->Fill(true, 2); filled = true; }
                      if (fEfficiencyTargetExt && hitExt && !filled) { hEfficiency->Fill(true, 3); filled = true; }
                      if (fEfficiencyTargetTc1 && hitTc1 && !filled) { hEfficiency->Fill(true, 4); filled = true; }
                      if (fEfficiencyTargetTc2 && hitTc2 && !filled) { hEfficiency->Fill(true, 5); filled = true; }
                    } else {
                      if (fEfficiencyTargetBh1 && !filled) { hEfficiency->Fill(false, 0); }
                      if (fEfficiencyTargetBh2 && !filled) { hEfficiency->Fill(false, 1); }
                   // if (fEfficiencyTargetHod && !filled) { hEfficiency->Fill(false, 2); }
                      if (fEfficiencyTargetExt && !filled) { hEfficiency->Fill(false, 3); }
                      if (fEfficiencyTargetTc1 && !filled) { hEfficiency->Fill(false, 4); }
                      if (fEfficiencyTargetTc2 && !filled) { hEfficiency->Fill(false, 5); }
                      break;
                    }
                  }

                  for (auto&& pair : hitmscounts) {
                    auto& tag = pair.first;
                    auto& ch  = tag.first;

                    hCoinExtTdcInSync->Fill(hitdtdcs[tag]);
                    hCoinExtMountain ->Fill(hitdtdcs[tag], mrSyncTime / msec);

                    fHitChannels      [fNofHits] = ch;
                    fHitMrSyncCounts  [fNofHits] = pair.second;
                    fHitTdcFromMrSyncs[fNofHits] = hitdtdcs[tag];
                    fHitTots          [fNofHits] = hittots [tag];
                    fNofHits++;

                    if (fNofHits == ExtinctionDetector::NofChannels) {
                      std::cerr << "[warning] too many multiple hit, mr sync count = " << fMrSyncCount << std::endl;
                      break;
                    }
                  }

                  fBlurWidth = coinTdcWidth - fCoinWidth + 1;
                  if (fCoinTree) {
                    fCoinTree->Fill();
                  }
                }
              }

              if (mscountSelection <= fMrSyncCount) {
                DrawTmpTimeline(0, 0);
                gPad->WaitPrimitive();
              }
            }

            // std::cout << "[debug] Shift timeline" << std::endl;
            // Shift timeline
            {
              fBh1Timeline.Clear();
              fBh2Timeline.Clear();
              fTc1Timeline.Clear();
              fTc2Timeline.Clear();
            }
            for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
              fExtTimeline[ch].Clear();
            }
          }

          // std::cout << "[debug] check end of spill" << std::endl;
          // Check end of spill
          if (reader->IsSpillEnded()) {
            reader->ClearLastSpill();
            // std::cout << "[debug] Throw away first mr sync" << std::endl;
            // Throw away first mr sync
            reader->Read(tdcDataInMrSync);
            tdcDataInMrSync.clear();

            fDate        = reader->GetDate();
            fEMCount     = reader->GetEMCount();
          }

          // std::cout << "[debug] check end of file" << std::endl;
          // Check end of file
          if (reader->IsFileEnded()) {
            break;
          }
        }
      }

      const clock_t stopClock = clock();
      std::cout << "time: " << (double)(stopClock - startClock) / CLOCKS_PER_SEC << " sec\n";

      return 0;
    }

    void TimelineCoincidence::DrawTmpTimeline(Int_t dtdc, Int_t range) {
      std::cout << "DrawTmpTimeline: MR Sync Count = " << fMrSyncCount << std::endl;

      // Timeline
      fTimelineCanvas->cd();

      hTmpBh1Timeline    ->Reset();
   // hTmpBh2Timeline    ->Reset();
      hTmpTc1Timeline    ->Reset();
      hTmpTc2Timeline    ->Reset();
      hTmpExtTimeline_Any->Reset();
      hTmpSumTimeline    ->Reset();

      auto setTimeline =
        [&] (Timeline_t& tl, TH1I* tmpTl) {
          const Int_t nbins = tmpTl->GetNbinsX();
          Int_t* px_any = tmpTl->GetArray();
          for (Int_t bin = 1; bin < nbins; ++bin) {
            px_any[bin] = (px_any[bin] || tl.HasHit(bin - 1));
          }
        };

      setTimeline(fBh1Timeline, hTmpBh1Timeline);
   // setTimeline(fBh2Timeline, hTmpBh2Timeline);
      setTimeline(fTc1Timeline, hTmpTc1Timeline);
      setTimeline(fTc2Timeline, hTmpTc2Timeline);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        setTimeline(fExtTimeline[ch], hTmpExtTimeline_Any);
      }

      hTmpSumTimeline->Add(hTmpBh1Timeline    );
   // hTmpSumTimeline->Add(hTmpBh2Timeline    );
   // hTmpSumTimeline->Add(hTmpHodTimeline    );
      hTmpSumTimeline->Add(hTmpExtTimeline_Any);
      hTmpSumTimeline->Add(hTmpTc1Timeline    );
      hTmpSumTimeline->Add(hTmpTc2Timeline    );

      hTmpBh1Timeline    ->Draw();
   // hTmpBh2Timeline    ->Draw("same");
   // hTmpHodTimeline    ->Draw("same");
      hTmpExtTimeline_Any->Draw("same");
      hTmpTc1Timeline    ->Draw("same");
      hTmpTc2Timeline    ->Draw("same");

      hTmpSumTimeline    ->Draw("same");

      if (dtdc || range) {
        hTmpBh1Timeline->GetXaxis()->SetRangeUser(dtdc - range + 1, dtdc + range + 1);
      } else {
        hTmpBh1Timeline->GetXaxis()->SetRangeUser(0, hTmpBh1Timeline->GetNbinsX() / 1.2);
      }

      hTmpTimelineLegend->Draw();

      fTimelineCanvas->Modified();
      fTimelineCanvas->Update();

      // Hitmap
      fHitMapCanvas->cd();

      fExtHitMap->Reset();
      if (dtdc) {
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          if (fExtTimeline[ch].HasHit(dtdc)) {
            ExtinctionDetector::Fill(fExtHitMap, ch, false);
          }
        }
      } else {
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          for (std::size_t i = 0, n = fExtTimeline[ch].Size(); i < n; ++i) {
            if (fExtTimeline[ch].HasHit(i)) {
              ExtinctionDetector::Fill(fExtHitMap, ch, false);
            }
          }
        }
      }

      fExtHitMap->SetMaximum(1.0);
      fExtHitMap    ->Draw("col");
      fExtBorderLine->Draw();

      fHitMapCanvas->Modified();
      fHitMapCanvas->Update();
    }
    
  }

}

#endif
