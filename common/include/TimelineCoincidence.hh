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
#include "TParameter.h"
#include "TEfficiency.h"

#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "Spill.hh"

#include "Linq.hh"
#include "String.hh"
#include "ObjectHelper.hh"
#include "ScopeSubstituter.hh"
#include "MargedReader.hh"

namespace Extinction {

  namespace Analyzer {

    class Timeline_t {
    private:
      std::size_t           fSize;
      Long64_t*             pTdcs = nullptr;

    public:
      Int_t                 fChannel;
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
        const Int_t last = std::min((std::size_t)(x + width), fSize);
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
      ITdcDataProvider*          fProvider               = nullptr;

      Timeline_t                 fBh1Timeline;
      Timeline_t                 fBh2Timeline;
      Timeline_t                 fTc1Timeline;
      Timeline_t                 fTc2Timeline;
      std::vector<Timeline_t>    fExtTimeline;

      TH1I*                      hTmpBh1Timeline         = nullptr;
      TH1I*                      hTmpBh2Timeline         = nullptr;
      TH1I*                      hTmpTc1Timeline         = nullptr;
      TH1I*                      hTmpTc2Timeline         = nullptr;
      std::vector<TH1I*>         hTmpExtTimeline;
      TH1I*                      hTmpExtTimeline_Any     = nullptr;
      TH1I*                      hTmpSumTimeline         = nullptr;
      TLegend*                   hTmpTimelineLegend      = nullptr;

      TH1D*                      hCoinTlTdcInSync        = nullptr;
      TH2D*                      hCoinTlMountain         = nullptr;
      TH1D*                      hCoinExtTdcInSync       = nullptr;
      TH2D*                      hCoinExtMountain        = nullptr;

      TH1*                       hEfficiencyFrame        = nullptr;
      TEfficiency*               hEfficiency             = nullptr;

      TFile*                     fCoinFile               = nullptr;
      TTree*                     fCoinTree               = nullptr;

      TFile*                     fSpillFile              = nullptr;
      TTree*                     fSpillTree              = nullptr;

      TCanvas*                   fTimelineCanvas         = nullptr;
      TCanvas*                   fHitMapCanvas           = nullptr;

      TH2*                       fExtHitMap              = nullptr;
      TList*                     fExtBorderLine          = nullptr;

      Bool_t                     fCoincidenceTargetBh1;
      Bool_t                     fCoincidenceTargetBh2;
      Bool_t                     fCoincidenceTargetHod;
      Bool_t                     fCoincidenceTargetExt;
      Bool_t                     fCoincidenceTargetTc1;
      Bool_t                     fCoincidenceTargetTc2;

      Bool_t                     fEfficiencyTargetBh1;
      Bool_t                     fEfficiencyTargetBh2;
      Bool_t                     fEfficiencyTargetHod;
      Bool_t                     fEfficiencyTargetExt;
      Bool_t                     fEfficiencyTargetTc1;
      Bool_t                     fEfficiencyTargetTc2;

      Long64_t                   fCoinTdcWidth;

      // Coincidence data
      ULong64_t                  fDate;
      Int_t                      fEMCount;
      UInt_t                     fCoinTarget;
      Int_t                      fMrSyncCount;
      Long64_t                   fTdcFromMrSync;
      Int_t                      fBlurWidth;
      Int_t                      fCoinWidth;
      Int_t                      fNofHits;
      Int_t                      fHitChannels      [Detectors::NofChannels];
   // Int_t                      fHitMrSyncCounts  [Detectors::NofChannels];
      Long64_t                   fHitTdcFromMrSyncs[Detectors::NofChannels];
      Int_t                      fHitTots          [Detectors::NofChannels];

      Int_t                      fBunchEdgeMargin;
      CoinSpillData              fSpillData;

      TGraph*                    fBunchRange             = nullptr;

      void InitializeCoinData() {
        fDate        = 0;
        fEMCount     = 0;
        fCoinTarget  = 0;
        fMrSyncCount = 0;
        ClearCoinEvent();
      }

      void ClearCoinEvent() {
        fTdcFromMrSync = 0;
        fBlurWidth     = 0;
        fCoinWidth     = 0;
        fNofHits       = 0;
        std::memset(fHitChannels      , 0, sizeof(fHitChannels      ));
     // std::memset(fHitMrSyncCounts  , 0, sizeof(fHitMrSyncCounts  ));
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
        std::cout << "SetCoinTimeWidth ... " << width / nsec << " nsec";
        fCoinTdcWidth = TMath::Nint(width / fProvider->GetTimePerTdc());
        std::cout << " (" << fCoinTdcWidth << " count)" << std::endl;
      }
      inline Long64_t      GetCoinTdcWidth() const { return fCoinTdcWidth; }

      inline void          SetBunchEdgeMargin(Double_t margin) {
        std::cout << "SetBunchEdgeMargin ... " << margin / nsec << " nsec";
        fBunchEdgeMargin = TMath::Nint(margin / fProvider->GetTimePerTdc());
        std::cout << " (" << fBunchEdgeMargin << " count)" << std::endl;
      }
      inline Long64_t      GetBunchEdgeMargin() const { return fBunchEdgeMargin; }

      Int_t                ReadPlots(const std::string& ifilename);
      void                 InitializePlots(const PlotsProfiles& profile);
      void                 InitializeCoinTree(const std::string& filename, const std::string& treename = "ctree");
      void                 InitializeSpillSummary(const std::string& filename, const std::string& treename = "cspill");
      void                 DrawPlots(const std::string& ofilename);
      void                 WritePlots(const std::string& ofilename);
      void                 WriteCoinTree();
      void                 WriteSpillSummary();
      void                 WriteBunchProfile(const std::string& ofilename);

      Int_t                GenerateEfficiency(MargedReader*                  reader,
                                              BoardMap_t<ITdcDataProvider*>  providers,
                                              const BoardMap_t<std::string>& ifilenames,
                                              const std::string&             itreename);

      Int_t                GeneratePlots(MargedReader* reader,
                                         Bool_t        drawCoinTimeline = false,
                                         Int_t         mscountSelection = -1);

      void                 CalcBunchProfile();

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

      fBunchRange = new TGraph();
      fBunchRange->SetFillStyle(1001);
      fBunchRange->SetFillColorAlpha(kYellow, 0.5);
    }

    TimelineCoincidence::~TimelineCoincidence() {
    }

    Int_t TimelineCoincidence::ReadPlots(const std::string& ifilename) {
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
        hCoinTlTdcInSync  = dynamic_cast<TH1D*>(getFromFile("hCoinTlTdcInSync" ));
        hCoinTlMountain   = dynamic_cast<TH2D*>(getFromFile("hCoinTlMountain"  ));
        hCoinExtTdcInSync = dynamic_cast<TH1D*>(getFromFile("hCoinExtTdcInSync"));
        hCoinExtMountain  = dynamic_cast<TH2D*>(getFromFile("hCoinExtMountain" ));
        hCoinTlTdcInSync ->SetDirectory(nullptr);
        hCoinTlMountain  ->SetDirectory(nullptr);
        hCoinExtTdcInSync->SetDirectory(nullptr);
        hCoinExtMountain ->SetDirectory(nullptr);
      }

      {
        hEfficiencyFrame = dynamic_cast<TH1*>(getFromFile("hEfficiencyFrame"));
        hEfficiencyFrame->SetDirectory(nullptr);

        hEfficiency = dynamic_cast<TEfficiency*>(getFromFile("hEfficiency"));
        hEfficiency->SetDirectory(nullptr);
      }

      file->Close();

      return 0;
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
      fBh1Timeline.Resize(xbinsTmpTimeline); fBh1Timeline.fChannel = 0 + BeamlineHodoscope::GlobalChannelOffset;
      fBh2Timeline.Resize(xbinsTmpTimeline); fBh2Timeline.fChannel = 1 + BeamlineHodoscope::GlobalChannelOffset;
      fTc1Timeline.Resize(xbinsTmpTimeline); fTc1Timeline.fChannel = 0 + TimingCounter    ::GlobalChannelOffset;
      fTc2Timeline.Resize(xbinsTmpTimeline); fTc2Timeline.fChannel = 1 + TimingCounter    ::GlobalChannelOffset;
      fExtTimeline.resize(ExtinctionDetector::NofChannels);
      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
        fExtTimeline[ch].Resize(xbinsTmpTimeline);
        fExtTimeline[ch].fChannel = ch + ExtinctionDetector::GlobalChannelOffset;
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

      fCoinTree = new TTree(treename.data(), "Coincidence results of residuals");

      fCoinTree->Branch("date"    , &fDate             , "date"           "/l");
      fCoinTree->Branch("emcount" , &fEMCount          , "emcount"        "/I");
      fCoinTree->Branch("ctarget" , &fCoinTarget       , "ctarget"        "/i");
      fCoinTree->Branch("mscount" , &fMrSyncCount      , "mscount"        "/I");
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

    void TimelineCoincidence::InitializeSpillSummary(const std::string& filename, const std::string& treename) {
      std::cout << "Initialize spill tree" << std::endl;

      fSpillFile = new TFile(filename.data(), "RECREATE");
      if (!fSpillFile->IsOpen()) {
        std::cout << "[error] spill file is not opened, " << filename << std::endl;
        return;
      }

      fSpillTree = new TTree(treename.data(), "Spill summary");
      fSpillData.CreateBranch(fSpillTree);
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

        const Double_t miny = TMath::Power(10., gPad->GetUymin());
        const Double_t maxy = TMath::Power(10., gPad->GetUymax());
        fBunchRange->Set(0);
        fBunchRange->SetMinimum(miny);
        fBunchRange->SetMaximum(maxy);
        for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
          if (fSpillData.BunchMinEdges[bunch] && fSpillData.BunchMaxEdges[bunch]) {
            fBunchRange->SetPoint(fBunchRange->GetN(), fSpillData.BunchMinEdges[bunch], miny);
            fBunchRange->SetPoint(fBunchRange->GetN(), fSpillData.BunchMinEdges[bunch], maxy);
            fBunchRange->SetPoint(fBunchRange->GetN(), fSpillData.BunchMaxEdges[bunch], maxy);
            fBunchRange->SetPoint(fBunchRange->GetN(), fSpillData.BunchMaxEdges[bunch], miny);
          }
        }

        if (fBunchRange->GetN()) {
          fBunchRange->Draw("AF");
          gPad->Update();
          fBunchRange->GetXaxis()->SetLimits(hCoinTlTdcInSync->GetXaxis()->GetXmin(), hCoinTlTdcInSync->GetXaxis()->GetXmax());
          hCoinTlTdcInSync->Draw("histsame");
          gPad->Print(ofilename.data());
        }
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

      {
        Tron::ObjectHelper::WriteValue<Long64_t>(fSpillData.Date.Convert(), "Time");
      }

      {
        hCoinTlTdcInSync ->Write();
        hCoinTlMountain  ->Write();
        hCoinExtTdcInSync->Write();
        hCoinExtMountain ->Write();
      }

      {
        hEfficiencyFrame ->Write();
        hEfficiency      ->Write();
      }

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

    void TimelineCoincidence::WriteSpillSummary() {
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

    void TimelineCoincidence::WriteBunchProfile(const std::string &ofilename) {
      std::cout << "Write bunch profile" << std::endl;
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cerr << "[error] output file is not opened, " << ofilename << std::endl;
        return;
      }

      for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
        ofile << bunch << "\t"
              << Form("%23.15e", fSpillData.BunchCenters [bunch]) << "\t"
              << Form("%23.15e", fSpillData.BunchWidths  [bunch]) << "\t"
              << Form("%lld", fSpillData.BunchMinEdges[bunch]) << "\t"
              << Form("%lld", fSpillData.BunchMaxEdges[bunch]) << std::endl;
      }

      ofile.close();

      std::cerr << "Info in <TimelineCoincidence::WriteBunchProfile>: dat file " << ofilename << " has been created" << std::endl;
    }

    void TimelineCoincidence::ClearLastSpill(Bool_t clearHists) {
      fSpillData.Clear();
      
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
      fSpillData.SetDate(fDate);
      fSpillData.EMCount = fEMCount;

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
                  fBh1Timeline.FillHit(data, fCoinTdcWidth);
                } else {
                  fBh2Timeline.FillHit(data, fCoinTdcWidth);
                }

              } else if (ExtinctionDetector::Contains(data.Channel)) {
                // std::cout << "[debug] fill ext timeline" << std::endl;
                const std::size_t extCh = ExtinctionDetector::GetChannel(data.Channel);
                fExtTimeline[extCh].FillHit(data, fCoinTdcWidth);

              } else if (TimingCounter::Contains(data.Channel)) {
                // std::cout << "[debug] fill tc timeline" << std::endl;
                const std::size_t tcCh = TimingCounter::GetChannel(data.Channel);
                if (tcCh == 0) {
                  fTc1Timeline.FillHit(data, fCoinTdcWidth);
                } else {
                  fTc2Timeline.FillHit(data, fCoinTdcWidth);
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
                  ++fSpillData.CoinCount;

                  if (drawCoinTimeline) {
                    DrawTmpTimeline(dtdc, 2 * fCoinTdcWidth);
                    gPad->WaitPrimitive();
                  }

                  // std::cout << "[debug] skip dead time" << std::endl;
                  // Skip dead time
                  fCoinWidth = 0;
                  Bool_t filled = false;
                  std::map<std::pair<Long64_t/*tdc*/, Int_t/*ch*/>, Long64_t> hitmscounts;
                  std::map<std::pair<Long64_t/*tdc*/, Int_t/*ch*/>, Long64_t> hitdtdcs;
                  std::map<std::pair<Long64_t/*tdc*/, Int_t/*ch*/>, Long64_t> hittots;
                  auto fillCoincidenceData =
                    [&] (Timeline_t& tl) {
                      if (tl.HasHit(dtdc)) {
                        const std::pair<Int_t, Long64_t> tag { tl.fTdcFromMrSyncs[dtdc], tl.fChannel };
                        hitmscounts.emplace(tag, tl.fMrSyncCounts  [dtdc]);
                        hitdtdcs   .emplace(tag, tl.fTdcFromMrSyncs[dtdc]);
                        hittots    .emplace(tag, tl.fTots          [dtdc]);
                      }
                    };
                  for (; dtdc < xmax; ++dtdc) {
                    if (isCoincident(dtdc)) {
                      ++fCoinWidth;

                      // std::cout << "[debug] fill coincidence data" << std::endl;
                      for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
                        fillCoincidenceData(fExtTimeline[ch]);
                      }
                      fillCoincidenceData(fTc1Timeline);
                      fillCoincidenceData(fTc2Timeline);
                      fillCoincidenceData(fBh1Timeline);
                      fillCoincidenceData(fBh2Timeline);

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

                  fBlurWidth = fTdcFromMrSync - hitdtdcs.begin()->first.first;

                  for (auto&& pair : hitmscounts) {
                    auto& tag = pair.first;
                    auto& gch = tag.second;

                    if        (ExtinctionDetector::Contains(gch)) {
                      hCoinExtTdcInSync->Fill(hitdtdcs[tag]);
                      hCoinExtMountain ->Fill(hitdtdcs[tag], mrSyncTime / msec);
                      ++fSpillData.ExtEntries;
                    } else if (BeamlineHodoscope::Contains(gch)) {
                      const Int_t bhch = BeamlineHodoscope::GetChannel(gch);
                      if (bhch == 0) {
                        ++fSpillData.Bh1Entries;
                      } else {
                        ++fSpillData.Bh2Entries;
                      }
                    } else if (TimingCounter::Contains(gch)) {
                      const Int_t tcch = TimingCounter::GetChannel(gch);
                      if (tcch == 0) {
                        ++fSpillData.Tc1Entries;
                      } else {
                        ++fSpillData.Tc2Entries;
                      }
                    }
                  }

                  for (auto&& pair : hitmscounts) {
                    auto& tag = pair.first;
                    auto& ch  = tag.second;

                    fHitChannels      [fNofHits] = ch;
                 // fHitMrSyncCounts  [fNofHits] = pair.second;
                    fHitTdcFromMrSyncs[fNofHits] = hitdtdcs[tag];
                    fHitTots          [fNofHits] = hittots [tag];
                    fNofHits++;

                    if (fNofHits == Detectors::NofChannels) {
                      std::cerr << "[warning] coincidence of too many multiple hit" << std::endl
                                << "  # of hits    \t" << hitmscounts.size() << std::endl
                                << "  begin of dtdc\t" << dtdc - fCoinWidth  << std::endl
                                << "  end   of dtdc\t" << dtdc - 1           << std::endl
                                << "  coin width   \t" << fCoinWidth         << std::endl
                                << "  mr sync count\t" << fMrSyncCount       << std::endl;
                      Int_t ipair = 0;
                      for (auto&& pair2 : hitmscounts) {
                        std::cerr << (ipair++ == 0 ?  "  dtdc ch      \t" : "               \t") << std::setw(8) << pair2.first.first << "\t" << std::setw(3) << pair2.first.second << std::endl;
                      }
                      break;
                    }
                  }

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
            std::cout << "[info] end of spill" << std::endl;
            // Calc spill summary
            CalcBunchProfile();

            // Fill spill summary
            if (fSpillTree) {
              fSpillTree->Fill();
            }

            ClearLastSpill(false);

            reader->ClearLastSpill();
            // std::cout << "[debug] Throw away first mr sync" << std::endl;
            // Throw away first mr sync
            reader->Read(tdcDataInMrSync);
            tdcDataInMrSync.clear();

            fDate        = reader->GetDate();
            fEMCount     = reader->GetEMCount();
            fSpillData.SetDate(fDate);
            fSpillData.EMCount = fEMCount;
          }

          // std::cout << "[debug] check end of file" << std::endl;
          // Check end of file
          if (reader->IsFileEnded()) {
            std::cout << "[info] end of file" << std::endl;
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

      if (!fTimelineCanvas || !fHitMapCanvas) {
        std::cout << "[error] timelines is not initialized" << std::endl;
        return;
      }
      
      // Timeline
      fTimelineCanvas->cd();

      hTmpBh1Timeline    ->Reset();
   // hTmpBh2Timeline    ->Reset();
      hTmpTc1Timeline    ->Reset();
      hTmpTc2Timeline    ->Reset();
      hTmpExtTimeline_Any->Reset();
      hTmpSumTimeline    ->Reset();

      hTmpBh1Timeline->SetTitle(Form("%s, Coincidence Timeline, MS %d", fProvider->GetName().data(), fMrSyncCount));

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

    void TimelineCoincidence::CalcBunchProfile() {
      std::cout << "_____ Bunch Profile _____" << std::endl;
      {
        TH1D* hTdcInSync = hCoinTlTdcInSync;

        for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
          fSpillData.BunchCenters [bunch] = 0.0;
          fSpillData.BunchWidths  [bunch] = 0.0;
          fSpillData.BunchMinEdges[bunch] = 0.0;
          fSpillData.BunchMaxEdges[bunch] = 0.0;
        }

        Int_t xbin = 1, nbinsx = hTdcInSync->GetNbinsX();
        const Double_t ymax = hTdcInSync->GetBinContent(hTdcInSync->GetMaximumBin());
        for (std::size_t bunch = 0; bunch < kNofBunches; ++bunch) {
          Int_t xbin1 = xbin, xbin2 = xbin;
          // std::cout << "... search start edge" << std::endl;
          for (; xbin <= nbinsx && hTdcInSync->GetBinContent(xbin) < ymax * 0.10; ++xbin) {
            // std::cout << xbin << "\t" << hTdcInSync->GetBinContent(xbin) << "\t" << ymax * 0.10 << std::endl;
            xbin1 = xbin;
          }
          // std::cout << "... search end edge" << std::endl;
          for (; xbin <= nbinsx && hTdcInSync->GetBinContent(xbin) > ymax * 0.01; ++xbin) {
            // std::cout << xbin << "\t" << hTdcInSync->GetBinContent(xbin) << "\t" << ymax * 0.01 << std::endl;
            xbin2 = xbin;
          }
          // std::cout << "xbin1 = " << xbin1 << "\txbin2 = " << xbin2 << std::endl;

          if (xbin2 <= xbin1) {
            std::cout << "[warning] " << (bunch == 0 ? "1st" :
                                          bunch == 1 ? "2nd" :
                                          bunch == 2 ? "3rd" :
                                          Form("%ldth", bunch + 1)) << " bunch was not found" << std::endl;
            break;
          }

          // Get bunch center
          const Long64_t maxbin = TMath::LocMax(xbin2 - xbin1 + 1, hTdcInSync->GetArray() + xbin1);
          fSpillData.BunchCenters[bunch] = hTdcInSync->GetBinCenter(xbin1 + maxbin);

          // Get bunch min/max edge
          for (Int_t edge = xbin1; edge > 0; --edge) {
            if (hTdcInSync->GetBinContent(edge) == 0) {
              fSpillData.BunchMinEdges[bunch] = hTdcInSync->GetBinCenter(edge) - fBunchEdgeMargin;
              break;
            }
          }
          if (fSpillData.BunchMinEdges[bunch] == 0.0) {
            std::cout << "[warning] bunch min edge reach underflow" << std::endl;
          }

          for (Int_t& edge = xbin2; edge <= nbinsx; ++edge) {
            if (hTdcInSync->GetBinContent(edge) == 0) {
              fSpillData.BunchMaxEdges[bunch] = hTdcInSync->GetBinCenter(edge) + fBunchEdgeMargin;
              break;
            }
          }
          if (fSpillData.BunchMaxEdges[bunch] == 0.0) {
            std::cout << "[warning] bunch max edge reach overflow" << std::endl;
          }

          // Get bunch width
          if (fSpillData.BunchMinEdges[bunch] && fSpillData.BunchMaxEdges[bunch]) {
            fSpillData.BunchWidths[bunch] = fSpillData.BunchMaxEdges[bunch] - fSpillData.BunchMinEdges[bunch] + 1;
          }
        }
      }
    }

  }

}

#endif
