#ifndef Extinction_Fct_MonitorWindow_hh
#define Extinction_Fct_MonitorWindow_hh

#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

#include "Units.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "Fct.hh"
#include "AnaCoin.hh"

namespace Extinction {

  namespace Fct {

    class MonitorWindow {
    private:
      const std::size_t kHistLimit     = 10000;
      const Double_t    kLastThreshold = 100.0 * nsec;
      const Double_t    kCoinTimeWidth =  25.0 * nsec;
      const Long64_t    kSpillLimit    = 1 * 60 * 60 / 5;

    private:
      using CoinOffset = Analyzer::AnaTimeOffset::CoinOffset;

      FctData       fProvider;
      std::map<std::size_t, std::map<std::size_t, Analyzer::AnaCoin::CoinInfo>> fCoinInfo;
      std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Bool_t>>    fContains;

      TApplication* fApplication;

      TCanvas*      fCanvas;

      TPad*         fPadBh1Tdc;
      TPad*         fPadBh2Tdc;
      TPad*         fPadExtTdc;
      TPad*         fPadHodTdc;
      TPad*         fPadTc1Tdc;
      TPad*         fPadTc2Tdc;

      TPad*         fPadHodHitMap;
      TPad*         fPadHodHitCnt;
      TPad*         fPadExtHitMap;
      TPad*         fPadExtHitCnt;

      TPad*         fPadMountain;
      TPad*         fPadTdcSync;
      TPad*         fPadHit;

      TH1*          hExtTdcInSpill_Any;
      TH1*          hHodTdcInSpill_Any;
      TH1**         hTcTdcInSpill;
      TH1**         hBhTdcInSpill;
      TH2*          hExtHitMap;
      TList*        lExtBorderLine;
      // TH1*          hExtEntryByCh;
      TH1*          hExtEntryByChBottom;
      TH1*          hExtEntryByChCenter1;
      TH1*          hExtEntryByChCenter2;
      TH1*          hExtEntryByChTop;
      TH2*          hHodHitMap;
      TList*        lHodBorderLine;
      TH1*          hHodEntryByCh;
      TH2*          hExtMountain_Any;  
      TH1*          hExtTdcInSync_Any;
      TGraphErrors* gHitInSpill;

      Long64_t      fSpillCount = 0;
      Long64_t      fCoinCount  = 0;

      std::vector<TdcData> fLastExtData;
      std::vector<TdcData> fLastHodData;
      std::vector<TdcData> fLastTcData;
      std::vector<TdcData> fLastBhData;
      std::vector<TdcData> fLastMrSyncData;

      Bool_t fIsTerminated = false;

    public:
      MonitorWindow();
      ~MonitorWindow();

      Int_t                LoadOffset(const std::string& ffilename);
      void                 InitializeWindow();
      void                 InitializePlots();

      void                 DrawPlots();

      Int_t                UpdatePlots(const std::string& ifilename);

      inline void          Run() { fApplication->Run(true); }

      inline void          Terminate() {
        fIsTerminated = true;
        gSystem->ExitLoop();
      }

    private:
      void                 ClearLastSpill();
      void                 FillCoincidence(TdcData& tdcData);
      std::vector<TdcData> CollectCoinExtData(const TdcData& tdcData, std::size_t i);
      std::size_t          RemoveOldTdc(std::vector<TdcData>* lastData, Double_t time);
      std::size_t          RemoveOldMrSync(std::vector<TdcData>* lastData, Int_t globalChannel);
    };

    MonitorWindow::MonitorWindow() {
      fApplication = new TApplication("monitor", nullptr, nullptr);
    }

    MonitorWindow::~MonitorWindow() {
      delete fApplication;
      fApplication = nullptr;
    }

    Int_t MonitorWindow::LoadOffset(const std::string& ffilename) {
      std::cout << "Load offset" << std::endl;
      for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
        Analyzer::AnaCoin::CoinInfo info;
        std::ifstream ffile(Form(ffilename.data(), extCh));
        if (ffile) {
          while (info.Read(ffile)) {
            fCoinInfo[extCh][info.Index] = info;
          }
        }
      }

      for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
        for (std::size_t index = 0; index < CoinOffset::N; ++index) {
          fContains[extCh][index] = false;
        }
      }
      for (auto&& info : fCoinInfo) {
        const std::size_t extCh = info.first;
        for (auto&& subInfo : info.second) {
          const std::size_t index = subInfo.first;
          fContains[extCh][index] = true;
        }
      }

      return fCoinInfo.size();
    }

    void MonitorWindow::InitializeWindow() {
      std::cout << "Initialize window" << std::endl;
      fCanvas = new TCanvas("c1", (fProvider.GetName() + " | Semi-online Monitor").data(), 1600, 1200);

      std::vector<TPad*> pads;
      const Double_t seg = 1.0 / 6.0;

      pads.push_back(fPadBh1Tdc    = new TPad("fPadBh1Tdc"   , "", 0.0 * seg, 5.0 * seg, 1.0 * seg, 6.0 * seg));
      pads.push_back(fPadBh2Tdc    = new TPad("fPadBh2Tdc"   , "", 0.0 * seg, 4.0 * seg, 1.0 * seg, 5.0 * seg));
      pads.push_back(fPadExtTdc    = new TPad("fPadExtTdc"   , "", 0.0 * seg, 3.0 * seg, 1.0 * seg, 4.0 * seg));
      pads.push_back(fPadHodTdc    = new TPad("fPadHodTdc"   , "", 0.0 * seg, 2.0 * seg, 1.0 * seg, 3.0 * seg));
      pads.push_back(fPadTc1Tdc    = new TPad("fPadTc1Tdc"   , "", 0.0 * seg, 1.0 * seg, 1.0 * seg, 2.0 * seg));
      pads.push_back(fPadTc2Tdc    = new TPad("fPadTc2Tdc"   , "", 0.0 * seg, 0.0 * seg, 1.0 * seg, 1.0 * seg));

      pads.push_back(fPadHodHitMap = new TPad("fPadHodHitMap", "", 1.0 * seg, 4.5 * seg, 3.0 * seg, 6.0 * seg));
      pads.push_back(fPadHodHitCnt = new TPad("fPadHodHitCnt", "", 1.0 * seg, 3.0 * seg, 3.0 * seg, 4.5 * seg));
      pads.push_back(fPadExtHitMap = new TPad("fPadExtHitMap", "", 1.0 * seg, 1.5 * seg, 3.0 * seg, 3.0 * seg));
      pads.push_back(fPadExtHitCnt = new TPad("fPadExtHitCnt", "", 1.0 * seg, 0.0 * seg, 3.0 * seg, 1.5 * seg));

      pads.push_back(fPadMountain  = new TPad("fPadMountain" , "", 3.0 * seg, 4.0 * seg, 6.0 * seg, 6.0 * seg));
      pads.push_back(fPadTdcSync   = new TPad("fPadTdcSync"  , "", 3.0 * seg, 2.0 * seg, 6.0 * seg, 4.0 * seg));
      pads.push_back(fPadHit       = new TPad("fPadHit"      , "", 3.0 * seg, 0.0 * seg, 6.0 * seg, 2.0 * seg));

      fPadBh1Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadBh1Tdc   ->SetGrid(); fPadBh1Tdc   ->SetLogy();
      fPadBh2Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadBh2Tdc   ->SetGrid(); fPadBh2Tdc   ->SetLogy();
      fPadExtTdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadExtTdc   ->SetGrid(); fPadExtTdc   ->SetLogy();
      fPadHodTdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHodTdc   ->SetGrid(); fPadHodTdc   ->SetLogy();
      fPadTc1Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTc1Tdc   ->SetGrid(); fPadTc1Tdc   ->SetLogy();
      fPadTc2Tdc   ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTc2Tdc   ->SetGrid(); fPadTc2Tdc   ->SetLogy();
      fPadHodHitMap->SetMargin(0.06, 0.06, 0.10, 0.10);
      fPadHodHitCnt->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHodHitCnt->SetGrid(); fPadHodHitCnt->SetLogy();
      fPadExtHitMap->SetMargin(0.06, 0.06, 0.10, 0.10);
      fPadExtHitCnt->SetMargin(0.06, 0.06, 0.10, 0.10); fPadExtHitCnt->SetGrid(); fPadExtHitCnt->SetLogy();
      fPadMountain ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadMountain ->SetGrid(); 
      fPadTdcSync  ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadTdcSync  ->SetGrid(); fPadTdcSync  ->SetLogy();
      fPadHit      ->SetMargin(0.06, 0.06, 0.10, 0.10); fPadHit      ->SetGrid(); 

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
      const Int_t    xbinsInSpill =   50;

      const Double_t  xminInSync  = (Int_t)(- 5400 * nsec / timePerTdc) - 0.5; // [count]
      const Int_t    xbinsInSync0 = (Int_t)(+16200 * nsec / timePerTdc) / 20; // [count/20]
      const Int_t    xbinsInSync  = xbinsInSync0 < 200 ? xbinsInSync0 : xbinsInSync0 / 4;
      const Double_t  xmaxInSync  = xbinsInSync0 < 200 ? xminInSync + xbinsInSync * 20 : xminInSync + xbinsInSync * 80; // [count]

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
      // hExtEntryByCh = new TH1D("hExtEntryByCh",
      //                          Form("%s, Extinction Detector Entries by Channel;"
      //                               "Channel;"
      //                               "", tdcName.data()),
      //                          ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5);
      // hExtEntryByCh->SetStats(false);

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
        // hExtEntryByCh->Draw();
        hExtEntryByChBottom ->Draw();
        hExtEntryByChCenter1->Draw("same");
        hExtEntryByChCenter2->Draw("same");
        hExtEntryByChTop    ->Draw("same");
        hExtEntryByChBottom->SetMaximum(2.0 * hExtHitMap->GetBinContent(hExtHitMap->GetMaximumBin()));
        hExtEntryByChBottom->GetXaxis()->SetLabelSize(0.05);
        hExtEntryByChBottom->GetYaxis()->SetLabelSize(0.05);
        hExtEntryByChBottom->GetXaxis()->SetTitleSize(0.05);
        hExtEntryByChBottom->GetYaxis()->SetTitleSize(0.05);
        hExtEntryByChBottom->GetXaxis()->SetTitleOffset(0.6);
        hExtEntryByChBottom->GetYaxis()->SetTitleOffset(0.5);
      }
        
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

      fCanvas->Modified();
      fCanvas->Update();
    }

    void MonitorWindow::ClearLastSpill() {
      fCoinCount = 0;

      fLastExtData   .clear();
      fLastHodData   .clear();
      fLastTcData    .clear();
      fLastBhData    .clear();
      fLastMrSyncData.clear();

      hBhTdcInSpill[0]  ->Reset();
      hBhTdcInSpill[1]  ->Reset();
      hHodTdcInSpill_Any->Reset();
      hExtTdcInSpill_Any->Reset();
      hTcTdcInSpill[0]  ->Reset();
      hTcTdcInSpill[1]  ->Reset();
      hHodHitMap        ->Reset();
      hHodEntryByCh     ->Reset();
      hExtHitMap        ->Reset();
      // hExtEntryByCh     ->Reset();
      hExtMountain_Any  ->Reset();
      hExtTdcInSync_Any ->Reset();
    }

    void MonitorWindow::FillCoincidence(TdcData& extData) {
      const std::size_t extCh      = ExtinctionDetector::GetChannel(extData.Channel);
      const Long64_t    tdc        = extData.Tdc;
      const Double_t    time       = extData.Time;
      const Double_t    timePerTdc = fProvider.GetTimePerTdc();

      Bool_t coincidence[CoinOffset::N];
      std::vector<std::size_t> coinHodChs;
      for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
        const std::size_t i = bhCh + CoinOffset::BH;
        if (fContains[extCh][i]) {
          coincidence[i] = false;
          for (auto&& lastData : fLastBhData) {
            const Double_t dt    = lastData.Time - time;
            const Double_t mean  = fCoinInfo[extCh][i].FitMean * timePerTdc;
            if (TMath::Abs(dt - mean) < kCoinTimeWidth) {
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
            const Double_t dt    = lastData.Time - time;
            const Double_t mean  = fCoinInfo[extCh][i].FitMean * timePerTdc;
            if (TMath::Abs(dt - mean) < kCoinTimeWidth) {
              coincidence[i] = true;
              coinHodChs.push_back(Hodoscope::GetChannel(lastData.Channel));
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
            const Double_t dt    = lastData.Time - time;
            const Double_t mean  = fCoinInfo[extCh][i].FitMean * timePerTdc;
            if (TMath::Abs(dt - mean) < kCoinTimeWidth) {
              coincidence[i] = true;
              break;
            }
          }
        } else {
          coincidence[i] = true;
        }
      }

      if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
        for (std::size_t n = fLastMrSyncData.size(), i = n; i != 0; --i) {
          auto lastData = fLastMrSyncData[i];
          if (lastData.Channel == extData.MrSyncChannel) {
            auto syncTdc = lastData.Tdc;
            ++fCoinCount;
            hExtTdcInSync_Any->Fill(tdc - syncTdc);
            hExtMountain_Any ->Fill(tdc - syncTdc, time / msec);
            break;
          }
        }
      }
    }

    std::vector<TdcData> MonitorWindow::CollectCoinExtData(const TdcData& tdcData, std::size_t i) {
      const Double_t time       = tdcData.Time;
      const Double_t timePerTdc = fProvider.GetTimePerTdc();

      std::vector<TdcData> coinExtData;
      for (auto&& lastData : fLastExtData) {
        const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
        if (fContains[extCh][i]) {
          const Double_t dt    = time - lastData.Time;
          const Double_t mean  = fCoinInfo[extCh][i].FitMean * timePerTdc;
          if (TMath::Abs(dt - mean) < kCoinTimeWidth) {
            coinExtData.push_back(lastData);
          }
        }
      }
      return coinExtData;
    }

    std::size_t MonitorWindow::RemoveOldTdc(std::vector<TdcData>* lastData, Double_t time) {
      for (std::size_t i = 0, n = lastData->size(); i < n; ++i) {
        if (TMath::Abs(time - lastData->at(0).Time) > kLastThreshold) {
          lastData->erase(lastData->begin());
        } else {
          break;
        }
      }
      return lastData->size();
    }

    std::size_t MonitorWindow::RemoveOldMrSync(std::vector<TdcData>* lastData, Int_t globalChannel) {
      Bool_t found = false;
      for (std::size_t n = lastData->size(), i = n; i != 0; --i) {
        if (lastData->at(i - 1).Channel == globalChannel) {
          if (found) {
            lastData->erase(lastData->begin() + i - 1);
          } else {
            // Save last one
            found = true;
          }
        }
      }
      return lastData->size();
    }

    Int_t MonitorWindow::UpdatePlots(const std::string& ifilename) {
      using CoinOffset = Analyzer::AnaTimeOffset::CoinOffset;

      std::cout << "Initialize decoder" << std::endl;
      Decoder decoder;
      Packet_t packet;

      std::cout << "Initialzie history" << std::endl;
      ClearLastSpill();

      std::cout << "Open file" << std::endl;
      std::ifstream ifile(ifilename, std::ios::binary);
      if (!ifile) {
        std::cout << "[error] input file is not opened, " << ifilename << std::endl;
        return 1;
      }

      {
        std::size_t count = 0UL;
        for (; decoder.Read(ifile, &packet); ++count) {
          if (fIsTerminated) {
            return 0;
          }

          if (count % 1000000UL == 0) {
            std::cout << ">> " << count << std::endl;
          }

          if (decoder.Data.Type == DataType::Header) {
            std::cout << "[info] detect header" << std::endl;

          } else if (decoder.Data.Type == DataType::GateStart) {
            std::cout << "[info] begin of spill " << decoder.Data.Spill << std::endl;
            ClearLastSpill();

          } else if (decoder.Data.Type == DataType::GateEnd) {
            std::cout << "[info] end of spill " << decoder.Data.Spill << std::endl;

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

            DrawPlots();

          } else {
            std::vector<TdcData> tdcData = decoder.Data.GetTdcData();

            for (auto&& data : tdcData) {
              const Int_t    globalChannel = data.Channel;
              const Double_t time          = data.Time;

              if (ExtinctionDetector::Contains(globalChannel)) {
                const Int_t ch = ExtinctionDetector::GetChannel(globalChannel);

                // hExtEntryByCh     ->Fill(ch);
                ExtinctionDetector::Fill(hExtHitMap, ch);
                hExtTdcInSpill_Any->Fill(time / msec);

                FillCoincidence(data);

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

                std::vector<TdcData> coinExtData = CollectCoinExtData(data, CoinOffset::Hod);
                fLastHodData.push_back(data);
                for (auto&& extData : coinExtData) {
                  FillCoincidence(extData);
                }

                if (RemoveOldTdc(&fLastHodData, time) > kHistLimit) {
                  std::cerr << "[error] size of fLastHodData reaches " << kHistLimit << std::endl;
                  return 1;
                }

              } else if (TimingCounter::Contains(globalChannel)) {
                const Int_t ch = TimingCounter::GetChannel(globalChannel);

                hTcTdcInSpill[ch]->Fill(time / msec);

                std::vector<TdcData> coinExtData = CollectCoinExtData(data, ch + CoinOffset::TC);
                fLastTcData.push_back(data);
                for (auto&& extData : coinExtData) {
                  FillCoincidence(extData);
                }

                if (RemoveOldTdc(&fLastTcData, time) > kHistLimit) {
                  std::cerr << "[error] size of fLastTcData reaches " << kHistLimit << std::endl;
                  return 1;
                }

              } else if (BeamlineHodoscope::Contains(globalChannel)) {
                const Int_t ch = BeamlineHodoscope::GetChannel(globalChannel);

                hBhTdcInSpill[ch]->Fill(time / msec);

                std::vector<TdcData> coinExtData = CollectCoinExtData(data, ch + CoinOffset::BH);
                fLastBhData.push_back(data);
                for (auto&& extData : coinExtData) {
                  FillCoincidence(extData);
                }

                if (RemoveOldTdc(&fLastBhData, time) > kHistLimit) {
                  std::cerr << "[error] size of fLastBhData reaches " << kHistLimit << std::endl;
                  return 1;
                }

              } else if (MrSync::Contains(globalChannel)) {
                if (RemoveOldMrSync(&fLastMrSyncData, globalChannel) > kHistLimit) {
                  std::cerr << "[error] size of fLastMrSyncData reaches " << kHistLimit << std::endl;
                  return 1;
                }
                fLastMrSyncData.push_back(data);
              }

            }
          }
        }
      }

      std::cout << "Close files" << std::endl;
      ifile.close();

      return 0;
    }

  }

}

#endif
