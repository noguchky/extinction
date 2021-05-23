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

#include "Math.hh"
#include "Linq.hh"
#include "String.hh"
#include "ObjectHelper.hh"
#include "ScopeSubstituter.hh"

namespace Extinction {

  namespace Analyzer {

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
      using CoinDiffs  = std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Double_t>>;

      struct CoinOffset {
        enum {
          BH  = 0,
          Hod = 0 + BeamlineHodoscope::NofChannels,
          TC  = 0 + BeamlineHodoscope::NofChannels + 1,
          N   = 0 + BeamlineHodoscope::NofChannels + 1 + TimingCounter::NofChannels,
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
      Bool_t                       fDebug                  = false;
      
      ITdcDataProvider*            fProvider               = nullptr;

      TH1D*                        hBh1Timeline            = nullptr;
      TH1D*                        hBh2Timeline            = nullptr;
      TH1D*                        hTc1Timeline            = nullptr;
      TH1D*                        hTc2Timeline            = nullptr;
      TH1D*                        hExtTimeline            = nullptr;

      TH1D*                        hTmpBh1Timeline         = nullptr;
      TH1D*                        hTmpBh2Timeline         = nullptr;
      TH1D*                        hTmpTc1Timeline         = nullptr;
      TH1D*                        hTmpTc2Timeline         = nullptr;
      std::vector<TH1D*>           hTmpExtTimeline;

      TH1D*                        hCoinTdcInSync          = nullptr;
      TH2D*                        hCoinMountain           = nullptr;

      TEfficiency*                 hEfficiency             = nullptr;
      TH1*                         hEfficiencyFrame        = nullptr;

      TFile*                       fSpillFile              = nullptr;
      TTree*                       fSpillTree              = nullptr;

      Long64_t                     fSpillCount             = 0;
      SpillData                    fSpillData;

      CoinDiffs                    fStdCoinDiffs;
      std::map<Int_t, Double_t>    fStdTimePerTdc;
      std::map<Int_t, Double_t>    fStdMrSyncInterval;
      Double_t                     fStdMrSyncIntervalAverage;
      std::map<Int_t, TdcData>     fLastMrSync;
      std::map<Int_t, TdcData>     fNextMrSync;
      std::map<Int_t, TdcData>     fNext2MrSync;

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
      
      Double_t                     fMrSyncOffset           =   0.0; // tdc count
      Double_t                     fCoinTimeWidth          =  50.0 * nsec;

    public:
      TimelineCoincidence(ITdcDataProvider* provider);
      ~TimelineCoincidence();

      void                 SetTimePerTdc(const std::map<Int_t, Double_t>& map);
      void                 SetMrSyncInterval(const std::map<Int_t, Double_t>& map);
      Int_t                LoadOffset(const std::string& ffilename);

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

      inline void          SetDebug(bool debug) { fDebug = debug; }

      inline void          SetMrSyncOffset(Double_t offset) {
        std::cout << "SetMrSyncOffset ... " << offset << std::endl;
        fMrSyncOffset = offset;
      }
      inline Double_t      GetMrSyncOffset() const { return fMrSyncOffset; }

      inline void          SetCoinTimeWidth(Double_t width) {
        std::cout << "SetCoinTimeWidth ... " << width / nsec << " nsec" << std::endl;
        fCoinTimeWidth = width;
      }
      inline Double_t      GetCoinTimeWidth() const { return fCoinTimeWidth; }

      void                 InitializePlots(const PlotsProfiles& profile);
      void                 InitializeSpillSummary(const std::string& filename, const std::string& treename = "spilltree");

      void                 DrawPlots(const std::string& ofilename);
      void                 WritePlots(const std::string& ofilename);
      void                 WriteSpillSummary();

      Int_t                GenerateEfficiency(std::map<Int_t, ITdcDataProvider*> providers,
                                              const std::map<Int_t, std::string>& ifilenames,
                                              const std::string& treename,
                                              const std::function<TDatime(const std::string&)>& parser = nullptr);

      Int_t                GeneratePlots(std::map<Int_t, ITdcDataProvider*> providers,
                                         const std::map<Int_t, std::string>& ifilenames,
                                         const std::string& treename,
                                         const std::function<TDatime(const std::string&)>& parser = nullptr);

    private:
      void                 ClearLastSpill(Bool_t clearHists);

      void                 DrawTmpTimeline(Int_t bin, Int_t range);

      template <typename T, typename V>
      inline void          FillToSeconds(std::map<T, V>* map, V value) {
        for (auto&& pair : *map) {
          pair.second = value;
        }
      }
    };

    TimelineCoincidence::TimelineCoincidence(ITdcDataProvider* provider)
      : fProvider(provider) {
      fEfficiencyTargetBh1 = false;
      fEfficiencyTargetBh2 = false;
      fEfficiencyTargetHod = false;
      fEfficiencyTargetExt = false;
      fEfficiencyTargetTc1 = false;
      fEfficiencyTargetTc2 = false;
    }

    TimelineCoincidence::~TimelineCoincidence() {
    }

    void TimelineCoincidence::SetTimePerTdc(const std::map<Int_t, Double_t>& map) {
      std::cout << "SetTimePerTdc" << std::endl;
      for (auto&& pair : map) {
        std::cout << pair.first << "\t" << pair.second << std::endl;
      }
      fStdTimePerTdc = map;
    }

    void TimelineCoincidence::SetMrSyncInterval(const std::map<Int_t, Double_t>& map) {
      std::cout << "SetMrSyncInterval" << std::endl;
      for (auto&& pair : map) {
        std::cout << pair.first << "\t" << pair.second << std::endl;
      }
      fStdMrSyncInterval        = map;
      fStdMrSyncIntervalAverage = Tron::Linq::From(map)
        .Select([](const std::pair<Int_t, Double_t>& pair) { return pair.second; })
        .Average();
    }

    Int_t TimelineCoincidence::LoadOffset(const std::string& ffilename) {
      std::cout << "Load offset" << std::endl;
      if (Tron::String::GetFormatSpecifiers(ffilename).empty()) {
        std::ifstream ffile(ffilename);
        if (ffile) {
          Int_t ch, index; Double_t mean;
          while (ffile >> ch >> index >> mean) {
            fStdCoinDiffs[ch][index] = mean;
            // std::cout << ch << "\t" << index << "\t" << mean << "\t" << fStdCoinDiffs[ch][index] / nsec << std::endl;
            // std::cout << "[" << ch << "]" << " [" << index << "] " << fStdCoinDiffs[ch][index] / nsec << " nsec" << std::endl;
          }
        } else {
        }
      } else {
        for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
          CoinInfo info;
          std::ifstream ffile(Form(ffilename.data(), extCh));
          if (ffile) {
            while (info.Read(ffile)) {
              fStdCoinDiffs[extCh][info.Index] = info.FitMean;
              // std::cout << extCh << "\t" << info.Index << "\t" << info.FitMean << std::endl;
              // std::cout << "[" << extCh << "]" << " [" << info.Index << "] " << fStdCoinDiffs[extCh][info.Index] / nsec << " nsec" << std::endl;
            }
          }
        }
      }

      return fStdCoinDiffs.size();
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

      const Int_t xbinsTmpTimeline = fStdMrSyncIntervalAverage * 2.5;
      const Int_t xminTmpTimeline  = 0;
      const Int_t xmaxTmpTimeline  = xbinsTmpTimeline;

      // Timeline
      hBh1Timeline  = new TH1D("hBh1Timeline",
                               Form("%s, BH1 Timeline;"
                                    "TDC [count]", tdcName.data()),
                               xbinsInSync, xminInSync, xmaxInSync);
      hBh2Timeline  = new TH1D("hBh2Timeline",
                               Form("%s, BH2 Timeline;"
                                    "TDC [count]", tdcName.data()),
                               xbinsInSync, xminInSync, xmaxInSync);
      hTc1Timeline  = new TH1D("hTc1Timeline",
                               Form("%s, TC1 Timeline;"
                                    "TDC [count]", tdcName.data()),
                               xbinsInSync, xminInSync, xmaxInSync);
      hTc2Timeline  = new TH1D("hTc2Timeline",
                               Form("%s, TC2 Timeline;"
                                    "TDC [count]", tdcName.data()),
                               xbinsInSync, xminInSync, xmaxInSync);
      hExtTimeline = new TH1D("hExtTimeline",
                              Form("%s, Ext Timeline;"
                                   "TDC [count]", tdcName.data()),
                              xbinsInSync, xminInSync, xmaxInSync);

      hTmpBh1Timeline  = new TH1D("hTmpBh1Timeline",
                                  Form("%s, BH1 Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpBh2Timeline  = new TH1D("hTmpBh2Timeline",
                                  Form("%s, BH2 Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpTc1Timeline  = new TH1D("hTmpTc1Timeline",
                                  Form("%s, TC1 Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpTc2Timeline  = new TH1D("hTmpTc2Timeline",
                                  Form("%s, TC2 Timeline;"
                                       "TDC [count]", tdcName.data()),
                                  xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      hTmpExtTimeline = std::vector<TH1D*>(MrSync::NofChannels);
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hTmpExtTimeline[ch] = new TH1D(Form("hTmpExtTimeline_%03lu", ch),
                                       Form("%s, Ext Timeline;"
                                            "TDC [count]", tdcName.data()),
                                       xbinsTmpTimeline, xminTmpTimeline, xmaxTmpTimeline);
      }

      // Extinction Detector TDC in sync (coincidence)
      hCoinTdcInSync = new TH1D("hCoinTdcInSync",
                                Form("%s, TDC in MR Sync (Coincidence);"
                                     "TDC [count]", tdcName.data()),
                                xbinsInSync, xminInSync, xmaxInSync);

      // Extinction Detector Mountain Plot (coincidence)
      hCoinMountain = new TH2D("hCoinMountain",
                               Form("%s, Mountain Plot (Coincidence);"
                                    "TDC [count];"
                                    "Time [ms]", tdcName.data()),
                               xbinsInSync / 2, xminInSync, xmaxInSync,
                               xbinsInSpill / 2, xminInSpill, xmaxInSpill);
      hCoinMountain->SetStats(false);

      hEfficiency = new TEfficiency("hEfficiency",
                                    "Efficiency",
                                    6, 0, 6);
    }

    void TimelineCoincidence::InitializeSpillSummary(const std::string& filename, const std::string& treename) {
      std::cout << "Initialize spill summary" << std::endl;

      fSpillFile = new TFile(filename.data(), "RECREATE");
      if (!fSpillFile->IsOpen()) {
        std::cout << "[error] spill summary file is not opened, " << filename << std::endl;
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

      std::cout << "hCoinTdcInSync" << std::endl;
      gPad->SetLogy(true);
      {
        hBh1Timeline->SetLineColor(51);
        hBh2Timeline->SetLineColor(61);
        hTc1Timeline->SetLineColor(71);
        hTc2Timeline->SetLineColor(81);
        hExtTimeline->SetLineColor(91);

        hBh1Timeline->Draw();
        hBh2Timeline->Draw("same");
        hTc1Timeline->Draw("same");
        hTc2Timeline->Draw("same");
        hExtTimeline->Draw("same");

        gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      std::cout << "hCoinTdcInSync" << std::endl;
      gPad->SetLogy(true);
      {
        hCoinTdcInSync->Draw("hist");
        hCoinTdcInSync->SetMinimum(0.2);
        gPad->Print(ofilename.data());

        // hCoinTdcInSync->Draw();
        // hCoinTdcInSync->SetMinimum(0.2);
        // gPad->Print(ofilename.data());
      }
      gPad->SetLogy(false);

      std::cout << "hCoinMountain" << std::endl;
      gPad->SetGrid(false, true);
      {
        hCoinMountain->Draw("colz");
        hCoinMountain->SetMinimum(0);
        gPad->Print(ofilename.data());
      }
      gPad->SetGrid(true, true);

      std::cout << "hEfficiency" << std::endl;
      {
        if (!hEfficiencyFrame) {
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
        }
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

      Tron::ObjectHelper::WriteValue<Long64_t>(fSpillData.Date.Convert(), "Time");
      hBh1Timeline  ->Write();
      hBh2Timeline  ->Write();
      hTc1Timeline  ->Write();
      hTc2Timeline  ->Write();
      hExtTimeline  ->Write();
      hCoinTdcInSync->Write();
      hCoinMountain ->Write();

      file->Close();
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

    void TimelineCoincidence::ClearLastSpill(Bool_t clearHists) {
      for (std::size_t ch = 0; ch <= MrSync::NofChannels; ++ch) {
        fLastMrSync [ch] = { };
        fNextMrSync [ch] = { };
        fNext2MrSync[ch] = { };
      }

      if (clearHists) {
        hBh1Timeline  ->Reset();
        hBh2Timeline  ->Reset();
        hTc1Timeline  ->Reset();
        hTc2Timeline  ->Reset();
        hExtTimeline  ->Reset();
        hCoinTdcInSync->Reset();
        hCoinMountain ->Reset();
      }
    }

    Int_t TimelineCoincidence::GenerateEfficiency(std::map<Int_t, ITdcDataProvider*> providers,
                                                  const std::map<int, std::string>& ifilenames,
                                                  const std::string& treename,
                                                  const std::function<TDatime(const std::string&)>& parser) {
      for (std::size_t i = 0; i < 6; ++i) {
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
            (fCoincidenceTargetTc2 && fEfficiencyTargetTc2)) {
          GeneratePlots(providers,
                        ifilenames,
                        treename,
                        parser);
        }
      }

      fEfficiencyTargetBh1 = false;
      fEfficiencyTargetBh2 = false;
      fEfficiencyTargetHod = false;
      fEfficiencyTargetExt = false;
      fEfficiencyTargetTc1 = false;
      fEfficiencyTargetTc2 = false;

      return 0;
    }
    
    Int_t TimelineCoincidence::GeneratePlots(std::map<Int_t, ITdcDataProvider*> providers,
                                             const std::map<int, std::string>& ifilenames,
                                             const std::string& treename,
                                             const std::function<TDatime(const std::string&)>& parser) {
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
      ClearLastSpill(true);

      std::cout << "Open file" << std::endl;
      std::map<Int_t, TFile*>   ifiles;
      std::map<Int_t, TTree*>   itrees;
      std::map<Int_t, Long64_t> entrieses;
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

        entrieses[board] = itrees[board]->GetEntries();
        std::cout << " + " << entrieses[board] << std::endl;

        providers[board]->SetBranchAddress(itrees[board]);
      }
      std::cout << " = " << Tron::Linq::From(entrieses).Sum([](std::pair<Int_t, Long64_t> p) { return p.second; }) << std::endl;
      const UInt_t averageTime = Tron::Linq::From(datimes)
        .Select([](TDatime datime) -> ULong64_t { return datime.Convert(); })
        .Average();
      fSpillData.SetDate(averageTime);

      {
        TdcData bhRef[2], tcRef[2];
        Int_t targetBoard = 0;
        std::map<Int_t, Int_t>                lastSpills;
        std::map<Int_t, Long64_t>             entries;
        std::map<Int_t, Bool_t>               spillEnded;
        std::map<Int_t, Bool_t>               fileEnded;
        for (auto&& pair : ifilenames) {
          targetBoard = pair.first;
          lastSpills       [pair.first] = -1;
          entries          [pair.first] = 0;
          spillEnded       [pair.first] = false;
          fileEnded        [pair.first] = false;
        }
        // std::cout << "targetBoard = " << targetBoard << std::endl; 

        Double_t* pBh1 = hTmpBh1Timeline->GetArray();
        Double_t* pBh2 = hTmpBh2Timeline->GetArray();
        Double_t* pTc1 = hTmpTc1Timeline->GetArray();
        Double_t* pTc2 = hTmpTc2Timeline->GetArray();
        std::vector<Double_t*> pExt(MrSync::NofChannels);
        for (std::size_t board = 0; board < MrSync::NofChannels; ++board ) {
          pExt[board] = hTmpExtTimeline[board]->GetArray();
        }

        const Int_t coinBinWidth = fCoinTimeWidth / fProvider->GetTimePerTdc();
        for (std::size_t count = 0UL;;) {

          // std::cout << "[debug] create timelines" << std::endl;
          // for (std::size_t i = 0; i < MrSync::NofChannels; ++i)
          //   std::cout << fLastMrSync[i].Tdc << "\t" << entries[i] << std::endl;
          // Create timeline
          for (auto&& pair : ifilenames) {
            targetBoard = pair.first;
            ITdcDataProvider* provider   = providers  [targetBoard];
            TTree*            itree      = itrees     [targetBoard];
            Long64_t&         entry      = entries    [targetBoard];
            Int_t&            lastSpill  = lastSpills [targetBoard];
            const Long64_t    ientries   = entrieses  [targetBoard];
            for (; !fNext2MrSync[targetBoard].Time; ++count, ++entry) {
              if (count % 1000000UL == 0) {
                std::cout << ">> " << count << std::endl;
              }

              // std::cout << "[debug] GetEntry" << std::endl;

              if (entry >= ientries) {
                std::cout << "[info] detect file end @ " << targetBoard << std::endl;
                fileEnded[targetBoard] = true;
                break;

              } else if (!itree->GetEntry(entry)) {
                std::cout << "[info] detect file end @ " << targetBoard << " (TTree::GetEntry)" << std::endl;
                fileEnded[targetBoard] = true;
                break;

              }

              // std::cout << "[debug] Data Check" << std::endl;

              if (!provider->IsData()) {
                continue;

              } else if (lastSpill != -1 && lastSpill != provider->GetSpill()) {
                std::cout << "[info] detect spill end @ " << targetBoard << std::endl;
                lastSpill               = provider->GetSpill();
                spillEnded[targetBoard] = true;
                break;

              } else {
                // std::cout << "[info] detect data @ " << targetBoard << std::endl;
                if (lastSpill != -1 && fSpillData.EMCount != provider->GetEMCount()) {
                  std::cout << "[error] conflict EMCount, " << fSpillData.EMCount << " <--> " << provider->GetEMCount() << std::endl;
                  exit(1);
                }
                lastSpill          = provider->GetSpill();
                fSpillData.EMCount = provider->GetEMCount();
                auto tdcData       = provider->GetTdcData(targetBoard);

                // std::cout << "[debug] data process" << std::endl;

                for (auto&& data : tdcData) {
                  if (MrSync::Contains(data.Channel)) {
                    // std::cout << "[debug] get mrsync" << std::endl;

                    data.Tdc += fMrSyncOffset;
                    if        (fNextMrSync[targetBoard].Tdc) {
                      // std::cout << "[debug] get next next mrsync" << std::endl;
                      fNext2MrSync[targetBoard] = data;
                    } else if (fLastMrSync[targetBoard].Tdc) {
                      // std::cout << "[debug] get next mrsync" << std::endl;
                      fNextMrSync [targetBoard] = data;
                    } else {
                      // std::cout << "[debug] get first mrsync" << std::endl;
                      fLastMrSync [targetBoard] = data;
                    }

                  } else if (BeamlineHodoscope::Contains  (data.Channel)) {
                    // std::cout << "[debug] fill bh timeline" << std::endl;

                    const std::size_t bhCh = BeamlineHodoscope::GetChannel(data.Channel);
                    const std::size_t tdcI = bhCh + CoinOffset::BH;
                    data.Tdc -= fStdCoinDiffs[40][tdcI];
                    if (auto syncTdc = fLastMrSync[targetBoard].Tdc) {
                      // (bhCh ? hTmpBh2Timeline : hTmpBh1Timeline)->Fill(data.Tdc - syncTdc + bin);
                      TH1D* hist = bhCh ? hTmpBh2Timeline : hTmpBh1Timeline;
                      const Long64_t tdc0 = data.Tdc - syncTdc;
                      for (Int_t dtdc = 0; dtdc <= coinBinWidth; ++dtdc) {
                        hist->Fill(tdc0 + dtdc);
                      }
                    }
                    bhRef[bhCh] = data;

                  } else if (TimingCounter    ::Contains  (data.Channel)) {
                    // std::cout << "[debug] fill tc timeline" << std::endl;

                    const std::size_t tcCh = TimingCounter::GetChannel(data.Channel);
                    const std::size_t tdcI = tcCh + CoinOffset::TC;
                    data.Tdc -= fStdCoinDiffs[40][tdcI];
                    if (auto syncTdc = fLastMrSync[targetBoard].Tdc) {
                      // (tcCh ? hTmpTc2Timeline : hTmpTc1Timeline)->Fill(data.Tdc - syncTdc);
                      TH1D* hist = tcCh ? hTmpTc2Timeline : hTmpTc1Timeline;
                      const Long64_t tdc0 = data.Tdc - syncTdc;
                      for (Int_t dtdc = 0; dtdc <= coinBinWidth; ++dtdc) {
                        hist->Fill(tdc0 + dtdc);
                      }
                    }
                    tcRef[tcCh] = data;

                  } else if (ExtinctionDetector::Contains(data.Channel)) {
                    // std::cout << "[debug] fill ext timeline" << std::endl;

                    const std::size_t extCh = ExtinctionDetector::GetChannel(data.Channel);
                    // (tdc - sync) - (bh - sync) - coinDiffs[ch][sync] = (40 - sync) - (bh - sync) - coinDiffs[40][sync]
                    // tdc + coinDiffs[40][sync] - coinDiffs[ch][sync] = 40
                    data.Tdc -= fStdCoinDiffs[extCh][CoinOffset::BH] - fStdCoinDiffs[40][CoinOffset::BH];
                    if (auto syncTdc = fLastMrSync[targetBoard].Tdc) {
                      // hTmpExtTimeline[targetBoard]->Fill(data.Tdc - syncTdc);
                      TH1D* hist = hTmpExtTimeline[targetBoard];
                      const Long64_t tdc0 = data.Tdc - syncTdc;
                      for (Int_t dtdc = 0; dtdc <= coinBinWidth; ++dtdc) {
                        hist->Fill(tdc0 + dtdc);
                      }
                    }

                  } else {
                    // std::cout << "[debug] skip others" << std::endl;
                    continue;
                  }
                }
              }
            }
          }

          {
            const Int_t syncEdge = fNextMrSync[0].Tdc - fLastMrSync[0].Tdc;

            Int_t sumBh1 = 0, sumBh2 = 0, sumTc1 = 0, sumTc2 = 0, sumExt = 0;
            auto isCoincident =
              [&] (Int_t bin) {
                {
                  sumExt = 0;
                  /*                   */ sumBh1  = pBh1[bin];
                  /*                   */ sumBh2  = pBh2[bin];
                  for (auto&& p : pExt) { sumExt += p   [bin];
                  /*                   */ sumTc1  = pTc1[bin];
                  /*                   */ sumTc2  = pTc2[bin]; }
                }

                return
                  (sumBh1 || !fCoincidenceTargetBh1 || fEfficiencyTargetBh1) &&
                  (sumBh2 || !fCoincidenceTargetBh2 || fEfficiencyTargetBh2) &&
               // (sumHod || !fCoincidenceTargetHod || fEfficiencyTargetHod) &&
                  (sumExt || !fCoincidenceTargetExt || fEfficiencyTargetExt) &&
                  (sumTc1 || !fCoincidenceTargetTc1 || fEfficiencyTargetTc1) &&
                  (sumTc2 || !fCoincidenceTargetTc2 || fEfficiencyTargetTc2);
              };

            for (Int_t bin = 1; bin <= syncEdge; ++bin) {

              if (isCoincident(bin)) {
                hCoinTdcInSync->Fill(bin - 1);
                hCoinMountain ->Fill(bin - 1, fLastMrSync[0].Time / msec);

                if (fDebug) {
                  DrawTmpTimeline(bin, 25 * coinBinWidth);
                }

                // Skip dead time
                Bool_t filled = false;
                for (; bin <= syncEdge; ++bin) {
                  if (isCoincident(bin)) {
                    if (fEfficiencyTargetBh1 && sumBh1 && !filled) { hEfficiency->Fill(true, 0); filled = true; }
                    if (fEfficiencyTargetBh2 && sumBh2 && !filled) { hEfficiency->Fill(true, 1); filled = true; }
                 // if (fEfficiencyTargetHod && sumHod && !filled) { hEfficiency->Fill(true, 2); filled = true; }
                    if (fEfficiencyTargetExt && sumExt && !filled) { hEfficiency->Fill(true, 3); filled = true; }
                    if (fEfficiencyTargetTc1 && sumTc1 && !filled) { hEfficiency->Fill(true, 4); filled = true; }
                    if (fEfficiencyTargetTc2 && sumTc2 && !filled) { hEfficiency->Fill(true, 5); filled = true; }
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
              }
            }
          }

          // // std::cout << "[debug] coincidence" << std::endl;
          // // Coincidence
          // {
          //   const Int_t syncEdge = fNextMrSync[0].Tdc - fLastMrSync[0].Tdc;

          //   Int_t bin1 = 1, bin2 = coinBinWidth;

          //   Int_t sumBh1 = 0, sumBh2 = 0, sumTc1 = 0, sumTc2 = 0, sumExt = 0;
          //   for (Int_t bin = bin1; bin <= bin2; ++bin) {
          //     /*                   */ sumBh1 += pBh1[bin];
          //     /*                   */ sumBh2 += pBh2[bin];
          //     /*                   */ sumTc1 += pTc1[bin];
          //     /*                   */ sumTc2 += pTc2[bin];
          //     for (auto&& p : pExt) { sumExt += p   [bin]; }
          //   }

          //   for (++bin1, ++bin2; bin1 <= syncEdge; ++bin1, ++bin2) {

          //     // if (         sumTc && sumExt) {
          //     //   hEfficiency->Fill(sumBh , 0);
          //     // }
          //     // if (sumBh          && sumExt) {
          //     //   hEfficiency->Fill(sumTc , 1);
          //     // }
          //     // if (sumBh && sumTc          ) {
          //     //   hEfficiency->Fill(sumExt, 2);
          //     // }

          //     if ((sumBh1 || !fCoincidenceTargetBh1) &&
          //         (sumBh2 || !fCoincidenceTargetBh2) &&
          //         (sumExt || !fCoincidenceTargetExt) &&
          //         (sumTc1 || !fCoincidenceTargetTc1) &&
          //         (sumTc2 || !fCoincidenceTargetTc2)) {
          //       hCoinTdcInSync->Fill(bin1 - 1);
          //       hCoinMountain ->Fill(bin1 - 1, fLastMrSync[0].Time / msec);

          //       if (!fCoincidenceTargetBh1) { hEfficiency->Fill(sumBh1, 0); }
          //       if (!fCoincidenceTargetBh2) { hEfficiency->Fill(sumBh2, 1); }
          //    // if (!fCoincidenceTargetHod) { hEfficiency->Fill(sumHod, 2); }
          //       if (!fCoincidenceTargetExt) { hEfficiency->Fill(sumTc2, 3); }
          //       if (!fCoincidenceTargetTc1) { hEfficiency->Fill(sumExt, 4); }
          //       if (!fCoincidenceTargetTc2) { hEfficiency->Fill(sumTc1, 5); }

          //       // Skip dead time
          //       bin1 += coinBinWidth;
          //       bin2 += coinBinWidth;
          //       if (bin1 > syncEdge) {
          //         break;
          //       }

          //       sumBh1 = sumBh2 = sumTc1 = sumTc2 = sumExt = 0;
          //       for (Int_t bin = bin1; bin <= bin2; ++bin) {
          //         /*                   */ sumBh1 += pBh1[bin];
          //         /*                   */ sumBh2 += pBh2[bin];
          //         for (auto&& p : pExt) { sumExt += p   [bin];
          //         /*                   */ sumTc1 += pTc1[bin];
          //         /*                   */ sumTc2 += pTc2[bin]; }
          //       }

          //     } else {
          //       /*                   */ sumBh1 += pBh1[bin2] - pBh1[bin1 - 1];
          //       /*                   */ sumBh2 += pBh2[bin2] - pBh2[bin1 - 1];
          //       for (auto&& p : pExt) { sumExt += p   [bin2] - p   [bin1 - 1];
          //       /*                   */ sumTc1 += pTc1[bin2] - pTc1[bin1 - 1];
          //       /*                   */ sumTc2 += pTc2[bin2] - pTc2[bin1 - 1]; }
          //     }

          //     // Data check
          //     if (sumBh1 < 0 || sumBh2 < 0 || sumTc1 < 0 || sumTc2 < 0 || sumExt < 0) {
          //       const Int_t binmin = TMath::Max(bin1 - 2 * coinBinWidth, 1);
          //       const Int_t binmax = TMath::Min(bin2 + 2 * coinBinWidth, syncEdge);

          //       std::cout << "sumBh1 = " << sumBh1 << std::endl;
          //       std::cout << "BH1 : ";
          //       for (Int_t bin = binmin; bin < binmax; ++bin) {
          //         std::cout << (pBh1[bin] ? Form("%d", (Int_t)pBh1[bin]) : ".");
          //       }
          //       std::cout << std::endl;

          //       std::cout << "sumBh2 = " << sumBh2 << std::endl;
          //       std::cout << "BH2 : ";
          //       for (Int_t bin = binmin; bin < binmax; ++bin) {
          //         std::cout << (pBh2[bin] ? Form("%d", (Int_t)pBh2[bin]) : ".");
          //       }
          //       std::cout << std::endl;

          //       std::cout << "sumTc1 = " << sumTc1 << std::endl;
          //       std::cout << "TC1 : ";
          //       for (Int_t bin = binmin; bin < binmax; ++bin) {
          //         std::cout << (pTc1[bin] ? Form("%d", (Int_t)pTc1[bin]) : ".");
          //       }
          //       std::cout << std::endl;

          //       std::cout << "sumTc2 = " << sumTc2 << std::endl;
          //       std::cout << "TC2 : ";
          //       for (Int_t bin = binmin; bin < binmax; ++bin) {
          //         std::cout << (pTc2[bin] ? Form("%d", (Int_t)pTc2[bin]) : ".");
          //       }
          //       std::cout << std::endl;

          //       std::cout << "sumExt = " << sumExt << std::endl;
          //       for (std::size_t board = 0; board < MrSync::NofChannels; ++board) {                  
          //         std::cout << "Ext : ";
          //         for (Int_t bin = binmin; bin < binmax; ++bin) {
          //           std::cout << (pExt[board][bin] ? Form("%d", (Int_t)pExt[board][bin]) : ".");
          //         }
          //         std::cout << std::endl;
          //       }

          //       // while (std::getchar() != '\n');
          //     }
          //   }
          // }

          // std::cout << "[debug] Shift timeline" << std::endl;
          // Shift timeline
          auto shiftTimeline =
            [&](Int_t board, TH1D* tmpTimeline, Double_t* p, TH1D* timeline) {
              const Int_t dTdc = fNextMrSync[board].Tdc - fLastMrSync[board].Tdc;
              const Int_t nbins = tmpTimeline->GetNbinsX();
              for (Int_t bin = 1, binmax = nbins - dTdc; bin <= binmax; ++bin) {
                p[bin] = p[bin + dTdc];
                if (p[bin] && bin - 1 < dTdc) {
                  timeline->Fill(bin - 1);
                }
              }
              for (Int_t bin = nbins - dTdc; bin <= nbins; ++bin) {
                p[bin] = 0;
              }
            };

          if (bhRef[0].Tdc) {
            shiftTimeline(bhRef[0].Board, hTmpBh1Timeline, pBh1, hBh1Timeline);
          }
          if (bhRef[1].Tdc) {
            shiftTimeline(bhRef[1].Board, hTmpBh2Timeline, pBh2, hBh2Timeline);
          }
          if (tcRef[0].Tdc) {
            shiftTimeline(tcRef[0].Board, hTmpTc1Timeline, pTc1, hTc1Timeline);
          }
          if (tcRef[0].Tdc) {
            shiftTimeline(tcRef[1].Board, hTmpTc2Timeline, pTc2, hTc2Timeline);
          }
          for (std::size_t board = 0; board < MrSync::NofChannels; ++board) {
            shiftTimeline(board, hTmpExtTimeline[board], pExt[board], hExtTimeline);
          }

          // Shift MrSync
          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            // const Long64_t preLastMrSyncTdc  = fLastMrSync [ch].Tdc;
            // const Long64_t preNextMrSyncTdc  = fNextMrSync [ch].Tdc;
            // const Long64_t preNext2MrSyncTdc = fNext2MrSync[ch].Tdc;
            fLastMrSync [ch] = fNextMrSync [ch];
            fNextMrSync [ch] = fNext2MrSync[ch];
            fNext2MrSync[ch] = { };
            // std::cout << "LastMrSync [" << ch << "]  " << preLastMrSyncTdc  << "\t->\t" << fLastMrSync [ch].Tdc << std::endl;
            // std::cout << "NextMrSync [" << ch << "]  " << preNextMrSyncTdc  << "\t->\t" << fNextMrSync [ch].Tdc << std::endl;
            // std::cout << "Next2MrSync[" << ch << "]  " << preNext2MrSyncTdc << "\t->\t" << fNext2MrSync[ch].Tdc << std::endl;
          }
          // while (std::getchar() != '\n');

          // std::cout << "[debug] Clear" << std::endl;
          // Clear
          if (spillEnded[MrSync::NofChannels - 1]) {
            std::cout << "[info] spill ended" << std::endl;
            // ClearLastSpill(true);
            // FillToSeconds(&spillEnded , false);
            break;
          }
          if (fileEnded[MrSync::NofChannels - 1]) {
            std::cout << "[info] file ended" << std::endl;
            break;
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

    void TimelineCoincidence::DrawTmpTimeline(Int_t bin, Int_t range) {
      hTmpBh1Timeline->SetLineColor(51); // pink
      hTmpBh2Timeline->SetLineColor(60); // blue
      // hTmpHodTimeline->SetLineColor();
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hTmpExtTimeline[ch]->SetLineColor(70); // cyan
      }
      hTmpTc1Timeline->SetLineColor(91); // green
      hTmpTc2Timeline->SetLineColor(99); // orange

      hTmpBh1Timeline->SetFillColorAlpha(51, 0.1);
      hTmpBh2Timeline->SetFillColorAlpha(60, 0.1);
      // hTmpHodTimeline->SetLineColor();
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hTmpExtTimeline[ch]->SetFillColorAlpha(70, 0.1);
      }
      hTmpTc1Timeline->SetFillColorAlpha(91, 0.1);
      hTmpTc2Timeline->SetFillColorAlpha(99, 0.1);

      hTmpBh1Timeline->SetFillStyle(3345);
      hTmpBh2Timeline->SetFillStyle(3454);
      // hTmpHodTimeline->SetLineColor();
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
        hTmpExtTimeline[ch]->SetFillStyle(3001);
      }
      hTmpTc1Timeline->SetFillStyle(3305);
      hTmpTc2Timeline->SetFillStyle(3395);

      hTmpBh1Timeline->Draw();
      hTmpBh2Timeline->Draw("same");
      // hTmpHodTimeline->Draw("same");
      for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) 
        hTmpExtTimeline[ch]->Draw("same");
      hTmpTc1Timeline->Draw("same");
      hTmpTc2Timeline->Draw("same");

      hTmpBh1Timeline->GetXaxis()->SetRangeUser(bin - range, bin + range);
                
      gPad->Modified();
      gPad->Update();
      gPad->WaitPrimitive();
    }
    
  }

}

#endif
