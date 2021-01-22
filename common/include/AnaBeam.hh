#ifndef Extinction_AnaBeam_hh
#define Extinction_AnaBeam_hh

#include <iostream>
#include <vector>
#include "Units.hh"
#include "Utility.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

namespace Extinction {

  namespace Analyzer {

    namespace AnaBeam {

      struct CoinOffset {
        enum {
              BH  = 0,
              Hod = 0 + BeamlineHodoscope::NofChannels,
              TC  = 0 + BeamlineHodoscope::NofChannels + Hodoscope::NofChannels,
              N   = 0 + BeamlineHodoscope::NofChannels + Hodoscope::NofChannels + TimingCounter::NofChannels,
        };
      };

      struct Results_t {
        // Extinction detector hit count
        TH1* hExtEntryByCh;

        // Extinction detector hit map
        TH2* hExtHitMap;

        // Hodoscope hit count
        TH1* hHodEntryByCh;

        // Hodoscope hit map
        TH2* hHodHitMap;

        // Extinction Detector TDC in spill
        TH1** hExtTdcInSpill; // [ExtinctionDetector::NofChannels]
        TH1*  hExtTdcInSpill_Any;

        // Extinction Detector TDC in sync
        TH1** hExtTdcInSync; // [ExtinctionDetector::NofChannels]
        TH1*  hExtTdcInSync_Any;

        // Extinction Detector Mountain Plot
        TH2** hExtMountain; // [ExtinctionDetector::NofChannels]
        TH2*  hExtMountain_Any;

        // Extinction Detector TDC difference
        TH2** hExtTdcDiff; // [ExtinctionDetector::NofChannels]
        TH2** hExtOscillo; // [ExtinctionDetector::NofChannels]

        // Hodoscope TDC in spill
        TH1** hHodTdcInSpill; // [Hodoscope::NofChannels]
        TH1*  hHodTdcInSpill_Any;

        // Hodoscope TDC in sync
        TH1** hHodTdcInSync; // [Hodoscope::NofChannels]
        TH1*  hHodTdcInSync_Any;

        // Hodoscope Mountain Plot
        TH2** hHodMountain; // [Hodoscope::NofChannels]
        TH2*  hHodMountain_Any;

        // Hodoscope TDC difference
        TH2** hHodTdcDiff; // [Hodoscope::NofChannels]
        TH2** hHodOscillo; // [Hodoscope::NofChannels]

        // Timing Counter TDC in spill
        TH1** hTcTdcInSpill; // [TimingCounter::NofChannels]

        // Timing Counter TDC in sync
        TH1** hTcTdcInSync; // [TimingCounter::NofChannels]

        // Timing Counter Mountain Plot
        TH2** hTcMountain; // [TimingCounter::NofChannels]

        // Timing Counter TDC difference
        TH1* hTcTdcDiff;
        TH1* hTcOscillo;

        // Beamline Hodoscope TDC in spill
        TH1** hBhTdcInSpill; // [BeamlineHodoscope::NofChannels]

        // Beamline Hodoscope TDC in sync
        TH1** hBhTdcInSync; // [BeamlineHodoscope::NofChannels]

        // Beamline Hodoscope Mountain Plot
        TH2** hBhMountain; // [BeamlineHodoscope::NofChannels]

        // Beamline Hodoscope TDC difference
        TH1* hBhTdcDiff;
        TH1* hBhOscillo;

        // MR Sync TDC in spill
        TH1** hMrSyncTdcInSpill; // [MrSync::NofChannels]
        TH1** hMrSyncTdcInterval1; // [MrSync::NofChannels]
        TH1** hMrSyncTdcInterval2; // [MrSync::NofChannels]
        TH2** hMrSyncTdcInterval3; // [MrSync::NofChannels]

        // Extinction Detector coincidence
        TH2** hExtTdcCoin; // [ExtinctionDetector::NofChannels]

        void Print(const std::string& ofilename) {
          TCanvas::MakeDefCanvas();
          gPad->SetGrid();

          gPad->Print((ofilename + "[").data());

          std::cout << "// Extinction Detector Hit" << std::endl;
          gPad->SetLogy(true);
          hExtEntryByCh->Draw();
          Utility::ResizeStats(hExtEntryByCh);
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          gPad->SetGrid(false, false);
          hExtHitMap->Draw("colz");
          ExtinctionDetector::CreateBorderLine()->Draw();
          gPad->Print(ofilename.data());
          gPad->SetGrid(true);
          
          std::cout << "// Hodoscope Hit" << std::endl;
          gPad->SetLogy(true);
          hHodEntryByCh->Draw();
          Utility::ResizeStats(hHodEntryByCh);
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          gPad->SetGrid(false, false);
          hHodHitMap->Draw("colz");
          Hodoscope::CreateBorderLine()->Draw();
          gPad->Print(ofilename.data());
          gPad->SetGrid(true);

          std::cout << "// Extinciton Detector TDC" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcInSpill[ch]->GetEntries()) {
              hExtTdcInSpill[ch]->Draw();
              Utility::ResizeStats(hExtTdcInSpill[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          hExtTdcInSpill_Any->Draw();
          Utility::ResizeStats(hExtTdcInSpill_Any);
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcInSync[ch]->GetEntries()) {
              hExtTdcInSync[ch]->Draw();
              Utility::ResizeStats(hExtTdcInSync[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          hExtTdcInSync_Any->Draw();
          Utility::ResizeStats(hExtTdcInSync_Any);
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          // gPad->SetLogz(true);
          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtMountain[ch]->GetEntries()) {
              hExtMountain[ch]->Draw("colz");
              gPad->Print(ofilename.data());
            }
          }
          // gPad->SetLogz(false);

          // gPad->SetLogz(true);
          hExtMountain_Any->Draw("colz");
          gPad->Print(ofilename.data());
          // gPad->SetLogz(false);

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcDiff[ch]->GetEntries()) {
              hExtTdcDiff[ch]->Draw("colz");
              Utility::ResizeStats(hExtTdcDiff[ch]);
              gPad->Print(ofilename.data());
            }
          }

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtOscillo[ch]->GetEntries()) {
              hExtOscillo[ch]->Draw("colz");
              Utility::ResizeStats(hExtOscillo[ch]);
              gPad->Print(ofilename.data());
            }
          }

          std::cout << "// Hodoscope TDC" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodTdcInSpill[ch]->GetEntries()) {
              hHodTdcInSpill[ch]->Draw();
              Utility::ResizeStats(hHodTdcInSpill[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          hHodTdcInSpill_Any->Draw();
          Utility::ResizeStats(hHodTdcInSpill_Any);
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodTdcInSync[ch]->GetEntries()) {
              hHodTdcInSync[ch]->Draw();
              Utility::ResizeStats(hHodTdcInSync[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          hHodTdcInSync_Any->Draw();
          Utility::ResizeStats(hHodTdcInSync_Any);
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          // gPad->SetLogz(true);
          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodMountain[ch]->GetEntries()) {
              hHodMountain[ch]->Draw("colz");
              gPad->Print(ofilename.data());
            }
          }
          // gPad->SetLogz(false);

          // gPad->SetLogz(true);
          hHodMountain_Any->Draw("colz");
          gPad->Print(ofilename.data());
          // gPad->SetLogz(false);

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodTdcDiff[ch]->GetEntries()) {
              hHodTdcDiff[ch]->Draw("colz");
              Utility::ResizeStats(hHodTdcDiff[ch]);
              gPad->Print(ofilename.data());
            }
          }

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodOscillo[ch]->GetEntries()) {
              hHodOscillo[ch]->Draw("colz");
              Utility::ResizeStats(hHodOscillo[ch]);
              gPad->Print(ofilename.data());
            }
          }

          std::cout << "// Timing Counter TDC" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (hTcTdcInSpill[ch]->GetEntries()) {
              hTcTdcInSpill[ch]->Draw();
              Utility::ResizeStats(hTcTdcInSpill[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (hTcTdcInSync[ch]->GetEntries()) {
              hTcTdcInSync[ch]->Draw();
              Utility::ResizeStats(hTcTdcInSync[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          if (hTcTdcDiff->GetEntries()) {
            hTcTdcDiff->Draw();
            Utility::ResizeStats(hTcTdcDiff);
            gPad->Print(ofilename.data());
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          if (hTcOscillo->GetEntries()) {
            hTcOscillo->Draw();
            Utility::ResizeStats(hTcOscillo);
            gPad->Print(ofilename.data());
          }
          gPad->SetLogy(false);

          // gPad->SetLogz(true);
          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (hTcMountain[ch]->GetEntries()) {
              hTcMountain[ch]->Draw("colz");
              gPad->Print(ofilename.data());
            }
          }
          // gPad->SetLogz(false);

          std::cout << "// Beamline Hodosocpe TDC" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (hBhTdcInSpill[ch]->GetEntries()) {
              hBhTdcInSpill[ch]->Draw();
              Utility::ResizeStats(hBhTdcInSpill[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (hBhTdcInSync[ch]->GetEntries()) {
              hBhTdcInSync[ch]->Draw();
              Utility::ResizeStats(hBhTdcInSync[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          // gPad->SetLogz(true);
          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (hBhMountain[ch]->GetEntries()) {
              hBhMountain[ch]->Draw("colz");
              gPad->Print(ofilename.data());
            }
          }
          // gPad->SetLogz(false);

          gPad->SetLogy(true);
          if (hBhTdcDiff->GetEntries()) {
            hBhTdcDiff->Draw();
            Utility::ResizeStats(hBhTdcDiff);
            gPad->Print(ofilename.data());
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          if (hBhOscillo->GetEntries()) {
            hBhOscillo->Draw();
            Utility::ResizeStats(hBhOscillo);
            gPad->Print(ofilename.data());
          }
          gPad->SetLogy(false);

          std::cout << "// MR Sync TDC" << std::endl;
          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInSpill[ch]->GetEntries()) {
              hMrSyncTdcInSpill[ch]->Draw();
              Utility::ResizeStats(hMrSyncTdcInSpill[ch]);
              gPad->Print(ofilename.data());
            }
          }

          std::cout << "// MR Sync Interval" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInterval1[ch]->GetEntries()) {
              hMrSyncTdcInterval1[ch]->Draw();
              Utility::ResizeStats(hMrSyncTdcInterval1[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInterval2[ch]->GetEntries()) {
              hMrSyncTdcInterval2[ch]->Draw();
              Utility::ResizeStats(hMrSyncTdcInterval2[ch]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInterval3[ch]->GetEntries()) {
              hMrSyncTdcInterval3[ch]->Draw("colz");
              gPad->Print(ofilename.data());
            }
          }

          std::cout << "// MPPC Coincidence" << std::endl;
          // gPad->SetLogz(true);
          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcCoin[ch]->GetEntries()) {
              hExtTdcCoin[ch]->Draw("colz");
              gPad->Print(ofilename.data());
            }
          }
          // gPad->SetLogz(false);

          gPad->Print((ofilename + "]").data());
        }

        void Write(TFile* file) {
          file->cd();

          hExtEntryByCh->Write();

          hExtHitMap->Write();
          
          hHodEntryByCh->Write();

          hHodHitMap->Write();

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcInSpill[ch]->GetEntries()) {
              hExtTdcInSpill[ch]->Write();
            }
          }

          hExtTdcInSpill_Any->Write();

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcInSync[ch]->GetEntries()) {
              hExtTdcInSync[ch]->Write();
            }
          }

          hExtTdcInSync_Any->Write();

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtMountain[ch]->GetEntries()) {
              hExtMountain[ch]->Write();
            }
          }

          hExtMountain_Any->Write();

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcDiff[ch]->GetEntries()) {
              hExtTdcDiff[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtOscillo[ch]->GetEntries()) {
              hExtOscillo[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodTdcInSpill[ch]->GetEntries()) {
              hHodTdcInSpill[ch]->Write();
            }
          }

          hHodTdcInSpill_Any->Write();

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodTdcInSync[ch]->GetEntries()) {
              hHodTdcInSync[ch]->Write();
            }
          }

          hHodTdcInSync_Any->Write();

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodMountain[ch]->GetEntries()) {
              hHodMountain[ch]->Write();
            }
          }

          hHodMountain_Any->Write();

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodTdcDiff[ch]->GetEntries()) {
              hHodTdcDiff[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (hHodOscillo[ch]->GetEntries()) {
              hHodOscillo[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (hTcTdcInSpill[ch]->GetEntries()) {
              hTcTdcInSpill[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (hTcTdcInSync[ch]->GetEntries()) {
              hTcTdcInSync[ch]->Write();
            }
          }

          if (hTcTdcDiff->GetEntries()) {
            hTcTdcDiff->Write();
          }

          if (hTcOscillo->GetEntries()) {
            hTcOscillo->Write();
          }

          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (hTcMountain[ch]->GetEntries()) {
              hTcMountain[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (hBhTdcInSpill[ch]->GetEntries()) {
              hBhTdcInSpill[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (hBhTdcInSync[ch]->GetEntries()) {
              hBhTdcInSync[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (hBhMountain[ch]->GetEntries()) {
              hBhMountain[ch]->Write();
            }
          }

          if (hBhTdcDiff->GetEntries()) {
            hBhTdcDiff->Write();
          }

          if (hBhOscillo->GetEntries()) {
            hBhOscillo->Write();
          }

          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInSpill[ch]->GetEntries()) {
              hMrSyncTdcInSpill[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInterval1[ch]->GetEntries()) {
              hMrSyncTdcInterval1[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInterval2[ch]->GetEntries()) {
              hMrSyncTdcInterval2[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
            if (hMrSyncTdcInterval3[ch]->GetEntries()) {
              hMrSyncTdcInterval3[ch]->Write();
            }
          }

          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (hExtTdcCoin[ch]->GetEntries()) {
              hExtTdcCoin[ch]->Write();
            }
          }
        }
      };

      Results_t Execute(TTree* tree,
                        ITdcDataProvider* provider) {
        std::cout << "=== Get Tree Information" << std::endl;
        const Long64_t entries = tree->GetEntries();

        std::cout << "=== Get TDC Information" << std::endl;
        const std::string tdcName = provider->GetName();
        const Double_t    timePerTdc = provider->GetTimePerTdc();

        const Double_t  xminInSpill =    0; // [msec]
        const Double_t  xmaxInSpill = 3000; // [msec]
        const Int_t    xbinsInSpill =  300;

        // const Double_t  xminInSync  = (Int_t)(-700 * nsec / timePerTdc) - 0.5; // [count]
        // const Int_t    xbinsInSync  = (Int_t)(6600 * nsec / timePerTdc) / 10; // [count/10]
        const Double_t  xminInSync  = (Int_t)(-1.2 * 5220 * nsec / timePerTdc) - 0.5; // [count]
        const Int_t    xbinsInSync  = (Int_t)(+3.4 * 5220 * nsec / timePerTdc) / 10; // [count/10]
        const Double_t  xmaxInSync  = xminInSync + xbinsInSync * 10; // [count]

        const Double_t  xminInDiff  = (Int_t)(-100 * nsec / timePerTdc) - 0.5; // [count]
        const Int_t    xbinsInDiff  = (Int_t)( 300 * nsec / timePerTdc); // [count]
        const Double_t  xmaxInDiff  = xminInDiff + xbinsInDiff; // [count]

        const Double_t meanMrSync   = (5220 * nsec / timePerTdc); // [count]
        const Int_t    xbinsInt1    = 104;
        const Double_t  xminInt1    = - 0.2; // [MrSync]
        const Double_t  xmaxInt1    = +10.2; // [MrSync]
        const Int_t    xbinsInt2    = (Int_t)(375 * nsec / timePerTdc); // [count]
        const Double_t  xminInt2    = (Int_t)meanMrSync - xbinsInt2 / 2 - 0.5; // [count]
        const Double_t  xmaxInt2    = xminInt2 + xbinsInt2; // [count]

        std::cout << "=== Initialzie History" << std::endl;
        const Double_t lastThreshold = 100 * nsec;
        std::vector<Extinction::TdcData> lastExtData;
        std::vector<Extinction::TdcData> lastHodData;
        std::vector<Extinction::TdcData> lastTcData;
        std::vector<Extinction::TdcData> lastBhData;
        std::vector<Extinction::TdcData> lastMrSyncData;
        std::map<std::size_t, Extinction::TdcData> recExtData;
        std::map<std::size_t, Extinction::TdcData> recHodData;
        std::map<std::size_t, Extinction::TdcData> recTcData;
        std::map<std::size_t, Extinction::TdcData> recBhData;

        std::cout << "=== Set Branch Address" << std::endl;
        provider->SetBranchAddress(tree);

        // Extinction detector hit count
        TH1* hExtEntryByCh = new TH1D("hExtEntryByCh",
                                      Form("%s, Extinction Detector Entries by Channel;"
                                           "Channel;"
                                           "Entries", tdcName.data()),
                                      ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5);

        // Extinction detector hit map
        TH2* hExtHitMap = Extinction::ExtinctionDetector::CreateHitMap("hExtHitMap");

        // Hodoscope hit count
        TH1* hHodEntryByCh = new TH1D("hHodEntryByCh",
                                      Form("%s, Hodoscope Entries by Channel;"
                                           "Channel;"
                                           "Entries", tdcName.data()),
                                      Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5);

        // Hodoscope hit map
        TH2* hHodHitMap = Extinction::Hodoscope::CreateHitMap("hHodHitMap");

        // Extinction Detector TDC in spill
        TH1** hExtTdcInSpill = new TH1*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcInSpill[ch] = new TH1D(Form("hExtTdcInSpill_%03ld", ch),
                                        Form("%s, Extinction Detector TDC in Spill @ ch%ld;"
                                             "Time [ms]", tdcName.data(), ch),
                                        xbinsInSpill, xminInSpill, xmaxInSpill);
        }
        TH1* hExtTdcInSpill_Any = new TH1D("hExtTdcInSpill_Any",
                                           Form("%s, Extinction Detector TDC in Spill;"
                                                "Time [ms]", tdcName.data()),
                                           xbinsInSpill, xminInSpill, xmaxInSpill);

        // Extinction Detector TDC in sync
        TH1** hExtTdcInSync = new TH1*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcInSync[ch] = new TH1D(Form("hExtTdcInSync_%03ld", ch),
                                       Form("%s, Extinction Detector TDC in MR Sync @ ch%ld;"
                                            "TDC [count]", tdcName.data(), ch),
                                       xbinsInSync, xminInSync, xmaxInSync);
        }
        TH1* hExtTdcInSync_Any = new TH1D("hExtTdcInSync_Any",
                                          Form("%s, Extinction Detector TDC in MR Sync;"
                                               "TDC [count]", tdcName.data()),
                                          xbinsInSync, xminInSync, xmaxInSync);

        // Extinction Detector Mountain Plot
        TH2** hExtMountain = new TH2*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtMountain[ch] = new TH2D(Form("hExtMountain_%03ld", ch),
                                      Form("%s, Extinction Detector Mountain Plot @ ch%ld;"
                                           "TDC [count];"
                                           "Time [ms]", tdcName.data(), ch),
                                      xbinsInSync / 2, xminInSync, xmaxInSync,
                                      xbinsInSpill / 2, xminInSpill, xmaxInSpill);
          hExtMountain[ch]->SetStats(false);
        }
        TH2* hExtMountain_Any = new TH2D("hExtMountain_Any",
                                         Form("%s, Extinction Detector Mountain Plot;"
                                              "TDC [count];"
                                              "Time [ms]", tdcName.data()),
                                         xbinsInSync / 2, xminInSync, xmaxInSync,
                                         xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hExtMountain_Any->SetStats(false);

        // Extinction Detector TDC difference
        TH2** hExtTdcDiff = new TH2*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcDiff[ch] = new TH2D(Form("hExtTdcDiff%ld", ch),
                                     Form("%s, Extinction Detector TDC Difference @ ch%ld;"
                                          "Channel;"
                                          "TDC [count]", tdcName.data(), ch),
                                     ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5,
                                     xbinsInDiff, xminInDiff, xmaxInDiff);
          hExtTdcDiff[ch]->SetStats(false);
        }

        TH2** hExtOscillo = new TH2*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtOscillo[ch] = new TH2D(Form("hExtOscillo%ld", ch),
                                     Form("%s, Extinction Detector TDC Difference @ ch%ld;"
                                          "Channel;"
                                          "TDC [count]", tdcName.data(), ch),
                                     ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5,
                                     xbinsInDiff, xminInDiff, xmaxInDiff);
          hExtOscillo[ch]->SetStats(false);
        }

        // Hodoscope TDC in spill
        TH1** hHodTdcInSpill = new TH1*[Hodoscope::NofChannels];
        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodTdcInSpill[ch] = new TH1D(Form("hHodTdcInSpill_%03ld", ch),
                                        Form("%s, Hodoscope TDC in Spill @ ch%ld;"
                                             "Time [ms]", tdcName.data(), ch),
                                        xbinsInSpill, xminInSpill, xmaxInSpill);
        }
        TH1* hHodTdcInSpill_Any = new TH1D("hHodTdcInSpill_Any",
                                           Form("%s, Hodoscope TDC in Spill;"
                                                "Time [ms]", tdcName.data()),
                                           xbinsInSpill, xminInSpill, xmaxInSpill);

        // Hodoscope TDC in sync
        TH1** hHodTdcInSync = new TH1*[Hodoscope::NofChannels];
        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodTdcInSync[ch] = new TH1D(Form("hHodTdcInSync_%03ld", ch),
                                       Form("%s, Hodoscope TDC in MR Sync @ ch%ld;"
                                            "TDC [count]", tdcName.data(), ch),
                                       xbinsInSync, xminInSync, xmaxInSync);
        }
        TH1* hHodTdcInSync_Any = new TH1D("hHodTdcInSync_Any",
                                          Form("%s, Hodoscope TDC in MR Sync;"
                                               "TDC [count]", tdcName.data()),
                                          xbinsInSync, xminInSync, xmaxInSync);

        // Hodoscope Mountain Plot
        TH2** hHodMountain = new TH2*[Hodoscope::NofChannels];
        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodMountain[ch] = new TH2D(Form("hHodMountain_%03ld", ch),
                                      Form("%s, Hodoscope Mountain Plot @ ch%ld;"
                                           "TDC [count];"
                                           "Time [ms]", tdcName.data(), ch),
                                      xbinsInSync / 2, xminInSync, xmaxInSync,
                                      xbinsInSpill / 2, xminInSpill, xmaxInSpill);
          hHodMountain[ch]->SetStats(false);
        }
        TH2* hHodMountain_Any = new TH2D("hHodMountain_Any",
                                         Form("%s, Hodoscope Mountain Plot;"
                                              "TDC [count];"
                                              "Time [ms]", tdcName.data()),
                                         xbinsInSync / 2, xminInSync, xmaxInSync,
                                         xbinsInSpill / 2, xminInSpill, xmaxInSpill);
        hHodMountain_Any->SetStats(false);

        // Hodoscope TDC difference
        TH2** hHodTdcDiff = new TH2*[Hodoscope::NofChannels];
        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodTdcDiff[ch] = new TH2D(Form("hHodTdcDiff_%03ld", ch),
                                     Form("%s, Hodoscope TDC Difference @ ch%ld;"
                                          "Channel;TDC [count]", tdcName.data(), ch),
                                     Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5,
                                     xbinsInDiff, xminInDiff, xmaxInDiff);
          hHodTdcDiff[ch]->SetStats(false);
        }

        TH2** hHodOscillo = new TH2*[Hodoscope::NofChannels];
        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          hHodOscillo[ch] = new TH2D(Form("hHodOscillo_%03ld", ch),
                                     Form("%s, Hodoscope TDC Difference @ ch%ld;"
                                          "Channel;TDC [count]", tdcName.data(), ch),
                                     Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5,
                                     xbinsInDiff, xminInDiff, xmaxInDiff);
          hHodOscillo[ch]->SetStats(false);
        }

        // Timing Counter TDC in spill
        TH1** hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcTdcInSpill[ch] = new TH1D(Form("hTcTdcInSpill_%03ld", ch),
                                       Form("%s, Timing Counter %ld TDC in Spill;"
                                            "Time [ms]", tdcName.data(), ch + 1),
                                       xbinsInSpill, xminInSpill, xmaxInSpill);
        }

        // Timing Counter TDC in sync
        TH1** hTcTdcInSync = new TH1*[TimingCounter::NofChannels];
        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcTdcInSync[ch] = new TH1D(Form("hTcTdcInSync_%03ld", ch),
                                      Form("%s, Timing Counter %ld TDC in MR Sync;"
                                           "TDC [Count]", tdcName.data(), ch + 1),
                                      xbinsInSync, xminInSync, xmaxInSync);
        }

        // Timing Counter Mountain Plot
        TH2** hTcMountain = new TH2*[TimingCounter::NofChannels];
        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          hTcMountain[ch] = new TH2D(Form("hTcMountain_%03ld", ch),
                                     Form("%s, Timing Counter Mountain Plot @ ch%ld;"
                                          "TDC [count];"
                                          "Time [ms]", tdcName.data(), ch),
                                     xbinsInSync / 2, xminInSync, xmaxInSync,
                                     xbinsInSpill / 2, xminInSpill, xmaxInSpill);
          hTcMountain[ch]->SetStats(false);
        }

        // Timing Counter TDC difference
        TH1* hTcTdcDiff = new TH1D("hTcTdcDiff",
                                   Form("%s, Timing Counter TDC Difference;"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInDiff, xminInDiff, xmaxInDiff);

        TH1* hTcOscillo = new TH1D("hTcOscillo",
                                   Form("%s, Timing Counter TDC Difference;"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInDiff, xminInDiff, xmaxInDiff);

        // Beamline Hodoscope TDC in spill
        TH1** hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhTdcInSpill[ch] = new TH1D(Form("hBhTdcInSpill_%03ld", ch),
                                       Form("%s, BH%ld TDC in Spill;"
                                            "Time [ms]", tdcName.data(), ch + 1),
                                       xbinsInSpill, xminInSpill, xmaxInSpill);
        }

        // Beamline Hodoscope TDC in sync
        TH1** hBhTdcInSync = new TH1*[BeamlineHodoscope::NofChannels];
        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhTdcInSync[ch] = new TH1D(Form("hBhTdcInSync_%03ld", ch),
                                      Form("%s, BH%ld TDC in MR Sync;"
                                           "TDC [Count]", tdcName.data(), ch + 1),
                                      xbinsInSync, xminInSync, xmaxInSync);
        }

        // Beamline Hodoscope Mountain Plot
        TH2** hBhMountain = new TH2*[BeamlineHodoscope::NofChannels];
        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          hBhMountain[ch] = new TH2D(Form("hBhMountain_%03ld", ch),
                                     Form("%s, BH%ld Mountain Plot;"
                                          "TDC [count];"
                                          "Time [ms]", tdcName.data(), ch + 1),
                                     xbinsInSync / 2, xminInSync, xmaxInSync,
                                     xbinsInSpill / 2, xminInSpill, xmaxInSpill);
          hBhMountain[ch]->SetStats(false);
        }

        // Beamline Hodoscope TDC difference
        TH1* hBhTdcDiff = new TH1D("hBhTdcDiff",
                                   Form("%s, BH TDC Difference;"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInDiff, xminInDiff, xmaxInDiff);

        TH1* hBhOscillo = new TH1D("hBhOscillo",
                                   Form("%s, BH TDC Difference;"
                                        "TDC [count]", tdcName.data()),
                                   xbinsInDiff, xminInDiff, xmaxInDiff);

        // MR Sync TDC in spill
        TH1** hMrSyncTdcInSpill = new TH1*[MrSync::NofChannels];
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncTdcInSpill[ch] = new TH1D(Form("hMrSyncTdcInSpill_%03ld", ch),
                                           Form("%s, MR Sync TDC @ board%ld;"
                                                "Time [ms];"
                                                "Entries", tdcName.data(), ch + 1),
                                           xbinsInSpill, xminInSpill, xmaxInSpill);
        }

        // MR Sync interval
        TH1** hMrSyncTdcInterval1 = new TH1*[MrSync::NofChannels];
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncTdcInterval1[ch] = new TH1D(Form("hMrSyncTdcInterval1_%03ld", ch),
                                             Form("%s, MR Sync Interval @ board%ld;"
                                                  "Time [MR Sync];"
                                                  "Entries", tdcName.data(), ch + 1),
                                             xbinsInt1, xminInt1, xmaxInt1);
        }

        TH1** hMrSyncTdcInterval2 = new TH1*[MrSync::NofChannels];
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncTdcInterval2[ch] = new TH1D(Form("hMrSyncTdcInterval2_%03ld", ch),
                                             Form("%s, MR Sync Interval @ board%ld;"
                                                  "TDC [count];"
                                                  "Entries", tdcName.data(), ch + 1),
                                             xbinsInt2, xminInt2, xmaxInt2);
        }

        TH2** hMrSyncTdcInterval3 = new TH2*[MrSync::NofChannels];
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hMrSyncTdcInterval3[ch] = new TH2D(Form("hMrSyncTdcInterval3_%03ld", ch),
                                             Form("%s, MR Sync Interval @ board%ld;"
                                                  "TDC [count];"
                                                  "Time [ms];"
                                                  "Entries", tdcName.data(), ch + 1),
                                             xbinsInt2, xminInt2, xmaxInt2,
                                             xbinsInSpill, xminInSpill, xmaxInSpill);
        }

        // Extinction Detector coincidence
        TH2** hExtTdcCoin = new TH2*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          hExtTdcCoin[ch] = new TH2D(Form("hExtTdcCoin_%03ld", ch),
                                     Form("%s, Extinction Detector Coincidence @ ch%ld;"
                                          "Channel;"
                                          "TDC [count]", tdcName.data(), ch),
                                     CoinOffset::N, 0 - 0.5, CoinOffset::N - 0.5,
                                     xbinsInDiff, xminInDiff, xmaxInDiff);
          hExtTdcCoin[ch]->SetStats(false);
        }

        std::cout << "=== Get Entry " << std::endl;
        Int_t targetSpill = -1;
        Int_t   lastSpill = -1;
        for (Long64_t entry = 0; entry < entries; ++entry) {
          if (entry % 1000000 == 0) {
            std::cout << ">> " << entry // << std::endl;
                      << " ("
                      << lastExtData.size() << " "
                      << lastHodData.size() << " "
                      << lastTcData.size() << " "
                      << lastBhData.size() << " "
                      << lastMrSyncData.size() << ")"
                      << std::endl;
          }

          tree->GetEntry(entry);

          if (!provider->IsData()) {
            continue;
          }

          if (provider->GetSpill() != lastSpill) {
            lastExtData.clear();
            lastHodData.clear();
            lastTcData.clear();
            lastBhData.clear();
            lastMrSyncData.clear();
          }

          std::vector<TdcData> tdcData = provider->GetTdcData();

          for (auto&& data : tdcData) {
            const Int_t    globalChannel = data.Channel;
            const Long64_t tdc           = data.Tdc;
            const Double_t time          = data.Time;

            if (ExtinctionDetector::Contains(globalChannel)) {
              const Int_t ch = ExtinctionDetector::GetChannel(globalChannel);

              hExtEntryByCh     ->Fill(ch);
              ExtinctionDetector::Fill(hExtHitMap, ch);
              hExtTdcInSpill[ch]->Fill(time / Extinction::msec);
              hExtTdcInSpill_Any->Fill(time / Extinction::msec);
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                auto lastData = lastMrSyncData[i];
                if (lastData.Channel == data.MrSyncChannel) {
                  auto syncTdc = lastData.Tdc;
                  hExtTdcInSync[ch]->Fill(tdc - syncTdc);
                  hExtTdcInSync_Any->Fill(tdc - syncTdc);
                  hExtMountain [ch]->Fill(tdc - syncTdc, time / Extinction::msec);
                  hExtMountain_Any ->Fill(tdc - syncTdc, time / Extinction::msec);
                  break;
                }                
              }

              Bool_t containsThis = false; 
              for (auto&& lastData : lastExtData) {
                auto lastCh  = ExtinctionDetector::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                hExtTdcDiff[lastCh]->Fill(    ch,     tdc - lastTdc);
                hExtTdcDiff[    ch]->Fill(lastCh, lastTdc -     tdc);
                if (lastData.Channel == globalChannel) {
                  containsThis = true;
                }
              }
              for (auto&& lastPair : recExtData) {
                auto lastData = lastPair.second;
                auto lastCh   = ExtinctionDetector::GetChannel(lastData.Channel);
                auto lastTdc  = lastData.Tdc;
                  hExtOscillo[lastCh]->Fill(    ch,     tdc - lastTdc);
                  hExtOscillo[    ch]->Fill(lastCh, lastTdc -     tdc);
              }
              if (containsThis == false) {
                recExtData[globalChannel] = data;
              }
              
              for (auto&& lastData : lastHodData) {
                auto lastCh  = Hodoscope::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                hExtTdcCoin[ch]->Fill(lastCh + CoinOffset::Hod, lastTdc - tdc);
              }
              for (auto&& lastData : lastTcData) {
                auto lastCh  = TimingCounter::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                hExtTdcCoin[ch]->Fill(lastCh + CoinOffset::TC, lastTdc - tdc);
              }
              for (auto&& lastData : lastBhData) {
                auto lastCh  = BeamlineHodoscope::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                hExtTdcCoin[ch]->Fill(lastCh + CoinOffset::BH, lastTdc - tdc);
              }

              for (std::size_t i = 0, n = lastExtData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastExtData[0].Time) > lastThreshold) {
                  lastExtData.erase(lastExtData.begin());
                } else {
                  break;
                }
              }
              lastExtData.push_back(data);
              if (lastExtData.size() > 10000) {
                std::cerr << "[error] size of lastExtData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (Hodoscope::Contains(globalChannel)) {
              const Int_t ch = Hodoscope::GetChannel(globalChannel);

              hHodEntryByCh     ->Fill(ch);
              Hodoscope         ::Fill(hHodHitMap, ch);
              hHodTdcInSpill[ch]->Fill(time / Extinction::msec);
              hHodTdcInSpill_Any->Fill(time / Extinction::msec);
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                auto lastData = lastMrSyncData[i];
                if (lastData.Channel == data.MrSyncChannel) {
                  auto syncTdc = lastData.Tdc;
                  hHodTdcInSync[ch]->Fill(tdc - syncTdc);
                  hHodTdcInSync_Any->Fill(tdc - syncTdc);
                  hHodMountain [ch]->Fill(tdc - syncTdc, time / Extinction::msec);
                  hHodMountain_Any ->Fill(tdc - syncTdc, time / Extinction::msec);
                  break;
                }
              }

              Bool_t containsThis = false;
              for (auto&& lastData : lastHodData) {
                auto lastCh  = Hodoscope::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                hHodTdcDiff[lastCh]->Fill(    ch,     tdc - lastTdc);
                hHodTdcDiff[    ch]->Fill(lastCh, lastTdc -     tdc);
                if (lastData.Channel == globalChannel) {
                  containsThis = true;
                }
              }
              for (auto&& lastPair : recHodData) {
                auto lastData = lastPair.second;
                auto lastCh   = Hodoscope::GetChannel(lastData.Channel);
                auto lastTdc  = lastData.Tdc;
                  hHodOscillo[lastCh]->Fill(    ch,     tdc - lastTdc);
                  hHodOscillo[    ch]->Fill(lastCh, lastTdc -     tdc);
              }
              if (containsThis == false) {
                recHodData[globalChannel] = data;
              }

              for (auto&& lastData : lastExtData) {
                auto lastCh  = ExtinctionDetector::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                hExtTdcCoin[lastCh]->Fill(ch + CoinOffset::Hod, tdc - lastTdc);
              }

              for (std::size_t i = 0, n = lastHodData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastHodData[0].Time) > lastThreshold) {
                  lastHodData.erase(lastHodData.begin());
                } else {
                  break;
                }
              }
              lastHodData.push_back(data);
              if (lastHodData.size() > 10000) {
                std::cerr << "[error] size of lastHodData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (TimingCounter::Contains(globalChannel)) {
              const Int_t ch = TimingCounter::GetChannel(globalChannel);

              hTcTdcInSpill[ch]->Fill(time / Extinction::msec);
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                auto lastData = lastMrSyncData[i];
                if (lastData.Channel == data.MrSyncChannel) {
                  auto syncTdc = lastData.Tdc;
                  hTcTdcInSync[ch]->Fill(tdc - syncTdc);
                  hTcMountain [ch]->Fill(tdc - syncTdc, time / Extinction::msec);
                  break;
                }
              }

              Bool_t containsThis = false;
              for (auto&& lastData : lastTcData) {
                auto lastCh  = TimingCounter::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                if        (ch > lastCh) {
                  hTcTdcDiff->Fill(lastTdc -     tdc);
                } else if (ch < lastCh) {
                  hTcTdcDiff->Fill(    tdc - lastTdc);
                }
                if (lastData.Channel == globalChannel) {
                  containsThis = true;
                }
              }
              for (auto&& lastPair : recTcData) {
                auto lastData = lastPair.second;
                auto lastCh   = TimingCounter::GetChannel(lastData.Channel);
                auto lastTdc  = lastData.Tdc;
                if        (ch > lastCh) {
                  hTcOscillo->Fill(lastTdc -     tdc);
                } else if (ch < lastCh) {
                  hTcOscillo->Fill(    tdc - lastTdc);
                }
              }
              if (containsThis == false) {
                recTcData[globalChannel] = data;
              }

              for (auto&& lastData : lastExtData) {
                auto lastCh   = ExtinctionDetector::GetChannel(lastData.Channel);
                auto lastTdc  = lastData.Tdc;
                hExtTdcCoin[lastCh]->Fill(ch + CoinOffset::TC, tdc - lastTdc);
              }

              for (std::size_t i = 0, n = lastTcData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastTcData[0].Time) > lastThreshold) {
                  lastTcData.erase(lastTcData.begin());
                } else {
                  break;
                }
              }
              lastTcData.push_back(data);
              if (lastTcData.size() > 10000) {
                std::cerr << "[error] size of lastTcData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (BeamlineHodoscope::Contains(globalChannel)) {
              const Int_t ch = BeamlineHodoscope::GetChannel(globalChannel);

              hBhTdcInSpill[ch]->Fill(time / Extinction::msec);
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                auto lastData = lastMrSyncData[i];
                if (lastData.Channel == data.MrSyncChannel) {
                  auto syncTdc = lastData.Tdc;
                  hBhTdcInSync[ch]->Fill(tdc - syncTdc);
                  hBhMountain [ch]->Fill(tdc - syncTdc, time / Extinction::msec);
                  break;
                }
              }

              Bool_t containsThis = false;
              for (auto&& lastData : lastBhData) {
                auto lastCh  = BeamlineHodoscope::GetChannel(lastData.Channel);
                auto lastTdc = lastData.Tdc;
                if        (ch > lastCh) {
                  hBhTdcDiff->Fill(lastTdc -     tdc);
                } else if (ch < lastCh) {
                  hBhTdcDiff->Fill(    tdc - lastTdc);
                }
                if (lastData.Channel == globalChannel) {
                  containsThis = true;
                }
              }
              for (auto&& lastPair : recBhData) {
                auto lastData = lastPair.second;
                auto lastCh   = BeamlineHodoscope::GetChannel(lastData.Channel);
                auto lastTdc  = lastData.Tdc;
                if        (ch > lastCh) {
                  hBhOscillo->Fill(lastTdc -     tdc);
                } else if (ch < lastCh) {
                  hBhOscillo->Fill(    tdc - lastTdc);
                }
              }
              if (containsThis == false) {
                recBhData[globalChannel] = data;
              }

              for (auto&& lastData : lastExtData) {
                auto lastCh   = ExtinctionDetector::GetChannel(lastData.Channel);
                auto lastTdc  = lastData.Tdc;
                hExtTdcCoin[lastCh]->Fill(ch + CoinOffset::BH, tdc - lastTdc);
              }

              for (std::size_t i = 0, n = lastBhData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastBhData[0].Time) > lastThreshold) {
                  lastBhData.erase(lastBhData.begin());
                } else {
                  break;
                }
              }
              lastBhData.push_back(data);
              if (lastBhData.size() > 10000) {
                std::cerr << "[error] size of lastBhData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (MrSync::Contains(globalChannel)) {
              const Int_t ch = MrSync::GetChannel(globalChannel);

              if (targetSpill >= 0) {
                if (data.Spill == targetSpill) {
                  hMrSyncTdcInSpill[ch]->Fill(time / Extinction::msec);
                } else if (data.Spill > targetSpill && !hMrSyncTdcInSpill[ch]->GetEntries()) {
                  targetSpill = data.Spill;
                  hMrSyncTdcInSpill[ch]->Fill(time / Extinction::msec);
                }
              } else {
                targetSpill = data.Spill + 1;
              }
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                if (lastMrSyncData[i - 1].Channel == globalChannel) {
                  const Long64_t dtdc = tdc - lastMrSyncData[i - 1].Tdc;
                  hMrSyncTdcInterval1[ch]->Fill(dtdc / meanMrSync);
                  hMrSyncTdcInterval2[ch]->Fill(dtdc);
                  hMrSyncTdcInterval3[ch]->Fill(dtdc, time / Extinction::msec);
                  break;
                }
              }

              Bool_t found = false;
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                if (lastMrSyncData[i - 1].Channel == globalChannel) {
                  if (found) {
                    lastMrSyncData.erase(lastMrSyncData.begin() + i - 1);
                  } else {
                    // Save last one
                    found = true;
                  }
                }
              }
              lastMrSyncData.push_back(data);
              if (lastMrSyncData.size() > 10000) {
                std::cerr << "[error] size of lastMrSyncData reaches 10000" << std::endl;
                exit(1);
              }
            }

          }

          lastSpill = provider->GetSpill();
        }

        return
          {
           hExtEntryByCh,
           hExtHitMap,
           hHodEntryByCh,
           hHodHitMap,
           hExtTdcInSpill,
           hExtTdcInSpill_Any,
           hExtTdcInSync,
           hExtTdcInSync_Any,
           hExtMountain,
           hExtMountain_Any,
           hExtTdcDiff,
           hExtOscillo,
           hHodTdcInSpill,
           hHodTdcInSpill_Any,
           hHodTdcInSync,
           hHodTdcInSync_Any,
           hHodMountain,
           hHodMountain_Any,
           hHodTdcDiff,
           hHodOscillo,
           hTcTdcInSpill,
           hTcTdcInSync,
           hTcMountain,
           hTcTdcDiff,
           hTcOscillo,
           hBhTdcInSpill,
           hBhTdcInSync,
           hBhMountain,
           hBhTdcDiff,
           hBhOscillo,
           hMrSyncTdcInSpill,
           hMrSyncTdcInterval1,
           hMrSyncTdcInterval2,
           hMrSyncTdcInterval3,
           hExtTdcCoin,
          };
      }

    }

  }

}

#endif
