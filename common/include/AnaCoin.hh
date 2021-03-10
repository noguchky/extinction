#ifndef Extinction_AnaCoin_hh
#define Extinction_AnaCoin_hh

#include <iostream>
#include <vector>
#include "Units.hh"
#include "Utility.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "AnaBeam.hh"
#include "AnaTimeOffset.hh"

namespace Extinction {

  namespace Analyzer {

    namespace AnaCoin {

      using CoinOffset = AnaTimeOffset::CoinOffset;
      using CoinInfo   = AnaTimeOffset::CoinInfo;

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

        // Extinction Detector Mountain Plot by Spill 
        std::vector<TH2**> hExtMountains; // [NofSpills][ExtinctionDetector::NofChannels]

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
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          gPad->SetGrid(false, false);
          hHodHitMap->Draw("colz");
          Hodoscope::CreateBorderLine()->Draw();
          gPad->Print(ofilename.data());
          gPad->SetGrid(true);

          std::cout << "// Extinction Detector TDC" << std::endl;
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

          gPad->Print((ofilename + "]").data());
        }

        void PrintMtns(const std::string& ofilename) {
          gPad->Print((ofilename + "[").data());

          for (auto&& mtns : hExtMountains) {
            for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
              if (mtns[ch]->GetEntries()) {
                mtns[ch]->Draw("colz");
                gPad->Print(ofilename.data());
              }
            }
          }

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
        }
      };

      Results_t Execute(TTree* tree,
                        ITdcDataProvider* provider,
                        std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, CoinInfo>> coinInfo) {
        std::cout << "=== Get Tree Information" << std::endl;
        const Long64_t entries = tree->GetEntries();

        std::cout << "=== Get TDC Information" << std::endl;
        const std::string tdcName = provider->GetName();
        const Double_t    timePerTdc = provider->GetTimePerTdc();

        const Double_t  xminInSpill =    0; // [msec]
        const Double_t  xmaxInSpill = 3000; // [msec]
        const Int_t    xbinsInSpill =  300;

        const Double_t  xminInSync  = (Int_t)(-700 * nsec / timePerTdc) - 0.5; // [count]
        const Int_t    xbinsInSync  = (Int_t)(6600 * nsec / timePerTdc) / 10; // [count/10]
        const Double_t  xmaxInSync  = xminInSync + xbinsInSync * 10; // [count]

        std::cout << "=== Initialzie Coincidence Information" << std::endl;
        std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Bool_t>> contains;
        for (std::size_t extCh = 0; extCh < ExtinctionDetector::NofChannels; ++extCh) {
          for (std::size_t index = 0; index < CoinOffset::N; ++index) {
            contains[extCh][index] = false;
          }
        }
        for (auto&& info : coinInfo) {
          const std::size_t extCh = info.first;
          for (auto&& subInfo : info.second) {
            const std::size_t index = subInfo.first;
            contains[extCh][index] = true;
          }
        }

        std::cout << "=== Initialzie History" << std::endl;
        const Double_t lastThreshold = 100 * nsec;
        std::vector<Extinction::TdcData> lastExtData;
        std::vector<Extinction::TdcData> lastHodData;
        std::vector<Extinction::TdcData> lastTcData;
        std::vector<Extinction::TdcData> lastBhData;
        std::vector<Extinction::TdcData> lastMrSyncData;

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
          hExtTdcInSpill[ch] = new TH1D(Form("hExtTdcInSpill%ld", ch),
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
          hExtTdcInSync[ch] = new TH1D(Form("hExtTdcInSync%ld", ch),
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
          hExtMountain[ch] = new TH2D(Form("hExtMountain%ld", ch),
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

        std::vector<TH2**> hExtMountains;
        
        std::cout << "=== Create Filler" << std::endl;
        auto fillCoin =
          [&](TdcData extData) {
            const std::size_t extCh = ExtinctionDetector::GetChannel(extData.Channel);
            const Long64_t    tdc   = extData.Tdc;
            const Double_t    time  = extData.Time;
            Bool_t coincidence[CoinOffset::N];
            std::vector<std::size_t> coinHodChs;

            for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
              const std::size_t i = bhCh + CoinOffset::BH;
              coincidence[i] = !contains[extCh][i];
            }
            for (auto&& lastData : lastBhData) {
              const std::size_t bhCh = BeamlineHodoscope::GetChannel(lastData.Channel);
              const std::size_t i = bhCh + CoinOffset::BH;
              if (!coincidence[i]) {
                const Double_t dt   = lastData.Time - time;
                const Double_t mean = coinInfo[extCh][i].FitMean * timePerTdc;
                const Double_t sigma = 25.0 * nsec;
                if (std::abs(dt - mean) < sigma) {
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
              coincidence[i] = !contains[extCh][i];
            }
            for (auto&& lastData : lastHodData) {
              if (lastData.Channel != Hodoscope::GlobalChannelOffset) {
                continue;
              } // TBD: for 20201217
              const std::size_t hodCh = 0;
              const std::size_t i = hodCh + CoinOffset::Hod;
              if (!coincidence[i]) {
                const Double_t dt   = lastData.Time - time;
                const Double_t mean = coinInfo[extCh][i].FitMean * timePerTdc;
                const Double_t sigma = 25.0 * nsec;
                if (std::abs(dt - mean) < sigma) {
                  coincidence[i] = true;
                  coinHodChs.push_back(Hodoscope::GetChannel(lastData.Channel));
                  break;
                }
              }
            }

            for (std::size_t tcCh = 0; tcCh < TimingCounter::NofChannels; ++tcCh) {
              const std::size_t i = tcCh + CoinOffset::TC;
              coincidence[i] = !contains[extCh][i];
            }
            for (auto&& lastData : lastTcData) {
              const std::size_t tcCh = TimingCounter::GetChannel(lastData.Channel);
              const std::size_t i = tcCh + CoinOffset::TC;
              if (!coincidence[i]) {
                const Double_t dt   = lastData.Time - time;
                const Double_t mean = coinInfo[extCh][i].FitMean * timePerTdc;
                const Double_t sigma = 25.0 * nsec;
                if (std::abs(dt - mean) < sigma) {
                  coincidence[i] = true;
                  if (std::all_of(coincidence + CoinOffset::TC, coincidence + CoinOffset::TC + TimingCounter::NofChannels, [](Bool_t b) { return b; })) {
                    break;
                  }
                }
              }
            }

            if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
              // std::cout << "contains: " << std::endl;
              // for (auto&& pair1 : contains) {
              //   for (auto&& pair2 : pair1.second) {
              //     std::cout << pair1.first << "\t"
              //               << pair2.first << "\t"
              //               << (pair2.second ? "yes" : "no") << "\t"
              //               << (pair2.second ? coinInfo[pair1.first][pair2.first].FitMean / nsec : 0.0) << "\t"
              //               << std::endl;
              //   }
              // }
              // std::cout << "extData: " << std::endl;
              // std::cout << extData.Time / nsec << std::endl;
              // std::cout << "bhData: " << std::endl;
              // for (auto&& lastData : lastBhData) {
              //   std::cout << lastData.Channel << "\t"
              //             << (lastData.Time - extData.Time) / nsec << std::endl;
              // }
              // std::cout << "hodData: " << std::endl;
              // for (auto&& lastData : lastHodData) {
              //   std::cout << lastData.Channel << "\t"
              //             << (lastData.Time - extData.Time) / nsec << std::endl;
              // }
              // std::cout << "tcData: " << std::endl;
              // for (auto&& lastData : lastTcData) {
              //   std::cout << lastData.Channel << "\t"
              //             << (lastData.Time - extData.Time) / nsec << std::endl;
              // }
              // std::cout << "mrSyncData: " << std::endl;
              // for (auto&& lastData : lastMrSyncData) {
              //   std::cout << lastData.Channel << "\t"
              //             << (lastData.Time - extData.Time) / nsec << std::endl;
              // }
              // std::cout << "> " << std::endl;
              // while (std::getchar() != '\n');

              hExtEntryByCh     ->Fill(extCh);
              ExtinctionDetector::Fill(hExtHitMap, extCh);
              for (auto&& hodCh : coinHodChs) {
                hHodEntryByCh->Fill(hodCh);
                Hodoscope    ::Fill(hHodHitMap, hodCh);
              }
              hExtTdcInSpill[extCh]->Fill(time / Extinction::msec);
              hExtTdcInSpill_Any   ->Fill(time / Extinction::msec);
              for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
                auto lastData = lastMrSyncData[i];
                if (lastData.Channel == extData.MrSyncChannel) {
                  auto syncTdc = lastData.Tdc;
                  hExtTdcInSync[extCh]->Fill(tdc - syncTdc);
                  hExtTdcInSync_Any   ->Fill(tdc - syncTdc);
                  hExtMountain [extCh]->Fill(tdc - syncTdc, time / Extinction::msec);
                  hExtMountain_Any    ->Fill(tdc - syncTdc, time / Extinction::msec);
                  if (hExtMountains.size()) {
                    hExtMountains.back()[extCh]->Fill(tdc - syncTdc, time / Extinction::msec);
                  }
                  break;
                }
              }
            }
          };

        std::cout << "=== Get Entry " << std::endl;
        Int_t lastSpill = -1;
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

            hExtMountains.push_back(new TH2*[ExtinctionDetector::NofChannels]);
            for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
              hExtMountains.back()[ch] = new TH2D(Form("hExtMountains%lu_%ld", hExtMountains.size(), ch),
                                                  Form("%s, Extinction Detector Mountain Plot @ spill %d, ch%ld;"
                                                       "TDC [count];"
                                                       "Time [ms]", tdcName.data(), provider->GetSpill(), ch),
                                                  xbinsInSync / 2, xminInSync, xmaxInSync,
                                                  xbinsInSpill / 2, xminInSpill, xmaxInSpill);
              hExtMountains.back()[ch]->SetStats(false);
            }
          }

          std::vector<TdcData> tdcData = provider->GetTdcData();

          for (auto&& data : tdcData) {
            const Int_t    globalChannel = data.Channel;
            // const Long64_t tdc           = data.Tdc;
            const Double_t time          = data.Time;
            
            if (ExtinctionDetector::Contains(globalChannel)) {
              fillCoin(data);

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
              std::vector<TdcData> coinExtData;
              {
                const std::size_t hodCh = 0;
                const std::size_t i = hodCh + CoinOffset::Hod;
                if (globalChannel != Hodoscope::GlobalChannelOffset) { // TBD: for 20201217
                  // Nothing to do
                } else {
                  for (auto&& lastData : lastExtData) {
                    const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
                    if (contains[extCh][i]) {
                      const Double_t dt    = time - lastData.Time;
                      const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
                      // const Double_t sigma = 3.0 * coinInfo[extCh][i].FitSigma * timePerTdc;
                      const Double_t sigma = 25.0 * nsec;
                      if (TMath::Abs(dt - mean) < sigma) {
                        coinExtData.push_back(lastData);
                      }
                    }
                  }
                }
              }

              lastHodData.push_back(data);
              for (auto&& extData : coinExtData) {
                fillCoin(extData);
              }

              for (std::size_t i = 0, n = lastHodData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastHodData[0].Time) > lastThreshold) {
                  lastHodData.erase(lastHodData.begin());
                } else {
                  break;
                }
              }
              if (lastHodData.size() > 10000) {
                std::cerr << "[error] size of lastHodData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (TimingCounter::Contains(globalChannel)) {
              const Int_t ch = TimingCounter::GetChannel(globalChannel);

              std::vector<TdcData> coinExtData;
              {
                const std::size_t i = ch + CoinOffset::TC;
                for (auto&& lastData : lastExtData) {
                  const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
                  if (contains[extCh][i]) {
                    const Double_t dt    = time - lastData.Time;
                    const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
                    // const Double_t sigma = 3.0 * coinInfo[extCh][i].FitSigma * timePerTdc;
                    const Double_t sigma = 25.0 * nsec;
                    if (TMath::Abs(dt - mean) < sigma) {
                      coinExtData.push_back(lastData);
                    }
                  }
                }
              }

              lastTcData.push_back(data);
              for (auto&& extData : coinExtData) {
                fillCoin(extData);
              }

              for (std::size_t i = 0, n = lastTcData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastTcData[0].Time) > lastThreshold) {
                  lastTcData.erase(lastTcData.begin());
                } else {
                  break;
                }
              }
              if (lastTcData.size() > 10000) {
                std::cerr << "[error] size of lastTcData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (BeamlineHodoscope::Contains(globalChannel)) {
              const Int_t ch = BeamlineHodoscope::GetChannel(globalChannel);

              std::vector<TdcData> coinExtData;
              {
                const std::size_t i = ch + CoinOffset::BH;
                for (auto&& lastData : lastExtData) {
                  const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
                  if (contains[extCh][i]) {
                    const Double_t dt    = time - lastData.Time;
                    const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
                    // const Double_t sigma = 3.0 * coinInfo[extCh][i].FitSigma * timePerTdc;
                    const Double_t sigma = 25.0 * nsec;
                    if (TMath::Abs(dt - mean) < sigma) {
                      coinExtData.push_back(lastData);
                    }
                  }
                }
              }

              lastBhData.push_back(data);
              for (auto&& extData : coinExtData) {
                fillCoin(extData);
              }

              for (std::size_t i = 0, n = lastBhData.size(); i < n; ++i) {
                if (TMath::Abs(time - lastBhData[0].Time) > lastThreshold) {
                  lastBhData.erase(lastBhData.begin());
                } else {
                  break;
                }
              }
              if (lastBhData.size() > 10000) {
                std::cerr << "[error] size of lastBhData reaches 10000" << std::endl;
                exit(1);
              }

            } else if (MrSync::Contains(globalChannel)) {
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
           hExtMountains,
          };
      }

    }

  }

}

#endif
