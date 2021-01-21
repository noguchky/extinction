#include <iostream>
#include <fstream>

#include "Units.hh"
#include "Detector.hh"
#include "Kc705.hh"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

#include "ArgReader.hh"
#include "AnaCoin.hh"

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"     ,                    "A rawdata filename");
  args->AddArg<std::string>("Offset"    ,                    "A offset filename format");
  args->AddOpt             ("Verbose"   , 'v', "verbose"   , "Output filled values");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename       = args->GetValue("Input");
  const std::string ffilename       = args->GetValue("Offset");

  TApplication* app = new TApplication("app", nullptr, nullptr);

  std::cout << "=== Open Input File" << std::endl;
  std::ifstream ifile(ifilename, std::ios::binary);
  if (!ifile) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  std::cout << "=== Set Style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1);

  std::cout << "=== Set Decoder" << std::endl;
  Extinction::Kc705::Decoder decoder;
  Extinction::Kc705::Kc705Data data;
  Extinction::Kc705::Packet_t packet;

  std::cout << "=== Get Offset" << std::endl;
  std::map<std::size_t, std::map<std::size_t, Extinction::Analyzer::AnaCoin::CoinInfo>> coinInfo;
  for (std::size_t extCh = 0; extCh < Extinction::ExtinctionDetector::NofChannels; ++extCh) {
    Extinction::Analyzer::AnaCoin::CoinInfo info;
    std::ifstream ffile(Form(ffilename.data(), extCh));
    if (ffile) {
      while (info.Read(ffile)) {
        coinInfo[extCh][info.Index] = info;
      }
    }
  }
  if (coinInfo.empty()) {
    std::cout << "[error] offset file was not found" << std::endl;
    exit(1);
  }

  TCanvas* canvas = new TCanvas();
  std::vector<TPad*> pads;

  const Double_t seg = 1.0 / 6.0;

  TPad* padBh1Tdc = new TPad("PadBh1Tdc", "", 0 * seg, 5 * seg, 1 * seg, 6 * seg);
  pads.push_back(padBh1Tdc);
  TPad* padBh2Tdc = new TPad("PadBh2Tdc", "", 0 * seg, 4 * seg, 1 * seg, 5 * seg);
  pads.push_back(padBh2Tdc);
  TPad* padExtTdc = new TPad("PadExtTdc", "", 0 * seg, 3 * seg, 1 * seg, 4 * seg);
  pads.push_back(padExtTdc);
  TPad* padHodTdc = new TPad("PadHodTdc", "", 0 * seg, 2 * seg, 1 * seg, 3 * seg);
  pads.push_back(padHodTdc);
  TPad* padTc1Tdc = new TPad("PadTc1Tdc", "", 0 * seg, 1 * seg, 1 * seg, 2 * seg);
  pads.push_back(padTc1Tdc);
  TPad* padTc2Tdc = new TPad("PadTc2Tdc", "", 0 * seg, 0 * seg, 1 * seg, 1 * seg);
  pads.push_back(padTc2Tdc);

  TPad* padHodHitMap = new TPad("PadHodHitMap", "", 1 * seg, 5 * seg, 3 * seg, 6 * seg);
  pads.push_back(padHodHitMap);
  TPad* padHodHitCnt = new TPad("PadHodHitCnt", "", 1 * seg, 4 * seg, 3 * seg, 5 * seg);
  pads.push_back(padHodHitCnt);
  TPad* padExtHitMap = new TPad("PadExtHitMap", "", 1 * seg, 2 * seg, 3 * seg, 4 * seg);
  pads.push_back(padExtHitMap);
  TPad* padExtHitCnt = new TPad("PadExtHitCnt", "", 1 * seg, 0 * seg, 3 * seg, 2 * seg);
  pads.push_back(padExtHitCnt);

  TPad* padMountain = new TPad("PadMountain", "", 3 * seg, 4 * seg, 6 * seg, 6 * seg);
  pads.push_back(padMountain);
  TPad* padTdcSync  = new TPad("PadTdcSync" , "", 3 * seg, 2 * seg, 6 * seg, 4 * seg);
  pads.push_back(padTdcSync);
  TPad* padHit      = new TPad("padHit"     , "", 3 * seg, 1 * seg, 6 * seg, 2 * seg);
  pads.push_back(padHit);
  TPad* padTotalHit = new TPad("padTotalHit", "", 3 * seg, 0 * seg, 6 * seg, 1 * seg);
  pads.push_back(padTotalHit);

  padBh1Tdc->SetGrid(); padBh1Tdc->SetLogy();  
  padBh2Tdc->SetGrid(); padBh2Tdc->SetLogy(); 
  padExtTdc->SetGrid(); padExtTdc->SetLogy(); 
  padHodTdc->SetGrid(); padHodTdc->SetLogy(); 
  padTc1Tdc->SetGrid(); padTc1Tdc->SetLogy(); 
  padTc2Tdc->SetGrid(); padTc2Tdc->SetLogy(); 

  padHodHitCnt->SetGrid(); padHodHitCnt->SetLogy();
  padExtHitCnt->SetGrid(); padExtHitCnt->SetLogy();

  padMountain->SetGrid();
  padTdcSync ->SetGrid(); padTdcSync ->SetLogy();
  padHit     ->SetGrid(); padHit     ->SetLogy();
  padTotalHit->SetGrid(); padTotalHit->SetLogy();

  for (Int_t ipad = 0, npad = pads.size(); ipad < npad; ++ipad) {
    canvas->cd();
    pads[ipad]->Draw();
    pads[ipad]->SetNumber(ipad + 1);
  }

  {
    using namespace Extinction;
    using CoinOffset = Analyzer::AnaTimeOffset::CoinOffset;

    std::cout << "=== Get TDC Information" << std::endl;
    const std::string tdcName = data.GetName();
    const Double_t    timePerTdc = data.GetTimePerTdc();

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

    std::cout << "=== Create hists" << std::endl;
    // Extinction Detector TDC in spill
    TH1* hExtTdcInSpill_Any = new TH1D("hExtTdcInSpill_Any",
                                       Form("%s, Extinction Detector TDC in Spill;"
                                            "Time [ms]", tdcName.data()),
                                       xbinsInSpill, xminInSpill, xmaxInSpill);

    // Hodoscope TDC in spill
    TH1* hHodTdcInSpill_Any = new TH1D("hHodTdcInSpill_Any",
                                       Form("%s, Hodoscope TDC in Spill;"
                                            "Time [ms]", tdcName.data()),
                                       xbinsInSpill, xminInSpill, xmaxInSpill);

    // Timing Counter TDC in spill
    TH1** hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
    for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
      hTcTdcInSpill[ch] = new TH1D(Form("hTcTdcInSpill_%03ld", ch),
                                   Form("%s, Timing Counter %ld TDC in Spill;"
                                        "Time [ms]", tdcName.data(), ch + 1),
                                   xbinsInSpill, xminInSpill, xmaxInSpill);
    }

    // Beamline Hodoscope TDC in spill
    TH1** hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
    for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
      hBhTdcInSpill[ch] = new TH1D(Form("hBhTdcInSpill_%03ld", ch),
                                   Form("%s, BH%ld TDC in Spill;"
                                        "Time [ms]", tdcName.data(), ch + 1),
                                   xbinsInSpill, xminInSpill, xmaxInSpill);
    }

    // Extinction detector hit map
    TH2* hExtHitMap = Extinction::ExtinctionDetector::CreateHitMap("hExtHitMap");
    TList* lExtBorderLine = Extinction::ExtinctionDetector::CreateBorderLine();

    // Extinction detector hit count
    TH1* hExtEntryByCh = new TH1D("hExtEntryByCh",
                                  Form("%s, Extinction Detector Entries by Channel;"
                                       "Channel;"
                                       "Entries", tdcName.data()),
                                  ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5);

    // Hodoscope hit map
    TH2* hHodHitMap = Extinction::Hodoscope::CreateHitMap("hHodHitMap");
    TList* lHodBorderLine = Extinction::Hodoscope::CreateBorderLine();

    // Hodoscope hit count
    TH1* hHodEntryByCh = new TH1D("hHodEntryByCh",
                                  Form("%s, Hodoscope Entries by Channel;"
                                       "Channel;"
                                       "Entries", tdcName.data()),
                                  Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5);

    // Extinction Detector Mountain Plot (coincidence)
    TH2* hExtMountain_Any = new TH2D("hExtMountain_Any",
                                     Form("%s, Extinction Detector Mountain Plot;"
                                          "TDC [count];"
                                          "Time [ms]", tdcName.data()),
                                     xbinsInSync / 2, xminInSync, xmaxInSync,
                                     xbinsInSpill / 2, xminInSpill, xmaxInSpill);
    hExtMountain_Any->SetStats(false);

    // Extinction Detector TDC in sync (coincidence)
    TH1* hExtTdcInSync_Any = new TH1D("hExtTdcInSync_Any",
                                      Form("%s, Extinction Detector TDC in MR Sync;"
                                           "TDC [count]", tdcName.data()),
                                      xbinsInSync, xminInSync, xmaxInSync);

    // Hit in spill (coincidence)
    TGraphErrors* gHitInSpill = new TGraphErrors();
    gHitInSpill->SetNameTitle("gHitInSpill",
                              Form("%s, # of Hits in Spill;"
                                   "Spill", tdcName.data()));

    // Total hit in spill (coincidence)
    TGraphErrors* gTotalHitInSpill = new TGraphErrors();
    gTotalHitInSpill->SetNameTitle("gTotalHitInSpill",
                                   Form("%s, # of Total Hits;"
                                        "Spill", tdcName.data()));

    // // Extinction
    // TGraphErrors* gExtinction = new TGraphErrors();
    // gExtinction->SetNameTitle("gExtinction",
    //                           Form("%s, Extinction;"
    //                                "Spill", tdcName.data()));

    std::cout << "=== Create Filler" << std::endl;
    Long64_t coinCount = 0;
    auto fillCoin =
      [&](TdcData extData) {
        const std::size_t extCh = ExtinctionDetector::GetChannel(extData.Channel);
        const Long64_t    tdc   = extData.Tdc;
        const Double_t    time  = extData.Time;
        Bool_t coincidence[CoinOffset::N];
        std::vector<std::size_t> coinHodChs;
        for (std::size_t bhCh = 0; bhCh < BeamlineHodoscope::NofChannels; ++bhCh) {
          const std::size_t i = bhCh + CoinOffset::BH;
          if (contains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : lastBhData) {
              const Double_t dt    = lastData.Time - time;
              const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
              const Double_t sigma = 25.0 * nsec;
              if (TMath::Abs(dt - mean) < sigma) {
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
          if (contains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : lastHodData) {
              if (lastData.Channel != Hodoscope::GlobalChannelOffset) {
                continue;
              } // TBD: for 20201217
              const Double_t dt    = lastData.Time - time;
              const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
              const Double_t sigma = 25.0 * nsec;
              if (TMath::Abs(dt - mean) < sigma) {
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
          if (contains[extCh][i]) {
            coincidence[i] = false;
            for (auto&& lastData : lastTcData) {
              const Double_t dt    = lastData.Time - time;
              const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
              const Double_t sigma = 25.0 * nsec;
              if (TMath::Abs(dt - mean) < sigma) {
                coincidence[i] = true;
                break;
              }
            }
          } else {
            coincidence[i] = true;
          }
        }

        if (std::all_of(coincidence, coincidence + CoinOffset::N, [](Bool_t b) { return b; })) {
          for (std::size_t n = lastMrSyncData.size(), i = n; i != 0; --i) {
            auto lastData = lastMrSyncData[i];
            if (lastData.Channel == extData.MrSyncChannel) {
              auto syncTdc = lastData.Tdc;
              ++coinCount;
              hExtTdcInSync_Any->Fill(tdc - syncTdc);
              hExtMountain_Any ->Fill(tdc - syncTdc, time / msec);
              break;
            }
          }
        }
      };

    std::cout << "=== Get Entry " << std::endl;
    std::size_t count = 0UL;
    Long64_t totalCoinCount = 0;
    for (; decoder.Read(ifile, &packet); ++count) {
      if (count % 1000000UL == 0) {
        std::cout << ">> " << count << std::endl;
      }

      if (decoder.Data.Type == Extinction::Kc705::DataType::Header) {
        std::cout << "[info] begin of spill " << decoder.Data.Spill << std::endl;

      } else if (decoder.Data.Type == Extinction::Kc705::DataType::HeaderError) {
        // Nothing to do

      } else if (decoder.Data.Type == Extinction::Kc705::DataType::Footer) {
        std::cout << "[info] end of spill " << decoder.Data.Spill << std::endl;
        lastExtData.clear();
        lastHodData.clear();
        lastTcData.clear();
        lastBhData.clear();
        lastMrSyncData.clear();

        const Int_t nhit = gHitInSpill->GetN();
        gHitInSpill->SetPoint(nhit, nhit, coinCount);
        gHitInSpill->SetPointError(nhit, 0.0, TMath::Sqrt(coinCount));

        totalCoinCount += coinCount;
        gTotalHitInSpill->SetPoint(nhit, nhit, totalCoinCount);
        gTotalHitInSpill->SetPointError(nhit, 0.0, TMath::Sqrt(totalCoinCount));

        // const Double_t leakCount = ***;
        // gExtinction->SetPoint
        //   (nhit, nhit,
        //    totalCoinCount ?                        leakCount        / totalCoinCount : 0.0);
        // gExtinction->SetPointError
        //   (nhit, nhit,
        //    totalCoinCount ? TMath::Sqrt(TMath::Max(leakCount, 1.0)) / totalCoinCount : 0.0);


        // Draw Hists
        std::cout << "=== Create hists" << std::endl;
        padBh1Tdc->cd();
        hBhTdcInSpill[0]->Draw();

        padBh2Tdc->cd();
        hBhTdcInSpill[1]->Draw();

        padHodTdc->cd();
        hHodTdcInSpill_Any->Draw();

        padExtTdc->cd();
        hExtTdcInSpill_Any->Draw();

        padTc1Tdc->cd();
        hTcTdcInSpill[0]->Draw();

        padTc2Tdc->cd();
        hTcTdcInSpill[1]->Draw();

        
        padHodHitMap->cd();
        hHodHitMap->Draw("col");
        lHodBorderLine->Draw();

        padHodHitCnt->cd();
        hHodEntryByCh->Draw();

        padExtHitMap->cd();
        hExtHitMap->Draw("col");
        lExtBorderLine->Draw();

        padExtHitCnt->cd();
        hExtEntryByCh->Draw();

        
        padMountain->cd();
        hExtMountain_Any->Draw("colz");

        padTdcSync->cd();
        hExtTdcInSync_Any->Draw();

        padHit->cd();
        gHitInSpill->Draw("AP");

        padTotalHit->cd();
        gTotalHitInSpill->Draw("AP");

        // padExtinction->cd();
        // gExtinction->Draw("AP");

        canvas->Modified();
        canvas->Update();
        canvas->WaitPrimitive();

        // Crear hists
        hBhTdcInSpill[0]->Clear();
        hBhTdcInSpill[1]->Clear();
        hHodTdcInSpill_Any->Clear();
        hExtTdcInSpill_Any->Clear();
        hTcTdcInSpill[0]->Clear();
        hTcTdcInSpill[1]->Clear();
        hHodHitMap->Clear();
        hHodEntryByCh->Clear();
        hExtHitMap->Clear();
        hExtEntryByCh->Clear();
        hExtMountain_Any->Clear();
        hExtTdcInSync_Any->Clear();
        
      } else {
        std::vector<TdcData> tdcData = decoder.Data.GetTdcData();

        for (auto&& data : tdcData) {
          const Int_t    globalChannel = data.Channel;
          const Double_t time          = data.Time;

          if (ExtinctionDetector::Contains(globalChannel)) {
            const Int_t ch = ExtinctionDetector::GetChannel(globalChannel);

            hExtEntryByCh     ->Fill(ch);
            ExtinctionDetector::Fill(hExtHitMap, ch);
            hExtTdcInSpill_Any->Fill(time / Extinction::msec);

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
            const Int_t ch = Hodoscope::GetChannel(globalChannel);

            hHodEntryByCh     ->Fill(ch);
            Hodoscope         ::Fill(hHodHitMap, ch);
            hHodTdcInSpill_Any->Fill(time / Extinction::msec);

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

            hTcTdcInSpill[ch]->Fill(time / Extinction::msec);

            std::vector<TdcData> coinExtData;
            {
              const std::size_t i = ch + CoinOffset::TC;
              for (auto&& lastData : lastExtData) {
                const std::size_t extCh = ExtinctionDetector::GetChannel(lastData.Channel);
                if (contains[extCh][i]) {
                  const Double_t dt    = time - lastData.Time;
                  const Double_t mean  = coinInfo[extCh][i].FitMean * timePerTdc;
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

            hBhTdcInSpill[ch]->Fill(time / Extinction::msec);

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
      }
    }
  }

  std::cout << "=== Close Files" << std::endl;
  ifile.close();

  app->Run();

  return 0;
}
