#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <sys/time.h>
#include <csignal>
#include <cstring>
#include <ctime>
#include <cerrno>
#include <unistd.h>
#include <mutex>
#include <limits>
#include <sys/inotify.h>
#include <pthread.h>

#include "Units.hh"
#include "Detector.hh"
#include "Fct.hh"

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"

#include "ArgReader.hh"
#include "AnaCoin.hh"

namespace {
  bool processing = false;
  std::mutex filenameMutex;
  std::string path;
  std::string filename;

  bool exiting = false;

  int fd, wd;

  Extinction::Fct::FctData data;
  std::map<std::size_t, std::map<std::size_t, Extinction::Analyzer::AnaCoin::CoinInfo>> coinInfo;
  std::map<std::size_t/*extCh*/, std::map<std::size_t/*index*/, Bool_t>> contains;

  TCanvas* canvas;

  TPad* padBh1Tdc;
  TPad* padBh2Tdc;
  TPad* padExtTdc;
  TPad* padHodTdc;
  TPad* padTc1Tdc;
  TPad* padTc2Tdc;

  TPad* padHodHitMap;
  TPad* padHodHitCnt;
  TPad* padExtHitMap;
  TPad* padExtHitCnt;

  TPad* padMountain;
  TPad* padTdcSync;
  TPad* padHit;
  TPad* padTotalHit;

  Long64_t ispill = 0;
  Long64_t totalCoinCount = 0;

  TH1* hExtTdcInSpill_Any;
  TH1* hHodTdcInSpill_Any;
  TH1** hTcTdcInSpill;
  TH1** hBhTdcInSpill;
  TH2* hExtHitMap;
  TList* lExtBorderLine;
  TH1* hExtEntryByCh;
  TH2* hHodHitMap;
  TList* lHodBorderLine;
  TH1* hHodEntryByCh;
  TH2* hExtMountain_Any;  
  TH1* hExtTdcInSync_Any;
  TGraphErrors* gHitInSpill;
  TGraphErrors* gTotalHitInSpill;
  // TGraphErrors* gExtinction;
}

Int_t updatePlots(const std::string& ifilename) {
  using namespace Extinction;
  using CoinOffset = Analyzer::AnaTimeOffset::CoinOffset;

  if (exiting) {
    return 0;
  }

  std::cout << "=== Set Decoder" << std::endl;
  Extinction::Fct::Decoder decoder;
  Extinction::Fct::Packet_t packet;
  const Double_t timePerTdc = data.GetTimePerTdc();

  std::cout << "=== Initialzie History" << std::endl;
  const Double_t lastThreshold = 100 * Extinction::nsec;
  std::vector<Extinction::TdcData> lastExtData;
  std::vector<Extinction::TdcData> lastHodData;
  std::vector<Extinction::TdcData> lastTcData;
  std::vector<Extinction::TdcData> lastBhData;
  std::vector<Extinction::TdcData> lastMrSyncData;

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

  std::cout << "=== Open File" << std::endl;
  std::ifstream ifile(ifilename, std::ios::binary);
  if (!ifile) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  const Long64_t spilllimit = 200;
  {
    std::cout << "=== Get Entry " << std::endl;
    std::size_t count = 0UL;
    for (; decoder.Read(ifile, &packet); ++count) {
      if (exiting) {
        return 0;
      }
      
      if (count % 1000000UL == 0) {
        std::cout << ">> " << count << std::endl;
      }

      if (decoder.Data.Type == Extinction::Fct::DataType::Header) {
        std::cout << "[info] detect header" << std::endl;

      } else if (decoder.Data.Type == Extinction::Fct::DataType::GateStart) {
        std::cout << "[info] begin of spill " << decoder.Data.Spill << std::endl;

        // Crear hists
        hBhTdcInSpill[0]->Reset();
        hBhTdcInSpill[1]->Reset();
        hHodTdcInSpill_Any->Reset();
        hExtTdcInSpill_Any->Reset();
        hTcTdcInSpill[0]->Reset();
        hTcTdcInSpill[1]->Reset();
        hHodHitMap->Reset();
        hHodEntryByCh->Reset();
        hExtHitMap->Reset();
        hExtEntryByCh->Reset();
        hExtMountain_Any->Reset();
        hExtTdcInSync_Any->Reset();

      } else if (decoder.Data.Type == Extinction::Fct::DataType::GateEnd) {
        std::cout << "[info] end of spill " << decoder.Data.Spill << std::endl;
        lastExtData.clear();
        lastHodData.clear();
        lastTcData.clear();
        lastBhData.clear();
        lastMrSyncData.clear();

        totalCoinCount += coinCount;
        gHitInSpill->SetPoint(ispill % spilllimit, ispill, coinCount);
        gHitInSpill->SetPointError(ispill % spilllimit, 0.0, TMath::Sqrt(coinCount));
        gTotalHitInSpill->SetPoint(ispill % spilllimit, ispill, totalCoinCount);
        gTotalHitInSpill->SetPointError(ispill % spilllimit, 0.0, TMath::Sqrt(totalCoinCount));
        ++ispill;
        coinCount = 0;

        // const Double_t leakCount = ***;
        // gExtinction->SetPoint
        //   (nhit, nhit,
        //    totalCoinCount ?                        leakCount        / totalCoinCount : 0.0);
        // gExtinction->SetPointError
        //   (nhit, nhit,
        //    totalCoinCount ? TMath::Sqrt(TMath::Max(leakCount, 1.0)) / totalCoinCount : 0.0);

        std::cout << "=== Draw hists" << std::endl;
        Int_t padnumber = 0;
        // padBh1Tdc->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hBhTdcInSpill[0]->Draw();

        // padBh2Tdc->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hBhTdcInSpill[1]->Draw();

        // padHodTdc->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hHodTdcInSpill_Any->Draw();

        // padExtTdc->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hExtTdcInSpill_Any->Draw();

        // padTc1Tdc->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hTcTdcInSpill[0]->Draw();

        // padTc2Tdc->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hTcTdcInSpill[1]->Draw();

        
        // padHodHitMap->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hHodHitMap->Draw("col");
        lHodBorderLine->Draw();

        // padHodHitCnt->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hHodEntryByCh->Draw();

        // padExtHitMap->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hExtHitMap->Draw("col");
        lExtBorderLine->Draw();

        // padExtHitCnt->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hExtEntryByCh->Draw();


        // padMountain->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hExtMountain_Any->Draw("colz");

        // padTdcSync->cd();
        if (!canvas->cd(++padnumber)) return 0;
        hExtTdcInSync_Any->Draw();

        // padHit->cd();
        if (!canvas->cd(++padnumber)) return 0;
        gHitInSpill->Draw("AP");

        // padTotalHit->cd();
        if (!canvas->cd(++padnumber)) return 0;
        gTotalHitInSpill->Draw("AP");

        // padExtinction->cd();
        // gExtinction->Draw("AP");

        canvas->Modified();
        canvas->Update();
        // canvas->WaitPrimitive();

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
              // exit(1);
              return 1;
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
              // exit(1);
              return 1;
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
              // exit(1);
              return 1;
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
              // exit(1);
              return 1;
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
              // exit(1);
              return 1;
            }

          }

        }
      }
    }
  }

  std::cout << "=== Close Files" << std::endl;
  ifile.close();

  canvas->Modified();
  canvas->Update();
  // canvas->WaitPrimitive();

  return 0;
}

void SignalHandler(int) {
  static unsigned long msec_cnt = 0;
  if (exiting) {
    return;
  }

  msec_cnt++;
  if (gPad) {
    // if (!(msec_cnt % 50)) {
    //   printf("SignalHandler:%lu sec\n", (msec_cnt / 50));
    {

      // Check processing
      if (processing) {
        return;
      }
      processing = true;

      // Get Filename
      std::string ifilename;
      {
        std::lock_guard<std::mutex> lock(filenameMutex);
        ifilename = filename;
        filename = "";
      }
      if (ifilename.empty()) {
        processing = false;
        return;
      }

      try {
        // Update hist
        std::cout << "[info] start modify plots, " << ifilename << std::endl;
        updatePlots(path + "/" + ifilename);
        processing = false;
      } catch (...) {
        std::cerr << "[info] something error was detected" << std::endl;
        gSystem->ExitLoop();
      }

    }
  } else {
    std::cerr << "[info] canvas was closed" << std::endl;
    exiting = true;
    gSystem->ExitLoop();
  }
}

void* monitorWritten(void*) {
  int i, aux = 0, ret;
  constexpr int BufferSize = 65536;
  char buffer[BufferSize];

  while (!exiting) {
    // Read events
    ret = read(fd, buffer + aux, BufferSize - aux);
    if (ret == -1) {
      if (ret == -EINTR) {
        continue;
      }
      perror("read");
      exit(EXIT_FAILURE);
    }
    ret += aux;

    // Read size check
    if (ret < (int)sizeof(inotify_event)) {
      fprintf(stderr, "short of red bytes\n");
      exit(EXIT_FAILURE);
    }

    for (i = 0; i < ret;) {
      inotify_event* inotify_p;
      inotify_p = (inotify_event*)(buffer + i);

      if (ret < i + (int)offsetof(inotify_event, name)) {
        // Get tail at next read
        aux = ret - i;
        memmove(buffer, buffer + i, aux);
        break;
      }

      int size = sizeof(inotify_event) + inotify_p->len;
      if (ret < i + size) {
        // Get tail at next read
        aux = ret - i;
        memmove(buffer, buffer + i, aux);
        break;
      }

      // File opened for writing was closed (*).
      if (inotify_p->mask & IN_CLOSE_WRITE ||
          inotify_p->mask & IN_MOVED_TO) {
        std::cerr << "[info] file was closed \"" << inotify_p->name << "\"" << std::endl;
        std::lock_guard<std::mutex> lock(filenameMutex);
        filename = inotify_p->name;
      }
      i += size;
    }
  }

  // Remove monitor
  {
    std::cerr << "[info] remove file monitor" << std::endl;
    int ret = inotify_rm_watch(fd, wd);
    if (ret == -1) {
      perror("inotify_rm_watch");
      // exit(EXIT_FAILURE);
      return 0;
    }
  }

  close(fd);

  return 0;
}

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Directory" ,                    "A rawdata directory");
  args->AddArg<std::string>("Offset"    ,                    "A offset filename format");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string directory = args->GetValue("Directory");
  const std::string ffilename = args->GetValue("Offset");
  path = directory;

  std::cout << "=== Initialize Application" << std::endl;
  TApplication* app = new TApplication("monitor", nullptr, nullptr);

  std::cout << "=== Set Style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);

  std::cout << "=== Get Offset" << std::endl;
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

  std::cout << "=== Initialize Canvas" << std::endl;
  canvas = new TCanvas("c1", "FCT | Semi-online Monitor", 1600, 1200);
  std::vector<TPad*> pads;

  const Double_t seg = 1.0 / 6.0;

  padBh1Tdc = new TPad("PadBh1Tdc", "", 0 * seg, 5 * seg, 1 * seg, 6 * seg);
  pads.push_back(padBh1Tdc);
  padBh2Tdc = new TPad("PadBh2Tdc", "", 0 * seg, 4 * seg, 1 * seg, 5 * seg);
  pads.push_back(padBh2Tdc);
  padExtTdc = new TPad("PadExtTdc", "", 0 * seg, 3 * seg, 1 * seg, 4 * seg);
  pads.push_back(padExtTdc);
  padHodTdc = new TPad("PadHodTdc", "", 0 * seg, 2 * seg, 1 * seg, 3 * seg);
  pads.push_back(padHodTdc);
  padTc1Tdc = new TPad("PadTc1Tdc", "", 0 * seg, 1 * seg, 1 * seg, 2 * seg);
  pads.push_back(padTc1Tdc);
  padTc2Tdc = new TPad("PadTc2Tdc", "", 0 * seg, 0 * seg, 1 * seg, 1 * seg);
  pads.push_back(padTc2Tdc);

  padHodHitMap = new TPad("PadHodHitMap", "", 1 * seg, 5 * seg, 3 * seg, 6 * seg);
  pads.push_back(padHodHitMap);
  padHodHitCnt = new TPad("PadHodHitCnt", "", 1 * seg, 4 * seg, 3 * seg, 5 * seg);
  pads.push_back(padHodHitCnt);
  padExtHitMap = new TPad("PadExtHitMap", "", 1 * seg, 2 * seg, 3 * seg, 4 * seg);
  pads.push_back(padExtHitMap);
  padExtHitCnt = new TPad("PadExtHitCnt", "", 1 * seg, 0 * seg, 3 * seg, 2 * seg);
  pads.push_back(padExtHitCnt);

  padMountain = new TPad("PadMountain", "", 3 * seg, 4 * seg, 6 * seg, 6 * seg);
  pads.push_back(padMountain);
  padTdcSync  = new TPad("PadTdcSync" , "", 3 * seg, 2 * seg, 6 * seg, 4 * seg);
  pads.push_back(padTdcSync);
  padHit      = new TPad("padHit"     , "", 3 * seg, 1 * seg, 6 * seg, 2 * seg);
  pads.push_back(padHit);
  padTotalHit = new TPad("padTotalHit", "", 3 * seg, 0 * seg, 6 * seg, 1 * seg);
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
  padHit     ->SetGrid(); // padHit     ->SetLogy();
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
    const std::string tdcName    = data.GetName();
    const Double_t    timePerTdc = data.GetTimePerTdc();

    const Double_t  xminInSpill =    0; // [msec]
    const Double_t  xmaxInSpill = 3000; // [msec]
    const Int_t    xbinsInSpill =   50;

    // const Double_t  xminInSync  = (Int_t)(-700 * nsec / timePerTdc) - 0.5; // [count]
    // const Int_t    xbinsInSync0 = (Int_t)(6600 * nsec / timePerTdc) / 20; // [count/10]
    const Double_t  xminInSync  = (Int_t)(- 5400 * nsec / timePerTdc) - 0.5; // [count]
    const Int_t    xbinsInSync0 = (Int_t)(+16200 * nsec / timePerTdc) / 20; // [count/10]
    const Int_t    xbinsInSync  = xbinsInSync0 < 200 ? xbinsInSync0 : xbinsInSync0 / 4;
    const Double_t  xmaxInSync  = xbinsInSync0 < 200 ? xminInSync + xbinsInSync * 20 : xminInSync + xbinsInSync * 80; // [count]

    std::cout << "=== Initialzie Coincidence Information" << std::endl;
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

    std::cout << "=== Create hists" << std::endl;
    // Extinction Detector TDC in spill
    hExtTdcInSpill_Any = new TH1D("hExtTdcInSpill_Any",
                                  Form("%s, Extinction Detector TDC in Spill;"
                                       "Time [ms]", tdcName.data()),
                                  xbinsInSpill, xminInSpill, xmaxInSpill);

    // Hodoscope TDC in spill
    hHodTdcInSpill_Any = new TH1D("hHodTdcInSpill_Any",
                                  Form("%s, Hodoscope TDC in Spill;"
                                       "Time [ms]", tdcName.data()),
                                  xbinsInSpill, xminInSpill, xmaxInSpill);

    // Timing Counter TDC in spill
    hTcTdcInSpill = new TH1*[TimingCounter::NofChannels];
    for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
      hTcTdcInSpill[ch] = new TH1D(Form("hTcTdcInSpill_%03ld", ch),
                                   Form("%s, Timing Counter %ld TDC in Spill;"
                                        "Time [ms]", tdcName.data(), ch + 1),
                                   xbinsInSpill, xminInSpill, xmaxInSpill);
    }

    // Beamline Hodoscope TDC in spill
    hBhTdcInSpill = new TH1*[BeamlineHodoscope::NofChannels];
    for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
      hBhTdcInSpill[ch] = new TH1D(Form("hBhTdcInSpill_%03ld", ch),
                                   Form("%s, BH%ld TDC in Spill;"
                                        "Time [ms]", tdcName.data(), ch + 1),
                                   xbinsInSpill, xminInSpill, xmaxInSpill);
    }

    // Extinction detector hit map
    hExtHitMap     = Extinction::ExtinctionDetector::CreateHitMap("hExtHitMap");
    lExtBorderLine = Extinction::ExtinctionDetector::CreateBorderLine(kBlack, kSolid, 1);

    // Extinction detector hit count
    hExtEntryByCh = new TH1D("hExtEntryByCh",
                             Form("%s, Extinction Detector Entries by Channel;"
                                  "Channel;"
                                  "Entries", tdcName.data()),
                             ExtinctionDetector::NofChannels, 0.0 - 0.5, ExtinctionDetector::NofChannels - 0.5);

    // Hodoscope hit map
    hHodHitMap     = Extinction::Hodoscope::CreateHitMap("hHodHitMap");
    lHodBorderLine = Extinction::Hodoscope::CreateBorderLine(kBlack, kSolid, 1);

    // Hodoscope hit count
    hHodEntryByCh = new TH1D("hHodEntryByCh",
                             Form("%s, Hodoscope Entries by Channel;"
                                  "Channel;"
                                  "Entries", tdcName.data()),
                             Hodoscope::NofChannels, 0.0 - 0.5, Hodoscope::NofChannels - 0.5);

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
    gHitInSpill->SetMarkerStyle(kOpenSquare);

    // Total hit in spill (coincidence)
    gTotalHitInSpill = new TGraphErrors();
    gTotalHitInSpill->SetNameTitle("gTotalHitInSpill",
                                   Form("%s, # of Total Hits;"
                                        "Spill", tdcName.data()));
    gTotalHitInSpill->SetMarkerStyle(kOpenSquare);

    // // Extinction
    // gExtinction = new TGraphErrors();
    // gExtinction->SetNameTitle("gExtinction",
    //                           Form("%s, Extinction;"
    //                                "Spill", tdcName.data()));
  }

  {
    std::cout << "=== Draw hists" << std::endl;
    Int_t padnumber = 0;
    // padBh1Tdc->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hBhTdcInSpill[0]->Draw();

    // padBh2Tdc->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hBhTdcInSpill[1]->Draw();

    // padHodTdc->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hHodTdcInSpill_Any->Draw();

    // padExtTdc->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hExtTdcInSpill_Any->Draw();

    // padTc1Tdc->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hTcTdcInSpill[0]->Draw();

    // padTc2Tdc->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hTcTdcInSpill[1]->Draw();

        
    // padHodHitMap->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hHodHitMap->Draw("col");
    lHodBorderLine->Draw();

    // padHodHitCnt->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hHodEntryByCh->Draw();

    // padExtHitMap->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hExtHitMap->Draw("col");
    lExtBorderLine->Draw();

    // padExtHitCnt->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hExtEntryByCh->Draw();


    // padMountain->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hExtMountain_Any->Draw("colz");

    // padTdcSync->cd();
    if (!canvas->cd(++padnumber)) return 0;
    hExtTdcInSync_Any->Draw();

    // padHit->cd();
    if (!canvas->cd(++padnumber)) return 0;
    gHitInSpill->Draw("AP");

    // padTotalHit->cd();
    if (!canvas->cd(++padnumber)) return 0;
    gTotalHitInSpill->Draw("AP");

    // padExtinction->cd();
    // gExtinction->Draw("AP");

    canvas->Modified();
    canvas->Update();
    // canvas->WaitPrimitive();
  }

  std::cout << "Set signal handler" << std::endl;
  struct sigaction action;
  memset(&action, 0, sizeof(action));

  action.sa_handler = SignalHandler;
  action.sa_flags = SA_RESTART;
  sigemptyset(&action.sa_mask);
  if (sigaction(SIGALRM, &action, nullptr) < 0) {
    perror("sigaction error");
    exit(1);
  }

  std::cout << "Set intarval timer" << std::endl;
  itimerval timer;
  timer.it_value.tv_sec = 0;
  timer.it_value.tv_usec = 20000; // = 20 ms
  timer.it_interval.tv_sec = 0;
  timer.it_interval.tv_usec = 20000;
  if (setitimer(ITIMER_REAL, &timer, nullptr) < 0) {
    // perror("setitimer error");
    std::cerr << "[error] setitimer error" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Initialize inotify
  fd = inotify_init();
  if (fd == -1) {
    std::cerr << "[error] inotify initialize error" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Add monitor directory
  std::cerr << "[info] add file monitor" << std::endl;
  wd = inotify_add_watch(fd, directory.data(),  IN_ALL_EVENTS);
  if (wd < 0) {
    // perror("inotify_add_watch");
    std::cerr << "[error] directory was not opened, " << directory << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "Set write monitor" << std::endl;
  pthread_t pthread;
  if (pthread_create(&pthread, nullptr, &monitorWritten, nullptr)) {
    std::cerr << "[error] pthread_create error" << std::endl;
    exit(EXIT_FAILURE);
  }

  app->Run(true);
  delete app; app = nullptr;

  exit(0);
}
