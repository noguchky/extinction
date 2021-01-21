#ifndef Extinction_AnaTree_hh
#define Extinction_AnaTree_hh

#include <iostream>
#include <vector>
#include "Units.hh"
#include "Utility.hh"
#include "Detector.hh"
#include "Tdc.hh"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"

namespace Extinction {

  namespace Analyzer {

    namespace AnaTree {

      void Execute(TTree* tree,
                   ITdcDataProvider* provider,
                   Long64_t initEntry = 0) {
        std::cout << "=== Get Tree Information" << std::endl;
        const Long64_t entries = tree->GetEntries();

        std::cout << "=== Get TDC Information" << std::endl;
        const std::string tdcName = provider->GetName();
        const Double_t    timePerTdc = provider->GetTimePerTdc();

        const Double_t  xminInSync  = (Int_t)(-700 * nsec / timePerTdc) - 0.5; // [count]
        const Int_t    xbinsInSync  = (Int_t)(6600 * nsec / timePerTdc) / 10; // [count/10]
        const Double_t  xmaxInSync  = xminInSync + xbinsInSync * 10; // [count]

        const Int_t    xbinsInCh    = GlobalChannel::NofChannels;
        const Double_t  xminInCh    = 0.0 - 0.5;
        const Double_t  xmaxInCh    = xminInCh + xbinsInCh;

        std::cout << "=== Initialzie History" << std::endl;
        std::map<std::size_t, Extinction::TdcData> lastMrSyncData;

        std::cout << "=== Set Branch Address" << std::endl;
        provider->SetBranchAddress(tree);

        // TDC in MR Sync
        TH2* hTdcInSyncByCh[MrSync::NofChannels];
        for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
          hTdcInSyncByCh[ch] = new TH2D(Form("hEntryByChInSync_%03ld", ch),
                                        Form("%s, TDC in MR Sync by Channel;"
                                             "TDC [count];"
                                             "Channel", tdcName.data()),
                                        xbinsInSync, xminInSync, xmaxInSync,
                                        xbinsInCh  , xminInCh  , xmaxInCh  );
        }

        std::cout << "=== Get Entry " << std::endl;
        Int_t lastSpill = -1;
        for (Long64_t entry = initEntry; entry < entries; ++entry) {
          tree->GetEntry(entry);

          if (!provider->IsData()) {
            continue;
          }

          if (provider->GetSpill() != lastSpill) {
            lastMrSyncData.clear();
            for (std::size_t ch = 0; ch < MrSync::NofChannels; ++ch) {
              hTdcInSyncByCh[ch]->Reset();
            }
          }

          std::vector<TdcData> tdcData = provider->GetTdcData();

          for (auto&& data : tdcData) {
            const Int_t    globalChannel = data.Channel;
            const Long64_t tdc           = data.Tdc;

            for (auto&& lastPair : lastMrSyncData) {
              auto lastCh   = lastPair.first;
              auto lastData = lastPair.second;
              if (lastData.Channel == data.MrSyncChannel) {
                const Long64_t lastTdc = lastData.Tdc;
                hTdcInSyncByCh[lastCh]->Fill(tdc - lastTdc, globalChannel);
              }
            }

            if (MrSync::Contains(data.Channel)) {
              const Int_t ch = MrSync::GetChannel(data.Channel);
              if (//hTdcInSyncByCh[ch]->GetEntries() > 10 ||
                  hTdcInSyncByCh[ch]->Integral(1, xbinsInSync, 1, ExtinctionDetector::NofChannels) > 0) {
                hTdcInSyncByCh[ch]->SetTitle(Form("%s, TDC in MR Sync, entry = %lld", tdcName.data(), entry));
                hTdcInSyncByCh[ch]->Draw("colz");
                gPad->SetLogz(true);
                gPad->Update();
                gPad->WaitPrimitive();
              }
              lastMrSyncData[ch] = data;
              hTdcInSyncByCh[ch]->Reset();
            }
          }

          lastSpill = provider->GetSpill();
        }
      }

    }

  }

}

#endif
