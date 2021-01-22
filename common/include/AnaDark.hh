#ifndef Extinction_AnaDark_hh
#define Extinction_AnaDark_hh

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
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
#include "TGraph.h"

namespace Extinction {

  namespace Analyzer {

    namespace AnaDark {

      struct DarkInfo {
        Int_t       SpillCount;
        std::size_t GlobalChannel;
        Long64_t    DarkCount;
        Long64_t    TotalCount;

        void Write(std::ofstream& file) {
          file << std::setw( 2) << SpillCount    << "    "
               << std::setw( 3) << GlobalChannel << "    "
               << std::setw(12) << DarkCount     << "    "
               << std::setw(12) << TotalCount    << std::endl;
        }
        std::basic_istream<char>& Read(std::ifstream& file) {
          return file >> SpillCount
                      >> GlobalChannel
                      >> DarkCount
                      >> TotalCount;
        }
      };
      
      struct Results_t {
        // Extinction detector dark count
        TGraph** gExtDarkBySpill;
        TGraph* gExtDarkBySpill_Any;

        // Hodoscope dark count
        TGraph** gHodDarkBySpill;
        TGraph* gHodDarkBySpill_Any;

        // Beamline hodoscope dark count
        TGraph** gBhDarkBySpill;

        // Timing counter dark count
        TGraph** gTcDarkBySpill;

        void Print(const std::string& ofilename) {
          TCanvas::MakeDefCanvas();
          gPad->SetGrid();

          gPad->Print((ofilename + "[").data());

          std::cout << "// Extinction Detector Dark Count" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (gExtDarkBySpill[ch]->GetN()) {
              gExtDarkBySpill[ch]->Draw("AP");
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          gExtDarkBySpill_Any->Draw("AP");
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          std::cout << "// Hodoscope Dark Count" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (gHodDarkBySpill[ch]->GetN()) {
              gHodDarkBySpill[ch]->Draw("AP");
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->SetLogy(true);
          gHodDarkBySpill_Any->Draw("AP");
          gPad->Print(ofilename.data());
          gPad->SetLogy(false);

          std::cout << "// Beamline Hodoscope Dark Count" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (gBhDarkBySpill[ch]->GetN()) {
              gBhDarkBySpill[ch]->Draw("AP");
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          std::cout << "// Timing Counter Dark Count" << std::endl;
          gPad->SetLogy(true);
          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (gTcDarkBySpill[ch]->GetN()) {
              gTcDarkBySpill[ch]->Draw("AP");
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->Print((ofilename + "]").data());
        }

        void Write(TFile* file) {
          file->cd();

          // Extinction Detector
          for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
            if (gExtDarkBySpill[ch]->GetN()) {
              gExtDarkBySpill[ch]->Write();
            }
          }
          gExtDarkBySpill_Any->Write();

          // Hodoscope
          for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
            if (gHodDarkBySpill[ch]->GetN()) {
              gHodDarkBySpill[ch]->Write();
            }
          }
          gHodDarkBySpill_Any->Write();

          // Beamline Hodoscope
          for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
            if (gBhDarkBySpill[ch]->GetN()) {
              gBhDarkBySpill[ch]->Write();
            }
          }

          // Timing Counter
          for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
            if (gTcDarkBySpill[ch]->GetN()) {
              gTcDarkBySpill[ch]->Write();
            }
          }
        }
      };

      Results_t Execute(TTree* tree,
                        ITdcDataProvider* provider,
                        const std::string& ofilename) {
        std::ofstream ofile(ofilename);
        if (!ofile.is_open()) {
          std::cout << "[error] output file is not opened, " << ofilename  << std::endl;
          return
            {
             nullptr,
             nullptr,
             nullptr,
             nullptr,
             nullptr,
             nullptr,
            };
        }

        std::cout << "=== Get Tree Information" << std::endl;
        const Long64_t entries = tree->GetEntries();

        std::cout << "=== Get TDC Information" << std::endl;
        const std::string tdcName = provider->GetName();

        const Double_t tmin =   0 * msec;
        const Double_t tmax = 400 * msec;

        std::cout << "=== Set Branch Address" << std::endl;
        provider->SetBranchAddress(tree);

        // Extinction detector dark count
        TGraph** gExtDarkBySpill = new TGraph*[ExtinctionDetector::NofChannels];
        for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
          gExtDarkBySpill[ch] = new TGraph();
          gExtDarkBySpill[ch]->SetNameTitle(Form("gExtDarkBySpill_%03ld", ch),
                                            Form("%s, Extinction Detector Dark Count by Spill @ ch%ld;"
                                                 "Spill;"
                                                 "Entries", tdcName.data(), ch));
          gExtDarkBySpill[ch]->SetMarkerStyle(kOpenCircle);
          gExtDarkBySpill[ch]->SetMarkerColor(kBlack);
          gExtDarkBySpill[ch]->SetMarkerSize(1);
        }

        TGraph* gExtDarkBySpill_Any = new TGraph();
        gExtDarkBySpill_Any->SetNameTitle("gExtDarkBySpill_Any",
                                          Form("%s, Extinction Detector Dark Count by Spill;"
                                               "Spill;"
                                               "Entries", tdcName.data()));
        gExtDarkBySpill_Any->SetMarkerStyle(kOpenCircle);
        gExtDarkBySpill_Any->SetMarkerColor(kBlack);
        gExtDarkBySpill_Any->SetMarkerSize(1);

        // Hodoscope dark count
        TGraph** gHodDarkBySpill = new TGraph*[Hodoscope::NofChannels];
        for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
          gHodDarkBySpill[ch] = new TGraph();
          gHodDarkBySpill[ch]->SetNameTitle(Form("gExtDarkBySpill_%03ld", ch),
                                            Form("%s, Extinction Detector Dark Count by Spill @ ch%ld;"
                                                 "Spill;"
                                                 "Entries", tdcName.data(), ch));
          gHodDarkBySpill[ch]->SetMarkerStyle(kOpenCircle);
          gHodDarkBySpill[ch]->SetMarkerColor(kBlack);
          gHodDarkBySpill[ch]->SetMarkerSize(1);
        }

        TGraph* gHodDarkBySpill_Any = new TGraph();
        gHodDarkBySpill_Any->SetNameTitle("hHodDarkBySpill_Any",
                                          Form("%s, Hodoscope Dark Count by Spill;"
                                               "Spill;"
                                               "Entries", tdcName.data()));
        gHodDarkBySpill_Any->SetMarkerStyle(kOpenCircle);
        gHodDarkBySpill_Any->SetMarkerColor(kBlack);
        gHodDarkBySpill_Any->SetMarkerSize(1);

        // Beamline Hodoscope dark count
        TGraph** gBhDarkBySpill = new TGraph*[BeamlineHodoscope::NofChannels];
        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          gBhDarkBySpill[ch] = new TGraph();
          gBhDarkBySpill[ch]->SetNameTitle(Form("hBhDarkBySpill_%ld", ch),
                                           Form("%s, BH%ld Dark Count by Spill;"
                                                "Spill;"
                                                "Entries", tdcName.data(), ch + 1));
          gBhDarkBySpill[ch]->SetMarkerStyle(kOpenCircle);
          gBhDarkBySpill[ch]->SetMarkerColor(kBlack);
          gBhDarkBySpill[ch]->SetMarkerSize(1);
        }

        // Timing Counter dark count
        TGraph** gTcDarkBySpill = new TGraph*[TimingCounter::NofChannels];
        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          gTcDarkBySpill[ch] = new TGraph();
          gTcDarkBySpill[ch]->SetNameTitle(Form("hTcDarkBySpill_%ld", ch),
                                           Form("%s, TC%ld Dark Count by Spill;"
                                                "Spill;"
                                                "Entries", tdcName.data(), ch + 1));
          gTcDarkBySpill[ch]->SetMarkerStyle(kOpenCircle);
          gTcDarkBySpill[ch]->SetMarkerColor(kBlack);
          gTcDarkBySpill[ch]->SetMarkerSize(1);
        }

        std::cout << "=== Create history " << std::endl;
        std::map<std::size_t/*GlobalChannel*/, Long64_t/*DarkCount*/> darkCounts;
        std::map<std::size_t/*GlobalChannel*/, Long64_t/*DarkCount*/> totalCounts;

        Double_t extDarkCount[ExtinctionDetector::NofChannels] = { 0 };
        Double_t extDarkCount_Any = 0;
        Double_t hodDarkCount[Hodoscope::NofChannels] = { 0 };
        Double_t hodDarkCount_Any = 0;
        Double_t  bhDarkCount[BeamlineHodoscope::NofChannels] = { 0 };
        Double_t  tcDarkCount[TimingCounter::NofChannels] = { 0 };

        std::cout << "=== Get Entry " << std::endl;
        Int_t   lastSpill = -1;
        Int_t   spillCount = 0;
        for (Long64_t entry = 0; entry < entries; ++entry) {
          if (entry % 1000000 == 0) {
            std::cout << ">> " << entry << std::endl;
          }

          tree->GetEntry(entry);

          if (!provider->IsData()) {
            continue;
          }

          if (provider->GetSpill() != lastSpill) {
            for (auto&& pair : totalCounts) {
              DarkInfo {
                spillCount,
                pair.first,
                darkCounts[pair.first],
                totalCounts[pair.first],  
                }.Write(ofile);
            }

            if (lastSpill > 0) {
              ++spillCount;

              for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
                if (extDarkCount[ch]) {
                  gExtDarkBySpill[ch]->SetPoint(gExtDarkBySpill[ch]->GetN(), spillCount, extDarkCount[ch]);
                }
              }
              if (extDarkCount_Any) {
                gExtDarkBySpill_Any->SetPoint(gExtDarkBySpill_Any->GetN(), spillCount, extDarkCount_Any);
              }
              for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
                if (hodDarkCount[ch]) {
                  gHodDarkBySpill[ch]->SetPoint(gHodDarkBySpill[ch]->GetN(), spillCount, hodDarkCount[ch]);
                }
              }
              if (hodDarkCount_Any) {
                gHodDarkBySpill_Any->SetPoint(gHodDarkBySpill_Any->GetN(), spillCount, hodDarkCount_Any);
              }
              for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
                if (bhDarkCount[ch]) {
                  gBhDarkBySpill[ch]->SetPoint(gBhDarkBySpill[ch]->GetN(), spillCount, bhDarkCount[ch]);
                }
              }
              for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
                if (tcDarkCount[ch]) {
                  gTcDarkBySpill[ch]->SetPoint(gTcDarkBySpill[ch]->GetN(), spillCount, tcDarkCount[ch]);
                }
              }
            }

            darkCounts.clear();
            totalCounts.clear();

            for (std::size_t ch = 0; ch < ExtinctionDetector::NofChannels; ++ch) {
              extDarkCount[ch] = 0;
            }
            extDarkCount_Any = 0;
            for (std::size_t ch = 0; ch < Hodoscope::NofChannels; ++ch) {
              hodDarkCount[ch] = 0;
            }
            hodDarkCount_Any = 0;
            for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
              bhDarkCount[ch] = 0;
            }
            for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
              tcDarkCount[ch] = 0;
            }
          }

          std::vector<TdcData> tdcData = provider->GetTdcData();

          for (auto&& data : tdcData) {
            const Int_t    globalChannel = data.Channel;
            const Double_t time          = data.Time;

            if (tmin <= time && time <= tmax) {
              ++darkCounts[globalChannel];
            }
            ++totalCounts[globalChannel];

            if (ExtinctionDetector::Contains(globalChannel)) {
              const Int_t ch = ExtinctionDetector::GetChannel(globalChannel);

              if (tmin <= time && time <= tmax) {
                ++extDarkCount[ch];
                ++extDarkCount_Any;
              }

            } else if (Hodoscope::Contains(globalChannel)) {
              const Int_t ch = Hodoscope::GetChannel(globalChannel);

              if (tmin <= time && time <= tmax) {
                ++hodDarkCount[ch];
                ++hodDarkCount_Any;
              }

            } else if (BeamlineHodoscope::Contains(globalChannel)) {
              const Int_t ch = BeamlineHodoscope::GetChannel(globalChannel);

              if (tmin <= time && time <= tmax) {
                ++bhDarkCount[ch];
              }

            } else if (MrSync::Contains(globalChannel)) {
              const Int_t ch = MrSync::GetChannel(globalChannel);

              if (tmin <= time && time <= tmax) {
                ++tcDarkCount[ch];
              }

            }

          }

          lastSpill = provider->GetSpill();
        }

        return
          {
           gExtDarkBySpill,
           gExtDarkBySpill_Any,
           gHodDarkBySpill,
           gHodDarkBySpill_Any,
           gBhDarkBySpill,
           gTcDarkBySpill,
          };
      }

    }

  }

}

#endif
