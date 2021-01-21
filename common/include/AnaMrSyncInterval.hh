#ifndef Extinction_AnaMrSyncInterval_hh
#define Extinction_AnaMrSyncInterval_hh

#include <iostream>
#include <fstream>
#include <vector>
#include "Units.hh"
#include "AnaBeam.hh"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPaveStats.h"

namespace Extinction {

  namespace Analyzer {

    namespace AnaMrSyncInterval {

      struct IntervalInfo {
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

      struct Results_t {
        TH1** hMrSyncTdcInterval; // [MrSync::NofChannels]

        void Print(const std::string& ofilename) {
          TCanvas::MakeDefCanvas();
          gPad->SetGrid();

          gPad->Print((ofilename + "[").data());

          gPad->SetLogy(true);
          for (std::size_t i = 0; i < MrSync::NofChannels; ++i) {
            if (hMrSyncTdcInterval[i]->Integral()) {
              hMrSyncTdcInterval[i]->Draw(); 
              gPad->Update();
              if (TPaveStats* st = dynamic_cast<TPaveStats*>(hMrSyncTdcInterval[i]->FindObject("stats"))) {
                const Double_t x1 = st->GetX1();
                const Double_t x2 = st->GetX2();
                const Double_t y1 = st->GetY1();
                const Double_t y2 = st->GetY2();
                st->SetY1(y2 - (y2 - y1) * 0.6);
                st->SetX1(x2 - (x2 - x1) * 0.6);
                gPad->Modified();
                gPad->Update();
              }
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->Print((ofilename + "]").data());
        }

        void Write(TFile* file) {
          file->cd();

          for (std::size_t i = 0; i < MrSync::NofChannels; ++i) {
            if (hMrSyncTdcInterval[i]->Integral()) {
              hMrSyncTdcInterval[i]->Write(); 
            }
          }
        }
      };

      Results_t Execute(TH1** hMrSyncTdcInterval,
                        const std::string& ofilename) {
        if (!hMrSyncTdcInterval) {
          std::cout << "[error] invalid hist" << std::endl;
          return { nullptr };
        }

        TF1* fGauss = nullptr;
        if (!(fGauss = dynamic_cast<TF1*>(gROOT->FindObject("fGauss")))) {
          fGauss = new TF1("fGauss", "gaus(0)");
          fGauss->SetParNames("Constant", "Mean", "Sigma");
          fGauss->SetNpx(200);
        }

        std::ofstream ofile(ofilename);
        if (!ofile.is_open()) {
          std::cout << "[error] output file is not opened, " << ofilename  << std::endl;
          return Results_t { nullptr };
        }

        for (std::size_t i = 0; i < MrSync::NofChannels; ++i) {
          if (hMrSyncTdcInterval[i]->Integral()) {
            const Int_t    maxbin = hMrSyncTdcInterval[i]->GetMaximumBin();
            const Double_t maxx   = hMrSyncTdcInterval[i]->GetBinCenter(maxbin);
            const Double_t maxy   = hMrSyncTdcInterval[i]->GetBinContent(maxbin);
            const Double_t width  = hMrSyncTdcInterval[i]->GetXaxis()->GetXmax() - hMrSyncTdcInterval[i]->GetXaxis()->GetXmin();
            const Double_t sigma  = hMrSyncTdcInterval[i]->GetStdDev();

            // // [1] fixing fitting (sigma & offset)
            // fGauss->FixParameter(0, maxy);
            // fGauss->FixParameter(1, maxx);
            // fGauss->SetParameter(2, TMath::Max(0.5 / TMath::Sqrt(12.0), hMrSyncTdcInterval[i]->GetNbinsX() * 0.1));
            // fGauss->SetParLimits(2,            0.5 / TMath::Sqrt(12.0), hMrSyncTdcInterval[i]->GetNbinsX() * 0.3 );
            // hMrSyncTdcInterval[i]->Fit(fGauss, "Q", "goff");
            // const Double_t rSigma     = fGauss->GetParameter(2);
            // const Double_t rSigmaErr  = fGauss->GetParError(2);

            // // [2] fixing fitting (content & mean)
            // fGauss->SetParLimits(0, 0.5 * maxy, 5.0 * maxy);
            // fGauss->SetParLimits(1, maxx - 0.1 * width, maxx + 0.1 * width);
            // fGauss->FixParameter(2, fGauss->GetParameter(2));
            // hMrSyncTdcInterval[i]->Fit(fGauss, "Q", "goff");
            // const Double_t rConst    = fGauss->GetParameter(0);
            // const Double_t rConstErr = fGauss->GetParError(0);
            // const Double_t rMean     = fGauss->GetParameter(1);
            // const Double_t rMeanErr  = fGauss->GetParError(1);

            // // [3] free fitting
            // fGauss->SetParLimits(0, rConst  - 3.0 * rConstErr, rConst + 3.0 * rConstErr);
            // fGauss->SetParLimits(1, rMean   - 3.0 * rMeanErr , rMean  + 3.0 * rMeanErr );
            // fGauss->SetParLimits(2, rSigma  - 3.0 * rSigmaErr, rSigma + 3.0 * rSigmaErr);

            fGauss->SetParameters(maxy, maxx, sigma);
            fGauss->SetParLimits(0, 0.5 * maxy, 5.0 * maxy);
            fGauss->SetParLimits(1, maxx - 0.1 * width, maxx + 0.1 * width);
            fGauss->SetParLimits(2, 0.5 / TMath::Sqrt(12.0), width * 0.3);
            hMrSyncTdcInterval[i]->Fit(fGauss, "Q", "goff");

            IntervalInfo {
              i,
              maxx,
              fGauss->GetParameter(1),
              fGauss->GetParameter(2),
            } .Write(ofile);

            fGauss->ReleaseParameter(0);
            fGauss->ReleaseParameter(1);
            fGauss->ReleaseParameter(2);
          }
        }
        ofile.close();

        return Results_t { hMrSyncTdcInterval };
      }

    }

  }
  
}

#endif
