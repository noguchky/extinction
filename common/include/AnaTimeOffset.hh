#ifndef Extinction_AnaTimeOffset_hh
#define Extinction_AnaTimeOffset_hh

#include <iostream>
#include <fstream>
#include <vector>
#include "Units.hh"
#include "Utility.hh"
#include "AnaBeam.hh"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

namespace Extinction {

  namespace Analyzer {

    namespace AnaTimeOffset {

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

      struct Results_t {
        TH1** hExtTdcCoinY; // [AnaBeam::CoinOffset::N]

        void Print(const std::string& ofilename) {
          TCanvas::MakeDefCanvas();
          gPad->SetGrid();

          gPad->Print((ofilename + "[").data());

          gPad->SetLogy(true);
          for (std::size_t i = 0; i < CoinOffset::N; ++i) {
            if (hExtTdcCoinY[i]->Integral()) {
              hExtTdcCoinY[i]->Draw();
              Utility::ResizeStats(hExtTdcCoinY[i]);
              gPad->Print(ofilename.data());
            }
          }
          gPad->SetLogy(false);

          gPad->Print((ofilename + "]").data());
        }

        void Write(TFile* file) {
          file->cd();

          for (std::size_t i = 0; i < CoinOffset::N; ++i) {
            if (hExtTdcCoinY[i]->Integral()) {
              hExtTdcCoinY[i]->Write(); 
            }
          }
        }
      };

      Results_t Execute(TH2* hExtTdcCoin,
                        const std::string& ofilename) {
        if (!hExtTdcCoin) {
          std::cout << "[error] invalid hist" << std::endl;
          return { nullptr };
        }

        TF1* fGaussLinear = nullptr;
        if (!(fGaussLinear = dynamic_cast<TF1*>(gROOT->FindObject("fGaussLinear")))) {
          fGaussLinear = new TF1("fGaussLinear", "gaus(0)+pol0(3)");
          fGaussLinear->SetParNames("Constant", "Mean", "Sigma", "Offset");
          fGaussLinear->SetNpx(200);
        }

        TH1** hExtTdcCoinY = new TH1*[CoinOffset::N];

        std::ofstream ofile(ofilename);
        if (!ofile.is_open()) {
          std::cout << "[error] output file is not opened, " << ofilename  << std::endl;
          return Results_t { nullptr };
        }

        for (std::size_t ch = 0; ch < BeamlineHodoscope::NofChannels; ++ch) {
          const std::size_t i   = ch + AnaBeam::CoinOffset::BH;
          const Int_t       bin = ch + AnaBeam::CoinOffset::BH + 1;
          hExtTdcCoinY[i] = hExtTdcCoin->ProjectionY(Form("%s_BH%ld", hExtTdcCoin->GetName(), ch + 1), bin, bin);
        }

        {
          const std::size_t i    = AnaBeam::CoinOffset::Hod;
          const Int_t       bin1 = AnaBeam::CoinOffset::Hod + 1;
          // const Int_t       bin2 = AnaBeam::CoinOffset::Hod + Hodoscope::NofChannels;
          // hExtTdcCoinY[i] = hExtTdcCoin->ProjectionY(Form("%s_HodOr", hExtTdcCoin->GetName()), bin1, bin2);
          hExtTdcCoinY[i] = hExtTdcCoin->ProjectionY(Form("%s_HodOr", hExtTdcCoin->GetName()), bin1, bin1); // for 20201217 (bin1 = HodOr)
        }

        for (std::size_t ch = 0; ch < TimingCounter::NofChannels; ++ch) {
          const std::size_t i   = ch + AnaBeam::CoinOffset::Hod + 1;
          const Int_t       bin = ch + AnaBeam::CoinOffset::TC + 1;
          hExtTdcCoinY[i] = hExtTdcCoin->ProjectionY(Form("%s_TC%ld", hExtTdcCoin->GetName(), ch + 1), bin, bin);
        }

        for (std::size_t i = 0; i < CoinOffset::N; ++i) {
          if (hExtTdcCoinY[i]->Integral()) {
            const Int_t    maxbin = hExtTdcCoinY[i]->GetMaximumBin();
            const Double_t maxx   = hExtTdcCoinY[i]->GetBinCenter(maxbin);
            const Double_t maxy   = hExtTdcCoinY[i]->GetBinContent(maxbin);
            const Double_t width  = hExtTdcCoinY[i]->GetXaxis()->GetXmax() - hExtTdcCoinY[i]->GetXaxis()->GetXmin();
            const Double_t firsty = hExtTdcCoinY[i]->GetBinContent(1);
            const Double_t lasty  = hExtTdcCoinY[i]->GetBinContent(hExtTdcCoinY[i]->GetNbinsX());
            // [1] fixing fitting (sigma & offset)
            fGaussLinear->FixParameter(0, maxy);
            fGaussLinear->FixParameter(1, maxx);
            fGaussLinear->SetParameter(2, TMath::Max(0.5 / TMath::Sqrt(12.0), hExtTdcCoinY[i]->GetNbinsX() * 0.1));
            fGaussLinear->SetParLimits(2,            0.5 / TMath::Sqrt(12.0), hExtTdcCoinY[i]->GetNbinsX() * 0.3 );
            fGaussLinear->SetParameter(3, 0.5 * (firsty + lasty));
            fGaussLinear->SetParLimits(3, 0.0, maxy);
            hExtTdcCoinY[i]->Fit(fGaussLinear, "QI", "goff");
            const Double_t rSigma     = fGaussLinear->GetParameter(2);
            const Double_t rSigmaErr  = fGaussLinear->GetParError(2);
            const Double_t rOffset    = fGaussLinear->GetParameter(3);
            const Double_t rOffsetErr = fGaussLinear->GetParError(3);

            // [2] fixing fitting (content & mean)
            fGaussLinear->SetParLimits(0, 0.5 * maxy, 5.0 * maxy);
            fGaussLinear->SetParLimits(1, maxx - 0.1 * width, maxx + 0.1 * width);
            fGaussLinear->FixParameter(2, fGaussLinear->GetParameter(2));
            fGaussLinear->FixParameter(3, fGaussLinear->GetParameter(3));
            hExtTdcCoinY[i]->Fit(fGaussLinear, "QI", "goff");
            const Double_t rConst    = fGaussLinear->GetParameter(0);
            const Double_t rConstErr = fGaussLinear->GetParError(0);
            const Double_t rMean     = fGaussLinear->GetParameter(1);
            const Double_t rMeanErr  = fGaussLinear->GetParError(1);

            // [3] free fitting
            fGaussLinear->SetParLimits(0, rConst  - 3.0 * rConstErr , rConst  + 3.0 * rConstErr );
            fGaussLinear->SetParLimits(1, rMean   - 3.0 * rMeanErr  , rMean   + 3.0 * rMeanErr  );
            fGaussLinear->SetParLimits(2, rSigma  - 3.0 * rSigmaErr , rSigma  + 3.0 * rSigmaErr );
            fGaussLinear->SetParLimits(3, rOffset - 3.0 * rOffsetErr, rOffset + 3.0 * rOffsetErr);
            hExtTdcCoinY[i]->Fit(fGaussLinear, "QI", "goff");
            CoinInfo {
              i,
              maxx,
              fGaussLinear->GetParameter(1),
              fGaussLinear->GetParameter(2),
            } .Write(ofile);

            fGaussLinear->ReleaseParameter(0);
            fGaussLinear->ReleaseParameter(1);
            fGaussLinear->ReleaseParameter(2);
            fGaussLinear->ReleaseParameter(3);
          }
        }
        ofile.close();

        return Results_t { hExtTdcCoinY };
      }

    }

  }
  
}

#endif
