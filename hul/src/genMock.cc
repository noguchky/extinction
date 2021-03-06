#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TRandom.h"

#include "Linq.hh"

#include "Tdc.hh"
#include "Detector.hh"
#include "Hul.hh"

namespace {
  using namespace Extinction;
  namespace K18BR {
    namespace X {
      Double_t Mu     = -6.722 * cm;
      Double_t SigmaP = 14.37  * cm;
      Double_t SigmaN =  6.621 * cm;
    }
    namespace Y {
      Double_t Mu     =  0.0   * cm;
      Double_t Sigma  =  3.574 * cm;
    }
  }

  Double_t       extinction            =    3.0e-11;
  Int_t          eventMatchNumber      = 1054;
  Int_t          eventMatchParity;

  Double_t       dtEM                  =   10    * nsec;
  Double_t       dtBH1                 =  100    * nsec;
  Double_t       dtBH2                 =  120    * nsec;
  Double_t       dtHod                 =  110    * nsec;
  Double_t       dtTC1                 =  140    * nsec;
  Double_t       dtTC2                 =  140    * nsec;
  Double_t       sigmaT                =    1    * nsec;
  Double_t       bunchSigma            =   15.0  * nsec;

  const Double_t nofParticles          =    1.6e11;
  const Double_t daqTime               = 1.0 * 24 * 60 * 60 * sec;

  const Double_t cycle                 =    5.52 *  sec;
  const Double_t spillLength           =    0.5  *  sec;
  const Double_t mrSyncInterval        = 5257.67 * nsec;
  const Double_t bunchInterval         = 1168.37 * nsec;
  const Double_t bunchT0               =  200.0  * nsec;

  const Double_t dataLength            = 1.5 * sec;
  const Double_t extractT0             = 0.5 * sec;
  const Double_t nofMrSync             = dataLength / mrSyncInterval;
  const Double_t nofMrSyncInSpill      = spillLength / mrSyncInterval;

  const Double_t nofSpill              = daqTime / cycle;
  const Double_t nofParticlesPerSpill  = nofParticles / nofSpill;

  const Double_t nofBunchPerMrSync     = 4;
  const Double_t nofParticlesPerMrSync = nofParticlesPerSpill / nofMrSyncInSpill;
  const Double_t nofParticlesPerBunch  = nofParticlesPerMrSync/ nofBunchPerMrSync;
}

Int_t main(Int_t argc, Char_t** argv) {
  if (argc < 2) {
    std::cout << "Usage: ./genMock [ConfFilename]" << std::endl;
    return 1;
  }

  TApplication* app = new TApplication("app", nullptr, nullptr);

  const std::string confFilename = argv[1];

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  const auto boards = conf->GetValues<Int_t>("Boards");
  {
    using namespace K18BR::X;
    Mu     = conf->GetValue<Double_t>("K18BR.X.Mu");
    SigmaP = conf->GetValue<Double_t>("K18BR.X.SigmaP");
    SigmaN = conf->GetValue<Double_t>("K18BR.X.SigmaN");
  }
  {
    using namespace K18BR::Y;
    Mu     = conf->GetValue<Double_t>("K18BR.Y.Mu");
    Sigma  = conf->GetValue<Double_t>("K18BR.Y.Sigma");
  }
  {
    extinction = conf->GetValue<Double_t>("Extinction");
    eventMatchNumber = conf->GetValue<Int_t>("EventMatchNumber");
    eventMatchParity = 0;
    for (std::size_t i = 0; i < 16; ++i) {
      eventMatchParity ^= (eventMatchNumber >> i) & 0x1;
    }

    dtBH1  = conf->GetValue<Double_t>("Offset.BH1");
    dtBH2  = conf->GetValue<Double_t>("Offset.BH2");
    dtHod  = conf->GetValue<Double_t>("Offset.Hod");
    dtTC1  = conf->GetValue<Double_t>("Offset.TC1");
    dtTC2  = conf->GetValue<Double_t>("Offset.TC2");

    sigmaT = conf->GetValue<Double_t>("TimeResolution");

    bunchSigma = conf->GetValue<Double_t>("BunchSigma");
  }

  Extinction::Hul::ChannelMapWithBoard::Load(conf, boards);

  // auto printChannels =
  //   [](const std::string& title, const std::map<Int_t, std::map<Int_t, Int_t>> map) {
  //     std::cout << title << std::endl;
  //     for (auto&& pair1 : map) {
  //       for (auto&& pair2 : pair1.second) {
  //         std::cout << pair1.first << "\t"
  //                   << pair2.first << "\t"
  //                   << pair2.second << std::endl;
  //       }
  //     }
  //   };
  // printChannels("BeamlineHodoscope" , Extinction::Hul::ChannelMapWithBoard::Bh    );
  // printChannels("ExtinctionDetector", Extinction::Hul::ChannelMapWithBoard::Ext   );
  // printChannels("Hodoscope"         , Extinction::Hul::ChannelMapWithBoard::Hod   );
  // printChannels("TimingCounter"     , Extinction::Hul::ChannelMapWithBoard::Tc    );
  // printChannels("MRSync"            , Extinction::Hul::ChannelMapWithBoard::MrSync);
  // printChannels("EventMatch"        , Extinction::Hul::ChannelMapWithBoard::Evm   );
  // return 0;

  TF1* fGaussX = new TF1("fGaussX", "[0]*TMath::Exp(-0.5*TMath::Sq((x-[1])/((x>=[1])*[2]+(x<[1])*[3])))", -40, 40);
  fGaussX->SetParNames("Constant", "Mu", "Sigma+", "Sigma-");
  fGaussX->SetNpx(80 * 10);
  fGaussX->SetParameters(1.0,
                         K18BR::X::Mu / Extinction::cm,
                         K18BR::X::SigmaP / Extinction::cm,
                         K18BR::X::SigmaN / Extinction::cm);

  TF1* fGaussY = new TF1("fGaussY", "[0]*TMath::Exp(-0.5*TMath::Sq((x-[1])/[2]))", -32, 32);
  fGaussY->SetParNames("Constant", "Mu", "Sigma");
  fGaussY->SetNpx(64 * 10);
  fGaussY->SetParameters(1.0,
                         K18BR::Y::Mu / Extinction::cm,
                         K18BR::Y::Sigma / Extinction::cm);

  TH2* hMap = new TH2D("hMap", "Hit Map;x [cm];y [cm]", 80, -40, 40, 64, -32, 32);

  // TH2* hCnl = new TH2D("hCnl", "Cnl Map;x [cm];y [cm]", 80, -40, 40, 64, -32, 32);
  // for (Int_t xbin = 1, xbins = hCnl->GetNbinsX(); xbin <= xbins; ++xbin) {
  //   for (Int_t ybin = 1, ybins = hCnl->GetNbinsY(); ybin <= ybins; ++ybin) {
  //     const Double_t x = hCnl->GetXaxis()->GetBinCenter(xbin);
  //     const Double_t y = hCnl->GetYaxis()->GetBinCenter(ybin);
  //     const Int_t ch = Extinction::ExtinctionDetector::FindChannel(x * Extinction::cm, y * Extinction::cm);
  //     std::cout << "x = " << x << ", y  = " << y << "  ->  " << ch << std::endl;
  //     hCnl->SetBinContent(xbin, ybin, (Double_t)ch);
  //   }
  // }
  // hCnl->Draw("colz");
  // gPad->Update();

  Extinction::Hul::HulData data;
  const Double_t timePerTdc = data.GetTimePerTdc();

  std::map<Int_t, std::ofstream> ofile;
  for (auto&& board : boards) {
    ofile[board].open(Form("hulmock_%d.dat", board), std::ios::binary);
    if (!ofile[board]) {
      std::cout << "[error] file is not opened" << std::endl;
      return 1;
    }
  }

  Int_t ch = 0;
  std::map<Int_t, Long64_t> heartbeat;

  for (auto&& board : boards) {
    heartbeat[board] = 1;
    data.WriteSpillStart(ofile[board]);
  }

  std::map<Int_t, std::map<Long64_t, std::set<Int_t>>> hits;
  std::cout << "nofParticles          = " << nofParticles          << std::endl;
  std::cout << "nofSpill              = " << nofSpill              << std::endl;
  std::cout << "nofMrSync             = " << nofMrSync             << std::endl;
  std::cout << "nofParticlesPerSpill  = " << nofParticlesPerSpill  << std::endl;
  std::cout << "nofParticlesPerMrSync = " << nofParticlesPerMrSync << std::endl;
  std::cout << "nofParticlesPerBunch  = " << nofParticlesPerBunch  << std::endl;
  std::cout << "timeResolution        = " << sigmaT / Extinction::nsec << " nsec" << std::endl;
  std::cout << "time/tdc              = " << timePerTdc / Extinction::nsec << " nsec" << std::endl;

  Long64_t particleCnt = 0, coinParticleCnt = 0;
  for (Long64_t mrSync = 0; mrSync < nofMrSync; ++mrSync) {
    // if (mrSync % 1000 == 0) {
    //   std::cout << " --- " << mrSync << std::endl;
    // }
    const Double_t t0 = (mrSync + 0.2) * mrSyncInterval;
    // std::cout << (Long64_t)(t0 / nsec) << std::endl;
    for (auto&& board : boards) {
      hits[board].clear();
    }
    for (auto&& board : boards) {
      const Long64_t tdc = TMath::Max(0.0, (t0 + gRandom->Gaus() * sigmaT) / timePerTdc);
      // std::cout << tdc << std::endl;
      for (auto&& pair : Extinction::Hul::ChannelMapWithBoard::MrSync[board]) {
        hits[board][tdc].insert(pair.first);
      }
    }

    if (mrSync >= 192 &&
        ((mrSync ==  0 + 192 || mrSync ==  1 + 192) ||
         (mrSync <  18 + 192 && ((eventMatchNumber >> (mrSync - 1)) & 0x1)) ||
         (mrSync == 18 + 192 && eventMatchParity))) {
      for (ch = 0; ch < (Int_t)Extinction::EventMatch::GlobalChannelOffset; ++ch) {
        if (Extinction::Hul::ChannelMapWithBoard::Board.count(ch + Extinction::EventMatch::GlobalChannelOffset)) {
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::EventMatch::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t0 + dtEM + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Evm[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
      }
    }

    if (extractT0 < t0 && t0 < extractT0 + spillLength) {
      // Extracting
      for (Long64_t residual = 0, residuals = gRandom->Poisson(nofParticlesPerBunch * extinction); residual < residuals; ++residual) {
        std::cout << " !! reakage" << std::endl;
        const Int_t bunch = gRandom->Integer(nofBunchPerMrSync);
        const Double_t mu = bunchT0 + (bunch + 0.5) * bunchInterval;
        const Double_t t = t0 + mu + gRandom->Gaus() * bunchSigma * 2.0;
        const Double_t x = fGaussX->GetRandom() * cm;
        const Double_t y = fGaussY->GetRandom() * cm;
        hMap->Fill(x / cm, y / cm);
        // BH1
        {
          ch = 0;
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtBH1 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Bh[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // BH2
        {
          ch = 1;
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtBH2 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Bh[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // HOD
        {
          ch = TMath::Max(Extinction::Hodoscope::FindChannel(x / 2, y), 0LL);
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::Hodoscope::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtHod + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Hod[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // EXT
        {
          ch = TMath::Max(Extinction::ExtinctionDetector::FindChannel(x, y), 0LL);
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::ExtinctionDetector::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Ext[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // TC1
        {
          ch = 0;
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtTC1 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Tc[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // TC2
        {
          ch = 1;
          const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtTC2 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Tc[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
          ++coinParticleCnt;
        }
      }
      for (Long64_t bunch = 0; bunch < nofBunchPerMrSync; ++bunch) {
        for (Long64_t particle = 0, particles = gRandom->Poisson(nofParticlesPerBunch * 2.1); particle < particles; ++particle) {
          if (++particleCnt % 200000 == 0) {
            std::cout << ">> " << particleCnt << std::endl;
          }
          const Double_t mu = bunchT0 + bunch * bunchInterval;
          const Double_t t = t0 + mu + gRandom->Gaus() * bunchSigma;
          const Double_t x = fGaussX->GetRandom() * cm;
          const Double_t y = fGaussY->GetRandom() * cm;
          hMap->Fill(x / cm, y / cm);
          const Double_t reach = gRandom->Uniform(0.0, 2.1);
          // BH1
          {
            ch = 0;
            const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtBH1 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Bh[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // BH2
          if (reach < 1.05) {
            ch = 1;
            const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtBH2 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Bh[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // HOD
          if (reach < 1.03 &&
              (ch = Extinction::Hodoscope::FindChannel(x / 2, y)) >= 0) {
            const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::Hodoscope::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtHod + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Hod[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // EXT
          if (reach < 1.02 &&
              (ch = Extinction::ExtinctionDetector::FindChannel(x, y)) >= 0) {
            const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::ExtinctionDetector::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Ext[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // TC1
          if (reach < 1.01) {
            ch = 0;
            const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtTC1 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Tc[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // TC2
          if (reach < 1.00) {
            ch = 1;
            const Int_t board = Extinction::Hul::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtTC2 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Hul::ChannelMapWithBoard::Tc[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
            ++coinParticleCnt;
          }
        }
      }
    } else {
      // No beam
    }

    for (auto&& hitsInBoard : hits) {
      const Int_t board = hitsInBoard.first;
      for (auto&& hitInTdc : hitsInBoard.second) {
        const Long64_t tdc = hitInTdc.first;
        const Long64_t currentHeartbeat = tdc >> 19;
        for (; heartbeat[board] <= currentHeartbeat; ++heartbeat[board]) {
          data.Heartbeat = heartbeat[board];
          data.WriteHeartbeat(ofile[board]);
        }

        for (auto&& hitCh : hitInTdc.second) {
          data.Tdc     = tdc;
          data.Channel = hitCh;
          data.WriteData(ofile[board]);
        }
      }
    }
  }

  for (auto&& board : boards) {
    data.WriteSpillEnd(ofile[board]);
    ofile[board].close();
  }

  std::cout << "coinParticleCnt       = " << coinParticleCnt << std::endl;

  hMap->Draw("col");
  
  app->Run();
  
  return 0;
}
