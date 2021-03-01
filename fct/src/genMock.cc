#include "TF1.h"
#include "TRandom.h"

#include "Linq.hh"

#include "Tdc.hh"
#include "Detector.hh"
#include "Fct.hh"

namespace {
  using namespace Extinction;
  namespace K18BR {
    namespace X {
      const Double_t Mu     = -6.722 * cm;
      const Double_t SigmaP = 14.37  * cm;
      const Double_t SigmaN =  6.621 * cm;
    }
    namespace Y {
      const Double_t Mu     =  0.0   * cm;
      const Double_t Sigma  =  3.574 * cm;
    }
  }

  const Double_t extinction            = 3.0e-11;
  const Double_t nofParticles          = 1.6e11;
  const Double_t daqTime               = 24 * 60 * 60 * sec;

  const Double_t cycle                 =    5.52 *  sec;
  const Double_t spillLength           =    0.5  *  sec;
  const Double_t mrSyncInterval        = 5257.67 * nsec;
  const Double_t bunchInterval         = 1168.37 * nsec;
  const Double_t bunchSigma            =   15.0  * nsec;
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

  const Double_t dtBH1 = 100 * nsec;
  const Double_t dtBH2 = 120 * nsec;
  const Double_t dtHod = 110 * nsec;
  const Double_t dtTC1 = 140 * nsec;
  const Double_t dtTC2 = 140 * nsec;
  const Double_t sigmaT =  1 * nsec;
}

Int_t main(Int_t argc, Char_t** argv) {
  if (argc < 2) {
    std::cout << "Usage: ./genMock [ConfFilename]" << std::endl;
    return 1;
  }
  
  const std::string confFilename = argv[1];

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  const auto boards = conf->GetValues<Int_t>("Boards");

  Extinction::Fct::ChannelMapWithBoard::Load(conf, boards);

  TF1* fGaussX = new TF1("fGaussX", "[0]*TMath::Exp(-0.5*TMath::Sq((x-[1])/((x>=[1])*[2]+(x<[1])*[3])))", -40, 40);
  fGaussX->SetParNames("Constant", "Mu", "Sigma+", "Sigma-");
  fGaussX->SetNpx(64*10);
  fGaussX->SetParameters(1.0,
                         K18BR::X::Mu / Extinction::cm,
                         K18BR::X::SigmaP / Extinction::cm,
                         K18BR::X::SigmaN / Extinction::cm);

  TF1* fGaussY = new TF1("fGaussY", "[0]*TMath::Exp(-0.5*TMath::Sq((x-[1])/[2]))", -32, 32);
  fGaussY->SetParNames("Constant", "Mu", "Sigma");
  fGaussY->SetNpx(64*10);
  fGaussY->SetParameters(1.0,
                         K18BR::Y::Mu / Extinction::cm,
                         K18BR::Y::Sigma / Extinction::cm);

  Extinction::Fct::FctData data;
  const Double_t timePerTdc = data.GetTimePerTdc();

  std::map<Int_t, std::ofstream> ofile;
  for (auto&& board : boards) {
    ofile[board].open(Form("fctmock_%d.dat", board), std::ios::binary);
    if (!ofile[board]) {
      std::cout << "[error] file is not opened" << std::endl;
      return 1;
    }
  }

  Int_t ch = 0;
  std::map<Int_t, Long64_t> carry;

  for (auto&& board : boards) {
    carry[board] = 1;
    data.WriteGateStart(ofile[board]);
  }

  std::map<Int_t, std::map<Long64_t, std::set<Int_t>>> hits;
  std::cout << "nofParticles          = " << nofParticles          << std::endl;
  std::cout << "nofSpill              = " << nofSpill              << std::endl;
  std::cout << "nofMrSync             = " << nofMrSync             << std::endl;
  std::cout << "nofParticlesPerSpill  = " << nofParticlesPerSpill  << std::endl;
  std::cout << "nofParticlesPerMrSync = " << nofParticlesPerMrSync << std::endl;
  std::cout << "nofParticlesPerBunch  = " << nofParticlesPerBunch  << std::endl;
  Long64_t particleCnt = 0;
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
      for (auto&& pair : Extinction::Fct::ChannelMapWithBoard::MrSync[board]) {
        hits[board][tdc].insert(pair.first);
      }
    }
    if (extractT0 < t0 && t0 < extractT0 + spillLength) {
      // Extracting
      for (Long64_t residual = 0; residual < gRandom->Poisson(nofParticlesPerBunch * extinction); ++residual) {
        std::cout << " !! reakage" << std::endl;
        const Int_t bunch = gRandom->Integer(nofBunchPerMrSync);
        const Double_t mu = bunchT0 + (bunch + 0.5) * bunchInterval;
        const Double_t t = t0 + mu + gRandom->Gaus() * bunchSigma * 2.0;
        const Double_t x = fGaussX->GetRandom() * cm;
        const Double_t y = fGaussY->GetRandom() * cm;
        const Double_t reach = gRandom->Uniform(0.0, 2.1);
        // BH1
        {
          ch = 0;
          const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtBH1 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Bh[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // BH2
        if (reach < 1.05) {
          ch = 1;
          const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtBH2 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Bh[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // HOD1
        if (reach < 1.03 &&
            (ch = Extinction::Hodoscope::FindChannel(x, y)) >= 0) {
          const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::Hodoscope::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtHod + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Hod[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // EXT
        if (reach < 1.02 &&
            (ch = Extinction::ExtinctionDetector::FindChannel(x, y)) >= 0) {
          const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::ExtinctionDetector::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + gRandom->Gaus() * sigmaT) / timePerTdc);

          const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Ext[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // TC1
        if (reach < 1.01) {
          ch = 0;
          const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtTC1 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Tc[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
        // TC2
        if (reach < 1.00) {
          ch = 1;
          const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
          const Long64_t tdc = TMath::Max(0.0, (t + dtTC2 + gRandom->Gaus() * sigmaT) / timePerTdc);
          const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Tc[board])
            .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
            .FirstOrDefault()
            .first;
          hits[board][tdc].insert(rawCh);
        }
      }
      for (Long64_t bunch = 0; bunch < nofBunchPerMrSync; ++bunch) {
        for (Long64_t particle = 0; particle < nofParticlesPerBunch * 2.1; ++particle) {
          if (++particleCnt % 200000 == 0) {
            std::cout << ">> " << particleCnt << std::endl;
          }
          const Double_t mu = bunchT0 + bunch * bunchInterval;
          const Double_t t = t0 + mu + gRandom->Gaus() * bunchSigma;
          const Double_t x = fGaussX->GetRandom() * cm;
          const Double_t y = fGaussY->GetRandom() * cm;
          const Double_t reach = gRandom->Uniform(0.0, 2.1);
          // BH1
          {
            ch = 0;
            const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtBH1 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Bh[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // BH2
          if (reach < 1.05) {
            ch = 1;
            const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::BeamlineHodoscope::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtBH2 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Bh[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // HOD1
          if (reach < 1.03 &&
              (ch = Extinction::Hodoscope::FindChannel(x, y)) >= 0) {
            const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::Hodoscope::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtHod + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Hod[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // EXT
          if (reach < 1.02 &&
              (ch = Extinction::ExtinctionDetector::FindChannel(x, y)) >= 0) {
            const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::ExtinctionDetector::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + gRandom->Gaus() * sigmaT) / timePerTdc);
            // if (ch == 50) {
            //   std::cout << (Long64_t)(t/nsec) << "\t" << tdc << std::endl;
            // }
            const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Ext[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // TC1
          if (reach < 1.01) {
            ch = 0;
            const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtTC1 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Tc[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
          }
          // TC2
          if (reach < 1.00) {
            ch = 1;
            const Int_t board = Extinction::Fct::ChannelMapWithBoard::Board[ch + Extinction::TimingCounter::GlobalChannelOffset];
            const Long64_t tdc = TMath::Max(0.0, (t + dtTC2 + gRandom->Gaus() * sigmaT) / timePerTdc);
            const Long64_t rawCh = Tron::Linq::From(Extinction::Fct::ChannelMapWithBoard::Tc[board])
              .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == ch; })
              .FirstOrDefault()
              .first;
            hits[board][tdc].insert(rawCh);
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
        const Long64_t currentCarry = tdc >> 24;
        // if (board == 0) {
        //   std::cout << Form("%30llx", tdc) << "\t" << currentCarry << std::endl;
        // }
        for (; carry[board] <= currentCarry; ++carry[board]) {
          data.Carry = carry[board];
          data.WriteCarry(ofile[board]);
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
    data.WriteGateEnd(ofile[board]);
  }
  
  return 0;
}