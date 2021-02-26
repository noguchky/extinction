#include <iostream>
#include <fstream>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

#include "Math.hh"
#include "ProcessLine.hh"

namespace TEKAFG {
  const std::size_t MinHeight = 0x0000;
  const std::size_t MaxHeight = 0x3ffe;
  namespace NIM {
    const std::size_t True  = MinHeight;
    const std::size_t False = MaxHeight;
  }
  namespace TTL {
    const std::size_t True  = MaxHeight;
    const std::size_t False = MinHeight;
  }

  std::vector<UShort_t> CreateHeader(const std::vector<UShort_t> data) {
    std::vector<UShort_t> header
      {
       0x5445u, 0x4b41u, 0x4647u, 0x3330u, 0x3030u, 0x0000u, 0x0000u, 0x0000u,
       0x0131u, 0xf0c2u,    
      };

    header.push_back(0x0000u);
    header.push_back(data.size());

    header.push_back(0x0000u);
    header.push_back(0x0001u);

    for (std::size_t i = 0; i < 206; ++i) {
      header.push_back(0x0000u);
    }

    for (std::size_t i = 0; i < 36; ++i) {
      header.push_back(0x0000u);
    }

    return header;
  }

  void WriteData(const std::string&           filename,
                 const std::vector<UShort_t>& header,
                 const std::vector<UShort_t>& data) {
    // Merge data
    std::vector<UShort_t> fileData;
    for (auto&& buff : header) {
      fileData.push_back(Tron::Math::RevertEndian(buff));
    }
    for (auto&& buff : data) {
      fileData.push_back(Tron::Math::RevertEndian(buff));
    }

    // Write data
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
      std::cout << "[error] file is not opened, " << filename << std::endl;
    } else {
      file.write((char*)fileData.data(), sizeof(UShort_t) * fileData.size());
    }
  }

  std::vector<UShort_t> Create8GeVNimMask(std::size_t dataSize, Double_t bunchFullWidth, std::size_t outOfBunchLevel, std::size_t inBunchLevel) {
    const Double_t keV = 1.0;
    const Double_t MeV = 1.0e3 * keV;
    const Double_t GeV = 1.0e3 * MeV;

    const Double_t protonMc2     = 938.2720813 * MeV;
    const Double_t kineticEnergy = 8.0 * GeV;
    const Double_t totalEnergy   = protonMc2 + kineticEnergy;
    const Double_t gamma         = totalEnergy / protonMc2;
    const Double_t beta          = TMath::Sqrt((gamma * gamma - 1.0) / (gamma * gamma));
    const Double_t velocity      = beta * 2.99792458e8; // meter/second

    const Double_t mrCircumference = 1567.5; // meter
    const Double_t mrSync          = mrCircumference / velocity * 1e9; // nsec
    const Double_t interval        = mrSync * 2.0 / 9.0;
    const Double_t timePerData     = mrSync / dataSize; // nsec/data

    std::cout << "MrSync         " << mrSync   << " nsec" << std::endl
              << "Bunch Interval " << interval << " nsec" << std::endl;

    std::vector<UShort_t> data;
    for (std::size_t x = 0, y = outOfBunchLevel, istep = 0; x < dataSize; ++x) {
      const Double_t t = x * timePerData;

      y = outOfBunchLevel;
      for (std::size_t bunch = 0; bunch < 4; ++bunch) {
        if (bunch * interval < t && t < bunch * interval + bunchFullWidth) {
          y = inBunchLevel;
          break;
        }
      }

      data.push_back(y);
    }

    return data;
  }

}

Int_t generateWaveform(const std::string filename = "TEKAFG00.TFW") {
  if (filename.empty()) {
    std::cout << "[error] invalid filename" << std::endl;
    return 1;
  }

  // Common Parameters
  const std::size_t bunchFullWidth = 300; // nsec

  // // Veto Parameters
  // const std::size_t    inBunchLevel = TEKAFG::NIM::True;
  // const std::size_t outOfBunchLevel = TEKAFG::NIM::False;

  // Enable Parameters
  const std::size_t    inBunchLevel = TEKAFG::NIM::False;
  const std::size_t outOfBunchLevel = TEKAFG::NIM::True;
  
  // Create waveform
  auto data = TEKAFG::Create8GeVNimMask(1000, 300, outOfBunchLevel, inBunchLevel);

  // Create header
  auto header = TEKAFG::CreateHeader(data);

  // Write data
  TEKAFG::WriteData(filename, header, data);

  // Draw waveform
  TCanvas::MakeDefCanvas();
  Tron::ProcessLine::From(new TGraph())
    .SetName(filename)
    .SetTitle(filename + ";;;")
    .SetLineStyle(kRed, kSolid, 2)
    .CallFor(data, [](TGraph* g, UShort_t value) {
                     g->SetPoint(g->GetN(), g->GetN(), value);
                   })
    .Draw("AL")
    .SetRangeX(0, data.size() - 1);

  return 0;
}

