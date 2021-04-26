#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDatime.h"
#include "TAxis.h"

#include "Linq.hh"
#include "String.hh"

class GL820DataReader {
protected:
  std::string              fTime;
  std::vector<std::string> fTemp;
public:
  GL820DataReader(std::size_t n) {
    fTemp.resize(n, "-200");
  }
  virtual void        ReadLine(const std::string& line) {
    std::vector<std::string> strs = Tron::String::Split(line, ",");
    if (strs.size() < 3 + fTemp.size()) {
      std::cout << "[warning] invalid line" << std::endl;
      fTime = "";
      fTemp.assign(fTemp.size(), "-200");
      return;
    }

    //      strs[0]; // No
    fTime = Tron::String::Replace(strs[1], "/", "-"); // Time (string)
    //      strs[2]  // ms
    for (std::size_t i = 0, n = fTemp.size(); i < n; ++i) {
      fTemp[i] = strs[3 + i];
    }
  }
  virtual std::string GetTime() const { return fTime; };
  virtual Double_t    GetTemperature(std::size_t i) {
    return Tron::String::Convert<Double_t>(fTemp[i]);
  }
  virtual std::vector<Double_t> GetTemperatures() {
    return Tron::Linq::From(fTemp)
      .Select([] (const std::string& temp) {
        return Tron::String::Convert<Double_t>(temp);
      })
      .ToVector();
  }
  const std::size_t size() const {
    return fTemp.size();
  }
};

struct GL820Data {
  GL820DataReader*               reader;
  const std::string              ifilename;
  const std::string              ofilename;
};

Int_t plotGL820() {
  std::vector<GL820Data> data
    {
      { new GL820DataReader(3),
        "./data_gl820/210423-175735_UG.CSV",
        "./data_gl820/210423-175735_UG.root" },
      { new GL820DataReader(3),
        "./data_gl820/210423-191135_UG.CSV",
        "./data_gl820/210423-191135_UG.root" },
    };

  const Double_t grad[3] = {
    0.987324,
    0.977127,
    0.977441,
  };
  const Double_t intercept[3] = {
    -5.29515,
    -5.27712,
    -5.09121,    
  };
  auto calibrate = [&](Double_t raw, std::size_t i){
    return raw * grad[i] + intercept[i];
  };

  for (auto&& datum : data) {
    TFile* ofile = new TFile(datum.ofilename.data(), "RECREATE");
    if (!ofile->IsOpen()) {
      std::cout << "[error] output file is not opened" << std::endl;
      return 0;
    }

    TTree* otree = new TTree("tree", "temperature");

    TDatime datime;
    otree->Branch("datime", &datime);

    Double_t raw[3] = { 0.0 };
    otree->Branch("raw", raw, "raw[3]/D");

    Double_t calib[3] = { 0.0 };
    otree->Branch("calib", calib, "calib[3]/D");

    std::vector<TGraph*> gTemps;
    gTemps.resize(datum.reader->size(), nullptr);
    for (std::size_t i = 0, n = gTemps.size(); i < n; ++i) {
      TGraph* gTemp = gTemps[i] = new TGraph();
      gTemp->SetName(Form("gTemp%c", 'A' + (char)i));
      gTemp->SetMinimum( 0);
      gTemp->SetMaximum(30);
      gTemp->SetMarkerStyle(kFullCircle);
      gTemp->SetMarkerColor(51 + (99 - 51) * i / n);
      gTemp->GetXaxis()->SetTimeDisplay(true);
      gTemp->GetXaxis()->SetTimeOffset(0);
      gTemp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gTemp->GetXaxis()->SetLabelOffset(0.03);
    }

    std::ifstream file(datum.ifilename);
    if (!file) {
      std::cout << "[warning] file is not opened, " << datum.ifilename << std::endl;
      continue;
    }

    // skip header
    std::string buff;
    for (Int_t i = 0; i < 14 + datum.reader->size(); ++i) {
      std::getline(file, buff);
    }

    // read data
    while (std::getline(file, buff)) {
      datum.reader->ReadLine(buff);
      datime.Set(datum.reader->GetTime().data());

      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        raw[i] = calib[i] = std::numeric_limits<Double_t>::min();
        const Double_t temp = datum.reader->GetTemperature(i);
        if (std::isfinite(raw[i])) {
          gTemps[i]->SetPoint(gTemps[i]->GetN(), datime.Convert(), temp);
          raw  [i] = temp;
          calib[i] = calibrate(raw[i], i);
        }
      }
      otree->Fill();
    }

    TCanvas::MakeDefCanvas();
    for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
      gTemps[i]->Draw(i ? "P" : "AP");
      gTemps[i]->Write();
    }

    otree->Write();
    ofile->Close();
  }

  return 0;
}
