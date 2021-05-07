#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TGraph.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TLegend.h"

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
  const std::string              pfilename;
  const std::string              datemin;
  const std::string              datemax;
};

Int_t plotGL820() {
  std::vector<GL820Data> data
    {
      // { new GL820DataReader(3),
      //   "./data_gl820/210423-175735_UG.CSV",
      //   "./data_gl820/210423-175735_UG.root",
      //   "./data_gl820/210423-175735_UG.pdf",
      //   "2021-04-23 00:00:00",
      //   "2021-05-24 00:00:00" },
      // { new GL820DataReader(3),
      //   "./data_gl820/210423-191135_UG.CSV",
      //   "./data_gl820/210423-191135_UG.root",
      //   "./data_gl820/210423-191135_UG.pdf",
      //   "2021-04-23 00:00:00",
      //   "2021-05-27 00:00:00" },
      { new GL820DataReader(3),
        "./data_gl820/210427-105914_UG.CSV",
        "./data_gl820/210427-105914_UG.root",
        "./data_gl820/210427-105914_UG.pdf",
        "2021-04-27 00:00:00",
        "2021-05-06 12:00:00" },
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

    std::vector<TGraph*> gRawTemps;
    gRawTemps.resize(datum.reader->size(), nullptr);
    std::vector<TGraph*> gCalibTemps;
    gCalibTemps.resize(datum.reader->size(), nullptr);
    std::vector<TH2*> hRawTemps;
    hRawTemps.resize(datum.reader->size(), nullptr);
    std::vector<TH2*> hCalibTemps;
    hCalibTemps.resize(datum.reader->size(), nullptr);
    TLegend* legend = new TLegend(0.93 - 0.15, 0.93 - 0.05 * datum.reader->size(), 0.93, 0.93);
    for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
      TGraph* gRawTemp = gRawTemps[i] = new TGraph();
      gRawTemp->SetName(Form("gRawTemp%c", 'A' + (char)i));
      gRawTemp->SetTitle("Raw Temperature;;Temperature [^{#circ}C]");
      gRawTemp->SetMinimum(20);
      gRawTemp->SetMaximum(40);
      // gRawTemp->SetMarkerStyle(kFullCircle);
      gRawTemp->SetMarkerStyle(kDot);
      gRawTemp->SetMarkerColor(51 + (99 - 51) * i / n);
      gRawTemp->GetXaxis()->SetTimeDisplay(true);
      gRawTemp->GetXaxis()->SetTimeOffset(0);
      gRawTemp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gRawTemp->GetXaxis()->SetLabelOffset(0.03);
      TGraph* gCalibTemp = gCalibTemps[i] = new TGraph();
      gCalibTemp->SetName(Form("gCalibTemp%c", 'A' + (char)i));
      gCalibTemp->SetTitle("Calibrated Temperature;;Temperature [^{#circ}C]");
      gCalibTemp->SetMinimum(16);
      gCalibTemp->SetMaximum(28);
      // gCalibTemp->SetMarkerStyle(kFullCircle);
      gCalibTemp->SetMarkerStyle(kDot);
      gCalibTemp->SetMarkerColor(51 + (99 - 51) * i / n);
      gCalibTemp->GetXaxis()->SetTimeDisplay(true);
      gCalibTemp->GetXaxis()->SetTimeOffset(0);
      gCalibTemp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      gCalibTemp->GetXaxis()->SetLabelOffset(0.03);

      legend->AddEntry(gRawTemp, Form("PT100-%c", (Char_t)('A' + i)), "P");

      TH2* hRawTemp = hRawTemps[i] =
        new TH2D(Form("hRawTemp%c", 'A' + (char)i),
                 Form("Raw Temperature (PT100-%c);;Temperature [^{#circ}C]", 'A' + (char)i),
                 // (TDatime(datum.datemax.data()).Convert() - TDatime(datum.datemin.data()).Convert()) / (3 * 60 * 60),
                 200,
                 TDatime(datum.datemin.data()).Convert(),
                 TDatime(datum.datemax.data()).Convert(),
                 (40 - 20) / 0.1, 20, 40);
      hRawTemp->GetXaxis()->SetTimeDisplay(true);
      hRawTemp->GetXaxis()->SetTimeOffset(0);
      hRawTemp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      hRawTemp->GetXaxis()->SetLabelOffset(0.03);
      TH2* hCalibTemp = hCalibTemps[i] =
        new TH2D(Form("hCalibTemp%c", 'A' + (char)i),
                 Form("Raw Temperature (PT100-%c);;Temperature [^{#circ}C]", 'A' + (char)i),
                 //(TDatime(datum.datemax.data()).Convert() - TDatime(datum.datemin.data()).Convert()) / (3 * 60 * 60),
                 200,
                 TDatime(datum.datemin.data()).Convert(),
                 TDatime(datum.datemax.data()).Convert(),
                 (28 - 22) / 0.1, 22, 28);
      hCalibTemp->GetXaxis()->SetTimeDisplay(true);
      hCalibTemp->GetXaxis()->SetTimeOffset(0);
      hCalibTemp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
      hCalibTemp->GetXaxis()->SetLabelOffset(0.03);
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
          raw  [i] = temp;
          calib[i] = calibrate(raw[i], i);
          gRawTemps  [i]->SetPoint(gRawTemps  [i]->GetN(), datime.Convert(), raw  [i]);
          gCalibTemps[i]->SetPoint(gCalibTemps[i]->GetN(), datime.Convert(), calib[i]);
          hRawTemps  [i]->Fill(datime.Convert(), raw  [i]);
          hCalibTemps[i]->Fill(datime.Convert(), calib[i]);
        }
      }
      otree->Fill();
    }

    if (gRawTemps.front()->GetN() < 1000) {
      TCanvas::MakeDefCanvas();
      gPad->SetGridy();
      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        gRawTemps[i]->Draw(i ? "P" : "AP");
        gRawTemps[i]->Write();
      }
      legend->Draw();

      gPad->Print((datum.pfilename + "[").data());
      gPad->Print(datum.pfilename.data());

      TCanvas::MakeDefCanvas();
      gPad->SetGridy();
      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        gCalibTemps[i]->Draw(i ? "P" : "AP");
        gCalibTemps[i]->Write();
      }
      legend->Draw();

      gPad->Print(datum.pfilename.data());
      gPad->Print((datum.pfilename + "]").data());

      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        hRawTemps[i]->Write();
      }

      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        hCalibTemps[i]->Write();
      }

    } else {
      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        gRawTemps[i]->Write();
      }

      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        gCalibTemps[i]->Write();
      }

      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        TCanvas::MakeDefCanvas();
        gPad->SetGridy();
        hRawTemps[i]->Draw();
        hRawTemps[i]->Write();

        if (i == 0) {
          gPad->Print((datum.pfilename + "[").data());
        }
        gPad->Print(datum.pfilename.data());
      }

      for (std::size_t i = 0, n = datum.reader->size(); i < n; ++i) {
        TCanvas::MakeDefCanvas();
        gPad->SetGridy();

        hCalibTemps[i]->Draw();
        hCalibTemps[i]->Write();

        gPad->Print(datum.pfilename.data());
      }

      gPad->Print((datum.pfilename + "]").data());
    }

    otree->Write();
    ofile->Close();
  }

  return 0;
}
