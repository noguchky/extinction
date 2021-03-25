#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDatime.h"
#include "TAxis.h"

#include "String.hh"

class JmaDataReader {
protected:
  std::string         fTime;
  std::string         fTemp;
public:
  virtual void        ReadLine(const std::string& line) = 0;
  virtual std::string GetTime() const { return fTime; };
  virtual Double_t    GetTemperature() {
    if (fTemp == "×") {
      return std::nan("");
    }
    return Tron::String::Convert<Double_t>(fTemp);
  }
};

class Mito : public JmaDataReader {
public:
  virtual void ReadLine(const std::string& line) override {
    std::stringstream ss(line);
    std::string buff;
    ss >> /*時分*/            fTime
       >> /*気圧(海面)[hPa]*/ buff
       >> /*気圧(平均)[hPa]*/ buff
       >> /*降水量[mm]*/      buff
       >> /*気温[C]*/         fTemp;
  }
};

class Hitachi : public JmaDataReader {
public:
  virtual void ReadLine(const std::string& line) override {
    std::stringstream ss(line);
    std::string buff;
    ss >> /*時分*/       fTime
       >> /*降水量[mm]*/ buff
       >> /*気温[C]*/    fTemp;
  }
};

struct JmaData {
  JmaDataReader*                 reader;
  const std::string              prefix;
  const std::string              suffix;
  const std::vector<std::string> dates;
  const std::string              ofilename;
};

Int_t plotJmaData() {
  const std::vector<std::string> march
    {
      "2021-03-10",
      "2021-03-11",
      "2021-03-12",
      "2021-03-13",
      "2021-03-14",
      "2021-03-15",
      "2021-03-16",
      "2021-03-17",
      "2021-03-18",
      "2021-03-19",
      "2021-03-20",
      "2021-03-21",
      "2021-03-22",
      "2021-03-23",
      "2021-03-24",
    };

  const std::vector<std::string> may
    {
      "2020-05-09",
      "2020-05-10",
      "2020-05-11",
      "2020-05-12",
      "2020-05-13",
      "2020-05-14",
      "2020-05-15",
      "2020-05-16",
      "2020-05-17",
      "2020-05-18",
      "2020-05-19",
      "2020-05-20",
      "2020-05-21",
      "2020-05-22",
      "2020-05-23",
      "2020-05-24",
      "2020-05-25",
      "2020-05-26",
      "2020-05-27",
      "2020-05-28",
      "2020-05-29",
      "2020-05-30",
      "2020-05-31",
    };

  std::vector<JmaData> data
    {
      { new Mito(),
        "./dataj/jma_mito_",
        ".txt",
        march,
        "./dataj/jma_mito_2021-03.root" },

      { new Hitachi(),
        "./dataj/jma_hitachi_",
        ".txt",
        march,
        "./dataj/jma_hitachi_2021-03.root" },

      { new Mito(),
        "./dataj/jma_mito_",
        ".txt",
        may,
        "./dataj/jma_mito_2020-05.root" },

      { new Hitachi(),
        "./dataj/jma_hitachi_",
        ".txt",
        may,
        "./dataj/jma_hitachi_2020-05.root" },
    };

  for (auto&& datum : data) {
    TFile* ofile = new TFile(datum.ofilename.data(), "RECREATE");
    if (!ofile->IsOpen()) {
      std::cout << "[error] output file is not opened" << std::endl;
      return 0;
    }

    TGraph* gTemp = new TGraph();
    gTemp->SetName("gTemp");
    gTemp->SetMinimum( 0);
    gTemp->SetMaximum(30);
    gTemp->GetXaxis()->SetTimeDisplay(true);
    gTemp->GetXaxis()->SetTimeOffset(0);

    std::string buff;
    TDatime datime;
    for (auto&& date : datum.dates) {
      std::ifstream file(datum.prefix + date + datum.suffix);
      if (!file) {
        std::cout << "[warning] file is not opened, " << datum.prefix + date + datum.suffix << std::endl;
        continue;
      }

      // skip header
      for (Int_t i = 0; i < 2; ++i) {
        std::getline(file, buff);
      }

      while (std::getline(file, buff)) {
        datum.reader->ReadLine(buff);

        if (std::isfinite(datum.reader->GetTemperature())) {
          datime.Set((date + " " + datum.reader->GetTime() + ":00").data());
          gTemp->SetPoint(gTemp->GetN(), datime.Convert(), datum.reader->GetTemperature());
        }
      }

    }

    TCanvas::MakeDefCanvas();
    gTemp->Draw("AP");

    gTemp->Write();
    ofile->Close();
  }

  return 0;
}
