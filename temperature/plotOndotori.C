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

class OndotoriDataReader {
protected:
  std::string         fTime;
  std::string         fTemp;
public:
  virtual void        ReadLine(const std::string& line) {
    std::vector<std::string> strs = Tron::String::Split(line, ",");
    if (strs.size() != 4) {
      fTime = "";
      fTemp = "-200";
      return;
    }
    
    for (auto&& str : strs) {
      str = Tron::String::Replace(str, "\"", "");
    }
    fTime = strs[0]; // Date/Time (string)
    //      strs[1]  // Date/Time (serial)
    fTemp = strs[2]; // Temperature
    //      strs[3]  // Humidity
  }
  virtual std::string GetTime() const { return fTime; };
  virtual Double_t    GetTemperature() {
    return Tron::String::Convert<Double_t>(fTemp);
  }
};

struct OndotoriData {
  OndotoriDataReader*            reader;
  const std::string              ifilename;
  const std::string              ofilename;
};

Int_t plotOndotori() {
  std::vector<OndotoriData> data
    {
      { new OndotoriDataReader(),
        "./data_ondotori/521455AD.1619164233.csv",
        "./data_ondotori/521455AD.1619164233.root" },
      { new OndotoriDataReader(),
        "./data_ondotori/521455AD.1619407372.csv",
        "./data_ondotori/521455AD.1619407372.root" },
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
    gTemp->SetMarkerStyle(kFullCircle);
    gTemp->GetXaxis()->SetTimeDisplay(true);
    gTemp->GetXaxis()->SetTimeOffset(0);
    gTemp->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
    gTemp->GetXaxis()->SetLabelOffset(0.03);

    std::ifstream file(datum.ifilename);
    if (!file) {
      std::cout << "[warning] file is not opened, " << datum.ifilename << std::endl;
      continue;
    }

    // skip header
    std::string buff;
    for (Int_t i = 0; i < 3; ++i) {
      std::getline(file, buff);
    }

    TDatime datime;
    while (std::getline(file, buff)) {
      datum.reader->ReadLine(buff);

      if (std::isfinite(datum.reader->GetTemperature())) {
        datime.Set(datum.reader->GetTime().data());
        gTemp->SetPoint(gTemp->GetN(), datime.Convert(), datum.reader->GetTemperature());
      }
    }

    TCanvas::MakeDefCanvas();
    gTemp->Draw("AP");

    gTemp->Write();
    ofile->Close();
  }

  return 0;
}
