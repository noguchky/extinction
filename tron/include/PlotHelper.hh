#ifndef Tron_PlotHelper_hh
#define Tron_PlotHelper_hh

#include <iostream>
#include <fstream>
#include <functional>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"

namespace Tron {

  namespace PlotHelper {

    using Func1_t     = std::function<Double_t(Double_t)>;
    using Func2_t     = std::function<Double_t(Double_t, Double_t)>;
    using Func3_t     = std::function<Double_t(Double_t, Double_t, Double_t)>;

    inline Bool_t ReadFile(std::vector<Double_t>& array,
                           const std::string& filename,
                           const std::string& format = "",
                           const Func1_t&     fx     = nullptr) {
      std::ifstream file(filename.data());
      if (!file) {
        std::cout << "ReadFile [error] file is not opened, "
                  << "\"" << filename << "\"" << std::endl;
        return false;
      }
      std::string buf;
      Double_t x;
      if (format.empty()) {
        while (file >> x) {
          array.push_back(fx ? fx(x) : x);
        }
      } else {
        while (std::getline(file, buf)) {
          if (std::sscanf(buf.data(), format.data(), &x)) {
            array.push_back(fx ? fx(x) : x);
          }
        }
      }
      if (array.empty()) {
        std::cout << "ReadFile [warning] array is empty, "
                  << "\"" << filename << "\"" << ", "
                  << "\"" << format   << "\"" << std::endl;
        return false;
      }
      return true;
    }

    inline Bool_t ReadFile(TH1*               hist,
                           const std::string& filename,
                           const std::string& format = "",
                           const Func1_t&     fx     = nullptr) {
      std::ifstream file(filename.data());
      if (!file) {
        std::cout << "ReadFile [error] file is not opened, "
                  << "\"" << filename << "\"" << std::endl;
        return false;
      }
      std::string buf;
      Double_t x;
      if (format.empty()) {
        while (file >> x) {
          hist->Fill(fx ? fx(x) : x);
        }
      } else {
        while (std::getline(file, buf)) {
          if (std::sscanf(buf.data(), format.data(), &x)) {
            hist->Fill(fx ? fx(x) : x);
          }
        }
      }
      if (!hist->GetEntries()) {
        std::cout << "ReadFile [warning] hist is empty, "
                  << "\"" << filename << "\"" << ", "
                  << "\"" << format   << "\"" << std::endl;
        return false;
      }
      return true;
    }

    inline Bool_t ReadFile(TH2*               hist,
                           const std::string& filename,
                           const std::string& format = "",
                           const Func2_t&     fx     = nullptr,
                           const Func2_t&     fy     = nullptr) {
      std::ifstream file(filename.data());
      if (!file) {
        std::cout << "ReadFile [error] file is not opened, "
                  << "\"" << filename << "\"" << std::endl;
        return false;
      }
      std::string buf;
      Double_t x, y;
      if (format.empty()) {
        while (file >> x >> y) {
          hist->Fill(fx ? fx(x, y) : x,
                     fy ? fy(x, y) : y);
        }
      } else {
        while (std::getline(file, buf)) {
          if (std::sscanf(buf.data(), format.data(), &x, &y)) {
            hist->Fill(fx ? fx(x, y) : x,
                       fy ? fy(x, y) : y);
          }
        }
      }
      if (!hist->GetEntries()) {
        std::cout << "ReadFile [warning] hist is empty, "
                  << "\"" << filename << "\"" << ", "
                  << "\"" << format   << "\"" << std::endl;
        return false;
      }
      return true;
    }

    inline Bool_t ReadFile(TGraph*            graph,
                           const std::string& filename,
                           const std::string& format = "",
                           const Func2_t&     fx     = nullptr,
                           const Func2_t&     fy     = nullptr) {
      std::ifstream file(filename.data());
      if (!file) {
        std::cout << "ReadFile [error] file is not opened, "
                  << "\"" << filename << "\"" << std::endl;
        return false;
      }
      std::string buf;
      Double_t x, y;
      if (format.empty()) {
        while (file >> x >> y) {
          graph->SetPoint(graph->GetN(),
                          fx ? fx(x, y) : x,
                          fy ? fy(x, y) : y);
        }
      } else {
        while (std::getline(file, buf)) {
          if (std::sscanf(buf.data(), format.data(), &x, &y)) {
            graph->SetPoint(graph->GetN(),
                            fx ? fx(x, y) : x,
                            fy ? fy(x, y) : y);
          }
        }
      }
      if (!graph->GetN()) {
        std::cout << "ReadFile [warning] graph is empty, "
                  << "\"" << filename << "\"" << ", "
                  << "\"" << format   << "\"" << std::endl;
        return false;
      }
      return true;
    }

    inline Bool_t ReadFile(TGraph2D*          graph,
                           const std::string& filename,
                           const std::string& format = "",
                           const Func3_t&     fx = nullptr,
                           const Func3_t&     fy = nullptr,
                           const Func3_t&     fz = nullptr) {
      std::ifstream file(filename.data());
      if (!file) {
        std::cout << "ReadFile [error] file is not opened, "
                  << "\"" << filename << "\"" << std::endl;
        return false;
      }
      std::string buf;
      Double_t x, y, z;
      if (format.empty()) {
        while (file >> x >> y >> z) {
          graph->SetPoint(graph->GetN(),
                          fx ? fx(x, y, z) : x,
                          fy ? fy(x, y, z) : y,
                          fz ? fz(x, y, z) : z);
        }
      } else {
        while (std::getline(file, buf)) {
          if (std::sscanf(buf.data(), format.data(), &x, &y, &z)) {
            graph->SetPoint(graph->GetN(),
                            fx ? fx(x, y, z) : x,
                            fy ? fy(x, y, z) : y,
                            fz ? fz(x, y, z) : z);
          }
        }
      }
      if (!graph->GetN()) {
        std::cout << "ReadFile [warning] graph is empty, "
                  << "\"" << filename << "\"" << ", "
                  << "\"" << format   << "\"" << std::endl;
        return false;
      }
      return true;
    }

    inline TH2* Projection(const TGraph* graph, const std::string& name, Int_t nbinsx = 50, Int_t nbinsy = 50) {
      const std::string title = Form("%s;%s;%s;%s", graph->GetTitle(), graph->GetXaxis()->GetTitle(), graph->GetYaxis()->GetTitle(), "Entries [A.U.]");
      const Double_t    xmin = graph->GetXaxis()->GetXmin();
      const Double_t    xmax = graph->GetXaxis()->GetXmax();
      const Double_t    ymin = graph->GetYaxis()->GetXmin();
      const Double_t    ymax = graph->GetYaxis()->GetXmax();
      TH2* hist = new TH2D(name.data(), title.data(), nbinsx, xmin, xmax, nbinsy, ymin, ymax);
      hist->FillN(graph->GetN(), graph->GetX(), graph->GetY(), nullptr);
      return hist;
    }

    inline TH1* ProjectionX(const TGraph* graph, const std::string& name, Int_t nbinsx = 50) {
      const std::string title = Form("%s;%s;%s", graph->GetTitle(), graph->GetXaxis()->GetTitle(), "Entries [A.U.]");
      const Double_t    xmin = graph->GetXaxis()->GetXmin();
      const Double_t    xmax = graph->GetXaxis()->GetXmax();
      TH1* hist = new TH1D(name.data(), title.data(), nbinsx, xmin, xmax);
      hist->FillN(graph->GetN(), graph->GetX(), nullptr);
      return hist;
    }

    inline TH1* ProjectionY(const TGraph* graph, const std::string& name, Int_t nbinsx = 50) {
      const std::string title = Form("%s;%s;%s", graph->GetTitle(), graph->GetYaxis()->GetTitle(), "Entries [A.U.]");
      const Double_t    xmin = graph->GetYaxis()->GetXmin();
      const Double_t    xmax = graph->GetYaxis()->GetXmax();
      TH1* hist = new TH1D(name.data(), title.data(), nbinsx, xmin, xmax);
      hist->FillN(graph->GetN(), graph->GetY(), nullptr);
      return hist;
    }

  }

}

#endif
