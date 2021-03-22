#ifndef Extinction_SpillData_hh
#define Extinction_SpillData_hh

#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TDatime.h"

#include "Units.hh"
#include "Detector.hh"

#include "Math.hh"

namespace Extinction {

  struct SpillData {
    static const std::size_t kNofBunches = 4;

    TDatime  Date;
    Int_t    Year;
    Int_t    Month;
    Int_t    Day;
    Int_t    Hour;
    Int_t    Minute;
    Int_t    Second;
    Int_t    EMCount;
    Long64_t Entries[GlobalChannel::NofChannels];
    Long64_t CoinCount;
    Long64_t Leakage;
    Long64_t InBunch;
    Double_t Extinction;
    Double_t ExtinctionErr;
    Double_t BunchCenters[kNofBunches];
    Double_t BunchWidths [kNofBunches];
    Double_t MrSyncInterval[MrSync::NofChannels];
    Double_t TimePerTdc    [MrSync::NofChannels];

    SpillData() {
      Clear();
    }

    void Clear() {
      SetDate(0);
      EMCount       = -1;
      CoinCount     = 0;
      Leakage       = 0;
      InBunch       = 0;
      Extinction    = 1.0;
      ExtinctionErr = 1.0;
      std::memset(Entries       , 0, sizeof(Entries       ));
      std::memset(BunchCenters  , 0, sizeof(BunchCenters  ));
      std::memset(BunchWidths   , 0, sizeof(BunchWidths   ));
      std::memset(MrSyncInterval, 0, sizeof(MrSyncInterval));
      std::memset(TimePerTdc    , 0, sizeof(TimePerTdc    ));
    }

    void SetDate(UInt_t date) {
      Date.Set(date);
      Year   = Date.GetYear();  
      Month  = Date.GetMonth(); 
      Day    = Date.GetDay();   
      Hour   = Date.GetHour();  
      Minute = Date.GetMinute();
      Second = Date.GetSecond();
    }

    void SetParticles(Long64_t leakage, Long64_t inBunch) {
      Leakage = leakage;
      InBunch = inBunch;
      if (InBunch) {
        Extinction = (Double_t)Leakage / (Double_t)InBunch;
        if (Leakage) {
          ExtinctionErr = Tron::Math::SqrtOfSumOfSquared(TMath::Sqrt(Leakage) / InBunch, Leakage / TMath::Power((Double_t)InBunch, 1.5));
        } else {
          ExtinctionErr = 1.0 / InBunch;
        }
      } else {
        Extinction = 1.0;
      }
    }

    void CreateBranch(TTree* tree) {
      // std::cout << "SpillData::CreateBranch()" << std::endl;
      tree->Branch("date"          , &Date                                                                        );
      tree->Branch("year"          , &Year          , Form("year"                "/I"                            ));
      tree->Branch("month"         , &Month         , Form("month"               "/I"                            ));
      tree->Branch("day"           , &Day           , Form("day"                 "/I"                            ));
      tree->Branch("hour"          , &Hour          , Form("hour"                "/I"                            ));
      tree->Branch("minute"        , &Minute        , Form("minute"              "/I"                            ));
      tree->Branch("second"        , &Second        , Form("second"              "/I"                            ));
      tree->Branch("emcount"       , &EMCount       , Form("emcount"             "/I"                            ));      
      tree->Branch("entries"       ,  Entries       , Form("entries"        "[%lu]/L", GlobalChannel::NofChannels));
      tree->Branch("coincount"     , &CoinCount     , Form("coincount"           "/L"                            ));
      tree->Branch("leakage"       , &Leakage       , Form("leakage"             "/L"                            ));
      tree->Branch("inbunch"       , &InBunch       , Form("inbunch"             "/L"                            ));
      tree->Branch("extinction"    , &Extinction    , Form("extinction"          "/D"                            ));
      tree->Branch("extinctionerr" , &ExtinctionErr , Form("extinctionerr"       "/D"                            ));
      tree->Branch("bunchcenter"   ,  BunchCenters  , Form("bunchcenter"    "[%lu]/D",                kNofBunches));
      tree->Branch("bunchwidth"    ,  BunchWidths   , Form("bunchwidth"     "[%lu]/D",                kNofBunches));
      tree->Branch("mrsyncinterval",  MrSyncInterval, Form("mrsyncinterval" "[%lu]/D", MrSync       ::NofChannels));
      tree->Branch("timepertdc"    ,  TimePerTdc    , Form("timepertdc"     "[%lu]/D", MrSync       ::NofChannels));
    }

    void SetBranchAddress(TTree* tree) {
      // std::cout << "SpillData::SetBranchAddress()" << std::endl;
      tree->SetBranchAddress("date"          , &Date          );
      tree->SetBranchAddress("year"          , &Year          );
      tree->SetBranchAddress("month"         , &Month         );
      tree->SetBranchAddress("day"           , &Day           );
      tree->SetBranchAddress("hour"          , &Hour          );
      tree->SetBranchAddress("minute"        , &Minute        );
      tree->SetBranchAddress("second"        , &Second        );
      tree->SetBranchAddress("emcount"       , &EMCount       );
      tree->SetBranchAddress("entries"       ,  Entries       );
      tree->SetBranchAddress("coincount"     , &CoinCount     );
      tree->SetBranchAddress("leakage"       , &Leakage       );
      tree->SetBranchAddress("inbunch"       , &InBunch       );
      tree->SetBranchAddress("extinction"    , &Extinction    );
      tree->SetBranchAddress("extinctionerr" , &ExtinctionErr );
      tree->SetBranchAddress("bunchcenter"   ,  BunchCenters  );
      tree->SetBranchAddress("bunchwidth"    ,  BunchWidths   );
      tree->SetBranchAddress("mrsyncinterval",  MrSyncInterval);
      tree->SetBranchAddress("timepertdc"    ,  TimePerTdc    );
    }

  };

  struct SpillPlotBase {
    virtual void FillData(SpillData*) = 0;
    virtual void Draw() = 0;
    Bool_t GridX = true, GridY = true, LogX = false, LogY = false, LogZ = false;
  };

  struct SpillHist1 : public SpillPlotBase {
    TH1* Hist;
    SpillHist1(const std::string& name, const std::string& title, Int_t xbins, Double_t xmin, Double_t xmax) {
      Hist = new TH1D(name.data(), title.data(), xbins, xmin, xmax);
    }

    std::function<Double_t(SpillData*)> Value;
    std::string GOption;

    virtual void FillData(SpillData* data) override {
      Double_t value;
      if (Value && std::isfinite(value = Value(data))) {
        Hist->Fill(value);
      }
    }

    virtual void Draw() override {
      Hist->Draw(GOption.data());
    }
  };

  struct SpillHist2 : public SpillPlotBase {
    TH2* Hist;
    SpillHist2(const std::string& name, const std::string& title, Int_t xbins, Double_t xmin, Double_t xmax, Int_t ybins, Double_t ymin, Double_t ymax) {
      Hist = new TH2D(name.data(), title.data(), xbins, xmin, xmax, ybins, ymin, ymax);
    }

    std::function<Double_t(SpillData*)> XValue;
    std::function<Double_t(SpillData*)> YValue;
    std::string GOption;

    virtual void FillData(SpillData* data) override {
      Double_t x, y;
      if (XValue && std::isfinite(x = XValue(data)) &&
          YValue && std::isfinite(y = YValue(data))) {
        Hist->Fill(x, y);
      }
    }

    virtual void Draw() override {
      Hist->Draw(GOption.data());
    }
  };

  struct SpillGraph : public SpillPlotBase {
    TGraph* Graph;
    SpillGraph(const std::string& name, const std::string& title) {
      Graph = new TGraph();
      Graph->SetNameTitle(name.data(), title.data());
    }

    std::function<Double_t(SpillData*)> XValue;
    std::function<Double_t(SpillData*)> YValue;
    std::string GOption;

    virtual void FillData(SpillData* data) override {
      Double_t x, y;
      if (XValue && std::isfinite(x = XValue(data)) &&
          YValue && std::isfinite(y = YValue(data))) {
        Graph->SetPoint(Graph->GetN(), x, y);
      }
    }

    virtual void Draw() override {
      Graph->Draw(GOption.data());
    }
  };

  struct SpillGraphErrors : public SpillPlotBase {
    TGraphErrors* Graph;
    SpillGraphErrors(const std::string& name, const std::string& title) {
      Graph = new TGraphErrors();
      Graph->SetNameTitle(name.data(), title.data());
    }

    std::function<Double_t(SpillData*)> XValue;
    std::function<Double_t(SpillData*)> YValue;
    std::function<Double_t(SpillData*)> ExValue;
    std::function<Double_t(SpillData*)> EyValue;
    std::string GOption;

    virtual void FillData(SpillData* data) override {
      Double_t x, y;
      if (XValue && std::isfinite(x = XValue(data)) &&
          YValue && std::isfinite(y = YValue(data))) {
        const Int_t np = Graph->GetN();
        Graph->SetPoint(np, x, y);
        const Double_t ex = ExValue ? ExValue(data) : 0.0;
        const Double_t ey = EyValue ? EyValue(data) : 0.0;
        Graph->SetPointError(np, std::isfinite(ex) ? ex : 0.0, std::isfinite(ey) ? ey : 0.0);
      }
    }

    virtual void Draw() override {
      Graph->Draw(GOption.data());
    }
  };

  class SpillPlotter {
  public:

  private:
    TTree*                      fTree;
    SpillData                   fData;
    std::vector<SpillPlotBase*> fPlots;

  public:
    SpillPlotter(TTree* tree)
      : fTree(tree) {
      fData.SetBranchAddress(fTree);
    }
    ~SpillPlotter() {
      fTree = nullptr;
    }

    void AddPlot(SpillPlotBase* plot) {
      fPlots.push_back(plot);
    }

    void FillData() {
      const Long64_t entries = fTree->GetEntries();
      for (Long64_t entry = 0; entry < entries; ++entry) {
        if (fTree->GetEntry(entry)) {
          for (auto&& plot : fPlots) {
            plot->FillData(&fData);
          }
        }
      }
    }

    void DrawPlots(const std::string& ofilename) {
      TCanvas::MakeDefCanvas();
      gPad->Print((ofilename + "[").data());

      for (auto&& plot : fPlots) {
        gPad->SetGrid(plot->GridX, plot->GridY);
        gPad->SetLogx(plot->LogX), gPad->SetLogy(plot->LogY), gPad->SetLogz(plot->LogZ);
        plot->Draw();
        gPad->Print(ofilename.data());
      }

      gPad->Print((ofilename + "]").data());
    }
  };

}

#endif
