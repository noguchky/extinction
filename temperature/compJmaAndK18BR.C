#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TDatime.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TMultiLayerPerceptron.h"
#include "TPaveStats.h"

#include "String.hh"
#include "Math.hh"
#include "LinearLeastSquareMethod.hh"
#include "Integrator.hh"
#include "ProgressBar.hh"

namespace {

  const Double_t stdTemperature = 20.0;
  
}

TGraph* GetGraphFromFile(const std::string& id, const std::string& filename, const std::string& graphname) {
  TFile* file = new TFile(filename.data());
  if (!file->IsOpen()) {
    std::cout << "[error] " << id << " file is not opened" << std::endl;
    exit(1);
  }

  TGraph* graph = dynamic_cast<TGraph*>(file->Get(graphname.data()));
  if (!graph) {
    std::cout << "[error] " << id << "graph is not found" << std::endl;
    exit(1);
  }

  return graph;
}

void RemoveOutlierPoints(TGraph* graph, Double_t ymin, Double_t ymax) {
  Double_t x, y;
  for (Int_t i = 0, n = graph->GetN(); i < n; ++i) {
    graph->GetPoint(i, x, y);
    if (y < ymin || ymax < y) {
      graph->RemovePoint(i);
      --i, --n;
    }
  }
}

Int_t compJmaAndK18BR() {
  // Get temperature graphs
  TGraph* gTemp_K18BR_March1  = GetGraphFromFile("K1.8BR(March-1)", "./data_easiroc/20210310_k18br.root" , "temp1");
  TGraph* gTemp_K18BR_March2  = GetGraphFromFile("K1.8BR(March-2)", "./data_easiroc/20210324_k18br.root" , "temp1");
  TGraph* gTemp_K18BR_March3  = GetGraphFromFile("K1.8BR(March-3)", "./data_easiroc/20210330_k18br.root" , "temp1");
  TGraph* gTemp_Hitachi_March = GetGraphFromFile("Hitachi(March)" , "./data_jma/jma_hitachi_2021-03.root", "gTemp");
  TGraph* gTemp_Mito_March    = GetGraphFromFile("Mito(March)"    , "./data_jma/jma_mito_" "2021-03.root", "gTemp");
  TGraph* gTemp_Hitachi_May   = GetGraphFromFile("Hitachi(May)"   , "./data_jma/jma_hitachi_2020-05.root", "gTemp");
  TGraph* gTemp_Mito_May      = GetGraphFromFile("Mito(May)"      , "./data_jma/jma_mito_" "2020-05.root", "gTemp");

  // Remove outlier points
  RemoveOutlierPoints(gTemp_K18BR_March1, 0, 40);
  RemoveOutlierPoints(gTemp_K18BR_March2, 0, 40);
  RemoveOutlierPoints(gTemp_K18BR_March3, 0, 40);

  // Set style
  for (auto&& pair : std::map<TGraph*, std::tuple<Int_t, Int_t>> {
     { gTemp_K18BR_March1  , { kGreen  + 1, kFullSquare       } },
     { gTemp_K18BR_March2  , { kSpring + 1, kFullSquare       } },
     { gTemp_K18BR_March3  , { kTeal   + 1, kFullSquare       } },
     { gTemp_Hitachi_March , { kBlue   + 1, kFullSquare /* kFullTriangleUp  */ } },
     { gTemp_Mito_March    , { kRed    + 1, kFullSquare /* kFullTriangleDown*/ } },
     { gTemp_Hitachi_May   , { kBlue   + 1, kFullSquare /* kFullTriangleUp  */ } },
     { gTemp_Mito_May      , { kRed    + 1, kFullSquare /* kFullTriangleDown*/ } },
    }) {
    TGraph* graph = pair.first;
    const Int_t color  = std::get<0>(pair.second);
    const Int_t marker = std::get<1>(pair.second);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.4);
    graph->SetLineColor(color);
    graph->SetLineWidth(4);
    graph->SetMinimum(- 5.0);
    graph->SetMaximum(+30.0);
    graph->SetTitle(";;Temperature [#circC]");
  }

  // Create difference
  TGraph* gDtemp_Hitachi_March = (TGraph*)gTemp_Hitachi_March->Clone();
  TGraph* gDtemp_Mito_March    = (TGraph*)gTemp_Mito_March   ->Clone();
  gDtemp_Hitachi_March->Set(0);
  gDtemp_Mito_March   ->Set(0);
  gDtemp_Hitachi_March->SetMinimum();
  gDtemp_Mito_March   ->SetMinimum();
  gDtemp_Hitachi_March->SetMaximum();
  gDtemp_Mito_March   ->SetMaximum();
  gDtemp_Hitachi_March->SetTitle(Form("Hitachi (March);;T_{Hitachi}-%g [#circC]", stdTemperature));
  gDtemp_Mito_March   ->SetTitle(Form("Mito " "(March);;T_{Mito" "}-%g [#circC]", stdTemperature));

  TGraph* gCtemp_Hitachi_March1 = new TGraph();
  TGraph* gCtemp_Hitachi_March2 = new TGraph();
  TGraph* gCtemp_Hitachi_March3 = new TGraph();
  TGraph* gCtemp_Mito_March1    = new TGraph();
  TGraph* gCtemp_Mito_March2    = new TGraph();
  TGraph* gCtemp_Mito_March3    = new TGraph();
  for (auto&& pair : std::map<TGraph*, std::tuple<std::string, Int_t, Int_t>> {
      { gCtemp_Hitachi_March1 , { "Hitachi" , kBlue   + 1, kFullSquare } },
      { gCtemp_Hitachi_March2 , { "Hitachi" , kAzure  + 1, kFullSquare } },
      { gCtemp_Hitachi_March3 , { "Hitachi" , kCyan   + 1, kFullSquare } },
      { gCtemp_Mito_March1    , { "Mito"    , kRed    + 1, kFullSquare } },
      { gCtemp_Mito_March2    , { "Mito"    , kOrange + 1, kFullSquare } },
      { gCtemp_Mito_March3    , { "Mito"    , kPink   + 1, kFullSquare } },
    }) {
    TGraph* graph = pair.first;
    const std::string place  = std::get<0>(pair.second);
    const Int_t       color  = std::get<1>(pair.second);
    const Int_t       marker = std::get<2>(pair.second);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.4);
    graph->SetLineColor(color);
    graph->SetLineWidth(4);
    graph->SetTitle(";;Temperature [#circC]");
    graph->SetTitle(Form("%s (March);T_{%s}-%g [#circC];T_{K1.8BR}-%g [#circC]", place.data(), place.data(),stdTemperature, stdTemperature));
  }

  const Double_t hour   = TDatime(2000, 1, 1, 1, 0, 0).Convert() - TDatime(2000, 1, 1, 0, 0, 0).Convert();
  const Double_t weight = 0.7;

  const Double_t tmin_March0 = Tron::Math::Greatest(TMath::MinElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MinElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()));
  const Double_t tmax_March0 = Tron::Math::Least   (TMath::MaxElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MaxElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()));
  const Double_t tmin_March1 = Tron::Math::Greatest(TMath::MinElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MinElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                    TMath::MinElement(gTemp_K18BR_March1 ->GetN(), gTemp_K18BR_March1 ->GetX()));
  const Double_t tmax_March1 = Tron::Math::Least   (TMath::MaxElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MaxElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                    TMath::MaxElement(gTemp_K18BR_March1 ->GetN(), gTemp_K18BR_March1 ->GetX()));
  const Double_t tmin_March2 = Tron::Math::Greatest(TMath::MinElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MinElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                    TMath::MinElement(gTemp_K18BR_March2 ->GetN(), gTemp_K18BR_March2 ->GetX()));
  const Double_t tmax_March2 = Tron::Math::Least   (TMath::MaxElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MaxElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                    TMath::MaxElement(gTemp_K18BR_March2 ->GetN(), gTemp_K18BR_March2 ->GetX()));
  const Double_t tmin_March3 = Tron::Math::Greatest(TMath::MinElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MinElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                    TMath::MinElement(gTemp_K18BR_March3 ->GetN(), gTemp_K18BR_March3 ->GetX()));
  const Double_t tmax_March3 = Tron::Math::Least   (TMath::MaxElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                    TMath::MaxElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                    TMath::MaxElement(gTemp_K18BR_March3 ->GetN(), gTemp_K18BR_March3 ->GetX()));

  for (Int_t i = 60, n = gDtemp_Hitachi_March->GetN(); i < n; ++i) {
    const Double_t t = gDtemp_Hitachi_March->GetX()[i];
    gDtemp_Hitachi_March->SetPoint(gDtemp_Hitachi_March->GetN(), t, gTemp_Hitachi_March->Eval(t) - stdTemperature);
  }
  for (Int_t i = 60, n = gDtemp_Mito_March->GetN(); i < n; ++i) {
    const Double_t t = gDtemp_Mito_March->GetX()[i];
    gDtemp_Mito_March   ->SetPoint(gDtemp_Mito_March   ->GetN(), t, gTemp_Mito_March   ->Eval(t) - stdTemperature);
  }

  for (Int_t i = 60, n = gTemp_K18BR_March1->GetN(); i < n; ++i) {
    const Double_t t = gTemp_K18BR_March1->GetX()[i];
    if (tmin_March1 <= t && t <= tmax_March1) {
      gCtemp_Hitachi_March1->SetPoint(gCtemp_Hitachi_March1->GetN(), (1.0 - weight) * (gTemp_Hitachi_March->Eval(t - 0.0 * hour) - stdTemperature) + weight * (gTemp_Hitachi_March->Eval(t - 2.0 * hour) - stdTemperature), gTemp_K18BR_March1->Eval(t) - stdTemperature);
      gCtemp_Mito_March1   ->SetPoint(gCtemp_Mito_March1   ->GetN(), (1.0 - weight) * (gTemp_Mito_March   ->Eval(t - 0.0 * hour) - stdTemperature) + weight * (gTemp_Mito_March   ->Eval(t - 2.0 * hour) - stdTemperature), gTemp_K18BR_March1->Eval(t) - stdTemperature);
    }
  }

  for (Int_t i = 60, n = gTemp_K18BR_March2->GetN(); i < n; ++i) {
    const Double_t t = gTemp_K18BR_March2->GetX()[i];
    if (tmin_March2 <= t && t <= tmax_March2) {
      gCtemp_Hitachi_March2->SetPoint(gCtemp_Hitachi_March2->GetN(), (1.0 - weight) * (gTemp_Hitachi_March->Eval(t - 0.0 * hour) - stdTemperature) + weight * (gTemp_Hitachi_March->Eval(t - 2.0 * hour) - stdTemperature), gTemp_K18BR_March2->Eval(t) - stdTemperature);
      gCtemp_Mito_March2   ->SetPoint(gCtemp_Mito_March2   ->GetN(), (1.0 - weight) * (gTemp_Mito_March   ->Eval(t - 0.0 * hour) - stdTemperature) + weight * (gTemp_Mito_March   ->Eval(t - 2.0 * hour) - stdTemperature), gTemp_K18BR_March2->Eval(t) - stdTemperature);
    }
  }

  for (Int_t i = 60, n = gTemp_K18BR_March3->GetN(); i < n; ++i) {
    const Double_t t = gTemp_K18BR_March3->GetX()[i];
    if (tmin_March3 <= t && t <= tmax_March3) {
      gCtemp_Hitachi_March3->SetPoint(gCtemp_Hitachi_March3->GetN(), (1.0 - weight) * (gTemp_Hitachi_March->Eval(t - 0.0 * hour) - stdTemperature) + weight * (gTemp_Hitachi_March->Eval(t - 2.0 * hour) - stdTemperature), gTemp_K18BR_March3->Eval(t) - stdTemperature);
      gCtemp_Mito_March3   ->SetPoint(gCtemp_Mito_March3   ->GetN(), (1.0 - weight) * (gTemp_Mito_March   ->Eval(t - 0.0 * hour) - stdTemperature) + weight * (gTemp_Mito_March   ->Eval(t - 2.0 * hour) - stdTemperature), gTemp_K18BR_March3->Eval(t) - stdTemperature);
    }
  }
  
  TF1* fline = new TF1("fline", "[0]*x + [1]");
  fline->SetParNames("grad", "intercept");
  fline->SetParameters(1.0, 0.0);

  gCtemp_Hitachi_March1->Fit(fline, "", "");
  gCtemp_Mito_March1   ->Fit(fline, "", "");

  const Double_t grad      = gCtemp_Mito_March1->GetFunction(fline->GetName())->GetParameter(0);
  const Double_t intercept = gCtemp_Mito_March1->GetFunction(fline->GetName())->GetParameter(1);
  auto estimate = [&](Double_t y, Double_t y0) {
    return (((1.0 - weight) * (y - stdTemperature) + weight * (y0 - stdTemperature)) * grad + intercept) + stdTemperature;
  };

  // Create stimated
  TGraph* gEstimatedTemp_K18BR_March1 = (TGraph*)gTemp_K18BR_March1->Clone();
  TGraph* gEstimatedTemp_K18BR_May    = (TGraph*)gTemp_K18BR_March1->Clone();
  gEstimatedTemp_K18BR_March1->Set(0);
  gEstimatedTemp_K18BR_May  ->Set(0);
  gEstimatedTemp_K18BR_March1->SetLineColor(kYellow + 1);
  gEstimatedTemp_K18BR_May  ->SetLineColor(kYellow + 1);

  for (Int_t i = 12, n = gTemp_Mito_March->GetN(); i < n; ++i) {
    gEstimatedTemp_K18BR_March1->SetPoint(gEstimatedTemp_K18BR_March1->GetN(),
                                          gTemp_Mito_March->GetX()[i],
                                          estimate(gTemp_Mito_March->GetY()[i - 0], gTemp_Mito_March->GetY()[i - 12]));
  }
  for (Int_t i = 12, n = gTemp_Mito_May->GetN(); i < n; ++i) {
    gEstimatedTemp_K18BR_May->SetPoint(gEstimatedTemp_K18BR_May->GetN(),
                                       gTemp_Mito_May->GetX()[i],
                                       estimate(gTemp_Mito_May->GetY()[i - 0], gTemp_Mito_May->GetY()[i - 12]));
  }

  const std::string ofilename = "./data/compJmaAndEasiroc.pdf";
  TCanvas::MakeDefCanvas();
  gPad->Print((ofilename + "[").data());
  
  // Draw graphs and hists
  {
    gPad->SetGrid();
    gTemp_Hitachi_March        ->Draw("AL");
    gTemp_Mito_March           ->Draw( "L");
    gTemp_K18BR_March1         ->Draw( "L");
    gTemp_K18BR_March2         ->Draw( "L");
    gTemp_K18BR_March3         ->Draw( "L");
    gEstimatedTemp_K18BR_March1->Draw( "L");
    // gTemp_Hitachi_March->GetXaxis()->SetLimits(TDatime(2021, 3, 10, 0, 0, 0).Convert(),
    //                                            TDatime(2021, 3, 20, 0, 0, 0).Convert());
    gTemp_Hitachi_March->GetXaxis()->SetLimits(TDatime(2021, 3, 10, 0, 0, 0).Convert(),
                                               TDatime(2021, 4, 20, 0, 0, 0).Convert());

    TLegend* legend = new TLegend(0.83, 0.97 - 0.04 * 6, 0.97, 0.97);
    legend->AddEntry(gTemp_K18BR_March1         , "K1.8BR(1)", "L");
    legend->AddEntry(gTemp_K18BR_March2         , "K1.8BR(2)", "L");
    legend->AddEntry(gTemp_K18BR_March3         , "K1.8BR(3)", "L");
    legend->AddEntry(gTemp_Hitachi_March        , "Hitachi"  , "L");
    legend->AddEntry(gTemp_Mito_March           , "Mito"     , "L");
    legend->AddEntry(gEstimatedTemp_K18BR_March1, "Estimated", "L");
    legend->Draw();

    gPad->Print(ofilename.data());
  }

  {
    TCanvas::MakeDefCanvas();
    gPad->SetGrid();
    gTemp_Hitachi_May       ->Draw("AL");
    gTemp_Mito_May          ->Draw( "L");
    gEstimatedTemp_K18BR_May->Draw( "L");

    TLegend* legend = new TLegend(0.83, 0.97 - 0.04 * 3, 0.97, 0.97);
    legend->AddEntry(gTemp_Hitachi_May       , "Hitachi"  , "L");
    legend->AddEntry(gTemp_Mito_May          , "Mito"     , "L");
    legend->AddEntry(gEstimatedTemp_K18BR_May, "Estimated", "L");
    legend->Draw();

    gPad->Print(ofilename.data());
  }

  const Double_t tempmin = TMath::MinElement(gEstimatedTemp_K18BR_May->GetN(), gEstimatedTemp_K18BR_May->GetY());
  const Double_t tempmax = TMath::MaxElement(gEstimatedTemp_K18BR_May->GetN(), gEstimatedTemp_K18BR_May->GetY());
  std::cout << "TempMin = " << tempmin << ", TempMax = " << tempmax << ", dT = " << tempmax - tempmin << std::endl;
  std::ofstream ofile("./data/compJmaAndEasiroc.txt");
  if (ofile) {
    ofile << "TempMin = " << tempmin << ", TempMax = " << tempmax << ", dT = " << tempmax - tempmin << std::endl;
  }

  // {
  //   TCanvas::MakeDefCanvas();
  //   gPad->SetGrid();
  //   TMultiGraph* mg = new TMultiGraph();
  //   mg->Add(gDtemp_Hitachi_March, "P");
  //   mg->Add(gDtemp_Mito_March   , "P");
  //   mg->Draw("A");

  //   TLegend* legend = new TLegend(0.78, 0.92 - 0.05 * 2, 0.92, 0.92);
  //   legend->AddEntry(gDtemp_Hitachi_March, "Hitachi"  , "L");
  //   legend->AddEntry(gDtemp_Mito_March   , "Mito"     , "L");
  //   legend->Draw();

  //   gPad->Print(ofilename.data());
  // }

  {
    TCanvas::MakeDefCanvas();
    gPad->SetGrid();
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gCtemp_Hitachi_March1, "P");
    mg->Add(gCtemp_Hitachi_March2, "P");
    mg->Add(gCtemp_Hitachi_March3, "P");
    mg->Draw("A");

    gPad->Print(ofilename.data());
  }

  {
    TCanvas::MakeDefCanvas();
    gPad->SetGrid();
    // TMultiGraph* mg = new TMultiGraph();
    // mg->Add(gCtemp_Mito_March1, "P");
    // mg->Add(gCtemp_Mito_March2, "P");
    // mg->Add(gCtemp_Mito_March3, "P");
    // mg->Draw("A");
    gCtemp_Mito_March1->Draw("AP");

    if (TPaveStats* st = dynamic_cast<TPaveStats*>(gCtemp_Mito_March1->FindObject("stats"))) {
      const Double_t dx = 0.97 - st->GetX2NDC();
      const Double_t dy = 0.97 - st->GetY2NDC();
      st->SetX2NDC(st->GetX2NDC() + dx);
      st->SetX1NDC(st->GetX1NDC() + dx);
      st->SetY2NDC(st->GetY2NDC() + dx);
      st->SetY1NDC(st->GetY1NDC() + dx);
    }

    gPad->Print(ofilename.data());
  }

  gPad->Print((ofilename + "]").data());

  return 0;
}
