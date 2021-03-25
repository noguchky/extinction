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

#include "String.hh"
#include "Math.hh"
#include "LinearLeastSquareMethod.hh"
#include "Integrator.hh"
#include "ProgressBar.hh"

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

Int_t compareTemp() {
  // Get temperature graphs
  TGraph* gTemp_K18BR_March   = GetGraphFromFile("K1.8BR(March)" , "./datae/20210310_k18br"   ".root", "temp1");
  TGraph* gTemp_Hitachi_March = GetGraphFromFile("Hitachi(March)", "./dataj/jma_hitachi_2021-03.root", "gTemp");
  TGraph* gTemp_Mito_March    = GetGraphFromFile("Mito(March)"   , "./dataj/jma_mito_" "2021-03.root", "gTemp");
  TGraph* gTemp_Hitachi_May   = GetGraphFromFile("Hitachi(May)"  , "./dataj/jma_hitachi_2020-05.root", "gTemp");
  TGraph* gTemp_Mito_May      = GetGraphFromFile("Mito(May)"     , "./dataj/jma_mito_" "2020-05.root", "gTemp");

  // Remove outlier points
  RemoveOutlierPoints(gTemp_K18BR_March, 0, 40);

  // Set style
  for (auto&& pair : std::map<TGraph*, std::tuple<Int_t, Int_t>> {
     { gTemp_K18BR_March   , { kGreen + 1, kFullSquare       } },
     { gTemp_Hitachi_March , { kBlue  + 1, kFullTriangleUp   } },
     { gTemp_Mito_March    , { kRed   + 1, kFullTriangleDown } },
     { gTemp_Hitachi_May   , { kBlue  + 1, kFullTriangleUp   } },
     { gTemp_Mito_May      , { kRed   + 1, kFullTriangleDown } },
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
  gDtemp_Hitachi_March->SetTitle("Hitachi (March);;T_{Hitachi}-20 [#circC]");
  gDtemp_Mito_March   ->SetTitle("Mito " "(March);;T_{Mito" "}-20 [#circC]");

  TGraph* gCtemp_Hitachi_March = (TGraph*)gTemp_Hitachi_March->Clone();
  TGraph* gCtemp_Mito_March    = (TGraph*)gTemp_Mito_May     ->Clone();
  gCtemp_Hitachi_March->Set(0);
  gCtemp_Mito_March   ->Set(0);
  gCtemp_Hitachi_March->SetMinimum();
  gCtemp_Mito_March   ->SetMinimum();
  gCtemp_Hitachi_March->SetMaximum();
  gCtemp_Mito_March   ->SetMaximum();
  gCtemp_Hitachi_March->GetXaxis()->SetTimeDisplay(false);
  gCtemp_Mito_March   ->GetXaxis()->SetTimeDisplay(false);
  gCtemp_Hitachi_March->SetTitle("Hitachi (March);T_{Hitachi}-20 [#circC];T_{K1.8BR}-20 [#circC]");
  gCtemp_Mito_March   ->SetTitle("Mito " "(March);T_{Mito" "}-20 [#circC];T_{K1.8BR}-20 [#circC]");

  const Double_t hour   = TDatime(2000, 1, 1, 1, 0, 0).Convert() - TDatime(2000, 1, 1, 0, 0, 0).Convert();
  const Double_t weight = 0.7;

  const Double_t tmin_March = Tron::Math::Greatest(TMath::MinElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                   TMath::MinElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                   TMath::MinElement(gTemp_K18BR_March  ->GetN(), gTemp_K18BR_March  ->GetX()));
  const Double_t tmax_March = Tron::Math::Least   (TMath::MaxElement(gTemp_Hitachi_March->GetN(), gTemp_Hitachi_March->GetX()),
                                                   TMath::MaxElement(gTemp_Mito_March   ->GetN(), gTemp_Mito_March   ->GetX()),
                                                   TMath::MaxElement(gTemp_K18BR_March  ->GetN(), gTemp_K18BR_March  ->GetX()));
  for (Int_t i = 60, n = gTemp_K18BR_March->GetN(); i < n; ++i) {
    const Double_t t = gTemp_K18BR_March->GetX()[i];
    if (tmin_March <= t && t <= tmax_March) {
      gDtemp_Hitachi_March->SetPoint(gDtemp_Hitachi_March->GetN(), t, gTemp_Hitachi_March->Eval(t) - 20.0);
      gDtemp_Mito_March   ->SetPoint(gDtemp_Mito_March   ->GetN(), t, gTemp_Mito_March   ->Eval(t) - 20.0);

      gCtemp_Hitachi_March->SetPoint(gCtemp_Hitachi_March->GetN(), (1.0 - weight) * (gTemp_Hitachi_March->Eval(t - 0.0 * hour) - 20.0) + weight * (gTemp_Hitachi_March->Eval(t - 2.0 * hour) - 20.0), gTemp_K18BR_March->Eval(t) - 20.0);
      gCtemp_Mito_March   ->SetPoint(gCtemp_Mito_March   ->GetN(), (1.0 - weight) * (gTemp_Mito_March   ->Eval(t - 0.0 * hour) - 20.0) + weight * (gTemp_Mito_March   ->Eval(t - 2.0 * hour) - 20.0), gTemp_K18BR_March->Eval(t) - 20.0);
    }
  }

  TF1* fline = new TF1("fline", "[0]*x + [1]");
  fline->SetParNames("grad", "intercept");
  fline->SetParameters(1.0, 0.0);

  gCtemp_Hitachi_March->Fit(fline, "", "");
  gCtemp_Mito_March   ->Fit(fline, "", "");

  const Double_t grad      = gCtemp_Mito_March->GetFunction(fline->GetName())->GetParameter(0);
  const Double_t intercept = gCtemp_Mito_March->GetFunction(fline->GetName())->GetParameter(1);
  auto estimate = [&](Double_t y, Double_t y0) {
    return (((1.0 - weight) * (y - 20.0) + weight * (y0 - 20.0)) * grad + intercept) + 20.0;
  };

  // Create stimated
  TGraph* gEstimatedTemp_K18BR_March = (TGraph*)gTemp_K18BR_March->Clone();
  TGraph* gEstimatedTemp_K18BR_May   = (TGraph*)gTemp_K18BR_March->Clone();
  gEstimatedTemp_K18BR_March->Set(0);
  gEstimatedTemp_K18BR_May  ->Set(0);
  gEstimatedTemp_K18BR_March->SetLineColor(kYellow + 1);
  gEstimatedTemp_K18BR_May  ->SetLineColor(kYellow + 1);

  for (Int_t i = 12, n = gTemp_Mito_March->GetN(); i < n; ++i) {
    gEstimatedTemp_K18BR_March->SetPoint(gEstimatedTemp_K18BR_March->GetN(),
                                         gTemp_Mito_March->GetX()[i],
                                         estimate(gTemp_Mito_March->GetY()[i - 0], gTemp_Mito_March->GetY()[i - 12]));
  }
  for (Int_t i = 12, n = gTemp_Mito_May->GetN(); i < n; ++i) {
    gEstimatedTemp_K18BR_May->SetPoint(gEstimatedTemp_K18BR_May->GetN(),
                                       gTemp_Mito_May->GetX()[i],
                                       estimate(gTemp_Mito_May->GetY()[i - 0], gTemp_Mito_May->GetY()[i - 12]));
  }

  const std::string ofilename = "./compareTemp.pdf";
  TCanvas::MakeDefCanvas();
  gPad->Print((ofilename + "[").data());
  
  // Draw graphs and hists
  {
    gPad->SetGrid();
    gTemp_Hitachi_March       ->Draw("AL");
    gTemp_Mito_March          ->Draw( "L");
    gTemp_K18BR_March         ->Draw( "L");
    gEstimatedTemp_K18BR_March->Draw( "L");

    TLegend* legend = new TLegend(0.78, 0.92 - 0.05 * 4, 0.92, 0.92);
    legend->AddEntry(gTemp_K18BR_March          , "K1.8BR"   , "L");
    legend->AddEntry(gTemp_Hitachi_March        , "Hitachi"  , "L");
    legend->AddEntry(gTemp_Mito_March           , "Mito"     , "L");
    legend->AddEntry(gEstimatedTemp_K18BR_March , "Estimated", "L");
    legend->Draw();

    gPad->Print(ofilename.data());
  }

  {
    TCanvas::MakeDefCanvas();
    gPad->SetGrid();
    gTemp_Hitachi_May       ->Draw("AL");
    gTemp_Mito_May          ->Draw( "L");
    gEstimatedTemp_K18BR_May->Draw( "L");

    TLegend* legend = new TLegend(0.78, 0.92 - 0.05 * 3, 0.92, 0.92);
    legend->AddEntry(gTemp_Hitachi_May       , "Hitachi"  , "L");
    legend->AddEntry(gTemp_Mito_May          , "Mito"     , "L");
    legend->AddEntry(gEstimatedTemp_K18BR_May, "Estimated", "L");
    legend->Draw();

    gPad->Print(ofilename.data());
  }

  const Double_t tempmin = TMath::MinElement(gEstimatedTemp_K18BR_May->GetN(), gEstimatedTemp_K18BR_May->GetY());
  const Double_t tempmax = TMath::MaxElement(gEstimatedTemp_K18BR_May->GetN(), gEstimatedTemp_K18BR_May->GetY());
  std::cout << "TempMin = " << tempmin << ", TempMax = " << tempmax << ", dT = " << tempmax - tempmin << std::endl;
  std::ofstream ofile("./compareTemp.txt");
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
    gCtemp_Hitachi_March->Draw("AP");

    gPad->Print(ofilename.data());
  }

  {
    TCanvas::MakeDefCanvas();
    gPad->SetGrid();
    gCtemp_Mito_March->Draw("AP");

    gPad->Print(ofilename.data());
  }

  gPad->Print((ofilename + "]").data());

  return 0;
}
