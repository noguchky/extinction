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
  if (!graph->GetN()) {
    std::cout << "[warning] all points were removed" << std::endl;
  }
}

Int_t compOndotoriAndGL820() {
  // Get temperature graphs
  TGraph* gTemp_Easiroc   = GetGraphFromFile("Easiroc" , "./data_easiroc/20210330_k18br.root"      , "temp1" );
  TGraph* gTemp_Ondotori  = GetGraphFromFile("Ondotori", "./data_ondotori/521455AD.1619407372.root", "gTemp" );
  TGraph* gTemp_GL820_A   = GetGraphFromFile("PT100-A" , "./data_gl820/210423-191135_UG.root"      , "gTempA");
  TGraph* gTemp_GL820_B   = GetGraphFromFile("PT100-B" , "./data_gl820/210423-191135_UG.root"      , "gTempB");
  TGraph* gTemp_GL820_C   = GetGraphFromFile("PT100-C" , "./data_gl820/210423-191135_UG.root"      , "gTempC");

  // Remove outlier points
  RemoveOutlierPoints(gTemp_Easiroc , 0, 40);
  RemoveOutlierPoints(gTemp_Ondotori, 0, 40);
  RemoveOutlierPoints(gTemp_GL820_A , 0, 40);
  RemoveOutlierPoints(gTemp_GL820_B , 0, 40);
  RemoveOutlierPoints(gTemp_GL820_C , 0, 40);

  // Set style
  for (auto&& pair : std::map<TGraph*, std::tuple<Int_t, Int_t>> {
     { gTemp_Easiroc , { kGreen  + 1, kFullSquare } },
     { gTemp_Ondotori, { kBlack     , kFullSquare } },
     { gTemp_GL820_A , { kRed    + 1, kFullSquare } },
     { gTemp_GL820_B , { kBlue   + 1, kFullSquare } },
     { gTemp_GL820_C , { kYellow + 1, kFullSquare } },
    }) {
    TGraph* graph = pair.first;
    const Int_t color  = std::get<0>(pair.second);
    const Int_t marker = std::get<1>(pair.second);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.4);
    graph->SetLineColor(color);
    graph->SetLineWidth(4);
    graph->SetMinimum(+10.0);
    graph->SetMaximum(+35.0);
    graph->SetTitle(";;Temperature [#circC]");
  }

  // Create difference
  TGraphErrors* gCtemp_GL820_A = new TGraphErrors();
  TGraphErrors* gCtemp_GL820_B = new TGraphErrors();
  TGraphErrors* gCtemp_GL820_C = new TGraphErrors();
  for (auto&& pair : std::map<TGraph*, std::tuple<std::string, std::string, Int_t, Int_t>> {
      { gCtemp_GL820_A , { "gCtemp_GL820_A", "PT100(A)" , kRed    + 1, kFullSquare } },
      { gCtemp_GL820_B , { "gCtemp_GL820_B", "PT100(B)" , kBlue   + 1, kFullSquare } },
      { gCtemp_GL820_C , { "gCtemp_GL820_C", "PT100(C)" , kYellow + 1, kFullSquare } },
    }) {
    TGraph* graph = pair.first;
    const std::string name   = std::get<0>(pair.second);
    const std::string label  = std::get<1>(pair.second);
    const Int_t       color  = std::get<2>(pair.second);
    const Int_t       marker = std::get<3>(pair.second);
    graph->SetName(name.data());
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.4);
    graph->SetLineColor(color);
    graph->SetLineWidth(4);
    graph->SetTitle(Form("%s;T_{%s} [#circC];T_{Ondotori} [#circC]", label.data(), label.data()));
  }

  const Double_t tmin = Tron::Math::Greatest(TMath::MinElement(gTemp_Ondotori->GetN(), gTemp_Ondotori->GetX()),
                                             TMath::MinElement(gTemp_GL820_A ->GetN(), gTemp_GL820_A ->GetX()),
                                             TMath::MinElement(gTemp_GL820_B ->GetN(), gTemp_GL820_B ->GetX()),
                                             TMath::MinElement(gTemp_GL820_C ->GetN(), gTemp_GL820_C ->GetX()));
  const Double_t tmax = Tron::Math::Least   (TMath::MaxElement(gTemp_Ondotori->GetN(), gTemp_Ondotori->GetX()),
                                             TMath::MaxElement(gTemp_GL820_A ->GetN(), gTemp_GL820_A ->GetX()),
                                             TMath::MaxElement(gTemp_GL820_B ->GetN(), gTemp_GL820_B ->GetX()),
                                             TMath::MaxElement(gTemp_GL820_C ->GetN(), gTemp_GL820_C ->GetX()));

  for (Int_t i = 0, n = gTemp_Ondotori->GetN(); i < n; ++i) {
    const Double_t t = gTemp_Ondotori->GetX()[i];
    if (tmin <= t && t <= tmax) {
      gCtemp_GL820_A->SetPoint(gCtemp_GL820_A->GetN(), gTemp_GL820_A->Eval(t), gTemp_Ondotori->Eval(t));
      gCtemp_GL820_B->SetPoint(gCtemp_GL820_B->GetN(), gTemp_GL820_B->Eval(t), gTemp_Ondotori->Eval(t));
      gCtemp_GL820_C->SetPoint(gCtemp_GL820_C->GetN(), gTemp_GL820_C->Eval(t), gTemp_Ondotori->Eval(t));
      gCtemp_GL820_A->SetPointError(gCtemp_GL820_A->GetN() - 1, 0.05, 0.05);
      gCtemp_GL820_B->SetPointError(gCtemp_GL820_B->GetN() - 1, 0.05, 0.05);
      gCtemp_GL820_C->SetPointError(gCtemp_GL820_C->GetN() - 1, 0.05, 0.05);
    }
  }

  TF1* fline = new TF1("fline", "[0]*x + [1]");
  fline->SetParNames("grad", "intercept");
  fline->SetParameters(1.0, 0.0);

  gCtemp_GL820_A->Fit(fline, "", "");
  gCtemp_GL820_B->Fit(fline, "", "");
  gCtemp_GL820_C->Fit(fline, "", "");

  const std::string ofilename = "./data/compOndotoriAndGL820.pdf";
  TCanvas::MakeDefCanvas();
  gPad->Print((ofilename + "[").data());
  
  // Draw graphs and hists
  {
    gPad->SetGrid();
    gTemp_Easiroc ->Draw("AP");
    gTemp_Ondotori->Draw( "P");
    gTemp_GL820_A ->Draw( "P");
    gTemp_GL820_B ->Draw( "P");
    gTemp_GL820_C ->Draw( "P");
    gTemp_Easiroc->GetXaxis()->SetLimits(TDatime(2021, 4, 24, 0, 0, 0).Convert(),
                                         TDatime(2021, 4, 27, 0, 0, 0).Convert());

    TLegend* legend = new TLegend(0.83, 0.97 - 0.04 * 5, 0.97, 0.97);
    legend->AddEntry(gTemp_Easiroc , "EASIROC" , "L");
    legend->AddEntry(gTemp_Ondotori, "Ondotori", "L");
    legend->AddEntry(gTemp_GL820_A , "PT100(A)", "L");
    legend->AddEntry(gTemp_GL820_B , "PT100(B)", "L");
    legend->AddEntry(gTemp_GL820_C , "PT100(C)", "L");
    legend->Draw();

    gPad->Print(ofilename.data());
  }

  {
    TCanvas::MakeDefCanvas();
    gPad->SetGrid();
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(";T_{PT100} [#circC];T_{Ondotori} [#circC]");
    mg->Add(gCtemp_GL820_A, "P");
    mg->Add(gCtemp_GL820_B, "P");
    mg->Add(gCtemp_GL820_C, "P");
    mg->Draw("A");

    gPad->Update();

    const Double_t st_left = 0.12;
    const Double_t st_margin = 0.02;
    Double_t st_a_bottom = 1;
    if (TPaveStats* st = dynamic_cast<TPaveStats*>(gCtemp_GL820_A->FindObject("stats"))) {
      st->SetLineColor(gCtemp_GL820_A->GetMarkerColor());
      st_a_bottom = st->GetY1NDC();
      const Double_t width = st->GetX2NDC() - st->GetX1NDC();
      st->SetX1NDC(st_left);
      st->SetX2NDC(st_left + width);
    }

    Double_t st_b_bottom = 1;
    if (TPaveStats* st = dynamic_cast<TPaveStats*>(gCtemp_GL820_B->FindObject("stats"))) {
      st->SetLineColor(gCtemp_GL820_B->GetMarkerColor());
      const Double_t height = st->GetY2NDC() - st->GetY1NDC();
      st->SetY1NDC(st_a_bottom - st_margin - height);
      st->SetY2NDC(st_a_bottom - st_margin);
      const Double_t width = st->GetX2NDC() - st->GetX1NDC();
      st->SetX1NDC(st_left);
      st->SetX2NDC(st_left + width);
      st_b_bottom = st->GetY1NDC();
    }

    if (TPaveStats* st = dynamic_cast<TPaveStats*>(gCtemp_GL820_C->FindObject("stats"))) {
      st->SetLineColor(gCtemp_GL820_C->GetMarkerColor());
      const Double_t height = st->GetY2NDC() - st->GetY1NDC();
      st->SetY1NDC(st_b_bottom - st_margin - height);
      st->SetY2NDC(st_b_bottom - st_margin);
      const Double_t width = st->GetX2NDC() - st->GetX1NDC();
      st->SetX1NDC(st_left);
      st->SetX2NDC(st_left + width);
    }

    gPad->Print(ofilename.data());
  }

  gPad->Print((ofilename + "]").data());

  const std::string dfilename = "./data/compOndotoriAndGL820.txt";
  std::ofstream dfile(dfilename);
  if (!dfile) {
    std::cout << "[warning] data file is not opened, " << dfilename << std::endl;
  } else {
    for (auto&& graph : std::vector<TGraph*>{
        gCtemp_GL820_A,
        gCtemp_GL820_B,
        gCtemp_GL820_C,
      }) {
      TF1* func = graph->GetFunction(fline->GetName());
      dfile << func->GetParameter(0) << "\t" << func->GetParameter(1) << std::endl;
    }
  }
  
  return 0;
}
