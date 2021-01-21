#include "PeakAnalyzer.hh"
#include "ExtremumSearcher.hh"
#include "RootFinder.hh"
#include "LinearLeastSquareMethod.hh"
#include "TSF1.hh"

double Tron::PeakAnalyzer::GetFWHM(const Func_t& func, double from, double to, int sign, double base, double dx) {
  static const auto& extremumSearcher = ExtremumSearcher::GoldenSectionSearchMethod::Get();
  ExtremumSearcher::Result_t extremum;
  if (sign >= 0) {
    extremum = extremumSearcher.FindMaximum(func, from, to, dx);
  } else {
    extremum = extremumSearcher.FindMinimum(func, from, to, dx);
  }

  const double threshold = (extremum.Y + base) * 0.5;
  auto funcMinusThreshold = [&](double x) {
    return func(x) - threshold;
  };

  static const auto& rootFinder = RootFinder::BisectionMethod::Get();
  const double xupper = rootFinder.Find(funcMinusThreshold, extremum.X, std::max(from, to), dx);
  const double xlower = rootFinder.Find(funcMinusThreshold, extremum.X, std::min(from, to), dx);
  
  return xupper - xlower;
}

double Tron::PeakAnalyzer::GetHWHM(const Func_t& func, double from, double to, int sign, double base, double dx) {
  static const auto& extremumSearcher = ExtremumSearcher::GoldenSectionSearchMethod::Get();
  ExtremumSearcher::Result_t extremum;
  if (sign >= 0) {
    extremum = extremumSearcher.FindMaximum(func, from, to, dx);
  } else {
    extremum = extremumSearcher.FindMinimum(func, from, to, dx);
  }

  const double threshold = (extremum.Y + base) * 0.5;
  auto funcMinusThreshold = [&](double x) {
    return func(x) - threshold;
  };

  static const auto& rootFinder = RootFinder::BisectionMethod::Get();
  double hwhm = 0.0;
  if (from <= to) {
    const double xupper = rootFinder.Find(funcMinusThreshold, extremum.X, std::max(from, to), dx);
    hwhm = xupper - extremum.X;
  } else {
    const double xlower = rootFinder.Find(funcMinusThreshold, extremum.X, std::min(from, to), dx);
    hwhm = extremum.X - xlower;
  }
  
  return hwhm;
}

Int_t Tron::PeakAnalyzer::FindMaximumBin(const TH1* h1, Double_t from, Double_t to) {
  const Int_t binx1 = h1->GetXaxis()->FindBin(from);
  const Int_t binx2 = h1->GetXaxis()->FindBin(to);
  Int_t   maxbin = binx1;
  Double_t maxval = h1->GetBinContent(maxbin);
  for (Int_t bin = binx1 + 1; bin <= binx2; ++bin) {
    const Double_t val = h1->GetBinContent(bin);
    if (val < maxval) {
      maxval = val;
      maxbin = bin;
    }
  }
  return maxbin;
}

Int_t Tron::PeakAnalyzer::FindMinimumBin(const TH1* h1, Double_t from, Double_t to) {
  const Int_t binx1 = h1->GetXaxis()->FindBin(from);
  const Int_t binx2 = h1->GetXaxis()->FindBin(to);
  Int_t   minbin = binx1;
  Double_t minval = h1->GetBinContent(minbin);
  for (Int_t bin = binx1 + 1; bin <= binx2; ++bin) {
    const Double_t val = h1->GetBinContent(bin);
    if (val < minval) {
      minval = val;
      minbin = bin;
    }
  }
  return minbin;
}

Double_t Tron::PeakAnalyzer::Find(const TH1* h1, Double_t from, Double_t to, Double_t y, Double_t dx, Double_t dy) {
  //--- Initialize
  std::string name = Form("%s_sf1", h1->GetName());
  const Double_t xmin = TMath::Min(from, to);
  const Double_t xmax = TMath::Max(from, to);
  TSF1* sf1 = new TSF1(name.data(), h1, xmin, xmax);
  const auto func = [&](Double_t x) {
    return sf1->Eval(x) - y;
  };

  //--- Find solution
  static const auto& rootFinder = RootFinder::BisectionMethod::Get();
  const Double_t solution = rootFinder.Find(func, from, to, dx, dy);

  //--- Finalize
  delete sf1; sf1 = nullptr;
  
  return solution;
}

Double_t Tron::PeakAnalyzer::Find(const TF1* f1, Double_t from, Double_t to, Double_t y, Double_t dx, Double_t dy) {
  //--- Initialize
  const auto func = [&](Double_t x) {
    return f1->Eval(x) - y;
  };

  //--- Find solution
  static const auto& rootFinder = RootFinder::BisectionMethod::Get();
  const Double_t solution = rootFinder.Find(func, from, to, dx, dy);
  
  return solution;
}

Double_t Tron::PeakAnalyzer::FindByLSM(const TH1* h1, Double_t from, Double_t to, Double_t y) {
  const Int_t bmin = h1->GetXaxis()->FindBin(from);
  const Int_t bmax = h1->GetXaxis()->FindBin(to);

  // x = ay + b
  std::vector<Double_t> xi, yi;
  for (Int_t bin = bmin; bin <= bmax; ++bin) {
    xi.push_back(h1->GetBinCenter(bin));
    yi.push_back(h1->GetBinContent(bin));
  }

  auto lsmResult = LinearLeastSquareMethod::Fit(yi.size(), yi.data(), xi.data());
  return lsmResult.Eval(y);
}

Double_t Tron::PeakAnalyzer::GetFWHM(const TH1* hist, Double_t from, Double_t to, Int_t sign, Double_t base) {
  const Int_t    maxbin  = sign >= 0 ? FindMaximumBin(hist, from, to) : FindMinimumBin(hist, from, to);
  const Double_t maxy    = hist->GetBinContent(maxbin);
  const Double_t maxx    = hist->GetBinCenter(maxbin);
  const Double_t halfmax = (maxy - base) * 0.5;
  const Double_t xupper  = Find(hist, maxx,   to, halfmax + base);
  const Double_t xlower  = Find(hist, maxx, from, halfmax + base);
  return xupper - xlower;
}

Double_t Tron::PeakAnalyzer::GetFWHMFast(const TH1* hist, Double_t from, Double_t to, Int_t sign, Double_t base) {
  const Int_t    nbins   = hist->GetNbinsX();
  const Int_t    maxbin  = sign >= 0 ? FindMaximumBin(hist, from, to) : FindMinimumBin(hist, from, to);
  const Double_t maxy    = hist->GetBinContent(maxbin);
  const Double_t halfmax = TMath::Abs(maxy - base) * 0.5;

  Int_t binlower = maxbin;
  while (TMath::Abs(hist->GetBinContent(binlower) - base) > halfmax && --binlower >     0);

  Int_t binupper = maxbin;
  while (TMath::Abs(hist->GetBinContent(binupper) - base) > halfmax && ++binupper < nbins);

  return hist->GetBinCenter(binupper) - hist->GetBinCenter(binlower);
}

Double_t Tron::PeakAnalyzer::GetHWHMFast(const TH1* hist, Double_t from, Double_t to, Int_t sign, Double_t base) {
  const Int_t    nbins   = hist->GetNbinsX();
  const Int_t    maxbin  = sign >= 0 ? FindMaximumBin(hist, from, to) : FindMinimumBin(hist, from, to);
  const Double_t maxy    = hist->GetBinContent(maxbin);
  const Double_t halfmax = TMath::Abs(maxy - base) * 0.5;

  Int_t bin = maxbin;
  if (from <= to) {
    while (TMath::Abs(hist->GetBinContent(bin) - base) > halfmax && --bin >     0);
  } else {
    while (TMath::Abs(hist->GetBinContent(bin) - base) > halfmax && ++bin < nbins);
  }
  
  return TMath::Abs(hist->GetBinCenter(maxbin) - hist->GetBinCenter(bin));
}
