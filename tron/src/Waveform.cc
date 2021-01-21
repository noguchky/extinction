#include "Waveform.hh"
#include "TMath.h"
#include "LinearLeastSquareMethod.hh"
#include "PeakAnalyzer.hh"
#include "Math.hh"

template <class THist>
Tron::Waveform<THist>::Waveform()
  : THist() {
}

template <class THist>
Tron::Waveform<THist>::Waveform(const Char_t* name, const Char_t* title, Int_t n)
  : THist(name, title, n, 0, n) {
}

template <class THist>
Tron::Waveform<THist>::Waveform(const Char_t* name, const Char_t* title, Int_t n, const YValue_t* yvalues, XValue_t xwidth)
  : THist(name, title, n, 0, n * xwidth) {
  std::copy(yvalues, yvalues + n, this->GetArray());
}

template <class THist>
Tron::Waveform<THist>::Waveform(const Char_t* name, const Char_t* title, Int_t n, Double_t xmin, Double_t xmax, std::initializer_list<Waveform*> waveforms)
  : THist(name, title, n, xmin, xmax) {
  for (Int_t bin = 0; bin <= n; ++bin) {
    const XValue_t x = this->GetBinCenter(bin);
    YValue_t sum = YValue_t();
    for (auto&& wf : waveforms) {
      sum += wf->GetValue(x);
    }
    this->GetArray()[bin] = sum / waveforms.size();
  }
}

template <class THist>
Tron::Waveform<THist>::Waveform(const THist& hist)
  : THist(hist) {
}

template <class THist>
Tron::Waveform<THist>::Waveform(const Waveform& wf)
  : THist(wf) {
}

template <class THist>
Tron::Waveform<THist>::~Waveform() {
}

template <class THist>
void Tron::Waveform<THist>::Draw(Option_t* option) {
  this->THist::Draw(option);
}

template <class THist>
void Tron::Waveform<THist>::SetYValue(const YValue_t* yvalues) {
  std::copy(yvalues, yvalues + this->GetNbinsX(), this->GetArray());
}

template <class THist>
void Tron::Waveform<THist>::SetYValue(const YValue_t* yvalues, Int_t cidx) {
  auto cyvalue = std::vector<YValue_t>(this->GetNbinsX());
  std::rotate_copy(yvalues, yvalues + cidx, yvalues + this->GetNbinsX(), cyvalue.begin());
  SetYValue(cyvalue.data());
}

template <class THist>
void Tron::Waveform<THist>::MovingAverage(Int_t mobility) {
  std::vector<YValue_t> average(this->GetNbinsX(), YValue_t());
  Math::MovingAverage(mobility, this->GetNbinsX(), this->GetArray(), average.data());
  std::copy(average.begin(), average.end(), this->GetArray());
}

template <class THist>
void Tron::Waveform<THist>::CorrectYPattern(const YValue_t* ypattern) {
  for (YValue_t *it = this->GetArray(), *end = it + this->GetNbinsX(), *base = (YValue_t*)ypattern; it != end; *it++ -= *base++);
}

template <class THist>
void Tron::Waveform<THist>::CorrectXPattern(const XValue_t* xpattern) {
  auto array = this->GetXaxis()->GetXbins()->fArray;
  for (XValue_t *it = array + 1, *end = it + this->GetNbinsX() + 1, *tIt = (XValue_t*)xpattern; it != end; *it = *(it - 1) + *tIt++, it++);
  this->GetXaxis()->Set(this->GetNbinsX(), array);
}

template <class THist>
void Tron::Waveform<THist>::CorrectYPattern(const YValue_t* ypattern, Int_t cidx) {
  auto cypattern = std::vector<YValue_t>(this->GetNbinsX());
  std::rotate_copy(ypattern, ypattern + cidx, ypattern + this->GetNbinsX(), cypattern.begin());
  CorrectYPattern(cypattern.data());
}

template <class THist>
void Tron::Waveform<THist>::CorrectXPattern(const XValue_t* xpattern, Int_t cidx) {
  auto cxpattern = std::vector<XValue_t>(this->GetNbinsX());
  std::rotate_copy(xpattern, xpattern + cidx, xpattern + this->GetNbinsX(), cxpattern.begin());
  CorrectXPattern(cxpattern.data());
}

template <class THist>
void Tron::Waveform<THist>::CorrectYOffset(YValue_t yoffset) {
  for (auto it = this->GetArray(), end = (it + this->GetNbinsX()); it != end; *it++ -= yoffset);
}

template <class THist>
void Tron::Waveform<THist>::CorrectXOffset(XValue_t xoffset) {
  auto array = this->GetXaxis()->GetXbins()->fArray;
  for (auto it = array, end = it + this->GetNbinsX() + 1; it != end; *it++ -= xoffset);
  this->GetXaxis()->Set(this->GetNbinsX(), array);
}

template <class THist>
typename Tron::Waveform<THist>::YValue_t Tron::Waveform<THist>::GetValue(XValue_t x) const {
  const Int_t bin = this->GetXaxis()->FindBin(x);
  
  const XValue_t x1 = this->GetXaxis()->GetBinCenter(bin - 1);
  const XValue_t x2 = this->GetXaxis()->GetBinCenter(bin);
  const XValue_t x3 = this->GetXaxis()->GetBinCenter(bin + 1);

  YValue_t* it = (YValue_t*)(this->GetArray() + bin - 1);
  const YValue_t y1 = *it++;
  const YValue_t y2 = *it++;
  const YValue_t y3 = *it++;

  return x < x2 ? Math::Interpolate(x, x1, y1, x2, y2) : Math::Interpolate(x, x2, y2, x3, y3);
}

template <class THist>
typename Tron::Waveform<THist>::YValue_t Tron::Waveform<THist>::GetError(XValue_t x) const {
  const Int_t bin = this->GetXaxis()->FindBin(x);
  return this->GetBinError(bin);
}

template <class THist>
Double_t Tron::Waveform<THist>::GetArea(XValue_t from, XValue_t to) const {
  const Int_t binx1 = this->GetXaxis()->FindBin(from);
  const Int_t binx2 = this->GetXaxis()->FindBin(to);
  return this->Integral(binx1, binx2, "w");
}

template <class THist>
typename Tron::Waveform<THist>::YValue_t Tron::Waveform<THist>::GetAverage(XValue_t from, XValue_t to) const {
  const Int_t binx1 = this->GetXaxis()->FindBin(from);
  const Int_t binx2 = this->GetXaxis()->FindBin(to);
  return this->Integral(binx1, binx2) / (binx2 - binx1 + 1);
}

template <class THist>
typename Tron::Waveform<THist>::YValue_t Tron::Waveform<THist>::GetRms(XValue_t from, XValue_t to) const {
  const Int_t binx1 = this->GetXaxis()->FindBin(from);
  const Int_t binx2 = this->GetXaxis()->FindBin(to);
  return TMath::RMS(binx2 - binx1 + 1, this->GetArray() + binx1);
}

template <class THist>
Double_t Tron::Waveform<THist>::GetSlope(XValue_t from, XValue_t to) const {
  const Int_t binx1 = this->GetXaxis()->FindBin(from);
  const Int_t binx2 = this->GetXaxis()->FindBin(to);

  // y = ax + b
  std::vector<Double_t> xi, yi;
  for (Int_t bin = binx1; bin <= binx2; ++bin) {
    xi.push_back(this->GetBinCenter (bin));
    yi.push_back(this->GetBinContent(bin));
  }

  auto lsmResult = LinearLeastSquareMethod::Fit(xi.size(), xi.data(), yi.data());
  return lsmResult.Slope;
}

template <class THist>
typename Tron::Waveform<THist>::YValue_t Tron::Waveform<THist>::GetMaximumY(XValue_t from, XValue_t to) const {
  const Int_t binx1 = this->GetXaxis()->FindBin(from);
  const Int_t binx2 = this->GetXaxis()->FindBin(to);
  return TMath::MaxElement(binx2 - binx1 + 1, this->GetArray() + binx1);
}

template <class THist>
typename Tron::Waveform<THist>::YValue_t Tron::Waveform<THist>::GetMinimumY(XValue_t from, XValue_t to) const {
  const Int_t binx1 = this->GetXaxis()->FindBin(from);
  const Int_t binx2 = this->GetXaxis()->FindBin(to);
  return TMath::MinElement(binx2 - binx1 + 1, this->GetArray() + binx1);
}

template <class THist>
typename Tron::Waveform<THist>::XValue_t Tron::Waveform<THist>::GetMaximumX(XValue_t from, XValue_t to) const {
  const Int_t maxbin = PeakAnalyzer::FindMaximumBin(this, from, to);
  return this->GetXaxis()->GetBinCenter(maxbin);
}

template <class THist>
typename Tron::Waveform<THist>::XValue_t Tron::Waveform<THist>::GetMinimumX(XValue_t from, XValue_t to) const {
  const Int_t minbin = PeakAnalyzer::FindMinimumBin(this, from, to);
  return this->GetXaxis()->GetBinCenter(minbin);
}

template <class THist>
typename Tron::Waveform<THist>::XValue_t Tron::Waveform<THist>::GetPosPeakWidth(XValue_t from, XValue_t to, YValue_t threshold) const {
  const XValue_t maxx = GetMaximumX(from, to);
  const XValue_t tlower = PeakAnalyzer::Find(this, maxx, TMath::Min(from, to), threshold);
  const XValue_t tupper = PeakAnalyzer::Find(this, maxx, TMath::Max(from, to), threshold);
  return tupper - tlower;
}

template <class THist>
typename Tron::Waveform<THist>::XValue_t Tron::Waveform<THist>::GetNegPeakWidth(XValue_t from, XValue_t to, YValue_t threshold) const {
  const XValue_t minx = GetMinimumX(from, to);
  const XValue_t tlower = PeakAnalyzer::Find(this, minx, TMath::Min(from, to), threshold);
  const XValue_t tupper = PeakAnalyzer::Find(this, minx, TMath::Max(from, to), threshold);
  return tupper - tlower;
}

template <class THist>
typename Tron::Waveform<THist>::XValue_t Tron::Waveform<THist>::FindCrossingX(XValue_t from, XValue_t to, YValue_t threshold,
                                                                              XValue_t dx, YValue_t dy) const {
  return PeakAnalyzer::Find(this, from, to, threshold, dx, dy);
}

template <class THist>
typename Tron::Waveform<THist>::XValue_t Tron::Waveform<THist>::FindCrossingXLSM(XValue_t from, XValue_t to, YValue_t threshold) const {
  return PeakAnalyzer::FindByLSM(this, from, to, threshold);
}

template <class THist>
typename Tron::Waveform<THist>* Tron::Waveform<THist>::Create(const Char_t* name,
                                                              const Char_t* title,
                                                              Int_t n,
                                                              Int_t cidx,
                                                              const YValue_t* yvalues,
                                                              const YValue_t* ypattern,
                                                              const XValue_t* xpattern) {
  auto wf = new Waveform(name, title, n);
  wf->SetYValue      (yvalues , cidx);
  wf->CorrectYPattern(ypattern, cidx);
  wf->CorrectXPattern(xpattern, cidx);
  return wf;
}

template class Tron::Waveform<TH1D>;
template class Tron::Waveform<TH1F>;
template class Tron::Waveform<TH1I>;
template class Tron::Waveform<TH1S>;
template class Tron::Waveform<TH1C>;
