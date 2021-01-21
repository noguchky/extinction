#ifndef Tron_Waveform_hh
#define Tron_Waveform_hh

#include <vector>
#include <type_traits>
#include "TH1.h"
#include "Types.hh"

namespace Tron {

  template <class THist>
  class Waveform : public THist {
  public:
    using XValue_t = Double_t;
    using YValue_t = typename std::remove_reference<decltype(*(THist().GetArray()))>::type;

    Waveform();
    Waveform(const Char_t* name, const Char_t* title, Int_t n);
    Waveform(const Char_t* name, const Char_t* title,
             Int_t n, const YValue_t* yvalues, XValue_t xwidth = 1.0);
    Waveform(const Char_t* name, const Char_t* title,
             Int_t n, Double_t xmin, Double_t xmax,
             std::initializer_list<Waveform*> waveforms); // Average
    Waveform(const THist& hist);
    Waveform(const Waveform& wf);

    virtual ~Waveform();

    virtual void Draw(Option_t* option = "") override;

    void     SetYValue(const YValue_t* yvalues);
    void     SetYValue(const YValue_t* yvalues, Int_t cidx);
    void     MovingAverage(Int_t mobility = 3);
    void     CorrectYPattern(const YValue_t* ypattern); // Pattern Correction
    void     CorrectYPattern(const YValue_t* ypattern, Int_t cidx);
    void     CorrectXPattern(const XValue_t* xpattern); // Timing Calibtation
    void     CorrectXPattern(const XValue_t* xpattern, Int_t cidx);
    void     CorrectYOffset(YValue_t xoffset);
    void     CorrectXOffset(XValue_t yoffset);

    YValue_t GetValue        (XValue_t x) const;
    YValue_t GetError        (XValue_t x) const;
    Double_t GetArea         (XValue_t from, XValue_t to) const;
    YValue_t GetAverage      (XValue_t from, XValue_t to) const;
    YValue_t GetRms          (XValue_t fron, XValue_t to) const;
    Double_t GetSlope        (XValue_t from, XValue_t to) const;
    YValue_t GetMaximumY     (XValue_t from, XValue_t to) const;
    YValue_t GetMinimumY     (XValue_t from, XValue_t to) const;
    XValue_t GetMaximumX     (XValue_t from, XValue_t to) const;
    XValue_t GetMinimumX     (XValue_t from, XValue_t to) const;
    XValue_t GetPosPeakWidth (XValue_t from, XValue_t to, YValue_t threshold = 0.0) const;
    XValue_t GetNegPeakWidth (XValue_t from, XValue_t to, YValue_t threshold = 0.0) const;
    XValue_t FindCrossingX   (XValue_t from, XValue_t to, YValue_t threshold,
                              XValue_t dx = 1.0, YValue_t dy = 1.0e-3) const;
    XValue_t FindCrossingXLSM(XValue_t from, XValue_t to, YValue_t threshold) const;

    static Waveform* Create(const Char_t* name,
                            const Char_t* title,
                            Int_t n,
                            Int_t cidx,
                            const YValue_t* yvalues,
                            const YValue_t* ypattern,
                            const XValue_t* xpattern);
  };

}

#endif
