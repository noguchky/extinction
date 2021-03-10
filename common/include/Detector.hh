#ifndef Extinction_Detector_hh
#define Extinction_Detector_hh

#include <vector>
#include "Units.hh"
#include "TH2.h"
#include "TAxis.h"
#include "TList.h"
#include "TLine.h"
#include "TMath.h"

namespace Extinction {

  // Global  Channel = Extinction Detector -> [0:131], Hodoscope -> [132:148], ...
  // (Local) Channel = Extinction Detector -> [0:131], Hodoscope -> [0:15], ...

  namespace ExtinctionDetector {
    constexpr std::size_t NofChannels         = 132U;
    constexpr std::size_t GlobalChannelOffset =   0U;

    struct ChannelType {
      enum {
            Center,
            TopBottom,
            LeftRight,
      };
    };
    constexpr std::size_t NofCenter    = 112U;
    constexpr std::size_t NofTopBottom =  16U;
    constexpr std::size_t NofLeftRight =   4U;
    struct ChannelOffset {
      enum {
            Center    = 0,
            TopBottom = Center + NofCenter,
            LeftRight = TopBottom + NofTopBottom,
            N         = LeftRight + NofLeftRight,
      };
    };

    inline Bool_t Contains(std::size_t globalChannel) {
      // return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
      return globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }

    inline Int_t GetType(std::size_t channel) {
      // if        (channel < ChannelOffset::Center   ) {
      //   return -1;
      // } else
      if (channel < ChannelOffset::TopBottom) {
        return ChannelType::Center;
      } else if (channel < ChannelOffset::LeftRight) {
        return ChannelType::TopBottom;
      } else if (channel < ChannelOffset::N        ) {
        return ChannelType::LeftRight;
      } else {
        return -1;
      }
    }

    inline Double_t GetWidth(Int_t type) {
      switch (type) {
      case ChannelType::Center:    return   8.0 * mm;
      case ChannelType::TopBottom: return  56.0 * mm;
      case ChannelType::LeftRight: return 100.0 * mm;
      }
      return 0.0;
    }

    inline Double_t GetHeight(Int_t type) {
      switch (type) {
      case ChannelType::Center:    return  60.0 * mm;
      case ChannelType::TopBottom: return 140.0 * mm;
      case ChannelType::LeftRight: return 200.0 * mm;
      }
      return 0.0;
    }

    inline Double_t GetArea(Int_t type) {
      return GetWidth(type) * GetHeight(type);
    }

    inline TH2* CreateHitMap(const std::string& name) {
      std::vector<Double_t> xbins, ybins;
      {
        Double_t x = -32.4;
        xbins.push_back(x);
        xbins.push_back(x += 10.0);
        for (std::size_t i = 0; i < 56; ++i) {
          xbins.push_back(x += 0.8);
        }
        xbins.push_back(x += 10.0);
      } {
        Double_t y = -20.0;
        ybins.push_back(y);
        ybins.push_back(y += 14.0);
        ybins.push_back(y +=  6.0);
        ybins.push_back(y +=  6.0);
        ybins.push_back(y += 14.0);
      }

      TH2* hist = new TH2D(name.data(),
                           "Extinction Detector;x [cm];y [cm];Entries/cm^{2}",
                           xbins.size() - 1, xbins.data(),
                           ybins.size() - 1, ybins.data());
      hist->SetStats(false);
      return hist;
    }

    inline TList* CreateBorderLine(Color_t color = kBlack, Style_t style = kSolid, Width_t width = 2) {
      TList* list = new TList();
      auto newLine =
        [&](Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
          TLine* line = new TLine(x1, y1, x2, y2);
          line->SetLineColor(color);
          line->SetLineStyle(style);
          line->SetLineWidth(width);
          return line;
        };

      const Double_t xmin = -32.4;
      const Double_t xmax = +32.4;
      const Double_t ymin = -20.0;
      const Double_t ymax = +20.0;

      list->Add(newLine(xmin       , ymin, xmax       , ymin));
      list->Add(newLine(xmin + 10.0, -6.0, xmax - 10.0, -6.0));
      list->Add(newLine(xmin       ,  0.0, xmax       ,  0.0));
      list->Add(newLine(xmin + 10.0, +6.0, xmax - 10.0, +6.0));
      list->Add(newLine(xmin       , ymax, xmax       , ymax));

      list->Add(newLine(xmin, ymin, xmin, ymax));
      list->Add(newLine(xmax, ymin, xmax, ymax));

      {
        Double_t x = -22.4;
        for (std::size_t i = 0; i <= 56; ++i) {
          list->Add(newLine(x, -6.0, x, +6.0));
          x += 0.8;
        }
      }

      {
        Double_t x = -22.4;
        for (std::size_t i = 0; i <= 8; ++i) {
          list->Add(newLine(x, -20.0, x, -6.0));
          list->Add(newLine(x, +20.0, x, +6.0));
          x += 5.6;
        }
      }

      return list;
    }

    inline Long64_t FindChannel(Double_t x, Double_t y) {
      if (TMath::Abs(x) > 8.0 * mm * 28 + 100.0 * mm ||
          TMath::Abs(y) >                 200.0 * mm) {
        // Outside
        return -1;
      } else if (TMath::Abs(x) >  8.0 * mm * 28) {
        // Left/Right
        const Int_t    ix    = (x <  0) ? 0 : 1;
        const Int_t    iy    = (y >= 0) ? 0 : 1;
        const Int_t    ch    = iy * 2 + ix;
        return ch + ChannelOffset::LeftRight + GlobalChannelOffset;
      } else if (TMath::Abs(y) > 60.0 * mm) {
        // Top/Bottom
        const Double_t dx    = (x + 56.0 * mm * 3.5);
        const Int_t    ix    = dx / (56.0 * mm);
        const Int_t    iy    = (y >= 0) ? 0 : 1;
        const Int_t    ch    = iy * 8 + ix;
        return ch + ChannelOffset::TopBottom + GlobalChannelOffset;
      } else {
        // Center
        const Double_t dx    = (x + 8.0 * mm * 28);
        const Int_t    board = dx / (8.0 * mm * 8);
        const Int_t    ix    = (dx - board * 8.0 * mm * 8) / (8.0 * mm);
        const Int_t    iy    = (y >= 0) ? 0 : 1;
        const Int_t    ch    = board * 16 + iy * 8 + ix;
        return ch + ChannelOffset::Center + GlobalChannelOffset;
      }
    }

    inline void Fill(TH2* hist, std::size_t channel) {
      const Int_t    type = GetType(channel);
      const Double_t area = GetArea(type) / cm2;
      switch (type) {
      case ChannelType::Center:
        {
          const std::size_t centerCh = channel - ChannelOffset::Center;
          const std::size_t boardId  = (centerCh / 16);
          const std::size_t boardCh  = (centerCh % 16);
          const std::size_t boardXi  = boardCh % 8;
          const std::size_t boardYi  = boardCh / 8;
          const Double_t    boardX   = 6.4 * (boardId - 3.0) + 0.8 * (boardXi - 3.5);
          const Double_t    boardY   = (boardYi ? -3.0 : +3.0);
          const Int_t       bin      = hist->FindBin(boardX, boardY);
          const Double_t    content  = hist->GetBinContent(bin);
          hist->SetBinContent(bin, content + 1.0 / area);
          hist->SetBinError(bin, TMath::Sqrt((content + 1.0 / area) / area));
        }
        break;

      case ChannelType::TopBottom:
        {
          const std::size_t topBottomCh = channel - ChannelOffset::TopBottom;
          const std::size_t boardXi     = (topBottomCh % 8);
          const std::size_t boardYi     = (topBottomCh / 8);
          const Double_t    boardX      = 5.6 * (boardXi - 3.5);
          const Double_t    boardY      = (boardYi ? -13.0 : +13.0);
          const Int_t       xbin1       = hist->GetXaxis()->FindBin(boardX - 2.7);
          const Int_t       xbin2       = hist->GetXaxis()->FindBin(boardX + 2.7);
          const Int_t       ybin        = hist->GetYaxis()->FindBin(boardY);
          const Double_t    content     = hist->GetBinContent(xbin1, ybin);
          for (Int_t xbin = xbin1; xbin <= xbin2; ++xbin) {
            hist->SetBinContent(xbin, ybin, content + 1.0 / area);
            hist->SetBinError(xbin, ybin, TMath::Sqrt((content + 1.0 / area) / area));
          }
        }
        break;

      case ChannelType::LeftRight:
        {
          const std::size_t leftRightCh = channel - ChannelOffset::LeftRight;
          const std::size_t boardXi     = (leftRightCh % 2);
          const std::size_t boardYi     = (leftRightCh / 2);
          const Double_t    boardX      = (boardXi ? +27.4 : -27.4);
          const Double_t    boardY1     = (boardYi ? - 3.0 : + 3.0);
          const Double_t    boardY2     = (boardYi ? -13.0 : +13.0);
          const Int_t       bin1        = hist->FindBin(boardX, boardY1);
          const Int_t       bin2        = hist->FindBin(boardX, boardY2);
          const Double_t    content     = hist->GetBinContent(bin1);
          hist->SetBinContent(bin1, content + 1.0 / area);
          hist->SetBinContent(bin2, content + 1.0 / area);
          hist->SetBinError(bin1, TMath::Sqrt((content + 1.0 / area) / area));
          hist->SetBinError(bin2, TMath::Sqrt((content + 1.0 / area) / area));
        }
        break;

      }
    }
  }

  namespace Hodoscope {
    constexpr std::size_t NofChannels         = 16U;
    constexpr std::size_t GlobalChannelOffset = ExtinctionDetector::GlobalChannelOffset + ExtinctionDetector::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }

    inline Double_t GetWidth() {
      return 4.0 * cm;
    }

    inline Double_t GetHeight() {
      return 24.0 * cm;
    }

    inline Double_t GetArea() {
      return GetWidth() * GetHeight();
    }

    inline TH2* CreateHitMap(const std::string& name) {
      std::vector<Double_t> xbins, ybins;
      {
        Double_t x = -32.4;
        xbins.push_back(x);
        xbins.push_back(x += 10.0);
        for (std::size_t i = 0; i < 56; ++i) {
          xbins.push_back(x += 0.8);
        }
        xbins.push_back(x += 10.0);
      } {
        Double_t y = -20.0;
        ybins.push_back(y);
        ybins.push_back(y += 14.0);
        ybins.push_back(y +=  6.0);
        ybins.push_back(y +=  6.0);
        ybins.push_back(y += 14.0);
      }

      TH2* hist = new TH2D(name.data(),
                           "Hodoscope;x [cm];y [cm];Entries/cm^{2}",
                           NofChannels, -32.0, +32.0,
                           1, -12.0, 12.0);
      hist->SetStats(false);
      return hist;
    }

    inline TList* CreateBorderLine(Color_t color = kBlack, Style_t style = kSolid, Width_t width = 2) {
      TList* list = new TList();
      auto newLine =
        [&](Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
          TLine* line = new TLine(x1, y1, x2, y2);
          line->SetLineColor(color);
          line->SetLineStyle(style);
          line->SetLineWidth(width);
          return line;
        };

      const Double_t xmin = -32.0;
      const Double_t xmax = +32.0;
      const Double_t ymin = -12.0;
      const Double_t ymax = +12.0;

      list->Add(newLine(xmin, ymin, xmax, ymin));
      list->Add(newLine(xmin, ymax, xmax, ymax));

      {
        Double_t x = xmin;
        for (std::size_t i = 0; i <= 16; ++i) {
          list->Add(newLine(x, ymin, x, ymax));
          x += 4.0;
        }
      }

      return list;
    }

    inline Long64_t FindChannel(Double_t x, Double_t y) {
      if (TMath::Abs(x) > 320.0 * mm ||
          TMath::Abs(y) > 120.0 * mm) {
        // Outside
        return -1;
      } else {
        // Inside
        const Double_t dx = (x + 320.0 * mm);
        const Int_t    ix = dx / (40.0 * mm);
        const Int_t    ch = ix;
        return ch;
      }
    }

    inline void Fill(TH2* hist, std::size_t channel) {
      const Double_t    area     = GetArea() / cm2;
      const Double_t    boardX   = 4.0 * (channel - 7.5);
      const Double_t    boardY   = 0.0;
      const Int_t       bin      = hist->FindBin(boardX, boardY);
      const Double_t    content  = hist->GetBinContent(bin);
      hist->SetBinContent(bin, content + 1.0 / area);
      hist->SetBinError(bin, TMath::Sqrt((content + 1.0 / area) / area));
    }
  }

  namespace TimingCounter {
    constexpr std::size_t NofChannels         = 2U;
    constexpr std::size_t GlobalChannelOffset = Hodoscope::GlobalChannelOffset + Hodoscope::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }
  }

  namespace BeamlineHodoscope {
    constexpr std::size_t NofChannels         = 2U;
    constexpr std::size_t GlobalChannelOffset = TimingCounter::GlobalChannelOffset + TimingCounter::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }
  }

  namespace MrRf {
    constexpr std::size_t NofChannels         = 1U;
    constexpr std::size_t GlobalChannelOffset = BeamlineHodoscope::GlobalChannelOffset + BeamlineHodoscope::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }
  }

  namespace MrP3 {
    constexpr std::size_t NofChannels         = 10U;
    constexpr std::size_t GlobalChannelOffset = MrRf::GlobalChannelOffset + MrRf::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }
  }

  namespace MrSync {
    constexpr std::size_t NofChannels         = 10U;
    constexpr std::size_t GlobalChannelOffset = MrP3::GlobalChannelOffset + MrP3::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }
  }

  namespace EventMatch {
    constexpr std::size_t NofChannels         = 1U;
    constexpr std::size_t GlobalChannelOffset = MrSync::GlobalChannelOffset + MrSync::NofChannels;

    inline Bool_t Contains(std::size_t globalChannel) {
      return GlobalChannelOffset <= globalChannel && globalChannel < GlobalChannelOffset + NofChannels;
    }

    inline Int_t GetChannel(std::size_t globalChannel) {
      if (Contains(globalChannel)) {
        return globalChannel - GlobalChannelOffset;
      }
      return -1;
    }
  }

  namespace GlobalChannel {
    constexpr std::size_t NofChannels = EventMatch::GlobalChannelOffset + EventMatch::NofChannels;
  }
  
}

#endif
