#ifndef Extinction_Tdc_hh
#define Extinction_Tdc_hh

#include <vector>
#include "TTree.h"
#include "TMath.h"
#include "Units.hh"

namespace Extinction {

  using Tag_t = std::tuple<Int_t/*MrSyncCount*/, Long64_t/*TdcFromMrSync*/, Int_t/*GlobalChannel*/>;

  class TdcData {
  public:
    // Spill information
    ULong64_t Date;
    Int_t     Spill;
    Int_t     EMCount;

    // Module information
    Int_t     Board;
    Double_t  TimePerTdc;

    // Event information
    Int_t     RawChannel;
    Int_t     Channel;
    Long64_t  Tdc;
    Double_t  Time;
    UChar_t   Tot;
    Int_t     MrSyncChannel;
    Long64_t  LastMrSyncCount;
    Long64_t  LastMrSyncTdc;
    Long64_t  NextMrSyncTdc;
    Long64_t  TdcFromMrSync;

    TdcData() {
      Clear();
    }
    TdcData(const TdcData& data)
      : Date           (data.Date           ),
        Spill          (data.Spill          ),
        EMCount        (data.EMCount        ),
        Board          (data.Board          ),
        TimePerTdc     (data.TimePerTdc     ),
        RawChannel     (data.RawChannel     ),
        Channel        (data.Channel        ),
        Tdc            (data.Tdc            ),
        Time           (data.Time           ),
        Tot            (data.Tot            ),
        MrSyncChannel  (data.MrSyncChannel  ),
        LastMrSyncCount(data.LastMrSyncCount),
        LastMrSyncTdc  (data.LastMrSyncTdc  ),
        NextMrSyncTdc  (data.NextMrSyncTdc  ),
        TdcFromMrSync  (data.TdcFromMrSync  ) {
    }
    virtual ~TdcData() {
    }

    TdcData& operator=(const TdcData& data) {
      Date            = data.Date;
      Spill           = data.Spill;
      EMCount         = data.EMCount;
      Board           = data.Board;
      TimePerTdc      = data.TimePerTdc;
      RawChannel      = data.RawChannel;
      Channel         = data.Channel;
      Tdc             = data.Tdc;
      Time            = data.Time;
      Tot             = data.Tot;
      MrSyncChannel   = data.MrSyncChannel;
      LastMrSyncCount = data.LastMrSyncCount;
      LastMrSyncTdc   = data.LastMrSyncTdc;
      NextMrSyncTdc   = data.NextMrSyncTdc;
      TdcFromMrSync   = data.TdcFromMrSync;
      return *this;
    }

    inline void Clear() {
      Date            =  0;
      Spill           = -1;
      EMCount         = -1;
      Board           = -1;
      TimePerTdc      =  1.0;
      RawChannel      = -1;
      Channel         = -1;
      Tdc             =  0;
      Time            =  0;
      Tot             =  0;
      MrSyncChannel   = -2;
      LastMrSyncCount =  0;
      LastMrSyncTdc   =  0;
      NextMrSyncTdc   =  0;
      TdcFromMrSync   =  0;
    }

    inline ULong64_t GetTdcTag(std::size_t mrcount, Long64_t mrtdc) const {
      // return (Tdc < 0 || Channel < 0) ? 0 : (Tdc * 1000ULL + Channel);
      return (Tdc < 0 || Channel < 0) ? 0 : ((mrcount * 10000000000ULL + (Tdc - mrtdc)) * 1000ULL + Channel);
    }

    inline void CreateBranch(TTree* tree) {
      // std::cout << "TdcData::CreateBranch()" << std::endl;
      tree->Branch("date"         , &Date           , "date"         "/l");
      tree->Branch("spill"        , &Spill          , "spill"        "/I");
      tree->Branch("emcount"      , &EMCount        , "emcount"      "/I");
      tree->Branch("board"        , &Board          , "board"        "/I");
      tree->Branch("timePerTdc"   , &TimePerTdc     , "timePerTdc"   "/D");
      tree->Branch("rawch"        , &RawChannel     , "rawch"        "/I");
      tree->Branch("ch"           , &Channel        , "ch"           "/I");
      tree->Branch("tdc"          , &Tdc            , "tdc"          "/L");
      tree->Branch("time"         , &Time           , "time"         "/D");
      tree->Branch("tot"          , &Tot            , "tot"          "/b");
      tree->Branch("msch"         , &MrSyncChannel  , "msch"         "/I");
      tree->Branch("mscount"      , &LastMrSyncCount, "mscount"      "/L");
      tree->Branch("mstdc1"       , &LastMrSyncTdc  , "mstdc1"       "/L");
      tree->Branch("mstdc2"       , &NextMrSyncTdc  , "mstdc2"       "/L");
      tree->Branch("dtdc"         , &TdcFromMrSync  , "dtdc"         "/L");
    }

    inline void SetBranchAddress(TTree* tree) {
      // std::cout << "TdcData::SetBranchAddress()" << std::endl;
      tree->SetBranchAddress("date"         , &Date           );
      tree->SetBranchAddress("spill"        , &Spill          );
      tree->SetBranchAddress("emcount"      , &EMCount        );
      tree->SetBranchAddress("board"        , &Board          );
      tree->SetBranchAddress("timePerTDc"   , &TimePerTdc     );
      tree->SetBranchAddress("rawch"        , &RawChannel     );
      tree->SetBranchAddress("ch"           , &Channel        );
      tree->SetBranchAddress("tdc"          , &Tdc            );
      tree->SetBranchAddress("time"         , &Time           );
      tree->SetBranchAddress("tot"          , &Tot            );
      tree->SetBranchAddress("mrch"         , &MrSyncChannel  );
      tree->SetBranchAddress("mrcount"      , &LastMrSyncCount);
      tree->SetBranchAddress("mrtdc1"       , &LastMrSyncTdc  );
      tree->SetBranchAddress("mrtdc2"       , &NextMrSyncTdc  );
      tree->SetBranchAddress("dtdc"         , &TdcFromMrSync  );
    }

    static Double_t GetTimeDifference(const TdcData& data1, const TdcData& data2) {
      const Int_t dmscount = data2.LastMrSyncCount - data1.LastMrSyncCount;
      Long64_t tdc1 = data1.TdcFromMrSync;
      Long64_t tdc2 = data2.TdcFromMrSync;
      if        (dmscount > 0) {
        if (data1.NextMrSyncTdc && data1.LastMrSyncTdc) {
          const Long64_t msinterval1 = data1.NextMrSyncTdc - data1.LastMrSyncTdc;
          tdc1 -= (+dmscount) * msinterval1;
        } else {
          tdc1 = data1.Tdc;
          tdc2 = data2.Tdc;
        }
      } else if (dmscount < 0) {
        if (data2.NextMrSyncTdc && data2.LastMrSyncTdc) {
          const Long64_t msinterval2 = data2.NextMrSyncTdc - data2.LastMrSyncTdc;
          tdc2 -= (-dmscount) * msinterval2;
        } else {
          tdc1 = data1.Tdc;
          tdc2 = data2.Tdc;
        }
      }
      const Double_t time1 = tdc1 * data1.TimePerTdc;
      const Double_t time2 = tdc2 * data2.TimePerTdc;
      return time1 - time2;
    }

    static Int_t DecodeEventMatchNumber(const std::vector<TdcData>& eventMatchData) {
      const std::size_t kHeaderSize     =  2;
      const std::size_t kEventMatchSize = 19;
      Int_t eventMatchNumber = -1;
      if (eventMatchData.size() < kHeaderSize) {
        std::cerr << "[warning] event match data is empty" << std::endl;
      } else {
        Int_t eventMatchBits[kEventMatchSize] = { 0 };
        // for (auto i : eventMatchBits) { std::cout << i; } std::cout << std::endl;
        // for (auto&& data : eventMatchData) { std::cout << data.Tdc << " "; } std::cout << std::endl;
        for (auto&& data : eventMatchData) {
          const Double_t norm = (Double_t)(data.Tdc - eventMatchData[0].Tdc) / (Double_t)(eventMatchData[1].Tdc - eventMatchData[0].Tdc);
          const std::size_t bit = TMath::Nint(norm);
          if (bit < kEventMatchSize) {
            // std::cout << bit << ", ";
            // std::cout << bit << ": " << data.Tdc << std::endl;
            eventMatchBits[bit] = 1;
          } else {
            // Maybe second event match signals
            break;
            // std::cerr << "[warning] invalid event match tdc" << std::endl;
            // std::cout << "(" << data.Tdc << " - " << eventMatchData[0].Tdc << ") / "
            //           << "(" << eventMatchData[1].Tdc << " - " << eventMatchData[0].Tdc << ")" << std::endl
            //           << "-> " << data.Tdc - eventMatchData[0].Tdc
            //           << " / " << eventMatchData[1].Tdc - eventMatchData[0].Tdc << std::endl
            //           << "-> " << norm << " -> " << bit << std::endl;
          }
        }
        // std::cout << std::endl;
        // for (auto i : eventMatchBits) { std::cout << i; } std::cout << std::endl;
        Int_t parityBit = 0;
        for (std::size_t bit = 0; bit < kEventMatchSize - 1; ++bit) {
          parityBit ^= eventMatchBits[bit];
        }
        // std::cout << parityBit << " <--> " << eventMatchBits[kEventMatchSize - 1] << std::endl;
        if (parityBit != eventMatchBits[kEventMatchSize - 1]) {
          std::cerr << "[warning] invalid parity bit" << std::endl;
          for (auto&& data : eventMatchData) {
            const Double_t norm = (Double_t)(data.Tdc - eventMatchData[0].Tdc) / (Double_t)(eventMatchData[1].Tdc - eventMatchData[0].Tdc);
            const std::size_t bit = TMath::Nint(norm);
            std::cerr << bit << ": " << data.Tdc << std::endl;
          }
        } else {
          eventMatchNumber = 0;
          for (std::size_t bit = kHeaderSize; bit < kEventMatchSize - 1; ++bit) {
            if (eventMatchBits[bit]) {
              eventMatchNumber += (0x1 << (bit - kHeaderSize));
            }
          }
        }
      }
      return eventMatchNumber;
    }
    
  };

  class ITdcDataProvider {
  public:
    virtual ~ITdcDataProvider() { }
    virtual Bool_t               IsData() const = 0;
    virtual Bool_t               IsFooter() const = 0;
    virtual std::string          GetName() const = 0;
    virtual void                 SetTimePerTdc(Double_t) = 0;
    virtual Double_t             GetTimePerTdc() const = 0;
    virtual ULong64_t            GetDate() const = 0;
    virtual Int_t                GetSpill() const = 0;
    virtual Int_t                GetEMCount() const = 0;
    virtual Double_t             GetTime() const = 0;
    virtual std::vector<TdcData> GetTdcData(Int_t) const = 0;
    virtual Int_t                FindBoard(Int_t) const = 0;
    virtual void                 SetBranchAddress(TTree*) = 0;
    virtual Int_t                DecodeEventMatchNumber(const std::vector<TdcData>&) const = 0;
  };

}

#endif
