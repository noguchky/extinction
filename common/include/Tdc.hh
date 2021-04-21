#ifndef Extinction_Tdc_hh
#define Extinction_Tdc_hh

#include <vector>
#include "TTree.h"
#include "Units.hh"

namespace Extinction {

  class TdcData {
  public:
    Int_t     Spill;
    Int_t     EMCount;
    Int_t     Channel;
    Long64_t  Tdc;
    Double_t  Time;
    Int_t     MrSyncChannel;
    Double_t  TimePerTdc;
    Int_t     Board;

    TdcData() {
      Clear();
    }
    TdcData(const TdcData& data)
      : Spill        (data.Spill        ),
        EMCount      (data.EMCount      ),
        Channel      (data.Channel      ),
        Tdc          (data.Tdc          ),
        Time         (data.Time         ),
        MrSyncChannel(data.MrSyncChannel),
        TimePerTdc   (data.TimePerTdc   ),
        Board        (data.Board        ) {
    }
    virtual ~TdcData() {
    }

    TdcData& operator=(const TdcData& data) {
      Spill         = data.Spill;
      EMCount       = data.EMCount;
      Channel       = data.Channel;
      Tdc           = data.Tdc;
      Time          = data.Time;
      MrSyncChannel = data.MrSyncChannel;
      TimePerTdc    = data.TimePerTdc;
      Board         = data.Board;
      return *this;
    }

    inline void Clear() {
      Spill         = -1;
      EMCount       = -1;
      Channel       = -1;
      Tdc           = 0;
      Time          = 0;
      MrSyncChannel = -2;
      TimePerTdc    = 1.0;
      Board         = 0;
    }

    inline ULong64_t GetTdcTag() const {
      // return (Tdc < 0 || Channel < 0) ? 0 : (Tdc * 1000ULL + Channel);
      return (Tdc < 0 || Channel < 0) ? 0 : (((ULong64_t)(Time / nsec)) * 1000ULL + Channel);
    }

    inline void CreateBranch(TTree* tree) {
      // std::cout << "TdcData::CreateBranch()" << std::endl;
      tree->Branch("spill"        , &Spill        , "spill"        "/I");
      tree->Branch("emcount"      , &EMCount      , "emcount"      "/I");
      tree->Branch("ch"           , &Channel      , "ch"           "/I");
      tree->Branch("tdc"          , &Tdc          , "tdc"          "/L");
      tree->Branch("time"         , &Time         , "time"         "/D");
      tree->Branch("mrSyncChannel", &MrSyncChannel, "mrSyncChannel""/I");
      tree->Branch("timePerTDc"   , &TimePerTdc   , "timePerTDc"   "/D");
      tree->Branch("board"        , &Board        , "board"        "/I");
    }

    inline void SetBranchAddress(TTree* tree) {
      // std::cout << "TdcData::SetBranchAddress()" << std::endl;
      tree->SetBranchAddress("spill"        , &Spill        );
      tree->SetBranchAddress("emcount"      , &EMCount      );
      tree->SetBranchAddress("ch"           , &Channel      );
      tree->SetBranchAddress("tdc"          , &Tdc          );
      tree->SetBranchAddress("time"         , &Time         );
      tree->SetBranchAddress("mrSyncChannel", &MrSyncChannel);
      tree->SetBranchAddress("timePerTDc"   , &TimePerTdc   );
      tree->SetBranchAddress("board"        , &Board        );
    }

  };

  class ITdcDataProvider {
  public:
    virtual ~ITdcDataProvider() { }
    virtual std::string          GetName() const = 0;
    virtual void                 SetTimePerTdc(Double_t) = 0;
    virtual Double_t             GetTimePerTdc() const = 0;
    virtual void                 SetBranchAddress(TTree*) = 0;
    virtual Bool_t               IsData() const = 0;
    virtual Bool_t               IsFooter() const = 0;
    virtual Int_t                GetSpill() const = 0;
    virtual Int_t                GetEMCount() const = 0;
    virtual Double_t             GetTime() const = 0;
    virtual std::vector<TdcData> GetTdcData() const = 0;
    virtual std::vector<TdcData> GetTdcData(Int_t) const = 0;
    virtual Int_t                DecodeEventMatchNumber(const std::vector<TdcData>&) const = 0;
  };

}

#endif
