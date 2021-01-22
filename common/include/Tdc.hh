#ifndef Extinction_Tdc_hh
#define Extinction_Tdc_hh

#include <vector>
#include "TTree.h"

namespace Extinction {

  class TdcData {
  public:
    Int_t    Spill;
    Int_t    Channel;
    Long64_t Tdc;
    Double_t Time;
    Int_t    MrSyncChannel;

    TdcData() {
      Clear();
    }
    TdcData(const TdcData& data)
      : Spill        (data.Spill        ),
        Channel      (data.Channel      ),
        Tdc          (data.Tdc          ),
        Time         (data.Time         ),
        MrSyncChannel(data.MrSyncChannel) {
    }
    virtual ~TdcData() {
    }

    TdcData& operator=(const TdcData& data) {
      Spill         = data.Spill;
      Channel       = data.Channel;
      Tdc           = data.Tdc;
      Time          = data.Time;
      MrSyncChannel = data.MrSyncChannel;
      return *this;
    }

    inline void Clear() {
      Spill         = -1;
      Channel       = -1;
      Tdc           = 0;
      Time          = 0;
      MrSyncChannel = -2;
    }

  };

  class ITdcDataProvider {
  public:
    virtual std::string          GetName() const = 0;
    virtual Double_t             GetTimePerTdc() const = 0;
    virtual void                 SetBranchAddress(TTree*) = 0;
    virtual Bool_t               IsData() const = 0;
    virtual Int_t                GetSpill() const = 0;
    virtual Double_t             GetTime() const = 0;
    virtual std::vector<TdcData> GetTdcData() const = 0;
  };

}

#endif
