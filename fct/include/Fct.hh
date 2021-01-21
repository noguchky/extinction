#ifndef Extinction_Fct_hh
#define Extinction_Fct_hh

#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "Units.hh"
#include "Tdc.hh"
#include "Detector.hh"

namespace Extinction {

  namespace Fct {

    const std::string Name = "FCT";

    namespace ChannelMap {
      // (*) key should be board no & raw channel
      const std::map<Int_t, Int_t> Ext
        {
         { 31, 48 +  1 - 1 },
         { 30, 48 +  3 - 1 },
         { 29, 48 +  4 - 1 },
         { 28, 48 +  6 - 1 },
         { 27, 48 +  7 - 1 },
         { 26, 48 + 12 - 1 },
         { 25, 48 +  8 - 1 },
         { 24, 48 +  2 - 1 },
         { 23, 48 +  9 - 1 },
         { 22, 48 + 10 - 1 },
         { 21, 48 + 11 - 1 },
         { 20, 48 + 16 - 1 },
         { 19, 48 +  5 - 1 },
         { 18, 48 + 13 - 1 },
         { 17, 48 + 14 - 1 },
         { 16, 48 + 15 - 1 },
        };
      const std::map<Int_t, Int_t> Hod
        {
         {  9, 11 },
         {  8, 10 },
         {  7,  9 },
         {  6,  8 },
         {  5,  7 },
         {  4,  6 },
         {  3,  5 },
         {  2,  4 },
         { 10,  0 }, // TBD: Hod OR
        };
      const std::map<Int_t, Int_t> Tc
        {
         // TC1
         // TC2
        };
      const std::map<Int_t, Int_t> Bh
        {
         { 12,  0 }, // BH1
         { 11,  1 }, // BH2
        };
      const std::map<Int_t, Int_t> MrSync
        {
         { 14, 0 },
        };
    }
 
    class FctData : public ITdcDataProvider {
    public:
      Int_t Spill;
      Int_t Channel;
      Int_t DTime;   // ? ... currently, all data are filled with 0
      Int_t Tdc;

      FctData() {
        Clear();
      }
      FctData(const FctData& data)
        : Spill    (data.Spill    ),
          Channel  (data.Channel  ),
          DTime    (data.DTime    ),
          Tdc      (data.Tdc      ) {
      }
      FctData& operator=(const FctData& data) {
        Spill     = data.Spill;
        Channel   = data.Channel;
        DTime     = data.DTime;
        Tdc       = data.Tdc;
        return *this;
      }

      inline void Clear() {
        Spill     = -1;
        Channel   = 0;
        DTime     = 0;
        Tdc       = 0;
      }

      inline virtual Double_t GetTimePerTdc() const override {
        return 7.5 * nsec;
      };
      
      inline virtual Double_t GetTime() const override {
        return Tdc * GetTimePerTdc();
      }

      inline void CreateBranch(TTree* tree) {
        // std::cout << "Fct::FctData::CreateBranch()" << std::endl;
        tree->Branch("CH"   , &Channel, "CH" "/I");
        tree->Branch("TDC"  , &Tdc    , "TDC""/I");
        tree->Branch("DTIME", &DTime  , "DTIME/I");
        tree->Branch("SPILL", &Spill  , "SPILL/I");
      }

      inline virtual void SetBranchAddress(TTree* tree) override {
        // std::cout << "Fct::FctData::SetBranchAddress()" << std::endl;
        tree->SetBranchAddress("CH"   , &Channel);
        tree->SetBranchAddress("TDC"  , &Tdc    );
        tree->SetBranchAddress("DTIME", &DTime  );
        tree->SetBranchAddress("SPILL", &Spill  );
      }

      inline void Show() const {
        std::cout << "Spill = " << Spill  << ", " 
                  << "Ch = " << Channel  << ", "
                  << "TDC = " << Tdc  << ", "
                  << "DTime = " << DTime << std::endl;
      }

      inline virtual std::string GetName() const override {
        return Name;
      }

      inline virtual Bool_t IsData() const override {
        return true;
      }

      inline virtual Int_t GetSpill() const override {
        return Spill;
      }

      inline virtual std::vector<TdcData> GetTdcData() const override {
        TdcData datum;
        datum.Spill         = Spill;
        datum.Tdc           = Tdc;
        datum.Time          = GetTime();
        if (ChannelMap::MrSync.find(14) != ChannelMap::MrSync.end()) {
          datum.MrSyncChannel = ChannelMap::MrSync.at(14) + MrSync::GlobalChannelOffset; // TBD
        }
        if        (ChannelMap::Ext   .find(Channel) != ChannelMap::Ext   .end()) {
          datum.Channel     = ChannelMap::Ext   .at(Channel) + ExtinctionDetector::GlobalChannelOffset;
        } else if (ChannelMap::Hod   .find(Channel) != ChannelMap::Hod   .end()) {
          datum.Channel     = ChannelMap::Hod   .at(Channel) + Hodoscope         ::GlobalChannelOffset;
        } else if (ChannelMap::Tc    .find(Channel) != ChannelMap::Tc    .end()) {
          datum.Channel     = ChannelMap::Tc    .at(Channel) + TimingCounter     ::GlobalChannelOffset;
        } else if (ChannelMap::Bh    .find(Channel) != ChannelMap::Bh    .end()) {
          datum.Channel     = ChannelMap::Bh    .at(Channel) + BeamlineHodoscope ::GlobalChannelOffset;
        } else if (ChannelMap::MrSync.find(Channel) != ChannelMap::MrSync.end()) {
          datum.Channel     = ChannelMap::MrSync.at(Channel) + MrSync            ::GlobalChannelOffset;
        }
        return { datum };
      }

    };

  }

}

#endif
