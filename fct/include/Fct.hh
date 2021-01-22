#ifndef Extinction_Fct_hh
#define Extinction_Fct_hh

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "TTree.h"
#include "Units.hh"
#include "Tdc.hh"
#include "Detector.hh"

namespace Extinction {

  namespace Fct {

    const std::string Name = "FCT";
    const std::size_t NofChannels = 32;

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

    struct DataType {
      enum {
            None,
            Header,
            GateStart,
            GateEnd,
            Carry, // ~ Heartbeat/Overflow
            Data,
            Error,
      };
    };
    
    // 4 byte
    using Packet_t = UInt_t;

    inline void ShowAsHex(Packet_t data) {
      printf("%08x \n", data);
    }

    inline Bool_t IsHeader(Packet_t data) {
      // return (data & 0xffffffff) == 0x12345678;
      return data == 0x12345678;
    }
    inline Bool_t IsGateStart(Packet_t data) {
      // return (data & 0xff0000ff) == 0xff0000aa;
      return data  == 0xffffaaaa;
    }
    inline Bool_t IsGateEnd(Packet_t data) {
      return (data & 0xff0000ff) == 0xff000055;
    }
    inline Bool_t IsCarry(Packet_t data) {
      return (data & 0xffff0000) == 0xffff0000;
    }
    inline Bool_t IsData(Packet_t data) {
      return
        ((data & 0xf0000000) == 0xc0000000) ||
        ((data & 0xf0000000) == 0xd0000000);
    }

    inline Int_t GetCarry(Packet_t data) {
      return data & 0xff;
    }
    inline Int_t GetChannel(Packet_t data) {
      return (data >> 24) & 0x1f;
    }
    inline Int_t GetTdc(Packet_t data) {
      return data & 0x00ffffff;
    }
    
    class FctData : public ITdcDataProvider {
    public:
      Int_t Type;
      Int_t Spill;
      Int_t Channel;
      Int_t DTime;   // ? ... currently, all data are filled with 0
      Int_t Tdc;
      Int_t Carry;

      Int_t  PreviousCarry[NofChannels];
      UInt_t PreviousTdc  [NofChannels];

      FctData() {
        Clear();
      }
      FctData(const FctData& data)
        : Type     (data.Type     ),
          Spill    (data.Spill    ),
          Channel  (data.Channel  ),
          DTime    (data.DTime    ),
          Tdc      (data.Tdc      ),
          Carry    (data.Carry    ) {
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = data.PreviousCarry[ch];
          PreviousTdc  [ch] = data.PreviousTdc  [ch];
        }
      }
      FctData& operator=(const FctData& data) {
        Type      = data.Type;
        Spill     = data.Spill;
        Channel   = data.Channel;
        DTime     = data.DTime;
        Tdc       = data.Tdc;
        Carry     = data.Carry;
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = data.PreviousCarry[ch];
          PreviousTdc  [ch] = data.PreviousTdc  [ch];
        }
        return *this;
      }

      inline void Clear() {
        Type      = DataType::None;
        Spill     = -1;
        Channel   = 0;
        DTime     = 0;
        Tdc       = 0;
        Carry     = 0;
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = 0;
          PreviousTdc  [ch] = 0;
        }
      }

      inline void SetDataAsHeader(Packet_t) {
        Type = DataType::Header;
      }

      inline void SetDataAsGateStart(Packet_t) {
        Type = DataType::GateStart;
        ++Spill;
        Channel = 0;
        DTime   = 0;
        Tdc     = 0;
        Carry   = 0;
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = 0;
          PreviousTdc  [ch] = 0;
        }
      }

      inline void SetDataAsGateEnd(Packet_t) {
        Type = DataType::GateEnd;
      }

      inline void SetDataAsCarry(Packet_t data) {
        Type  = DataType::Carry;
        Carry = GetCarry(data);
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = Carry - 1;
        }
      }

      inline void SetDataAsData(Packet_t data) {
        Type      = DataType::Data;
        Channel   = GetChannel(data);
        Int_t tdc = GetTdc(data);

        if (Channel < 0 || (Long64_t)NofChannels <= Channel) {
          Type = DataType::Error;
          return;
        }

        Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
        if ((Long64_t)Tdc < PreviousTdc[Channel]) {
          PreviousCarry[Channel] = Carry;
          Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
        }
        PreviousTdc[Channel] = Tdc;
      }

      std::basic_istream<char>& Read(std::basic_istream<char>& file,
                                     Packet_t* packet = nullptr) {
        Packet_t buff;
        auto& ret = file.read((char*)&buff, sizeof(Packet_t));
        if (!ret) {
          return ret;
        }

        if        (IsHeader(buff)) {
          SetDataAsHeader(buff);
        } else if (IsGateStart(buff)) {
          SetDataAsGateStart(buff);
        } else if (IsGateEnd(buff)) {
          SetDataAsGateEnd(buff);
        } else if (IsCarry(buff)) {
          SetDataAsCarry(buff);
        } else if (Fct::IsData(buff)) {
          SetDataAsData(buff);
        } else {
          Type = DataType::Error;
        }

        if (packet) {
          *packet = buff;
        }

        return ret;
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
        Type = DataType::Data;
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
        return Type == DataType::Data;
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

    class Decoder {
    public:
      FctData Data;
      TTree*  Tree;

      Decoder() : Tree(nullptr) {
      }

      void InitializeTree(const std::string& treename = "FCTTDC") {
        // std::cout << "Fct::Decoder::InitializeTree()" << std::endl;
        if (Tree) {
          delete Tree;
        }
        Tree = new TTree(treename.data(), treename.data());

        Data.CreateBranch(Tree);
      }

      std::size_t Decode(std::basic_istream<char>& file) {
        // std::cout << "Fct::Decoder::Decode()" << std::endl;
        std::size_t count = 0UL;
        if (!Tree) {
          std::cout << "[warning] tree has not initialized yet" << std::endl;
        } else {
          for (; Read(file); ++count) {
            // Data.Show();
            if (Data.Type == DataType::Data) {
              // std::cout << "TTree::Fill()" << std::endl;
              Tree->Fill();
            }
          }
        }
        return count;
      }

      std::size_t Decode(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
          std::cerr << "[error] file is not opened, " << filename << std::endl;
          return 0UL;
        }
        return Decode(file);
      }

      inline std::basic_istream<char>& Read(std::basic_istream<char>& file,
                                            Packet_t* packet = nullptr) {
        return Data.Read(file, packet);
      }
    };
    
  }

}

#endif
