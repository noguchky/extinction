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

#include "ConfReader.hh"
#include "String.hh"

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

    namespace ChannelMapWithBoard {
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Ext;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Hod;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Tc;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Bh;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> MrSync;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Evm;
      std::map<Int_t/*global*/, Int_t/*board*/>               Board;

      void Load(const Tron::ConfReader* conf, const std::vector<int>& boards) {
        for (auto&& board : boards) {
          const std::string key = Form("ChannelMap.%d", board);
          if (conf->Exists(key)) {
            const std::vector<std::string> tuples = conf->GetValues(key);
            for (auto&& tuple : tuples) {
              const std::vector<std::string> elems = Tron::String::Split(tuple, ",");
              if (elems.size() == 3) {
                const Int_t       raw      = Tron::String::Convert<Int_t>(elems[0]);
                const std::string detector =                              elems[1] ;
                const Int_t       channel  = Tron::String::Convert<Int_t>(elems[2]);
                if      (detector == "Ext"   ) {
                  Ext   [board][raw] = channel;
                  Board [channel + ExtinctionDetector::GlobalChannelOffset] = board;
                } else if (detector == "Hod"   ) {
                  Hod   [board][raw] = channel;
                  Board [channel + Hodoscope         ::GlobalChannelOffset] = board;
                } else if (detector == "Bh"    ) {
                  Bh    [board][raw] = channel;
                  Board [channel + BeamlineHodoscope ::GlobalChannelOffset] = board;
                } else if (detector == "Tc"    ) {
                  Tc    [board][raw] = channel;
                  Board [channel + TimingCounter     ::GlobalChannelOffset] = board;
             // } else if (detector == "MrP3"  ) {
             //   MrP3  [board][raw] = channel;
             //   Board [channel + MrP3              ::GlobalChannelOffset] = board;
             // } else if (detector == "MrRf"  ) {
             //   MrRf  [board][raw] = channel;
             //   Board [channel + MrRf              ::GlobalChannelOffset] = board;
                } else if (detector == "MrSync") {
                  MrSync[board][raw] = channel;
                  Board [channel + MrSync            ::GlobalChannelOffset] = board;
                } else if (detector == "Evm"   ) {
                  Evm[board][raw] = channel;
                  Board [channel + EventMatch        ::GlobalChannelOffset] = board;
                }
              }
            }
          }

        }
      }
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
      Int_t    Type;
      Int_t    Spill;
      Int_t    Channel;
      Int_t    DTime;   // ? ... currently, all data are filled with 0
      Int_t    Tdc;
      Int_t    Carry;

      Double_t TimePerTdc = 7.5 * nsec;

      Int_t    PreviousCarry[NofChannels];
      UInt_t   PreviousTdc  [NofChannels];

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

        // Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
        // if ((Long64_t)Tdc < PreviousTdc[Channel]) {
        //   PreviousCarry[Channel] = Carry;
        //   Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
        // }
        if (PreviousTdc[Channel] == 0 || tdc < 0xF00000) {
          Tdc = tdc + (PreviousCarry[Channel] = Carry) * 0x1000000;
        } else {
          Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
          if ((Long64_t)Tdc < PreviousTdc[Channel]) {
            PreviousCarry[Channel] = Carry;
            Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
          }
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

      inline virtual void     SetTimePerTdc(Double_t timePerTdc) override {
        TimePerTdc = timePerTdc;
      }
      inline virtual Double_t GetTimePerTdc() const override {
        return TimePerTdc;
      };
      
      inline virtual Double_t GetTime() const override {
        return Tdc * TimePerTdc;
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
        std::cout << "Spill = " << Spill   << ", " 
                  << "Ch = "    << Channel << ", "
                  << "TDC = "   << Tdc     << ", "
                  << "DTime = " << DTime   << std::endl;
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
        using namespace ChannelMap;
        TdcData datum;
        datum.Spill         = Spill;
        datum.Tdc           = Tdc;
        datum.Time          = GetTime();
        datum.TimePerTdc    = TimePerTdc;
        datum.Board         = 1;
        if (MrSync.begin() != MrSync.end()) {
          datum.MrSyncChannel = MrSync.begin()->second + MrSync::GlobalChannelOffset;
        }
        if        (Ext   .find(Channel) != Ext   .end()) {
          datum.Channel     = Ext   .at(Channel) + ExtinctionDetector::GlobalChannelOffset;
        } else if (Hod   .find(Channel) != Hod   .end()) {
          datum.Channel     = Hod   .at(Channel) + Hodoscope         ::GlobalChannelOffset;
        } else if (Tc    .find(Channel) != Tc    .end()) {
          datum.Channel     = Tc    .at(Channel) + TimingCounter     ::GlobalChannelOffset;
        } else if (Bh    .find(Channel) != Bh    .end()) {
          datum.Channel     = Bh    .at(Channel) + BeamlineHodoscope ::GlobalChannelOffset;
        } else if (MrSync.find(Channel) != MrSync.end()) {
          datum.Channel     = MrSync.at(Channel) + MrSync            ::GlobalChannelOffset;
     // } else if (Evm   .find(Channel) != Evm   .end()) {
     //   datum.Channel     = Evm   .at(Channel) + EventMatch        ::GlobalChannelOffset;
        }
        return { datum };
      }

      inline virtual std::vector<TdcData> GetTdcData(int board) const override {
        using namespace ChannelMapWithBoard;
        TdcData datum;
        datum.Spill         = Spill;
        datum.Tdc           = Tdc;
        datum.Time          = GetTime();
        datum.TimePerTdc    = TimePerTdc;
        datum.Board         = board;
        typename decltype(MrSync      )::const_iterator itr1;
        typename decltype(itr1->second)::const_iterator itr2;
        if        ((itr1 = MrSync      .find(board  )) != Ext         .end() &&
                   (itr2 = itr1->second.begin()      ) != itr1->second.end()) {
          datum.MrSyncChannel = itr2->second + MrSync::GlobalChannelOffset;
        }
        if        ((itr1 = Ext         .find(board  )) != Ext         .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + ExtinctionDetector::GlobalChannelOffset;
        } else if ((itr1 = Hod         .find(board  )) != Hod         .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + Hodoscope         ::GlobalChannelOffset;
        } else if ((itr1 = Tc          .find(board  )) != Tc          .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + TimingCounter     ::GlobalChannelOffset;
        } else if ((itr1 = Bh          .find(board  )) != Bh          .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + BeamlineHodoscope ::GlobalChannelOffset;
        } else if ((itr1 = MrSync      .find(board  )) != MrSync      .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + MrSync            ::GlobalChannelOffset;
        } else if ((itr1 = Evm         .find(board  )) != Evm         .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + EventMatch        ::GlobalChannelOffset;
        }
        return { datum };
      }

      inline void WriteHeader(std::ofstream& file) const {
        Packet_t data = 0x12345678;
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteGateStart(std::ofstream& file) const {
        Packet_t data = 0xffffaaaa;
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteData(std::ofstream& file) const {
        Packet_t data = 0xc0000000U + ((Channel & 0x1fU) << 24) + (Tdc & 0xffffffU);
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteGateEnd(std::ofstream& file) const {
        Packet_t data = 0xffff5555;
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteCarry(std::ofstream& file) const {
        Packet_t data = 0xffff0000 + (Carry & 0xff);
        file.write((Char_t*)&data, sizeof(Packet_t));
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
