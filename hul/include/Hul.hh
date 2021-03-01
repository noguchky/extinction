#ifndef Extinction_Hul_hh
#define Extinction_Hul_hh

#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "Units.hh"
#include "Tdc.hh"
#include "Detector.hh"

#include "ConfReader.hh"
#include "String.hh"

namespace Extinction {

  namespace Hul {

    const std::string Name = "HUL";

    namespace ChannelMap {
      // (*) key should be board no & raw channel
      const std::map<Int_t, Int_t> Ext
        {
         { 16 - 16 - 1, 48 +  1 - 1 },
         { 17 - 16 - 1, 48 +  3 - 1 },
         { 18 - 16 - 1, 48 +  4 - 1 },
         { 19 - 16 - 1, 48 +  6 - 1 },
         { 20 - 16 - 1, 48 +  7 - 1 },
         { 21 - 16 - 1, 48 + 12 - 1 },
         { 22 - 16 - 1, 48 +  8 - 1 },
         { 23 - 16 - 1, 48 +  2 - 1 },
         { 24 - 16 - 1, 48 +  9 - 1 },
         { 25 - 16 - 1, 48 + 10 - 1 },
         { 26 - 16 - 1, 48 + 11 - 1 },
         { 27 - 16 - 1, 48 + 16 - 1 },
         { 28 - 16 - 1, 48 +  5 - 1 },
         { 29 - 16 - 1, 48 + 13 - 1 },
         { 30 - 16 - 1, 48 + 14 - 1 },
         { 31 - 16 - 1, 48 + 15 - 1 },
        };
      const std::map<Int_t, Int_t> Hod
        {
         {  6 + 16 - 1, 11 },
         {  7 + 16 - 1, 10 },
         {  8 + 16 - 1,  9 },
         {  9 + 16 - 1,  8 },
         { 10 + 16 - 1,  7 },
         { 11 + 16 - 1,  6 },
         { 12 + 16 - 1,  5 },
         { 14 + 16 - 1,  4 },
         {  5 + 16 - 1,  0 }, // TBD: Hod OR
        };
      const std::map<Int_t, Int_t> Tc
        {
         // TC1
         // TC2
        };
      const std::map<Int_t, Int_t> Bh
        {
         {  3 + 16 - 1,  0 }, // BH1
         {  4 + 16 - 1,  1 }, // BH2
        };
      const std::map<Int_t, Int_t> MrSync
        {
         {  1 + 16 - 1, 0 },
        };
    }

    namespace ChannelMapWithBoard {
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Ext;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Hod;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Tc;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Bh;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> MrSync;
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
                }
              }
            }
          }
            
        }
      }
    }

    struct DataType {
      enum {
            None       = 0x0U,
            SpillStart = 0x1U,
            SpillEnd   = 0x4U,
            Data       = 0xDU,
            Error      = 0xEU,
            Heartbeat  = 0xFU,
      };
    };

    // 5 byte
    using Packet1_t = UInt_t;
    using Packet2_t = UChar_t;

    struct Packet_t {
      Packet1_t Packet1;
      Packet2_t Packet2;
    };

    inline void ShowAsHex(Packet1_t data1, Packet2_t data2) {
      printf("%08x %02x \n", data1, data2);
    }

    inline UInt_t GetHeader(Packet1_t data) {
      return (data & 0xFE000000U) >> 25;
    }
    inline UInt_t GetChannel(Packet1_t data) {
      return (data & 0x01F80000U) >> 19;
    }
    inline UInt_t GetTdc(Packet1_t data) {
      return (data & 0x0007FFFFU);
    }
    inline UInt_t GetHeartbeat(Packet1_t data) {
      return (data & 0x0000FFFFU);
    }
    inline UChar_t GetType(Packet2_t data) {
      return (data & 0xF0U) >> 4;
    }
    inline UChar_t GetFooter(Packet2_t data) {
      return (data & 0x0FU);
    }

    class HulData : public ITdcDataProvider {
    public:
      UChar_t  Type;
      Int_t    Spill;     // log_2(60 * 60 * 24) = 16.39 -> need more than 16 bit
      UShort_t Channel;
      UInt_t   Tdc;       // 19 bit
      UShort_t Heartbeat; // 16 bit
      // TDC (19 bit) + Heartbeat (16 bit) = 35 bit
      // 1.04 GHz sampling -> max 33 sec
      Double_t ClockFreq;

      HulData()
        : ClockFreq(1.04 * GHz) {
        Clear();
      }
      HulData(const HulData& data)
        : Type     (data.Type     ),
          Spill    (data.Spill    ),
          Channel  (data.Channel  ),
          Tdc      (data.Tdc      ),
          Heartbeat(data.Heartbeat),
          ClockFreq(data.ClockFreq) {
      }
      HulData& operator=(const HulData& data) {
        Type      = data.Type;
        Spill     = data.Spill;
        Channel   = data.Channel;
        Tdc       = data.Tdc;
        Heartbeat = data.Heartbeat;
        ClockFreq = data.ClockFreq;
        return *this;
      }

      inline void Clear() {
        Type      = DataType::None;
        Spill     = -1;
        Channel   = 0;
        Tdc       = 0;
        Heartbeat = 0;
      }

      std::basic_istream<char>& Read(std::basic_istream<char>& file,
                                     Packet_t* packet = nullptr) {
        Type    = DataType::None;
        Channel = 0;
        Tdc     = 0;

        Packet1_t buff1;
        std::basic_istream<char>& ret1 = file.read((char*)&buff1, sizeof(Packet1_t));
        if (!ret1) {
          return ret1;
        }

        Packet2_t buff2;
        std::basic_istream<char>& ret2 = file.read((char*)&buff2, sizeof(Packet2_t));
        if (!ret2) {
          return ret2;
        }

        Type = GetType(buff2);
        switch (Type) {
        case DataType::SpillStart:
          ++Spill;
          Heartbeat = 0;
          break;
        case DataType::SpillEnd:
          break;
        case DataType::Data:
          Channel = GetChannel(buff1);
          Tdc     = GetTdc(buff1);
          break;
        case DataType::Error:
          Heartbeat = GetHeartbeat(buff1);
          std::cerr << "error detected, spill = " << Spill
                    << ", heartbeat = " << Heartbeat << std::endl;
          break;
        case DataType::Heartbeat:
          Heartbeat = GetHeartbeat(buff1);
          break;
        }

        if (packet) {
          packet->Packet1 = buff1;
          packet->Packet2 = buff2;
        }

        return ret2;
      }

      inline Long64_t GetTdc2() const {
        return Tdc + 0x80000UL * Heartbeat;
      }

      inline virtual void SetTimePerTdc(Double_t timePerTdc) override {
        ClockFreq = 1.0 / timePerTdc;
      }

      inline virtual Double_t GetTimePerTdc() const override {
        return 1.0 / ClockFreq;
      }

      inline virtual Double_t GetTime() const override {
        return GetTdc2() * GetTimePerTdc();
      }

      inline void CreateBranch(TTree* tree) {
        // std::cout << "Hul::HulData::CreateBranch()" << std::endl;
        tree->Branch("type"     , &Type     , "type"   "/b");
        tree->Branch("spill"    , &Spill    , "spill"  "/I");
        tree->Branch("ch"       , &Channel  , "ch"     "/s");
        tree->Branch("tdc"      , &Tdc      , "tdc"    "/i");
        tree->Branch("heartbeat", &Heartbeat, "heartbeat/s");
        tree->SetAlias("tdc2", "tdc + 0x80000 * heartbeat");
      }

      inline virtual void SetBranchAddress(TTree* tree) override {
        // std::cout << "Hul::HulData::SetBranchAddress()" << std::endl;
        tree->SetBranchAddress("type"     , &Type     );
        tree->SetBranchAddress("spill"    , &Spill    );
        tree->SetBranchAddress("ch"       , &Channel  );
        tree->SetBranchAddress("tdc"      , &Tdc      );
        tree->SetBranchAddress("heartbeat", &Heartbeat);
      }

      inline void Show() const {
        std::string tname =
          Type == DataType::None       ? "None"       :
          Type == DataType::SpillStart ? "SpillStart" :
          Type == DataType::SpillEnd   ? "SpillEnd"   :
          Type == DataType::Data       ? "Data"       :
          Type == DataType::Error      ? "Error"      :
          Type == DataType::Heartbeat  ? "Heartbeat"  : "Other";
        std::cout << "Spill = " << Spill  << ", " 
                  << (UInt_t)Type << "(" << tname << ") " << ", "
                  << "Ch = " << Channel  << ", "
                  << "TDC = " << Tdc << " + 0x8000 * " << Heartbeat << std::endl;
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
        datum.TimePerTdc    = GetTimePerTdc();
        datum.Board         = 1;
        if (ChannelMap::MrSync.begin() != ChannelMap::MrSync.end()) {
          datum.MrSyncChannel = ChannelMap::MrSync.begin()->second + MrSync::GlobalChannelOffset;
        }
        if        (ChannelMap::Ext.find(Channel) != ChannelMap::Ext.end()) {
          datum.Channel     = ChannelMap::Ext.at(Channel) + ExtinctionDetector::GlobalChannelOffset;
        } else if (ChannelMap::Hod.find(Channel) != ChannelMap::Hod.end()) {
          datum.Channel     = ChannelMap::Hod.at(Channel) + Hodoscope::GlobalChannelOffset;
        } else if (ChannelMap::Tc.find(Channel) != ChannelMap::Tc.end()) {
          datum.Channel     = ChannelMap::Tc.at(Channel) + TimingCounter::GlobalChannelOffset;
        } else if (ChannelMap::Bh.find(Channel) != ChannelMap::Bh.end()) {
          datum.Channel     = ChannelMap::Bh.at(Channel) + BeamlineHodoscope::GlobalChannelOffset;
        } else if (ChannelMap::MrSync.find(Channel) != ChannelMap::MrSync.end()) {
          datum.Channel     = ChannelMap::MrSync.at(Channel) + MrSync::GlobalChannelOffset;
        }
        return { datum };
      }

      inline virtual std::vector<TdcData> GetTdcData(Int_t board) const override {
        using namespace ChannelMapWithBoard;
        TdcData datum;
        datum.Spill         = Spill;
        datum.Tdc           = Tdc;
        datum.Time          = GetTime();
        datum.TimePerTdc    = GetTimePerTdc();
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
        }
        return { datum };
      }

      inline void WriteSpillStart(std::ofstream& file) const {
        Packet1_t data1 = 0x0000ffff;
        Packet2_t data2 = (DataType::SpillStart << 4);
        file.write((Char_t*)&data1, sizeof(Packet1_t));
        file.write((Char_t*)&data2, sizeof(Packet2_t));
      }
      inline void WriteData(std::ofstream& file) const {
        Packet1_t data1 = ((Channel & 0x3F) << 19) + (Tdc & 0x7FFFF);
        Packet2_t data2 = (DataType::Data << 4);
        file.write((Char_t*)&data1, sizeof(Packet1_t));
        file.write((Char_t*)&data2, sizeof(Packet2_t));
      }
      inline void WriteSpillEnd(std::ofstream& file) const {
        Packet1_t data1 = 0x0;
        Packet2_t data2 = (DataType::SpillEnd << 4);
        file.write((Char_t*)&data1, sizeof(Packet1_t));
        file.write((Char_t*)&data2, sizeof(Packet2_t));
      }
      inline void WriteError(std::ofstream& file) const {
        Packet1_t data1 = (Heartbeat & 0xFFFF);
        Packet2_t data2 = (DataType::Error << 4);
        file.write((Char_t*)&data1, sizeof(Packet1_t));
        file.write((Char_t*)&data2, sizeof(Packet2_t));
      }
      inline void WriteHeartbeat(std::ofstream& file) const {
        Packet1_t data1 = (Heartbeat & 0xFFFF);
        Packet2_t data2 = (DataType::Heartbeat << 4);
        file.write((Char_t*)&data1, sizeof(Packet1_t));
        file.write((Char_t*)&data2, sizeof(Packet2_t));
      }

    };

    class Decoder {
    public:
      HulData Data;
      TTree*  Tree;

      Decoder() : Tree(nullptr) {
      }

      void InitializeTree(const std::string& treename = "tree") {
        // std::cout << "Hul::Decoder::InitializeTree()" << std::endl;
        if (Tree) {
          delete Tree;
        }
        Tree = new TTree(treename.data(), treename.data());

        Data.CreateBranch(Tree);
      }

      std::size_t Decode(std::basic_istream<char>& file) {
        // std::cout << "Hul::Decoder::Decode()" << std::endl;
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
