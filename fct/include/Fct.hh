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
#include "Linq.hh"

// #define FCT_FORMAT_VERSION 1 // initial version           ---           ~2021/05/11
#define FCT_FORMAT_VERSION 2    // w/ timestamp after header --- 2021/05/12~

namespace Extinction {

  namespace Fct {

    const std::string Name = "FCT";
    const std::size_t NofChannels = 32;

    namespace ChannelMapWithBoard {
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Ext;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Hod;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Tc;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Bh;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> MrSync;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Evm;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, Int_t>> Veto;
      std::map<Int_t/*global*/, Int_t/*board*/>               Board;

      void Load(const Tron::ConfReader* conf, const std::vector<Int_t>& boards) {
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
                if        (detector == "Ext"   ) {
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
                  Evm   [board][raw] = channel;
                  Board [channel + EventMatch        ::GlobalChannelOffset] = board;
                } else if (detector == "Veto"  ) {
                  Veto  [board][raw] = channel;
                  Board [channel + Veto              ::GlobalChannelOffset] = board;
                }
              }
            }
          }
        }

        auto checkDouplicate =
          [&](const std::string& name, std::map<Int_t, std::map<Int_t, Int_t>>& map, std::size_t nch) {
            for (auto&& board : boards) {
              for (std::size_t ch = 0; ch < nch; ++ch) {
                auto count = Tron::Linq::From(map[board])
                  .Where([&](const std::pair<Int_t, Int_t>& pair) { return pair.second == (Int_t)ch; })
                  .Count();
                if (count > 1) {
                  std::cerr << "[warning] channel duplicate at " << name << "'s ch" << ch << std::endl;
                }
              }
            }
          };
        checkDouplicate("ExtinctionDetector", Ext   , ExtinctionDetector::NofChannels);
        checkDouplicate("Hodoscope"         , Hod   , Hodoscope         ::NofChannels);
        checkDouplicate("BeamlineHodoscope" , Bh    , BeamlineHodoscope ::NofChannels);
        checkDouplicate("TimingCounter"     , Tc    , TimingCounter     ::NofChannels);
        checkDouplicate("MRSync"            , MrSync, MrSync            ::NofChannels);
        checkDouplicate("EventMatch"        , Evm   , EventMatch        ::NofChannels);
        checkDouplicate("Veto"              , Veto  , Veto              ::NofChannels);
      }
    }

    struct DataType {
      enum : UChar_t {
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

    namespace Packet {

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

    }

    class FctData : public ITdcDataProvider {
    public:
      ULong64_t Date;
      Int_t     Spill;
      Int_t     EMCount;
      UChar_t   Type;
      Int_t     Channel;
      Int_t     Tdc;
      Int_t     Carry;
      Int_t     MrSyncCount;
      Int_t     MrSyncTdc;
      Int_t     TdcFromMrSync;

      Double_t  TimePerTdc = 7.5 * nsec;

      Int_t     PreviousCarry[NofChannels];
      Int_t     PreviousTdc  [NofChannels];

      FctData() {
        Clear();
      }
      FctData(const FctData& data)
        : Date         (data.Date         ),
          Spill        (data.Spill        ),
          EMCount      (data.EMCount      ),
          Type         (data.Type         ),
          Channel      (data.Channel      ),
          Tdc          (data.Tdc          ),
          Carry        (data.Carry        ),
          MrSyncCount  (data.MrSyncCount  ),
          MrSyncTdc    (data.MrSyncTdc    ),
          TdcFromMrSync(data.TdcFromMrSync),
          TimePerTdc   (data.TimePerTdc   ) {
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = data.PreviousCarry[ch];
          PreviousTdc  [ch] = data.PreviousTdc  [ch];
        }
      }
      FctData& operator=(const FctData& data) {
        Date          = data.Date;
        Spill         = data.Spill;
        EMCount       = data.EMCount;
        Type          = data.Type;
        Channel       = data.Channel;
        Tdc           = data.Tdc;
        Carry         = data.Carry;
        MrSyncCount   = data.MrSyncCount;
        MrSyncTdc     = data.MrSyncTdc;
        TdcFromMrSync = data.TdcFromMrSync;
        TimePerTdc    = data.TimePerTdc;
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = data.PreviousCarry[ch];
          PreviousTdc  [ch] = data.PreviousTdc  [ch];
        }
        return *this;
      }

      inline void Clear() {
        Date          =  0;
        Spill         = -1;
        EMCount       = -1;
        Type          = DataType::None;
        Channel       =  0;
        Tdc           =  0;
        Carry         =  0;
        MrSyncCount   =  0;
        MrSyncTdc     =  0;
        TdcFromMrSync =  0;
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
        Carry = Packet::GetCarry(data);
        for (std::size_t ch = 0; ch < NofChannels; ++ch) {
          PreviousCarry[ch] = Carry - 1;
        }
      }

      inline void SetDataAsData(Packet_t data) {
        Type      = DataType::Data;
        Channel   = Packet::GetChannel(data);
        Int_t tdc = Packet::GetTdc(data);

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
          if (Tdc < PreviousTdc[Channel]) {
            Tdc = tdc + (PreviousCarry[Channel] = ++Carry) * 0x1000000;
            std::cerr << "[warning] Carry was skipped, "
                      << "Carry " << Carry << ", "
                      << "Tdc " << Form("%x", PreviousTdc[Channel]) << " -> " << Form("%x", Tdc) << std::endl;
          }
        } else {
          Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
          if (Tdc < PreviousTdc[Channel]) {
            PreviousCarry[Channel] = Carry;
            Tdc = tdc + PreviousCarry[Channel] * 0x1000000;
          }
        }
        PreviousTdc[Channel] = Tdc;
      }

      std::basic_istream<char>& Read(std::basic_istream<char>& file,
                                     Packet_t* packet = nullptr) {
        Packet_t buff;
        auto& ret1 = file.read((char*)&buff, sizeof(Packet_t));
        if (!ret1) {
          return ret1;
        }

        if        (Packet::IsHeader(buff)) {
          SetDataAsHeader(buff);
        } else if (Packet::IsGateStart(buff)) {
          SetDataAsGateStart(buff);
#if FCT_FORMAT_VERSION == 1
          // No time stamp
#else
          auto& ret2 = file.read((char*)&Date, sizeof(ULong64_t));
          if (!ret2) {
            return ret2;
          }
#endif
        } else if (Packet::IsGateEnd(buff)) {
          SetDataAsGateEnd(buff);
        } else if (Packet::IsCarry(buff)) {
          SetDataAsCarry(buff);
        } else if (Packet::IsData(buff)) {
          SetDataAsData(buff);
        } else {
          Type = DataType::Error;
        }

        if (packet) {
          *packet = buff;
        }

        return file;
      }

      inline virtual Bool_t IsData() const override {
        return Type == DataType::Data;
      }

      inline virtual Bool_t IsFooter() const override {
        return Type == DataType::GateEnd;
      }

      inline virtual std::string GetName() const override {
        return Name;
      }

      inline virtual void     SetTimePerTdc(Double_t timePerTdc) override {
        TimePerTdc = timePerTdc;
      }
      inline virtual Double_t GetTimePerTdc() const override {
        return TimePerTdc;
      };

      inline virtual ULong64_t GetDate() const override {
        return Date;
      }

      inline virtual Int_t GetSpill() const override {
        return Spill;
      }

      inline virtual Int_t GetEMCount() const override {
        return EMCount;
      }

      inline virtual Double_t GetTime() const override {
        return Tdc * TimePerTdc;
      }

      inline void Show() const {
        std::cout << "Spill = "   << Spill   << ", "
                  << "EMCount = " << EMCount << ", "
                  << "Ch = "      << Channel << ", "
                  << "TDC = "     << Tdc     << std::endl;
      }

      inline virtual std::vector<TdcData> GetTdcData(int board) const override {
        using namespace ChannelMapWithBoard;
        TdcData datum;
        datum.Date            = Date;
        datum.Spill           = Spill;
        datum.EMCount         = EMCount;
        datum.Board           = board;
        datum.TimePerTdc      = TimePerTdc;
        datum.RawChannel      = Channel;
        datum.Tdc             = Tdc;
        datum.Time            = GetTime();
        datum.Tot             = 0;
        datum.LastMrSyncCount = MrSyncCount;
        datum.LastMrSyncTdc   = MrSyncTdc;
        datum.NextMrSyncTdc   = 0;
        datum.TdcFromMrSync   = TdcFromMrSync;
        typename decltype(MrSync      )::const_iterator itr1;
        typename decltype(itr1->second)::const_iterator itr2;
        if        ((itr1 = MrSync      .find(board  )) != MrSync      .end() &&
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
        } else if ((itr1 = Veto        .find(board  )) != Veto        .end() &&
                   (itr2 = itr1->second.find(Channel)) != itr1->second.end()) {
          datum.Channel     = itr2->second + Veto              ::GlobalChannelOffset;
        }
        return { datum };
      }

      inline virtual Int_t FindBoard(Int_t globalChannel) const override {
        std::map<Int_t, Int_t>::const_iterator itr;
        if ((itr = ChannelMapWithBoard::Board.find(globalChannel)) != ChannelMapWithBoard::Board.end()) {
          return itr->second;
        }
        return -1;
      }

      inline void CreateBranch(TTree* tree) {
        // std::cout << "Fct::FctData::CreateBranch()" << std::endl;
        tree->Branch("date"   , &Date         , "date"   "/l");
        tree->Branch("spill"  , &Spill        , "spill"  "/I");
        tree->Branch("type"   , &Type         , "type"   "/b");
        tree->Branch("ch"     , &Channel      , "ch"     "/I");
        tree->Branch("tdc"    , &Tdc          , "tdc"    "/I");
     // tree->Branch("carry"  , &Carry        , "carry"  "/I");
     // tree->Branch("mscount", &MrSyncCount  , "mscount""/I");
     // tree->Branch("mstdc"  , &MrSyncTdc    , "mstdc"  "/I");
        tree->Branch("dtdc"   , &TdcFromMrSync, "dtdc"   "/I");
      }
      inline TBranch* AddEMBranch(TTree* tree) {
        return tree->Branch("emcount", &EMCount, "emcount/I");
      }

      inline virtual void SetBranchAddress(TTree* tree) override {
        // std::cout << "Fct::FctData::SetBranchAddress()" << std::endl;
        tree->SetBranchAddress("date"   , &Date         );
        tree->SetBranchAddress("spill"  , &Spill        );
        tree->SetBranchAddress("emcount", &EMCount      );
        tree->SetBranchAddress("type"   , &Type         );
        tree->SetBranchAddress("ch"     , &Channel      );
        tree->SetBranchAddress("tdc"    , &Tdc          );
     // tree->SetBranchAddress("carry"  , &Carry        );
     // tree->SetBranchAddress("mscount", &MrSyncCount  );
     // tree->SetBranchAddress("mstdc"  , &MrSyncTdc    );
        tree->SetBranchAddress("dtdc"   , &TdcFromMrSync);
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

      virtual Int_t DecodeEventMatchNumber(const std::vector<TdcData>& eventMatchData) const override {
        return TdcData::DecodeEventMatchNumber(eventMatchData);
      }
    };

    class Decoder {
    public:
      FctData Data;
      TTree*  Tree;

      Decoder() : Tree(nullptr) {
      }

      void InitializeTree(const std::string& treename = "tree") {
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
          std::cerr << "[warning] tree has not initialized yet" << std::endl;
        } else {
          for (; Read(file); ++count) {
            // Data.Show();
            if (Data.IsData() || Data.IsFooter()) {
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
