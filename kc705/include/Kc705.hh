#ifndef Extinction_Kc705_hh
#define Extinction_Kc705_hh

#include <iostream>
#include <fstream>
#include <cstring>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "Units.hh"
#include "Inet.hh"
#include "Tdc.hh"
#include "Detector.hh"

#include "ConfReader.hh"
#include "String.hh"

#define KC705_FORMAT_VERSION 2

namespace Extinction {

  namespace Kc705 {

    const std::string Name = "KC705";

    namespace ChannelMap {
      // (*) key should be board no & raw channel
      const std::map<Int_t/*MppcCh*/, Int_t> ExtMppc
        {
         { 32, 49 - 1 },
         { 33, 51 - 1 },
         { 42, 52 - 1 },
         { 43, 54 - 1 },
         { 52, 55 - 1 },
         { 53, 60 - 1 },
         { 54, 56 - 1 },
         { 55, 50 - 1 },
         { 56, 57 - 1 },
         { 57, 58 - 1 },
         { 58, 59 - 1 },
         { 59, 64 - 1 },
         { 60, 53 - 1 },
         { 61, 61 - 1 },
         { 62, 62 - 1 },
         { 63, 63 - 1 },
        };
      const std::map<Int_t/*SubCh*/, Int_t> ExtSub
        {
         // LeftRight (ch128-131)
        };
      const std::map<Int_t/*SubCh*/, Int_t> Hod
        {
         {  3, 11 },
         {  4, 10 },
         {  5,  9 },
         {  6,  8 },
         {  7,  7 },
         {  8,  6 },
         {  9,  5 },
         { 10,  4 },
         {  2,  0 }, // TBD: Hod OR
        };
      const std::map<Int_t/*SubCh*/, Int_t> Tc
        {
         // TC1
         // TC2
        };
      const std::map<Int_t/*SubCh*/, Int_t> Bh
        {
         { 0, 0 }, // BH1
         { 1, 1 }, // BH2
        };
      const std::map<Int_t/*BoardId*/, Int_t> MrSync
        {
         { 0, 0 },
         { 1, 1 },
        };
    }

    namespace ChannelMapWithBoard {
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, std::set<Int_t>>> ExtMppc;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, std::set<Int_t>>> ExtSub;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, std::set<Int_t>>> Hod;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, std::set<Int_t>>> Tc;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, std::set<Int_t>>> Bh;
      std::map<Int_t/*board*/, std::map<Int_t/*raw*/, std::set<Int_t>>> MrSync;
      std::map<Int_t/*global*/, std::set<Int_t/*board*/>> Board;

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
                if        (detector == "ExtMppc") {
                  ExtMppc[board][raw].insert(channel);
                  Board  [channel + ExtinctionDetector::GlobalChannelOffset].insert(board);
                } else if (detector == "ExtSub" ) {
                  ExtSub [board][raw].insert(channel);
                  Board  [channel + ExtinctionDetector::GlobalChannelOffset].insert(board);
                } else if (detector == "Hod"    ) {
                  Hod    [board][raw].insert(channel);
                  Board  [channel + Hodoscope         ::GlobalChannelOffset].insert(board);
                } else if (detector == "Bh"     ) {
                  Bh     [board][raw].insert(channel);
                  Board  [channel + BeamlineHodoscope ::GlobalChannelOffset].insert(board);
                } else if (detector == "Tc"     ) {
                  Tc     [board][raw].insert(channel);
                  Board  [channel + TimingCounter     ::GlobalChannelOffset].insert(board);
             // } else if (detector == "MrP3"   ) {
             //   MrP3   [board][raw].insert(channel);
             //   Board  [channel + MrP3              ::GlobalChannelOffset].insert(board);
             // } else if (detector == "MrRf"   ) {
             //   MrRf   [board][raw].insert(channel);
             //   Board  [channel + MrRf              ::GlobalChannelOffset].insert(board);
                } else if (detector == "MrSync" ) {
                  MrSync [board][raw].insert(channel);
                  Board  [channel + MrSync            ::GlobalChannelOffset].insert(board);
                }
              }
            }
          }
        }

        // Nothing to check
      }
    }
    
    struct DataType {
      enum {
            None,
            Header,
            Data,
            Footer,
            HeaderError,
            PacketLoss = -1,
      };
    };

    const std::size_t MppcNch = 64UL;
    const std::size_t SubNch  = 12UL;

    using Packet_t = UChar_t[13];
    inline void ShowAsHex(Packet_t data) {
      for (std::size_t i = 0, n = sizeof(Packet_t); i < n; ++i) {
        printf("%02x ", data[i]);
      }
      puts("");
    }
    inline Bool_t HasBit(Packet_t data, std::size_t bit) {
      std::size_t byte = bit / 8;
      return data[sizeof(Packet_t) - 1 - byte] & (1U << (bit % 8));
    }
    inline std::size_t GetBitCount(UChar_t data) {
      std::size_t bitCount = 0UL;
      for (std::size_t i = 0UL; i < 8UL; ++i) {
        if (data & (1U << i)) {
          ++bitCount;
        }
      }
      return bitCount;
    }

#if KC705_FORMAT_VERSION == 1

    /// Header
    // Format : Id[83:0], BoardId[3:0], Spill[15:0]
    //          II II II II II II II II II II IB SS SS
    //                                       ^*F
    //                                          ^FF-FF
    //     Id : AB B0 00 12 34 56 70 12 34 56 7_ __ __
    //         ^FF-FF-F*-**
    //         ^**-**-*F-FF-FF-FF-F*-**
    //                     ^**-**-*F-FF-FF-FF-F*-**
    inline UInt_t GetHeaderId1(Packet_t data) {
      return pntohl(data) >> 12;
    }
    inline UInt_t GetHeaderId2(Packet_t data) {
      return (pntohll(data    ) >> 12) & 0xFFFFFFFFU;
    }
    inline UChar_t GetHeaderId3(Packet_t data) {
      return (pntohll(data + 4) >> 12) & 0xFFFFFFFFU;
    }
    inline UChar_t GetBoardId(Packet_t data) {
      return pntohc(data + 10) & 0x0FU;
    }
    inline UShort_t GetSpill(Packet_t data) {
      return pntohs(data + 11);
    }
    inline Bool_t IsHeader(Packet_t data) {
      // return
      //   GetHeaderId1(data) ==    0xABB00 &&
      //   GetHeaderId2(data) == 0x01234567 &&
      //   GetHeaderId3(data) == 0x01234567;
      return
        (data[ 0]        ) == 0xABU &&
        (data[ 1]        ) == 0xB0U &&
        (data[ 2]        ) == 0x00U &&
        (data[ 3]        ) == 0x12U &&
        (data[ 4]        ) == 0x34U &&
        (data[ 5]        ) == 0x56U &&
        (data[ 6]        ) == 0x70U &&
        (data[ 7]        ) == 0x12U &&
        (data[ 8]        ) == 0x34U &&
        (data[ 9]        ) == 0x56U &&
        (data[10] & 0xF0U) == 0x70U;
    }
    /// Data
    // Format : MPPC[63:0], Sub[11:0], MR_SYNC, TDC[26:0]
    //          MM MM MM MM MM MM MM MM SS SM TT TT TT
    //         ^FF-FF-FF-FF-FF-FF-FF-FF
    //                                 ^FF-F*
    //                                    ^*8
    //                                    ^*7-FF-FF-FF
    inline ULong64_t GetMppcBit(Packet_t data) {
      return pntohll(data);
    }
    inline UShort_t GetSubBit(Packet_t data) {
      return pntohs(data + 8) >> 4;
    }
    inline Bool_t GetMrSync(Packet_t data) {
      return (pntohc(data + 9) >> 3) & 0x01U;
    }
    inline UInt_t GetTdc(Packet_t data) {
      return pntohl(data + 9) & 0x07FFFFFFU;
    }
    /// Footer
    // Format : EMCount[15:0], ????[23:0], Id[63:0]
    //          EE EE ?? ?? ?? II II II II II II II II
    //         ^FF-FF         ^FF-FF-FF-FF-FF-FF-FF-FF
    //     Id :                AA-AA-AA-AA-AA-AA-AA-AA
    inline UShort_t GetEMCount(Packet_t data) {
      return pntohs(data);
    }
    inline UShort_t GetWRCount(Packet_t) {
      return 0;
    }
    inline ULong64_t GetFooterId(Packet_t data) {
      return pntohll(data + 5);
    }
    inline Bool_t IsFooter(Packet_t data) {
      // return GetFooterId(data) == 0xAAAAAAAAAAAAAAAAULL;
      return
        data[ 5] == 0xAAU &&
        data[ 6] == 0xAAU &&
        data[ 7] == 0xAAU &&
        data[ 8] == 0xAAU &&
        data[ 9] == 0xAAU &&
        data[10] == 0xAAU &&
        data[11] == 0xAAU &&
        data[12] == 0xAAU;
    }

#else

    /// Header
    // Format : Id1[31:0], Spill[15:0], Zero[3:0], BoardId[3:0], Id2[47:0]
    //          II II II II SS SS 0B II II II II II II
    //         ^FF-FF-FF-FF
    //                     ^FF-FF
    //                           ^FF
    //                        ^**-**-FF-FF-FF-FF-FF-FF
    //     Id : 01-23-45-67          01-23-45-67-89-AB
    inline UInt_t GetHeaderId1(Packet_t data) {
      return pntohl(data);
    }
    inline UShort_t GetSpill(Packet_t data) {
      return pntohs(data + 4);
    }
    inline UChar_t GetBoardId(Packet_t data) {
      return pntohc(data + 6) & 0x0FU;
    }
    inline ULong64_t GetHeaderId2(Packet_t data) {
      return pntohll(data + 5) & 0x00FFFFFFU;
    }
    inline Bool_t IsHeader(Packet_t data) {
      return GetHeaderId1(data) == 0x01234567U;
    }
    /// Data
    // Format : MPPC[63:0], Sub[11:0], MR_SYNC, TDC[26:0]
    //          MM MM MM MM MM MM MM MM SS SM TT TT TT
    //         ^FF-FF-FF-FF-FF-FF-FF-FF
    //                                 ^FF-F*
    //                                    ^*8
    //                                    ^*7-FF-FF-FF
    inline ULong64_t GetMppcBit(Packet_t data) {
      return pntohll(data);
    }
    inline UShort_t GetSubBit(Packet_t data) {
      return pntohs(data + 8) >> 4;
    }
    inline Bool_t GetMrSync(Packet_t data) {
      return (pntohc(data + 9) >> 3) & 0x01U;
    }
    inline UInt_t GetTdc(Packet_t data) {
      return pntohl(data + 9) & 0x07FFFFFFU;
    }
    /// Footer
    // Format : EMCount[15:0], ????[23:0], Id[63:0]
    //          Id1[31:0], Spill[15:0], EMCount[15:0], WrCnt[31:0], Id2[1:0]
    //          II II II II SS SS EE EE WW WW WW WW II
    //         ^FF-FF-FF-FF
    //                     ^FF-FF
    //                           ^FF-FF
    //                                 ^FF-FF-FF-FF
    //                                             ^FF
    //     Id : AA-AA-AA-AA                         AB
    inline UInt_t GetFooterId1(Packet_t data) {
      return pntohl(data);
    }
    // inline UShort_t GetSpill(Packet_t data) {
    //   return pntohs(data + 4);
    // } // the same as header
    inline UShort_t GetEMCount(Packet_t data) {
      return pntohs(data + 6);
    }
    inline UShort_t GetWRCount(Packet_t data) {
      return pntohl(data + 8);
    }
    inline UInt_t GetFooterId2(Packet_t data) {
      return pntohc(data + 12);
    }
    inline Bool_t IsFooter(Packet_t data) {
      return GetFooterId1(data) == 0xAAAAAAAAULL;
    }

#endif

    class Kc705Data : public ITdcDataProvider {
    public:
      Char_t    Type;

      // Header
      UChar_t   BoardId;  //  2 bit
      Int_t     Spill;    //  8 bit

      // Data
      ULong64_t MppcBit;  // 64 bit: MPPCs of extinction detector
      UShort_t  SubBit;   // 12 bit: PMTs of detectors and MR RF
      Bool_t    MrSync;   //  1 bit
      UInt_t    Tdc;      // 27 bit

      // Footer
      UShort_t  EMCount;  // 16 bit
      UShort_t  WRCount;  // 16 bit

      // Additional
      UInt_t    Overflow; // TrueTdc = Tdc + 0x8000000 * Overflow
      Double_t  TimePerTdc = 5.0 * nsec;
      UChar_t   Mppcs[MppcNch];
      UChar_t   Subs[SubNch];

      Kc705Data() {
        Clear();
      }
      Kc705Data(const Kc705Data& data)
        : Type    (data.Type    ),
          BoardId (data.BoardId ),
          Spill   (data.Spill   ),
          MppcBit (data.MppcBit ),
          SubBit  (data.SubBit  ),
          MrSync  (data.MrSync  ),
          Tdc     (data.Tdc     ),
          EMCount (data.EMCount ),
          WRCount (data.WRCount ),
          Overflow(data.Overflow) {
        std::memcpy(Mppcs, data.Mppcs, MppcNch);
        std::memcpy( Subs, data. Subs,  SubNch);
      }
      Kc705Data& operator=(const Kc705Data& data) {
        Type     = data.Type;
        BoardId  = data.BoardId;
        Spill    = data.Spill;
        MppcBit  = data.MppcBit;
        SubBit   = data.SubBit;
        MrSync   = data.MrSync;
        Tdc      = data.Tdc;
        EMCount  = data.EMCount;
        WRCount  = data.WRCount;
        Overflow = data.Overflow;
        std::memcpy(Mppcs, data.Mppcs, MppcNch);
        std::memcpy( Subs, data. Subs,  SubNch);
        return *this;
      }

      inline void Clear() {
        Type     = DataType::None;
        BoardId  = 0;
        Spill    = -1;
        MppcBit  = 0;
        SubBit   = 0;
        MrSync   = 0;
        Tdc      = 0;
        EMCount  = 0;
        WRCount  = 0;
        Overflow = 0;
        std::memset(Mppcs, 0, MppcNch);
        std::memset( Subs, 0,  SubNch);
      }

#if KC705_FORMAT_VERSION == 1

      inline void SetDataAsHeader(Packet_t packet) {
        Type     = DataType::Header;
        BoardId  = Kc705::GetBoardId(packet);
        Spill    = Kc705::GetSpill(packet);
        MppcBit  = 0;
        SubBit   = 0;
        MrSync   = 0;
        Tdc      = 0;
        EMCount  = 0;
        WRCount  = 0;
        Overflow = 0;
        std::memset(Mppcs, 0, MppcNch);
        std::memset( Subs, 0,  SubNch);
      }

      inline void SetDataAsData(Packet_t packet) {
        const UInt_t lastTdc = Tdc;
        Type    = DataType::Data;
        MppcBit = Kc705::GetMppcBit(packet);
        SubBit  = Kc705::GetSubBit(packet);
        MrSync  = Kc705::GetMrSync(packet);
        Tdc     = Kc705::GetTdc(packet);
        if (Tdc < lastTdc) { ++Overflow; }
        for (std::size_t ch = 0; ch < MppcNch; ++ch) {
          Mppcs[ch] = IsMppcHit(ch);
        }
        for (std::size_t ch = 0; ch < SubNch; ++ch) {
          Subs[ch] = IsSubHit(ch);
        }
      }

      inline void SetDataAsFooter(Packet_t packet) {
        Type    = DataType::Footer;
        MppcBit = 0;
        SubBit  = 0;
        MrSync  = 0;
        Tdc     = 0;
        EMCount = Kc705::GetEMCount(packet);
        WRCount = Kc705::GetWRCount(packet);
        std::memset(Mppcs, 0, MppcNch);
        std::memset( Subs, 0,  SubNch);
      }

      inline void GetData(Packet_t& packet) const {
        packet[ 0] = ((MppcBit >> 56) & 0xFFU);
        packet[ 1] = ((MppcBit >> 48) & 0xFFU);
        packet[ 2] = ((MppcBit >> 40) & 0xFFU);
        packet[ 3] = ((MppcBit >> 32) & 0xFFU);
        packet[ 4] = ((MppcBit >> 24) & 0xFFU);
        packet[ 5] = ((MppcBit >> 16) & 0xFFU);
        packet[ 6] = ((MppcBit >>  8) & 0xFFU);
        packet[ 7] = ((MppcBit >>  0) & 0xFFU);
        packet[ 8] = ((SubBit  >>  4) & 0xFFU);
        packet[ 9] = ((SubBit  <<  4) & 0xF0U) + (MrSync * 0x08U) + ((Tdc >> 24) & 0x07U);
        packet[10] = ((Tdc     >> 16) & 0xFFU);
        packet[11] = ((Tdc     >>  8) & 0xFFU);
        packet[12] = ((Tdc     >>  0) & 0xFFU);
      }

#else

      inline void SetDataAsHeader(Packet_t packet) {
        Type     = DataType::Header;
        BoardId  = Kc705::GetBoardId(packet);
        Spill    = Kc705::GetSpill(packet);
        MppcBit  = 0;
        SubBit   = 0;
        MrSync   = 0;
        Tdc      = 0;
        EMCount  = 0;
        WRCount  = 0;
        Overflow = 0;
        std::memset(Mppcs, 0, MppcNch);
        std::memset( Subs, 0,  SubNch);
      }

      inline void SetDataAsData(Packet_t packet) {
        const UInt_t lastTdc = Tdc;
        Type    = DataType::Data;
        MppcBit = Kc705::GetMppcBit(packet);
        SubBit  = Kc705::GetSubBit(packet);
        MrSync  = Kc705::GetMrSync(packet);
        Tdc     = Kc705::GetTdc(packet);
        if (Tdc < lastTdc) { ++Overflow; }
        for (std::size_t ch = 0; ch < MppcNch; ++ch) {
          Mppcs[ch] = IsMppcHit(ch);
        }
        for (std::size_t ch = 0; ch < SubNch; ++ch) {
          Subs[ch] = IsSubHit(ch);
        }
      }

      inline void SetDataAsFooter(Packet_t packet) {
        Type    = DataType::Footer;
        const Int_t lastSpill = Spill;
        Spill   = Kc705::GetSpill(packet);
        if (Spill != lastSpill) {
          std::cerr << "[warning] spill inconsistent" << std::endl;
        }
        MppcBit = 0;
        SubBit  = 0;
        MrSync  = 0;
        Tdc     = 0;
        EMCount = Kc705::GetEMCount(packet);
        WRCount = Kc705::GetWRCount(packet);
        std::memset(Mppcs, 0, MppcNch);
        std::memset( Subs, 0,  SubNch);
      }

      inline void GetData(Packet_t& packet) const {
        packet[ 0] = ((MppcBit >> 56) & 0xFFU);
        packet[ 1] = ((MppcBit >> 48) & 0xFFU);
        packet[ 2] = ((MppcBit >> 40) & 0xFFU);
        packet[ 3] = ((MppcBit >> 32) & 0xFFU);
        packet[ 4] = ((MppcBit >> 24) & 0xFFU);
        packet[ 5] = ((MppcBit >> 16) & 0xFFU);
        packet[ 6] = ((MppcBit >>  8) & 0xFFU);
        packet[ 7] = ((MppcBit >>  0) & 0xFFU);
        packet[ 8] = ((SubBit  >>  4) & 0xFFU);
        packet[ 9] = ((SubBit  <<  4) & 0xF0U) + (MrSync * 0x08U) + ((Tdc >> 24) & 0x07U);
        packet[10] = ((Tdc     >> 16) & 0xFFU);
        packet[11] = ((Tdc     >>  8) & 0xFFU);
        packet[12] = ((Tdc     >>  0) & 0xFFU);
      }

#endif

      std::basic_istream<char>& ReadHeader(std::basic_istream<char>& file,
                                           Packet_t* packet = nullptr) {
        Packet_t buff;
        auto& ret = file.read((char*)buff, sizeof(Packet_t));
        if (!ret) {
          return ret;
        }

        if (Kc705::IsHeader(buff)) {
          SetDataAsHeader(buff);
        } else if (Type == DataType::HeaderError) {
          // Nothing to do
        } else {
          std::cout << "[error] invalid header, maybe "
                    << (Kc705::IsFooter(buff) ? "footer" : "data") << std::endl;
          Type = DataType::HeaderError;
        }

        if (packet) {
          std::memcpy(*packet, buff, sizeof(Packet_t));
        }

        return ret;
      }

      std::basic_istream<char>& ReadDataOrFooter(std::basic_istream<char>& file,
                                                 Packet_t* packet = nullptr) {
        Packet_t buff;
        auto& ret = file.read((char*)buff, sizeof(Packet_t));
        if (!ret) {
          return ret;
        }

        if (Kc705::IsFooter(buff)) {
          SetDataAsFooter(buff);
        } else {
          SetDataAsData(buff);
        }

        if (packet) {
          std::memcpy(*packet, buff, sizeof(Packet_t));
        }

        return ret;
      }

      inline Bool_t IsMppcHit(Int_t ch) const {
        return MppcBit & (0x1ULL << ch);
      }
      inline std::vector<Int_t> GetMppcHitChannels() const {
        std::vector<Int_t> channels;
        for (std::size_t ch = 0; ch < MppcNch; ++ch) {
          if (IsMppcHit(ch)) {
            channels.push_back(ch);
          }
        }
        return channels;
      }

      inline Bool_t IsSubHit(Int_t ch) const {
        return SubBit & (0x1U << ch);
      }
      inline std::vector<Int_t> GetSubHitChannels() const {
        std::vector<Int_t> channels;
        for (std::size_t ch = 0; ch < SubNch; ++ch) {
          if (IsSubHit(ch)) {
            channels.push_back(ch);
          }
        }
        return channels;
      }

      inline Long64_t GetTdc2() const {
        return Tdc + 0x8000000ULL * Overflow;
      }

      inline virtual void SetTimePerTdc(Double_t timePerTdc) override {
        TimePerTdc = timePerTdc;
      }

      inline virtual Double_t GetTimePerTdc() const override {
        return TimePerTdc;
      }

      inline virtual Double_t GetTime() const override {
        return GetTdc2() * GetTimePerTdc();
      }

      inline void CreateBranch(TTree* tree) {
        // std::cout << "Kc705::Kc705Data::CreateBranch()" << std::endl;
        tree->Branch("type"    , &Type    , "type"    "/B");
        tree->Branch("boardId" , &BoardId , "boardId" "/b");
        tree->Branch("spill"   , &Spill   , "spillC"  "/I");
        tree->Branch("tdc"     , &Tdc     , "tdc"     "/i");
        tree->Branch("mrSync"  , &MrSync  , "MrSync"  "/O");
        tree->Branch("mppc"    , &MppcBit , "mppc"    "/l");
        tree->Branch("sub"     , &SubBit  , "sub"     "/s");
        tree->Branch("overflow", &Overflow, "overflow""/i");
        tree->Branch("mppcs"   ,  Mppcs   , Form("mppcs[%lu]" "/b", MppcNch));
        tree->Branch("subs"    ,  Subs    , Form("subs [%lu]" "/b", SubNch));
        tree->SetAlias("tdc2", "tdc + 0x8000000 * overflow");
      }
      inline TBranch* AddEMBranch(TTree* tree) {
        return tree->Branch("emcount", &EMCount, "emcount/I");
      }

      inline virtual void SetBranchAddress(TTree* tree) override {
        // std::cout << "Kc705::Kc705Data::SetBranchAddress()" << std::endl;
        tree->SetBranchAddress("type"    , &Type    );
        tree->SetBranchAddress("boardId" , &BoardId );
        tree->SetBranchAddress("spill"   , &Spill   );
        tree->SetBranchAddress("emcount" , &EMCount );
        tree->SetBranchAddress("tdc"     , &Tdc     );
        tree->SetBranchAddress("mrSync"  , &MrSync  );
        tree->SetBranchAddress("mppc"    , &MppcBit );
        tree->SetBranchAddress("sub"     , &SubBit  );
        tree->SetBranchAddress("overflow", &Overflow);
        tree->SetBranchAddress("mppcs"   ,  Mppcs   );
        tree->SetBranchAddress("subs"    ,  Subs    );
      }

      inline void ShowAsHex() const {
        Packet_t data;
        GetData(data);
        for (std::size_t i = 0, n = sizeof(Packet_t); i < n; ++i) {
          printf("%02x ", data[i]);
        }
        printf(" (%02d) \n", (Int_t)Type);
      }
      inline void Show() const {
        std::cout << "Type = " << (Int_t)Type << ", "
                  << "Board = " << (UInt_t)BoardId  << ", " 
                  << "Spill = " << Spill  << ", "
                  << "EMCount = " << EMCount << ", "
                  << "Tdc = " << Tdc << ", "
                  << "Overflow = " << Overflow << ", "
                  << "MPPC =";
        std::vector<Int_t> mppcHit = GetMppcHitChannels();
        if (mppcHit.empty()) {
          std::cout << " None";
        } else {
          for (std::size_t i = 0, n = mppcHit.size(); i < n; ++i) {
            std::cout << " " << mppcHit[i];
          }
        }
        std::cout << ", Sub =";
        std::vector<Int_t> subHit = GetSubHitChannels();
        if (subHit.empty()) {
          std::cout << " None";
        } else {
          for (std::size_t i = 0, n = subHit.size(); i < n; ++i) {
            std::cout << " " << subHit[i];
          }
        }
        std::cout << ", MrSync = " << (MrSync ? "true" : "false") << std::endl;
      }

      inline virtual std::string GetName() const override {
        return Name;
      }

      inline virtual Bool_t IsData() const override {
        return Type == DataType::Data;
      }

      inline virtual Bool_t IsFooter() const override {
        return Type == DataType::Footer;
      }
      
      inline virtual Int_t GetSpill() const override {
        return Spill;
      }

      inline virtual Int_t GetEMCount() const override {
        return EMCount;
      }

      inline virtual std::vector<TdcData> GetTdcData() const override {
        using namespace ChannelMap;
        namespace CM = ChannelMap;
        std::vector<TdcData> data;
        const Long64_t tdc  = GetTdc2();
        const Double_t time = GetTime();
        for (auto&& mppcCh : GetMppcHitChannels()) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if (CM::MrSync.find(BoardId) != CM::MrSync.end()) {
            datum.MrSyncChannel = CM::MrSync.at(BoardId) + MrSync            ::GlobalChannelOffset;
          }
          if (ExtMppc.find(mppcCh) != ExtMppc.end()) {
            datum.Channel       = ExtMppc   .at(mppcCh ) + ExtinctionDetector::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        for (auto&& subCh : GetSubHitChannels()) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if (CM::MrSync.find(BoardId) != CM::MrSync.end()) {
            datum.MrSyncChannel = CM::MrSync.at(BoardId) + MrSync            ::GlobalChannelOffset;
          }
          if        (ExtSub.find(subCh) != ExtSub.end()) {
            datum.Channel       = ExtSub    .at(subCh  ) + ExtinctionDetector::GlobalChannelOffset;
          } else if (Hod   .find(subCh) != Hod   .end()) {
            datum.Channel       = Hod       .at(subCh  ) + Hodoscope         ::GlobalChannelOffset;
          } else if (Tc    .find(subCh) != Tc    .end()) {
            datum.Channel       = Tc        .at(subCh  ) + TimingCounter     ::GlobalChannelOffset;
          } else if (Bh    .find(subCh) != Bh    .end()) {
            datum.Channel       = Bh        .at(subCh  ) + BeamlineHodoscope ::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        if (MrSync) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if (CM::MrSync.find(BoardId) != CM::MrSync.end()) {
            datum.Channel       = CM::MrSync.at(BoardId) + MrSync            ::GlobalChannelOffset;
            datum.MrSyncChannel = datum.Channel;
          }
          data.push_back(datum);
        }
        return data;
      }

      inline virtual std::vector<TdcData> GetTdcData(Int_t board) const override {
        using namespace ChannelMapWithBoard;
        namespace CM = ChannelMapWithBoard;
        std::vector<TdcData> data;
        const Long64_t tdc  = GetTdc2();
        const Double_t time = GetTime();
        typename decltype(CM::MrSync  )::const_iterator itr1;
        typename decltype(itr1->second)::const_iterator itr2;
        for (auto&& mppcCh : GetMppcHitChannels()) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if        ((itr1 = CM::MrSync  .find(board )) != CM::MrSync  .end() &&
                     (itr2 = itr1->second.begin()     ) != itr1->second.end()) {
            datum.MrSyncChannel = *itr2->second.begin() + MrSync            ::GlobalChannelOffset;
          }
          if        ((itr1 = ExtMppc     .find(board )) != ExtMppc     .end() &&
                     (itr2 = itr1->second.find(mppcCh)) != itr1->second.end()) {
            datum.Channel       = *itr2->second.begin() + ExtinctionDetector::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        for (auto&& subCh : GetSubHitChannels()) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if        ((itr1 = CM::MrSync  .find(board)) != CM::MrSync  .end() &&
                     (itr2 = itr1->second.begin()    ) != itr1->second.end()) {
            datum.MrSyncChannel = *itr2->second.begin() + MrSync            ::GlobalChannelOffset;
          }
          if        ((itr1 = ExtSub      .find(board)) != ExtSub      .end() &&
                     (itr2 = itr1->second.find(subCh)) != itr1->second.end()) {
            datum.Channel       = *itr2->second.begin() + ExtinctionDetector::GlobalChannelOffset;
          } else if ((itr1 = Hod         .find(board)) != Hod         .end() &&
                     (itr2 = itr1->second.find(subCh)) != itr1->second.end()) {
            datum.Channel       = *itr2->second.begin() + Hodoscope         ::GlobalChannelOffset;
          } else if ((itr1 = Tc          .find(board)) != Tc          .end() &&
                     (itr2 = itr1->second.find(subCh)) != itr1->second.end()) {
            datum.Channel       = *itr2->second.begin() + TimingCounter     ::GlobalChannelOffset;
          } else if ((itr1 = Bh          .find(board)) != Bh          .end() &&
                     (itr2 = itr1->second.find(subCh)) != itr1->second.end()) {
            datum.Channel       = *itr2->second.begin() + BeamlineHodoscope ::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        if (MrSync) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if ((itr1 = CM::MrSync  .find(board)) != CM::MrSync  .end() &&
              (itr2 = itr1->second.begin()    ) != itr1->second.end()) {
            datum.Channel       = *itr2->second.begin() + MrSync            ::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        return data;
      }

#if KC705_FORMAT_VERSION == 1

      inline void WriteHeader(std::ofstream& file) const {
        Packet_t data = { 0 };
        data[ 0] = 0xAB;
        data[ 1] = 0xB0;
        data[ 2] = 0x00;
        data[ 3] = 0x12;
        data[ 4] = 0x34;
        data[ 5] = 0x56;
        data[ 6] = 0x70;
        data[ 7] = 0x12;
        data[ 8] = 0x34;
        data[ 9] = 0x56;
        data[10] = 0x70 + (BoardId & 0xF);
        data[11] = ((Spill >> 8) & 0xFF);
        data[12] = ( Spill       & 0xFF);
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteData(std::ofstream& file) const {
        Packet_t data = { 0 };
        data[ 0] = ((MppcBit >> 56) & 0xFFU);
        data[ 1] = ((MppcBit >> 48) & 0xFFU);
        data[ 2] = ((MppcBit >> 40) & 0xFFU);
        data[ 3] = ((MppcBit >> 32) & 0xFFU);
        data[ 4] = ((MppcBit >> 24) & 0xFFU);
        data[ 5] = ((MppcBit >> 16) & 0xFFU);
        data[ 6] = ((MppcBit >>  8) & 0xFFU);
        data[ 7] = ((MppcBit >>  0) & 0xFFU);
        data[ 8] = ((SubBit  >>  4) & 0xFFU);
        data[ 9] = ((SubBit  <<  4) & 0xF0U) + (MrSync * 0x08U) + ((Tdc >> 24) & 0x07U);
        data[10] = ((Tdc     >> 16) & 0xFFU);
        data[11] = ((Tdc     >>  8) & 0xFFU);
        data[12] = ((Tdc     >>  0) & 0xFFU);
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteFooter(std::ofstream& file) const {
        Packet_t data = { 0 };
        data[ 0] = ((EMCount >> 8) & 0xFF);
        data[ 1] = ( EMCount       & 0xFF);
        data[ 2] = 0x00;
        data[ 3] = 0x00;
        data[ 4] = 0x00;
        data[ 5] = 0xAA;
        data[ 6] = 0xAA;
        data[ 7] = 0xAA;
        data[ 8] = 0xAA;
        data[ 9] = 0xAA;
        data[10] = 0xAA;
        data[11] = 0xAA;
        data[12] = 0xAA;
        file.write((Char_t*)&data, sizeof(Packet_t));
      }

#else

      inline void WriteHeader(std::ofstream& file) const {
        Packet_t data = { 0 };
        data[ 0] = 0x01;
        data[ 1] = 0x23;
        data[ 2] = 0x45;
        data[ 3] = 0x67;
        data[ 4] = ((Spill >> 8) & 0xFF);
        data[ 5] = ((Spill     ) & 0xFF);
        data[ 6] = (BoardId & 0x0F);
        data[ 7] = 0x01;
        data[ 8] = 0x23;
        data[ 9] = 0x45;
        data[10] = 0x67;
        data[11] = 0x89;
        data[12] = 0xAB;
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteData(std::ofstream& file) const {
        Packet_t data = { 0 };
        data[ 0] = ((MppcBit >> 56) & 0xFFU);
        data[ 1] = ((MppcBit >> 48) & 0xFFU);
        data[ 2] = ((MppcBit >> 40) & 0xFFU);
        data[ 3] = ((MppcBit >> 32) & 0xFFU);
        data[ 4] = ((MppcBit >> 24) & 0xFFU);
        data[ 5] = ((MppcBit >> 16) & 0xFFU);
        data[ 6] = ((MppcBit >>  8) & 0xFFU);
        data[ 7] = ((MppcBit >>  0) & 0xFFU);
        data[ 8] = ((SubBit  >>  4) & 0xFFU);
        data[ 9] = ((SubBit  <<  4) & 0xF0U) + (MrSync * 0x08U) + ((Tdc >> 24) & 0x7U);
        data[10] = ((Tdc     >> 16) & 0xFFU);
        data[11] = ((Tdc     >>  8) & 0xFFU);
        data[12] = ((Tdc     >>  0) & 0xFFU);
        file.write((Char_t*)&data, sizeof(Packet_t));
      }
      inline void WriteFooter(std::ofstream& file) const {
        Packet_t data = { 0 };
        data[ 0] = 0xAA;
        data[ 1] = 0xAA;
        data[ 2] = 0xAA;
        data[ 3] = 0xAA;
        data[ 4] = ((Spill   >>  8) & 0xFF);
        data[ 5] = ( Spill          & 0xFF);
        data[ 6] = ((EMCount >>  8) & 0xFF);
        data[ 7] = ( EMCount        & 0xFF);
        data[ 8] = ((WRCount >> 24) & 0xFF);
        data[ 9] = ((WRCount >> 16) & 0xFF);
        data[10] = ((WRCount >>  8) & 0xFF);
        data[11] = ((WRCount      ) & 0xFF);
        data[12] = 0xAB;
        file.write((Char_t*)&data, sizeof(Packet_t));
      }

#endif

      virtual Int_t DecodeEventMatchNumber(const std::vector<TdcData>& /*data*/) const override {
        return EMCount;
      }

    };

    class Decoder {
    public:
      Kc705Data Data;
      TTree*    Tree;

      Decoder() : Tree(nullptr) {
      }

      void InitializeTree(const std::string& treename = "tree") {
        // std::cout << "Kc705::Decoder::InitializeTree()" << std::endl;
        if (Tree) {
          delete Tree;
        }
        Tree = new TTree(treename.data(), treename.data());

        Data.CreateBranch(Tree);
      }

      std::size_t Decode(std::basic_istream<char>& file) {
        // std::cout << "Kc705::Decoder::Decode()" << std::endl;
        std::size_t count = 0UL;
        if (!Tree) {
          std::cout << "[warning] tree has not initialized yet" << std::endl;
        } else {
          for (; Read(file); ++count) {
            // Data.Show();
            if (Data.Type == DataType::Data || Data.Type == DataType::Footer) {
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
        switch (Data.Type) {
        case DataType::None:
        case DataType::HeaderError:
        case DataType::Footer:
          return Data.ReadHeader(file, packet);
        case DataType::Header:
        case DataType::Data:
          return Data.ReadDataOrFooter(file, packet);
        default:
          std::cout << "[error] invalid data type" << std::endl;
          exit(1);
        }
      }

    };

  }

}

#endif
