#ifndef Extinction_Kc705_hh
#define Extinction_Kc705_hh

#include <iostream>
#include <fstream>
#include <cstring>
#include "TFile.h"
#include "TTree.h"
#include "Units.hh"
#include "Inet.hh"
#include "Tdc.hh"
#include "Detector.hh"

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
      }

#if KC705_FORMAT_VERSION == 1

      inline void SetDataAsHeader(Packet_t packet) {
        Type     = DataType::Header;
        BoardId  = GetBoardId(packet);
        Spill    = Extinction::Kc705::GetSpill(packet);
        MppcBit  = 0;
        SubBit   = 0;
        MrSync   = 0;
        Tdc      = 0;
        EMCount  = 0;
        WRCount  = 0;
        Overflow = 0;
      }

      inline void SetDataAsData(Packet_t packet) {
        const UInt_t lastTdc = Tdc;
        Type    = DataType::Data;
        MppcBit = GetMppcBit(packet);
        SubBit  = GetSubBit(packet);
        MrSync  = GetMrSync(packet);
        Tdc     = GetTdc(packet);
        if (Tdc < lastTdc) { ++Overflow; }
      }

      inline void SetDataAsFooter(Packet_t packet) {
        Type    = DataType::Footer;
        const Int_t lastSpill = Spill;
        Spill   = Extinction::::Kc705::GetSpill(packet);
        if (Spill != lastSpill) {
          std::cerr << "[warning] spill inconsistent" << std::endl;
        }
        MppcBit = 0;
        SubBit  = 0;
        MrSync  = 0;
        Tdc     = 0;
        EMCount = GetEMCount(packet);
        WRCount = GetWRCount(packet);
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
        packet[ 9] = ((SubBit  <<  4) & 0xF0U) + (MrSync * 0x08U) + ((Tdc >> 24) & 0x7FU);
        packet[10] = ((Tdc     >> 16) & 0xFFU);
        packet[11] = ((Tdc     >>  8) & 0xFFU);
        packet[12] = ((Tdc     >>  0) & 0xFFU);
      }

#else

      inline void SetDataAsHeader(Packet_t packet) {
        Type     = DataType::Header;
        BoardId  = GetBoardId(packet);
        Spill    = Extinction::Kc705::GetSpill(packet);
        MppcBit  = 0;
        SubBit   = 0;
        MrSync   = 0;
        Tdc      = 0;
        EMCount  = 0;
        Overflow = 0;
      }

      inline void SetDataAsData(Packet_t packet) {
        const UInt_t lastTdc = Tdc;
        Type    = DataType::Data;
        MppcBit = GetMppcBit(packet);
        SubBit  = GetSubBit(packet);
        MrSync  = GetMrSync(packet);
        Tdc     = GetTdc(packet);
        if (Tdc < lastTdc) { ++Overflow; }
      }

      inline void SetDataAsFooter(Packet_t packet) {
        Type    = DataType::Footer;
        MppcBit = 0;
        SubBit  = 0;
        MrSync  = 0;
        Tdc     = 0;
        EMCount = GetEMCount(packet);
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
        packet[ 9] = ((SubBit  <<  4) & 0xF0U) + (MrSync * 0x08U) + ((Tdc >> 24) & 0x7FU);
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

        if (IsHeader(buff)) {
          SetDataAsHeader(buff);
        } else if (Type == DataType::HeaderError) {
          // Nothing to do
        } else {
          std::cout << "[error] invalid header, maybe "
                    << (IsFooter(buff) ? "footer" : "data") << std::endl;
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

        if (IsFooter(buff)) {
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

      inline virtual Double_t GetTimePerTdc() const override {
        return 5.0 * nsec;
      }

      inline virtual Double_t GetTime() const override {
        return GetTdc2() * GetTimePerTdc();
      }

      inline void CreateBranch(TTree* tree) {
        // std::cout << "Kc705::Kc705Data::CreateBranch()" << std::endl;
        tree->Branch("type"    , &Type    , "type"    "/B");
        tree->Branch("boardId" , &BoardId , "boardId" "/b");
        tree->Branch("spill"   , &Spill   , "spillC"  "/I");
        tree->Branch("emCount" , &EMCount , "emCount" "/s");
        tree->Branch("tdc"     , &Tdc     , "tdc"     "/i");
        tree->Branch("mrSync"  , &MrSync  , "MrSync"  "/O");
        tree->Branch("mppc"    , &MppcBit , "mppc"    "/l");
        tree->Branch("sub"     , &SubBit  , "sub"     "/s");
        tree->Branch("overflow", &Overflow, "overflow""/i");
        tree->SetAlias("tdc2", "tdc + 0x8000000 * overflow");
      }

      inline virtual void SetBranchAddress(TTree* tree) override {
        // std::cout << "Kc705::Kc705Data::SetBranchAddress()" << std::endl;
        tree->SetBranchAddress("type"    , &Type    );
        tree->SetBranchAddress("boardId" , &BoardId );
        tree->SetBranchAddress("spill"   , &Spill   );
        tree->SetBranchAddress("emCount" , &EMCount );
        tree->SetBranchAddress("tdc"     , &Tdc     );
        tree->SetBranchAddress("mrSync"  , &MrSync  );
        tree->SetBranchAddress("mppc"    , &MppcBit );
        tree->SetBranchAddress("sub"     , &SubBit  );
        tree->SetBranchAddress("overflow", &Overflow);
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

      inline virtual Int_t GetSpill() const override {
        return Spill;
      }

      inline virtual std::vector<TdcData> GetTdcData() const override {
        std::vector<TdcData> data;
        const Long64_t tdc  = GetTdc2();
        const Double_t time = GetTime();
        for (auto&& mppcCh : GetMppcHitChannels()) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if (ChannelMap::MrSync.find(BoardId) != ChannelMap::MrSync.end()) {
            datum.MrSyncChannel = ChannelMap::MrSync.at(BoardId) + MrSync::GlobalChannelOffset;
          }
          if (ChannelMap::ExtMppc.find(mppcCh) != ChannelMap::ExtMppc.end()) {
            datum.Channel = ChannelMap::ExtMppc.at(mppcCh) + ExtinctionDetector::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        for (auto&& subCh : GetSubHitChannels()) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if (ChannelMap::MrSync.find(BoardId) != ChannelMap::MrSync.end()) {
            datum.MrSyncChannel = ChannelMap::MrSync.at(BoardId) + MrSync::GlobalChannelOffset;
          }
          if        (ChannelMap::ExtSub.find(subCh) != ChannelMap::ExtSub.end()) {
            datum.Channel = ChannelMap::ExtSub.at(subCh) + ExtinctionDetector::GlobalChannelOffset;
          } else if (ChannelMap::Hod   .find(subCh) != ChannelMap::Hod   .end()) {
            datum.Channel = ChannelMap::Hod   .at(subCh) + Hodoscope        ::GlobalChannelOffset;
          } else if (ChannelMap::Tc    .find(subCh) != ChannelMap::Tc    .end()) {
            datum.Channel = ChannelMap::Tc    .at(subCh) + TimingCounter    ::GlobalChannelOffset;
          } else if (ChannelMap::Bh    .find(subCh) != ChannelMap::Bh    .end()) {
            datum.Channel = ChannelMap::Bh    .at(subCh) + BeamlineHodoscope::GlobalChannelOffset;
          }
          data.push_back(datum);
        }
        if (MrSync) {
          TdcData datum;
          datum.Spill   = Spill;
          datum.Tdc     = tdc;
          datum.Time    = time;
          if (ChannelMap::MrSync.find(BoardId) != ChannelMap::MrSync.end()) {
            datum.Channel = ChannelMap::MrSync.at(BoardId) + MrSync::GlobalChannelOffset;
            datum.MrSyncChannel = datum.Channel;
          }
          data.push_back(datum);
        }
        return data;
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
