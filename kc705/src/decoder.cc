#include <iostream>
#include <fstream>
#include "Units.hh"
#include "Kc705.hh"
#include "TStyle.h"
#include "TFile.h"
#include "TH2.h"
#include "ArgReader.hh"

Int_t main(Int_t argc, Char_t** argv) {
  // if (argc < 2) {
  //   std::cout << "Usage: " << argv[0] << " [Input] [(optional)Output]" << std::endl;
  //   return 1;
  // }

  // const std::string ifilename(argv[1]);
  // const std::string ofilename(argc > 2 ? argv[2] : "");
  // const Bool_t      isVerbose = false;

  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Input"     ,                    "A rawdata filename");
  args->AddOpt<std::string>("Output"    , 'o', "output"    , "A output filename", "");
  args->AddOpt<std::size_t>("Size"      , 's', "size"      , "A byte size of rawdata", "13000000UL");
  args->AddOpt<Int_t>      ("PacketLoss", 'p', "packetloss", "Check packet loss", "0");
  args->AddOpt<Int_t>      ("Debug"     , 'd', "debug"     , "Output debug infos", "0"); // 0: none, 1: show info and packet history, 2: show check result, 3: show each check result
  args->AddOpt             ("Wait"      , 'w', "wait"      , "Wait enter key to see debug infos");
  args->AddOpt             ("Verbose"   , 'v', "verbose"   , "Output filled values");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string ifilename       = args->GetValue("Input");
  const std::string ofilename       = args->GetValue("Output");
  const std::size_t fileSize        = args->GetValue<std::size_t>("Size");
  const Int_t       checkPacketLoss = args->GetValue<Int_t>("PacketLoss");
  const Int_t       debugLevel      = args->GetValue<Int_t>("Debug");
  const Bool_t      shouldWait      = args->IsSet("Wait");
  const Bool_t      isVerbose       = args->IsSet("Verbose");

  std::string ofilename2;
  if (ofilename.empty()) {
    TString buff = ifilename.data();
    if (buff.EndsWith(".dat") ||
        buff.EndsWith(".txt")) {
      buff.Replace(buff.Length() - 4, 4, "");
    } else if (buff.EndsWith(".data")) {
      buff.Replace(buff.Length() - 5, 5, "");
    }
    ofilename2 = buff + ".root";
  } else {
    ofilename2 = ofilename;
  }

  std::cout << "=== Open Input File" << std::endl;
  std::ifstream ifile(ifilename, std::ios::binary);
  if (!ifile) {
    std::cout << "[error] input file is not opened, " << ifilename << std::endl;
    return 1;
  }

  std::cout << "=== Create Output File" << std::endl;
  std::cout << ofilename2 << std::endl;
  TFile* ofile = new TFile(ofilename2.data(), "RECREATE");
  if (!ofile->IsOpen()) {
    std::cout << "[error] output file is not opened, " << ofilename2 << std::endl;
    return 1;
  }

  std::cout << "=== Setup Packetloss Checker" << std::endl;
  std::vector<std::size_t> usedCh { 32, 33,42, 43, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 };
  std::vector<std::size_t> unusedCh;
  for (std::size_t ch = 0; ch < Extinction::Kc705::MppcNch; ++ch) {
    if (std::find(usedCh.begin(), usedCh.end(), ch) == usedCh.end()) {
      unusedCh.push_back(ch);
    }
  }

  Extinction::Kc705::Decoder decoder;
  Extinction::Kc705::Kc705Data lastData1, lastData2, lastData3, lastData4;

  Extinction::Kc705::Packet_t packet;
  const std::size_t bytes   = sizeof(packet);
  const std::size_t entries = fileSize / bytes;

  std::vector<UChar_t*> lastPackets(2U);
  std::size_t lastPacketsSize = 0;
  for (auto&& lastPacket : lastPackets) {
    lastPacket = new UChar_t[bytes];
    for (std::size_t i = 0; i < bytes; ++i) {
      lastPacket[i] = 0U;
    }
  }

  auto hasPacketLoss =
    [&] (Extinction::Kc705::Packet_t data) {
      for (auto&& ch : unusedCh) {
        const std::size_t bit = ch + 40;
        if (Extinction::Kc705::HasBit(data, bit)) {
          return true;
        }
      }
      return false;
    };

  auto hasNoHitBit =
    [&] (Extinction::Kc705::Packet_t data) {
      std::size_t cnt = 0;
      for (std::size_t bit = 27; bit < 40; ++bit) {
        if (Extinction::Kc705::HasBit(data, bit)) {
          ++cnt;
        }
      }
      for (auto&& ch : usedCh) {
        const std::size_t bit = ch + 40;
        if (Extinction::Kc705::HasBit(data, bit)) {
          ++cnt;
        }
      }
      return cnt == 0;
    };

  auto invalidTdc1 =
    [&] (Extinction::Kc705::Packet_t data) {
      if (lastData3.Type != Extinction::Kc705::DataType::None) {
          const Int_t tdc3 = lastData3.Tdc;
          const Int_t tdc2 = lastData2.Tdc;
          const Int_t tdc1 = lastData1.Tdc;
          const Int_t tdc0 = Extinction::Kc705::GetTdc(data);
          auto dTdc = [](Int_t tdc1, Int_t tdc2) { return tdc1 > tdc2 ? (tdc1 - tdc2) : (tdc1 + 0x8000000 - tdc2); };
          if (lastData1.Type == Extinction::Kc705::DataType::Data) {
            return 
              !(dTdc(tdc0, tdc1) < dTdc(tdc0, tdc2) &&
                dTdc(tdc0, tdc2) < dTdc(tdc0, tdc3) &&
                (tdc2 < tdc0 || ((tdc2 & 0x7000000) && dTdc(tdc0, tdc2) < 0x0800000)));
          } else {
            return 
              !(dTdc(tdc0, tdc2) < dTdc(tdc0, tdc3) &&
                (tdc2 < tdc0 || ((tdc2 & 0x7000000) && dTdc(tdc0, tdc2) < 0x0800000)));
          }
      }
      return false;
    };

  auto invalidTdc2 =
    [&] (Extinction::Kc705::Packet_t data) {
      if (lastData3.Type != Extinction::Kc705::DataType::None) {
          const Int_t tdc3 = lastData3.Tdc;
          const Int_t tdc2 = lastData2.Tdc;
          const Int_t tdc0 = Extinction::Kc705::GetTdc(data);
          auto dTdc = [](Int_t tdc1, Int_t tdc2) { return tdc1 > tdc2 ? (tdc1 - tdc2) : (tdc1 + 0x8000000 - tdc2); };
          return
            !(dTdc(tdc0, tdc2) < dTdc(tdc0, tdc3) &&
              (tdc2 < tdc0 || ((tdc2 & 0x7000000) && dTdc(tdc0, tdc2) < 0x0800000)));
      }
      return false;
    };

  auto goMore1 =
    [&] (Extinction::Kc705::Packet_t data) {
      if (lastData3.Type != Extinction::Kc705::DataType::None) {
        const std::size_t bitCnt11 = Extinction::Kc705::GetBitCount((UChar_t)( lastData1.Tdc       & 0xFF));
        const std::size_t bitCnt12 = Extinction::Kc705::GetBitCount((UChar_t)((lastData1.Tdc >> 8) & 0xFF));
        const std::size_t bitCnt01 = Extinction::Kc705::GetBitCount(               data[bytes - 1]        );
        const std::size_t bitCnt02 = Extinction::Kc705::GetBitCount(               data[bytes - 2]        );

        const Bool_t may1shift1 = bitCnt11            < 3; // packet loss at head or tail of lastData1
        const Bool_t may1shift2 = bitCnt11 + bitCnt12 < 3;
        const Bool_t may0shift1 = bitCnt01            < 3; // packet loss at head or tail of data
        const Bool_t may0shift2 = bitCnt01 + bitCnt02 < 3;

        if (may1shift1 || may1shift2 || may0shift1 || may0shift2) {
          Extinction::Kc705::Packet_t data1;
          lastData1.GetData(data1);
          const Int_t tdc3shift0 = lastData3.Tdc;
          const Int_t tdc2shift0 = lastData2.Tdc;
          const Int_t tdc1shift0 = lastData1.Tdc;
          const Int_t tdc1shift1 = Extinction::Kc705::GetTdc(data1 - 1);
          const Int_t tdc1shift2 = Extinction::Kc705::GetTdc(data1 - 2);
          const Int_t tdc0shift0 = Extinction::Kc705::GetTdc(data     );
          const Int_t tdc0shift1 = Extinction::Kc705::GetTdc(data  - 1);
          const Int_t tdc0shift2 = Extinction::Kc705::GetTdc(data  - 2);
          auto dTdc = [](Int_t tdc1, Int_t tdc2) { return tdc1 > tdc2 ? (tdc1 - tdc2) : (tdc1 + 0x8000000 - tdc2); };
          return
             // 1 byte packet loss at head of lastData1
             ((may1shift1 && may0shift1) &&
              dTdc(tdc1shift1, tdc3shift0) < dTdc(tdc1shift0, tdc3shift0) &&
              dTdc(tdc1shift1, tdc2shift0) < dTdc(tdc1shift0, tdc2shift0) &&
              dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
              dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift1, tdc3shift0) &&
              dTdc(tdc0shift1, tdc1shift1) < dTdc(tdc0shift1, tdc2shift0)) ||
             // 2 byte packet loss at head of lastData1
             ((may1shift2 && may0shift2) &&
              dTdc(tdc1shift2, tdc3shift0) < dTdc(tdc1shift0, tdc3shift0) &&
              dTdc(tdc1shift2, tdc2shift0) < dTdc(tdc1shift0, tdc2shift0) &&
              dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
              dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift2, tdc3shift0) &&
              dTdc(tdc0shift2, tdc1shift2) < dTdc(tdc0shift2, tdc2shift0)) ||
             // 1 byte packet loss at tail of lastData1
             ((may1shift1 && may0shift1) &&
              dTdc(tdc0shift1, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
              dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
              dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift1, tdc3shift0)) ||
             // 2 byte packet loss at tail of lastData1
             ((may1shift2 && may0shift2) &&
              dTdc(tdc0shift2, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
              dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
              dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift2, tdc3shift0)) ||
             // 1 byte packet loss at head of data
             ((may0shift1) &&
              dTdc(tdc0shift1, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
              dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
              dTdc(tdc0shift1, tdc1shift0) < dTdc(tdc0shift0, tdc1shift0) &&
              dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift1, tdc3shift0) &&
              dTdc(tdc0shift1, tdc1shift0) < dTdc(tdc0shift1, tdc2shift0)) ||
             // 2 byte packet loss at head of data
             ((may0shift2) &&
              dTdc(tdc0shift2, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
              dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
              dTdc(tdc0shift2, tdc1shift0) < dTdc(tdc0shift0, tdc1shift0) &&
              dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift2, tdc3shift0) &&
              dTdc(tdc0shift2, tdc1shift0) < dTdc(tdc0shift2, tdc2shift0));
        }
      }
      return false;
    };

  auto goMore2 =
    [&] (Extinction::Kc705::Packet_t data) {
      if (lastData3.Type != Extinction::Kc705::DataType::None) {
        const Int_t tdc3shift0 = lastData3.Tdc;
        const Int_t tdc2shift0 = lastData2.Tdc;
        const Int_t tdc1shift0 = lastData1.Tdc;
        const Int_t tdc0shift0 = Extinction::Kc705::GetTdc(data    );
        const Int_t tdc0shift1 = Extinction::Kc705::GetTdc(data - 1);
        const Int_t tdc0shift2 = Extinction::Kc705::GetTdc(data - 2);
        auto dTdc = [](Int_t tdc1, Int_t tdc2) { return tdc1 > tdc2 ? (tdc1 - tdc2) : (tdc1 + 0x8000000 - tdc2); };
        if (lastData1.Type == Extinction::Kc705::DataType::Data) {
          return
            (dTdc(tdc0shift1, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
             dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
             dTdc(tdc0shift1, tdc1shift0) < dTdc(tdc0shift0, tdc1shift0) &&
             dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift1, tdc3shift0)) ||
            (dTdc(tdc0shift2, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
             dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
             dTdc(tdc0shift2, tdc1shift0) < dTdc(tdc0shift0, tdc1shift0) &&
             dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift2, tdc3shift0));
        } else {
          return
            (dTdc(tdc0shift1, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
             dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
             dTdc(tdc0shift1, tdc2shift0) < dTdc(tdc0shift1, tdc3shift0)) ||
            (dTdc(tdc0shift2, tdc3shift0) < dTdc(tdc0shift0, tdc3shift0) &&
             dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift0, tdc2shift0) &&
             dTdc(tdc0shift2, tdc2shift0) < dTdc(tdc0shift2, tdc3shift0));
        }
      }
      return false;
    };

  auto goNext1 =
    [&] (Extinction::Kc705::Packet_t data) {
      if (lastData1.Type != Extinction::Kc705::DataType::None) {
        if (Extinction::Kc705::GetBitCount(data[0]) > 3 &&
            Extinction::Kc705::GetBitCount(lastData1.MppcBit >> 56) > 3) {
          return true;
        }
      }
      return false;
    };

  auto goNext2 =
    [&] (Extinction::Kc705::Packet_t data) {
      if (lastData3.Type != Extinction::Kc705::DataType::None) {
        const Int_t tdc0         = Extinction::Kc705::GetTdc(data);
        const Int_t tdc3shift0_1 =  lastData3.Tdc         & 0x7FFFF00;
     // const Int_t tdc3shift0_2 =  lastData3.Tdc         & 0x7FF0000;
        const Int_t tdc2shift0_1 =  lastData2.Tdc         & 0x7FFFF00;
     // const Int_t tdc2shift0_2 =  lastData2.Tdc         & 0x7FF0000;
        const Int_t tdc0shift0_1 = (          tdc0      ) & 0x7FFFF00;
     // const Int_t tdc0shift0_2 = (          tdc0      ) & 0x7FF0000;
        const Int_t tdc0shift1   = (          tdc0 <<  8) & 0x7FFFF00;
     // const Int_t tdc0shift2   = (          tdc0 << 16) & 0x7FF0000;
        auto dTdc = [](Int_t tdc1, Int_t tdc2) { return tdc1 >= tdc2 ? tdc1 - tdc2 : tdc1 + 0x8000000 - tdc2; };
        return
          (dTdc(tdc0shift1, tdc3shift0_1) <= dTdc(tdc0shift0_1, tdc3shift0_1) &&
           dTdc(tdc0shift1, tdc2shift0_1) <= dTdc(tdc0shift0_1, tdc2shift0_1))/* ||
          (dTdc(tdc0shift2, tdc3shift0_2) <= dTdc(tdc0shift0_2, tdc3shift0_2) &&
           dTdc(tdc0shift2, tdc2shift0_2) <= dTdc(tdc0shift0_2, tdc2shift0_2))*/;
      }
      return false;
    };

  std::cout << "=== Initialize Output Objects" << std::endl;
  gStyle->SetOptStat(111111);

  decoder.Tree = new TTree("tree", ifilename.data());
  lastData2.CreateBranch(decoder.Tree);

  TH2* hBitMap = new TH2D("hBitMap",
                          Form("%s;Entries;[byte]", ifilename.data()),
                          1000, 0, entries, bytes * 2, 0, bytes);

  std::size_t readSize = 0UL;
  std::size_t dataEntry = 0UL;
  auto fillDataAndShift =
    [&]() {
      if (lastData2.Type == Extinction::Kc705::DataType::Data) {
        if (isVerbose) printf("%ld   %016llx %03x %01x %07x \n",
                              dataEntry, lastData2.MppcBit, lastData2.SubBit, (Int_t)lastData2.MrSync, lastData2.Tdc);
        decoder.Tree->Fill();
        ++dataEntry;
      }

      lastData4 = lastData3;
      lastData3 = lastData2;
      lastData2 = lastData1;
      lastData1 = decoder.Data;
    };

  auto fillPacketLossAndShift =
    [&](Extinction::Kc705::Packet_t data) {
      const Int_t tmpType = lastData2.Type;
      lastData2.Type = Extinction::Kc705::DataType::PacketLoss;
      decoder.Tree->Fill();
      ++dataEntry;
      lastData2.Type = tmpType;
      lastData1    = lastData2;
      decoder.Data = lastData2;
      lastData1   .SetDataAsData(data);
      decoder.Data.SetDataAsData(data);
    };

  auto eraseNewestHistory =
    [&]() {
      // Erase newest history
      --lastPacketsSize;
    };

  auto fillAndEraseOldestHistory =
    [&]() {
      // Fill oldest history
      for (std::size_t bit = 0; bit < bytes * 8; ++bit) {
        if (Extinction::Kc705::HasBit(lastPackets.front(), bit)) {
          hBitMap->Fill(readSize / bytes, bit / 8.0);
        }
      }
      // Erase oldest history
      for (std::size_t i = 1, n = lastPackets.size(); i < n; ++i) {
        std::iter_swap(lastPackets.begin() + i - 1, lastPackets.begin() + i);
      }
      --lastPacketsSize;
    };

  // Show debug information function
  std::vector<UChar_t> thrownPacket;
  auto showPacketHistory =
    [&](Extinction::Kc705::Packet_t data) {
      lastData4.ShowAsHex();
      lastData3.ShowAsHex();
      lastData2.ShowAsHex();
      lastData1.ShowAsHex();
      if (!thrownPacket.empty()) {
        for (UChar_t& val : thrownPacket) {
          printf("%02x ", val);
        }
        std::cout << std::endl;
      }
      Extinction::Kc705::ShowAsHex(data);
    };

  auto showFirstCheckResult =
    [&](Extinction::Kc705::Packet_t data) {
      std::cout << "--- hasPacketLoss = " << hasPacketLoss(data) << std::endl
                << "--- hasNoHitBit   = " << hasNoHitBit  (data) << std::endl
                << "--- invalidTdc1   = " << invalidTdc1  (data) << std::endl
            //  << "--- invalidTdc2   = " << invalidTdc2  (data) << std::endl
                << "--- goMore1       = " << goMore1      (data) << std::endl
                << "--- goNext1       = " << goNext1      (data) << std::endl;
    };

  auto showSecondCheckResult =
    [&](Extinction::Kc705::Packet_t data) {
      std::cout << "--- hasPacketLoss = " << hasPacketLoss(data) << std::endl
                << "--- hasNoHitBit   = " << hasNoHitBit  (data) << std::endl
                << "--- invalidTdc2   = " << invalidTdc2  (data) << std::endl
                << "--- goMore1       = " << goMore1      (data) << std::endl
                << "--- goMore2       = " << goMore2      (data) << std::endl
                << "--- goNext2       = " << goNext2      (data) << std::endl;
    };

  auto waitEnter =
    [&]() {
      if (shouldWait) {
        std::cout << std::flush;
        while (std::getchar() != '\n');
      }
    };

  std::cout << "=== Decode" << std::endl;
  std::size_t count = 0UL;
  Char_t lastType = Extinction::Kc705::DataType::None;
  for (; decoder.Read(ifile, &packet); ++count) {
    readSize += bytes;
    if (count % 1000000UL == 0) {
      std::cout << ">> " << count << std::endl;
    }

    if (checkPacketLoss && decoder.Data.Type == Extinction::Kc705::DataType::HeaderError) {
      if (!Extinction::Kc705::IsFooter(packet)) {
        for (std::size_t i = 0; i < lastPacketsSize; ++i) {
          Extinction::Kc705::ShowAsHex(lastPackets[i]);
        }
        Extinction::Kc705::ShowAsHex(packet);
        decoder.Data.Type = Extinction::Kc705::DataType::Data;
      }
    }

    if (decoder.Data.Type == Extinction::Kc705::DataType::Header) {
      std::cout << "[info] begin of spill " << decoder.Data.Spill << std::endl;

    } else if (decoder.Data.Type == Extinction::Kc705::DataType::HeaderError) {
      if (checkPacketLoss && lastType != Extinction::Kc705::DataType::HeaderError) {
        for (std::size_t i = 0; i < lastPacketsSize; ++i) {
          Extinction::Kc705::ShowAsHex(lastPackets[i]);
        }
        Extinction::Kc705::ShowAsHex(packet);
      }

    } else if (decoder.Data.Type == Extinction::Kc705::DataType::Footer) {
      std::cout << "[info] end of spill " << decoder.Data.Spill << std::endl;
      fillDataAndShift();
      fillDataAndShift();
      lastData1.Clear();
      lastData2.Clear();
      lastData3.Clear();
      lastData4.Clear();

    } else {
      if (checkPacketLoss &&
          (hasPacketLoss(packet) ||
           hasNoHitBit  (packet) ||
           invalidTdc1  (packet) ||
        // invalidTdc2  (packet) ||
           goMore1      (packet) ||
           goNext1      (packet))) {
        thrownPacket.clear();
        if (debugLevel > 1) {
          std::cout << "[info] packet loss detected at entry = " << dataEntry << std::endl;
          showPacketHistory(packet);
          showFirstCheckResult(packet);
        }

        while (hasPacketLoss(packet) ||
               hasNoHitBit  (packet) ||
               invalidTdc2  (packet) ||
               goMore1      (packet) ||
               goMore2      (packet) ||
               goNext2      (packet)) {
          if (thrownPacket.size() == bytes) {
            break;
          }
          if (debugLevel > 2) {
            showPacketHistory(packet);
            showSecondCheckResult(packet);
          }

          // Throw first byte and read next 1 byte
          thrownPacket.push_back(packet[0]);
          for (std::size_t i = 1; i < bytes; ++i) {
            packet[i - 1] = packet[i];
          }
          ifile.read((char*)(packet + bytes - 1), 1);
          readSize += 1;
        }

        if (thrownPacket.size()) {
          if (debugLevel) {
            std::cout << "[info] packet loss detected at entry = " << dataEntry << ", throw " << thrownPacket.size() << " bytes" << std::endl;
            showPacketHistory(packet);
          }

          if (thrownPacket.size() == bytes) {
            if (!hasPacketLoss(thrownPacket.data()) &&
                !hasNoHitBit  (thrownPacket.data()) &&
                !invalidTdc1  (thrownPacket.data()) &&
             // !invalidTdc2  (thrownPacket.data()) &&
                 goMore1      (thrownPacket.data()) &&
                !goNext1      (thrownPacket.data())) {
              Extinction::Kc705::Kc705Data tmpData3 = lastData3;
              lastData3 = lastData2;
              lastData2 = lastData1;
              lastData1.SetDataAsData(thrownPacket.data());
              const Bool_t missChecked =
                !(hasPacketLoss(packet) ||
                  hasNoHitBit  (packet) ||
                  invalidTdc1  (packet) ||
               // invalidTdc2  (packet) ||
                  goMore1      (packet) ||
                  goNext1      (packet));
              lastData1 = lastData2;
              lastData2 = lastData3;
              lastData3 =  tmpData3;

              if (missChecked) {
                if (debugLevel) {
                  std::cout << "[info] miss checked as packet loss" << std::endl;
                }

                fillDataAndShift();

                decoder.Data.SetDataAsData(packet);
                fillDataAndShift();

                continue;
              } else if (debugLevel > 1) {
                showSecondCheckResult(thrownPacket.data());
                waitEnter();
              }
            } else if (!hasPacketLoss(thrownPacket.data()) &&
                        hasNoHitBit  (thrownPacket.data()) &&
                       !invalidTdc1  (thrownPacket.data()) &&
                    // !invalidTdc2  (thrownPacket.data()) &&
                       !goMore1      (thrownPacket.data()) &&
                       !goNext1      (thrownPacket.data())) {
              // TDC does not have any problems
              if (debugLevel) {
                std::cout << "[info] empty data record was detected" << std::endl;
              }
            } else if (debugLevel > 1) {
              showSecondCheckResult(thrownPacket.data());
              waitEnter();
            }
          }

          // Record packet loss
          fillPacketLossAndShift(packet);

          // Erase last history
          if (lastPacketsSize) {
            eraseNewestHistory();
          }

          continue;
        }
      }

      fillDataAndShift();
    }

    if (lastPacketsSize >= lastPackets.size()) {
      fillAndEraseOldestHistory();
    }

    lastType = decoder.Data.Type;
    std::memcpy(lastPackets[lastPacketsSize], packet, sizeof(packet));
    ++lastPacketsSize;
  }
  std::cout << "[info] # of data record = " << dataEntry << std::endl;

  while (lastPacketsSize) {
    fillAndEraseOldestHistory();
  }

  std::cout << "=== Write Objects" << std::endl;
  std::cout << decoder.Tree->GetName() << std::endl;
  decoder.Tree->Write();
  hBitMap->Write();

  std::cout << "=== Close Files" << std::endl;
  ifile.close();
  ofile->Close();

  return 0;
}
