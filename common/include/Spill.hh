#ifndef Extinction_SpillData_hh
#define Extinction_SpillData_hh

#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TDatime.h"

#include "Units.hh"
#include "Detector.hh"

#include "Math.hh"

namespace Extinction {
  static const std::size_t kNofBunches = 4;

  struct SpillData {
    TDatime Date;
    Int_t   Year;
    Int_t   Month;
    Int_t   Day;
    Int_t   Hour;
    Int_t   Minute;
    Int_t   Second;
    Int_t   EMCount;

    virtual void Clear() {
      SetDate(0);
      EMCount = -1;
    }

    void SetDate(UInt_t date) {
      Date.Set(date);
      Year   = Date.GetYear();
      Month  = Date.GetMonth();
      Day    = Date.GetDay();
      Hour   = Date.GetHour();
      Minute = Date.GetMinute();
      Second = Date.GetSecond();
    }

    virtual void CreateBranch(TTree* tree) {
      tree->Branch("date"   , &Date                   );
      tree->Branch("year"   , &Year   , "year"    "/I");
      tree->Branch("month"  , &Month  , "month"   "/I");
      tree->Branch("day"    , &Day    , "day"     "/I");
      tree->Branch("hour"   , &Hour   , "hour"    "/I");
      tree->Branch("minute" , &Minute , "minute"  "/I");
      tree->Branch("second" , &Second , "second"  "/I");
      tree->Branch("emcount", &EMCount, "emcount" "/I");
    }

    virtual void SetBranchAddress(TTree* tree) {
      tree->SetBranchAddress("date"   , &Date   );
      tree->SetBranchAddress("year"   , &Year   );
      tree->SetBranchAddress("month"  , &Month  );
      tree->SetBranchAddress("day"    , &Day    );
      tree->SetBranchAddress("hour"   , &Hour   );
      tree->SetBranchAddress("minute" , &Minute );
      tree->SetBranchAddress("second" , &Second );
      tree->SetBranchAddress("emcount", &EMCount);
    }

  };

  struct RawSpillData : SpillData {
    Long64_t Entries[GlobalChannel::NofChannels];
    Double_t MrSyncInterval[MrSync::NofChannels];

    RawSpillData() {
      Clear();
    }

    virtual void Clear() override {
      SpillData::Clear();
      std::memset(Entries       , 0, sizeof(Entries       ));
      std::memset(MrSyncInterval, 0, sizeof(MrSyncInterval));
    }

    virtual void CreateBranch(TTree* tree) override {
      // std::cout << "RawSpillData::CreateBranch()" << std::endl;
      SpillData::CreateBranch(tree);
      tree->Branch("entries"       ,  Entries       , Form("entries"        "[%lu]/L", GlobalChannel::NofChannels));
      tree->Branch("mrsyncinterval",  MrSyncInterval, Form("mrsyncinterval" "[%lu]/D", MrSync       ::NofChannels));
    }

    virtual void SetBranchAddress(TTree* tree) override {
      // std::cout << "RawSpillData::SetBranchAddress()" << std::endl;
      SpillData::SetBranchAddress(tree);
      tree->SetBranchAddress("entries"       ,  Entries       );
      tree->SetBranchAddress("mrsyncinterval",  MrSyncInterval);
    }

  };

  struct CoinSpillData : SpillData {
    Long64_t CoinCount;
    Long64_t ExtEntries;
    Long64_t Bh1Entries;
    Long64_t Bh2Entries;
    Long64_t Tc1Entries;
    Long64_t Tc2Entries;
    Double_t BunchCenters [kNofBunches];
    Double_t BunchWidths  [kNofBunches];
    Long64_t BunchMinEdges[kNofBunches];
    Long64_t BunchMaxEdges[kNofBunches];

    CoinSpillData() {
      Clear();
    }

    virtual void Clear() override {
      SpillData::Clear();
      CoinCount      = 0;
      ExtEntries     = 0;
      Bh1Entries     = 0;
      Bh2Entries     = 0;
      Tc1Entries     = 0;
      Tc2Entries     = 0;
      std::memset(BunchCenters , 0, sizeof(BunchCenters ));
      std::memset(BunchWidths  , 0, sizeof(BunchWidths  ));
      std::memset(BunchMinEdges, 0, sizeof(BunchMinEdges));
      std::memset(BunchMaxEdges, 0, sizeof(BunchMaxEdges));
    }

    virtual void CreateBranch(TTree* tree) override {
      // std::cout << "CoinSpillData::CreateBranch()" << std::endl;
      SpillData::CreateBranch(tree);
      tree->Branch("coincount"    , &CoinCount     , Form("coincount"           "/L"             ));
      tree->Branch("extentries"   , &ExtEntries    , Form("extentries"          "/L"             ));
      tree->Branch("bh1entries"   , &Bh1Entries    , Form("bh1entries"          "/L"             ));
      tree->Branch("bh2entries"   , &Bh2Entries    , Form("bh2entries"          "/L"             ));
      tree->Branch("tc1entries"   , &Tc1Entries    , Form("tc1entries"          "/L"             ));
      tree->Branch("tc2entries"   , &Tc2Entries    , Form("tc2entries"          "/L"             ));
      tree->Branch("bunchcenter"  ,  BunchCenters  , Form("bunchcenter"    "[%lu]/D", kNofBunches));
      tree->Branch("bunchwidth"   ,  BunchWidths   , Form("bunchwidth"     "[%lu]/D", kNofBunches));
      tree->Branch("bunchminedge" ,  BunchMinEdges , Form("bunchminedge"   "[%lu]/L", kNofBunches));
      tree->Branch("bunchmaxedge" ,  BunchMaxEdges , Form("bunchmaxedge"   "[%lu]/L", kNofBunches));
    }

    virtual void SetBranchAddress(TTree* tree) override {
      // std::cout << "CoinSpillData::SetBranchAddress()" << std::endl;
      SpillData::SetBranchAddress(tree);
      tree->SetBranchAddress("coincount"    , &CoinCount     );
      tree->SetBranchAddress("extentries"   , &ExtEntries    );
      tree->SetBranchAddress("bh1entries"   , &Bh1Entries    );
      tree->SetBranchAddress("bh2entries"   , &Bh2Entries    );
      tree->SetBranchAddress("tc1entries"   , &Tc1Entries    );
      tree->SetBranchAddress("tc2entries"   , &Tc2Entries    );
      tree->SetBranchAddress("bunchcenter"  ,  BunchCenters  );
      tree->SetBranchAddress("bunchwidth"   ,  BunchWidths   );
      tree->SetBranchAddress("bunchminedge" ,  BunchMinEdges );
      tree->SetBranchAddress("bunchmaxedge" ,  BunchMaxEdges );
    }

  };

  struct LeakSpillData : SpillData {
    Long64_t LeakageEntries;
    Long64_t InBunchEntries;
    Double_t Extinction;
    Double_t ExtinctionErr;

    LeakSpillData() {
      Clear();
    }

    virtual void Clear() override {
      SpillData::Clear();
      LeakageEntries = 0;
      InBunchEntries = 0;
      Extinction     = 1.0;
      ExtinctionErr  = 1.0;
    }

    void SetParticles(Long64_t leakage, Long64_t inBunch) {
      LeakageEntries = leakage;
      InBunchEntries = inBunch;
      if (InBunchEntries) {
        Extinction = (Double_t)LeakageEntries / (Double_t)InBunchEntries;
        if (LeakageEntries) {
          ExtinctionErr = Tron::Math::SqrtOfSumOfSquared(TMath::Sqrt(LeakageEntries) / InBunchEntries, LeakageEntries / TMath::Power((Double_t)InBunchEntries, 1.5));
        } else {
          ExtinctionErr = 1.0 / InBunchEntries;
        }
      } else {
        Extinction = 1.0;
      }
    }

    virtual void CreateBranch(TTree* tree) override {
      // std::cout << "LeakSpillData::CreateBranch()" << std::endl;
      SpillData::CreateBranch(tree);
      tree->Branch("leakage"      , &LeakageEntries, Form("leakage"             "/L"             ));
      tree->Branch("inbunch"      , &InBunchEntries, Form("inbunch"             "/L"             ));
      tree->Branch("extinction"   , &Extinction    , Form("extinction"          "/D"             ));
      tree->Branch("extinctionerr", &ExtinctionErr , Form("extinctionerr"       "/D"             ));
    }

    virtual void SetBranchAddress(TTree* tree) override {
      // std::cout << "LeakSpillData::SetBranchAddress()" << std::endl;
      SpillData::SetBranchAddress(tree);
      tree->SetBranchAddress("leakage"      , &LeakageEntries);
      tree->SetBranchAddress("inbunch"      , &InBunchEntries);
      tree->SetBranchAddress("extinction"   , &Extinction    );
      tree->SetBranchAddress("extinctionerr", &ExtinctionErr );
    }

  };

}

#endif
