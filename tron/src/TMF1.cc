#include "TMF1.hh"
#include "Linq.hh"

namespace {

  Int_t NparSum(Int_t nfun, const TF1** funcs) {
    return Tron::Linq::From(funcs, funcs + nfun)
      .Select([](const TF1* func) { return func->GetNpar(); })
      .Sum();
  }

  Int_t NparSum(const std::initializer_list<const TF1*>& funcs) {
    return Tron::Linq::From(funcs.begin(), funcs.end())
      .Select([](const TF1* func) { return func->GetNpar(); })
      .Sum();
  }
  
}

//==================================================
Tron::TMF1::TMF1()
  : TF1(),
    fFuncs(),
    fIparOffset() {
}

Tron::TMF1::TMF1(const Char_t* name, Int_t nfun, const TF1** funcs, Double_t xmin, Double_t xmax)
  : TF1(name, this, &TMF1::Func, xmin, xmax, NparSum(nfun, funcs), "", ""),
    fFuncs(nfun),
    fIparOffset(nfun + 1) {

  //--- Init par offsets
  fIparOffset[0] = 0;
  for (Int_t ifun = 1; ifun <= nfun; ++ifun) {
    const TF1* func = funcs[ifun - 1];
    fIparOffset[ifun] = fIparOffset[ifun - 1] + func->GetNpar();
  }

  //--- Init funcs
  SetLineColor(kRed);
  for (Int_t ifun = 0; ifun < nfun; ++ifun) {
    const TF1* func = funcs[ifun];
    fFuncs[ifun] = new TF1(*func);
    fFuncs[ifun]->SetLineColor(kBlue);
    fFuncs[ifun]->SetLineWidth(GetLineWidth() - 1);
    
    const Int_t npar = funcs[ifun]->GetNpar();
    for (Int_t ipar = 0; ipar < npar; ++ipar) {
      const Int_t thisIpar = GetIpar(ifun, ipar);
      SetParName(thisIpar, Form("f%d_%s", ifun, func->GetParName(ipar)));
      SetParameter(thisIpar, func->GetParameter(ipar));
    }
  }
  
  //--- Import parameters
  ExportParams();
}

Tron::TMF1::TMF1(const Char_t* name, const std::initializer_list<const TF1*>& funcs, Double_t xmin, Double_t xmax)
  : TF1(name, this, &TMF1::Func, xmin, xmax, NparSum(funcs), "", ""),
    fFuncs(funcs.size()),
    fIparOffset(funcs.size() + 1) {
  const Int_t nfun = funcs.size();

  //--- Init par offsets
  fIparOffset[0] = 0;
  for (Int_t ifun = 1; ifun <= nfun; ++ifun) {
    const TF1* func = *(funcs.begin() + ifun - 1);
    fIparOffset[ifun] = fIparOffset[ifun - 1] + func->GetNpar();
  }

  //--- Init funcs
  SetLineColor(kRed);
  for (Int_t ifun = 0; ifun < nfun; ++ifun) {
    const TF1* func = *(funcs.begin() + ifun);
    fFuncs[ifun] = new TF1(*func);
    fFuncs[ifun]->SetLineColor(kBlue);
    fFuncs[ifun]->SetLineWidth(GetLineWidth() - 1);

    const Int_t npar = func->GetNpar();
    for (Int_t ipar = 0; ipar < npar; ++ipar) {
      const Int_t thisIpar = GetIpar(ifun, ipar);
      SetParName(thisIpar, Form("f%d_%s", ifun, func->GetParName(ipar)));
      SetParameter(thisIpar, func->GetParameter(ipar));
    }
  }
  
  //--- Import parameters
  ExportParams();
}

Tron::TMF1::TMF1(const TMF1& mf1)
  : TF1(),
    fFuncs(),
    fIparOffset() {
  ((TF1&)mf1).Copy((TF1&)*this);

  const Int_t nfun = mf1.GetNfun();
  fFuncs.resize(nfun);
  fIparOffset.resize(nfun + 1);

  //--- Init par offsets
  fIparOffset[0] = 0;
  for (Int_t ifun = 1; ifun <= nfun; ++ifun) {
    const TF1* func = mf1.GetFunc(ifun - 1);
    fIparOffset[ifun] = fIparOffset[ifun - 1] + func->GetNpar();
  }

  //--- Init funcs
  for (Int_t ifun = 0; ifun < nfun; ++ifun) {
    const TF1* func = mf1.GetFunc(ifun);
    fFuncs[ifun] = new TF1(*func);

    const Int_t npar = func->GetNpar();
    for (Int_t ipar = 0; ipar < npar; ++ipar) {
      const Int_t thisIpar = GetIpar(ifun, ipar);
      SetParName(thisIpar, Form("f%d_%s", ifun, func->GetParName(ipar)));
      SetParameter(thisIpar, func->GetParameter(ipar));
    }
  }

  //--- Export parameters
  ExportParams();
}

Tron::TMF1::~TMF1() {
  for (auto&& func : fFuncs) {
    delete func;
    func = nullptr;
  }
}

//==================================================
Double_t Tron::TMF1::Func(Double_t* x, Double_t* p) {
  Double_t sum = 0.0;
  for (Int_t ifun = 0, nfun = GetNfun(); ifun < nfun; ++ifun) {
    TF1* func = fFuncs[ifun];
    if (func) {
      sum += func->EvalPar(x, p + fIparOffset[ifun]);
    }
  }
  return sum;
}

//==================================================
void Tron::TMF1::ExportParams() {
  Double_t xmin, xmax;
  GetRange(xmin, xmax);
  
  for (Int_t ifun = 0, nfun = GetNfun(); ifun < nfun; ++ifun) {
    TF1* func = fFuncs[ifun];
    if (!func) {
      continue;
    }
    
    //--- Export range, etc
    func->SetRange(xmin, xmax);
    func->SetNpx(GetNpx());
    
    for (Int_t ipar = 0, npar = func->GetNpar(); ipar < npar; ++ipar) {
      //--- Export params
      func->SetParameter(ipar, GetParameter(GetIpar(ifun, ipar)));
      
      //--- Export par limits
      Double_t parmin, parmax;
      GetParLimits(GetIpar(ifun, ipar), parmin, parmax);
      if (parmin != 0 || parmax != 0) {
        func->SetParLimits(ipar, parmin, parmax);
      }
    }

    func->Update();
  }
}

//==================================================
void Tron::TMF1::Draw(Option_t* option) {
  ExportParams();

  this->TF1::Draw(option);
  for (Int_t i = 0, n = fFuncs.size(); i < n; ++i) {
    TF1* func = fFuncs[i];
    if (func) {
      func->Draw("same");
    }
  }
}

//==================================================
void Tron::TMF1::Update() {
  this->TF1::Update();
  ExportParams();
}
