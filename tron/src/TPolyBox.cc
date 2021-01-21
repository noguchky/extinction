#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include "String.hh"
#include "TPolyBox.hh"

//--------------------------------------------------
Tron::TPolyBox2D::TPolyBox2D() {
  fList = new TList();
  Set(0,0, 0,0);
  Update();
}

Tron::TPolyBox2D::TPolyBox2D(const Double_t* x, const Double_t* y) {
  fList = new TList();
  Set(x, y);
  Update();
}

Tron::TPolyBox2D::TPolyBox2D(Double_t rx, Double_t ry, Double_t px, Double_t py) {
  fList = new TList();
  Set(rx, ry, px, py);
  Update();
}

Tron::TPolyBox2D::~TPolyBox2D() {
}

//--------------------------------------------------
void Tron::TPolyBox2D::Set(const Double_t* x, const Double_t* y) {
  for (Int_t i = 0; i < 8; ++i) {
    fX[i] = x[i];
    fY[i] = y[i];
  }
}

void Tron::TPolyBox2D::Set(Double_t rx, Double_t ry, Double_t px, Double_t py) {
  fX[0] = fX[3] = fX[4] = fX[7] = rx - px;
  fX[1] = fX[2] = fX[5] = fX[6] = rx + px;

  fY[0] = fY[1] = fY[2] = fY[3] = ry - py;
  fY[4] = fY[5] = fY[6] = fY[7] = ry + py;
}

//--------------------------------------------------
void Tron::TPolyBox2D::Update() {
  TPolyLine* surfaces[4] = {};
  const int n = sizeof(surfaces) / sizeof(*surfaces);

  for (Int_t i = 0; i < n;) {
    {
      // First process
      surfaces[i] = (TPolyLine*)fList->First();
      goto COMMON_PROCESS;
    }
    for (; i < n; ++i) {
      // Second or later processes
      surfaces[i] = (TPolyLine*)fList->After(surfaces[i - 1]);
    COMMON_PROCESS:
      if (!surfaces[i]) {
        fList->Add(surfaces[i] = new TPolyLine(4));
      }
    }
  }

  for (Int_t i = 0; i < n; ++i) {
    TPolyLine* surface = surfaces[0];
    for (Int_t j : { 0, 1, 5, 4 }) {
      const Int_t index = (j + i) % 4 + (j / 4) * 4;
      surface->SetPoint(j, fX[index], fY[index]);
    }
  }
}


//--------------------------------------------------
void Tron::TPolyBox2D::Draw(Option_t* option) {
  Update();
  fList->Draw(option);
}


//--------------------------------------------------
void Tron::TPolyBox2D::SetLineColor(Color_t color) {
  TIter it(fList);
  std::for_each(it.Begin(), it.End(), [&](TObject* obj) {
      auto surface = (TPolyLine*)obj;
      surface->SetLineColor(color);
    });
}

void Tron::TPolyBox2D::SetLineStyle(Style_t style) {
  TIter it(fList);
  std::for_each(it.Begin(), it.End(), [&](TObject* obj) {
      auto surface = (TPolyLine*)obj;
      surface->SetLineStyle(style);
    });
}

void Tron::TPolyBox2D::SetLineWidth(Width_t width) {
  TIter it(fList);
  std::for_each(it.Begin(), it.End(), [&](TObject* obj) {
      auto surface = (TPolyLine*)obj;
      surface->SetLineWidth(width);
    });
}


//--------------------------------------------------
Color_t Tron::TPolyBox2D::GetLineColor() {
  TPolyLine* surface = (TPolyLine*)fList->First();
  return surface ? surface->GetLineColor() : 0;
}

Style_t Tron::TPolyBox2D::GetLineStyle() {
  TPolyLine* surface = (TPolyLine*)fList->First();
  return surface ? surface->GetLineStyle() : 0;
}

Width_t Tron::TPolyBox2D::GetLineWidth() {
  TPolyLine* surface = (TPolyLine*)fList->First();
  return surface ? surface->GetLineWidth() : 0;
}


//==================================================
Tron::TPolyBox3D::TPolyBox3D() {
  fList = new TList();
  Set(0, 0, 0, 0, 0, 0);
  Update();
}

Tron::TPolyBox3D::TPolyBox3D(const Double_t* x, const Double_t* y, const Double_t* z) {
  fList = new TList();
  Set(x, y, z);
  Update();
}

Tron::TPolyBox3D::TPolyBox3D(Double_t rx, Double_t ry, Double_t rz,
                             Double_t px, Double_t py, Double_t pz) {
  fList = new TList();
  Set(rx, ry, rz, px, py, pz);
  Update();
}

Tron::TPolyBox3D::~TPolyBox3D() {
}


//--------------------------------------------------
void Tron::TPolyBox3D::Set(const Double_t* x, const Double_t* y, const Double_t* z) {
  for (Int_t i = 0; i < 8; ++i) {
    fX[i] = x[i];
    fY[i] = y[i];
    fZ[i] = z[i];
  }
}

void Tron::TPolyBox3D::Set(Double_t rx, Double_t ry, Double_t rz,
                           Double_t px, Double_t py, Double_t pz) {
  fX[0] = fX[3] = fX[4] = fX[7] = rx - px;
  fX[1] = fX[2] = fX[5] = fX[6] = rx + px;

  fY[0] = fY[1] = fY[2] = fY[3] = ry - py;
  fY[4] = fY[5] = fY[6] = fY[7] = ry + py;

  fZ[0] = fZ[1] = fZ[4] = fZ[5] = rz - pz;
  fZ[2] = fZ[3] = fZ[6] = fZ[7] = rz + pz;
}


//--------------------------------------------------
void Tron::TPolyBox3D::Update() {
  TPolyLine3D* surfaces[4] = {};

  Int_t i = 0, n = sizeof(surfaces) / sizeof(*surfaces);
  if (i < n) {
    {
      // First process
      surfaces[i] = (TPolyLine3D*)fList->First();
      goto COMMON_PROCESS;
    }
    for (; i < n; ++i) {
      // Second or later processes
      surfaces[i] = (TPolyLine3D*)fList->After(surfaces[i - 1]);
    COMMON_PROCESS:
      if (!surfaces[i]) {
        fList->Add(surfaces[i] = new TPolyLine3D(4));
      }
    }
  }

  for (Int_t i = 0; i < n; ++i) {
    TPolyLine3D* surface = surfaces[i];
    for (Int_t j : { 0, 1, 5, 4 }) {
      const Int_t k = (j + i) % 4 + (j / 4) * 4;
      surface->SetPoint(j, fX[k], fY[k], fZ[k]);
    }
  }
}


//--------------------------------------------------
void Tron::TPolyBox3D::Draw(const Char_t* option) {
  Update();
  fList->Draw(option);
}


//--------------------------------------------------
void Tron::TPolyBox3D::SetLineColor(Color_t color) {
  TIter it(fList);
  std::for_each(it.Begin(), it.End(), [&](TObject* obj) {
      auto surface = (TPolyLine*)obj;
      surface->SetLineColor(color);
    });
}

void Tron::TPolyBox3D::SetLineStyle(Style_t style) {
  TIter it(fList);
  std::for_each(it.Begin(), it.End(), [&](TObject* obj) {
      auto surface = (TPolyLine*)obj;
      surface->SetLineStyle(style);
    });
}

void Tron::TPolyBox3D::SetLineWidth(Width_t width) {
  TIter it(fList);
  std::for_each(it.Begin(), it.End(), [&](TObject* obj) {
      auto surface = (TPolyLine*)obj;
      surface->SetLineWidth(width);
    });
}


//--------------------------------------------------
Color_t Tron::TPolyBox3D::GetLineColor() {
  TPolyLine3D* surface = (TPolyLine3D*)fList->First();
  return surface ? surface->GetLineColor() : 0;
}

Style_t Tron::TPolyBox3D::GetLineStyle() {
  TPolyLine3D* surface = (TPolyLine3D*)fList->First();
  return surface ? surface->GetLineStyle() : 0;
}

Width_t Tron::TPolyBox3D::GetLineWidth() {
  TPolyLine3D* surface = (TPolyLine3D*)fList->First();
  return surface ? surface->GetLineWidth() : 0;
}


//--------------------------------------------------
Tron::TPolyBox2D* Tron::TPolyBox3D::Profile(Option_t* option) {
  std::string opt = String::ToLower(option);

  TPolyBox2D* profile = nullptr;
  if (opt.empty()) {
    profile = new TPolyBox2D(fX, fY);
  } else if (opt == "xy") {
    profile = new TPolyBox2D(fX, fY);
  } else if (opt == "yz") {
    profile = new TPolyBox2D(fY, fZ);
  } else if (opt == "zx") {
    profile = new TPolyBox2D(fY, fX);
  } else if (opt == "xz") {
    profile = new TPolyBox2D(fX, fZ);
  } else if (opt == "zy") {
    profile = new TPolyBox2D(fZ, fY);
  } else if (opt == "yx") {
    profile = new TPolyBox2D(fY, fX);
  } else {
    profile = new TPolyBox2D(fX, fY);
  }

  if (profile) {
    profile->SetLineColor(GetLineColor());
    profile->SetLineStyle(GetLineStyle());
    profile->SetLineWidth(GetLineWidth());
  }
  
  return profile;
}
