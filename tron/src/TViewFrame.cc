#include <string>
#include "TViewFrame.hh"
#include "String.hh"

Tron::TViewFrame2D* Tron::TViewFrame3D::Projection(const Char_t* name, Option_t* option) {
  struct {
    TAxis*   axis;
    Double_t min;
    Double_t max;
  } xaxis, yaxis, zaxis, *pxaxis, *pyaxis;
  
  xaxis.axis = this->GetXaxis();
  yaxis.axis = this->GetYaxis();
  zaxis.axis = this->GetZaxis();

  Double_t* px = this->GetX();
  Double_t* py = this->GetY();
  Double_t* pz = this->GetZ();
  xaxis.min = *px;
  yaxis.min = *py;
  zaxis.min = *pz;    
  xaxis.max = *(px + 1);
  yaxis.max = *(py + 1);
  zaxis.max = *(pz + 1);
  
  std::string opt = String::ToLower(option);
  if (opt.empty()) {
    pxaxis = &xaxis, pyaxis = &yaxis;
  } else if (opt == "xy") {
    pxaxis = &xaxis, pyaxis = &yaxis;
  } else if (opt == "xz") {
    pxaxis = &xaxis, pyaxis = &zaxis;
  } else if (opt == "yx") {
    pxaxis = &yaxis, pyaxis = &xaxis;
  } else if (opt == "yz") {
    pxaxis = &yaxis, pyaxis = &zaxis;
  } else if (opt == "zx") {
    pxaxis = &zaxis, pyaxis = &xaxis;
  } else if (opt == "zy") {
    pxaxis = &zaxis, pyaxis = &yaxis;
  } else {
    pxaxis = &xaxis, pyaxis = &yaxis;
  }
  
  auto newOne = new TViewFrame2D(name,
                                 this->GetTitle(),
                                 pxaxis->min, pxaxis->max,
                                 pyaxis->min, pyaxis->max);
  newOne->GetXaxis()->SetTitle(pxaxis->axis->GetTitle());
  newOne->GetYaxis()->SetTitle(pyaxis->axis->GetTitle());

  return newOne;
}
