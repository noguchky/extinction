#ifndef Tron_TPolyBox_hh
#define Tron_TPolyBox_hh

#include <TList.h>

namespace Tron {

  class TPolyBox2D {
  private:
    TList*   fList;
    Double_t fX[8];
    Double_t fY[8];
  
  public:
    TPolyBox2D();
    TPolyBox2D(const Double_t* x, const Double_t* y);
    TPolyBox2D(Double_t rx, Double_t ry, Double_t px, Double_t py);
    virtual ~TPolyBox2D();

    void    Set(const Double_t* x, const Double_t* y);
    void    Set(Double_t rx, Double_t ry, Double_t px, Double_t py);
    void    Update();
    void    Draw(Option_t* option = "");

    void    SetLineColor(Color_t color);
    void    SetLineStyle(Style_t style);
    void    SetLineWidth(Width_t width);

    Color_t GetLineColor();
    Style_t GetLineStyle();
    Width_t GetLineWidth();
  };

  class TPolyBox3D {
  private:
    TList*   fList;
    Double_t fX[8];
    Double_t fY[8];
    Double_t fZ[8];

  public:
    TPolyBox3D();
    TPolyBox3D(const Double_t* x, const Double_t* y, const Double_t* z);
    TPolyBox3D(Double_t rx, Double_t ry, Double_t rz,
               Double_t px, Double_t py, Double_t pz);
    virtual ~TPolyBox3D();

    void        Set(const Double_t* x, const Double_t* y, const Double_t* z);
    void        Set(Double_t rx, Double_t ry, Double_t rz,
                    Double_t px, Double_t py, Double_t pz);
    void        Update();
    void        Draw(Option_t* option = "");

    void        SetLineColor(Color_t color);
    void        SetLineStyle(Style_t style);
    void        SetLineWidth(Width_t width);

    Color_t     GetLineColor();
    Style_t     GetLineStyle();
    Width_t     GetLineWidth();

    TPolyBox2D* Profile(Option_t* option="xy");
  };

}

#endif
