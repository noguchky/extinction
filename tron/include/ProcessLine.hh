#ifndef Tron_ProcessLine_hh
#define Tron_ProcessLine_hh

#include "PlotHelper.hh"
#include "StyleHelper.hh"

namespace Tron {

  namespace ProcessLine {
    namespace PH = PlotHelper;
    namespace SH = StyleHelper;
  }
  namespace PL = ProcessLine;

  namespace ProcessLine {

    template <class Object_t>
    class Processor {
    private:
      Object_t* fObject;

      template <typename T>
      struct TypeToObject {
        using type = Object_t;
      };

      template <Bool_t V>
      using EnableIf = typename std::enable_if<V, std::nullptr_t>::type;

      template <typename T, typename S>
      using IsBaseOf = typename std::is_base_of<S, typename TypeToObject<T>::type>;

    public:
      Processor(Object_t*& object)
        : fObject(object) {
      }

      Object_t* operator->() const {
        return fObject;
      }

      explicit operator Object_t*&() {
        return fObject;
      }

      explicit operator const Object_t*&() const {
        return fObject;
      }

      inline Processor& Push(Object_t*& obj) {
        obj = fObject;
        return *this;
      }

      template <typename Predicate_t, typename... Args_t>
      inline Processor& Call(Predicate_t predicate, Args_t... args) {
        if (fObject) {
          predicate(fObject, args...);
        }
        return *this;
      }

      template <typename Iterable_t, typename Predicate_t>
      inline Processor& CallFor(Iterable_t values, Predicate_t predicate) {
        if (fObject) {
          for (auto&& value : values) {
            predicate(fObject, value);
          }
        }
        return *this;
      }

      template <typename Predicate_t, typename... Args_t>
      inline Processor& Execute(Predicate_t predicate, Args_t...  args) {
        if (fObject) {
          (fObject->*predicate)(args...);
        }
        return *this;
      }

      template <typename ChildObject_t,
                typename T = void,
                EnableIf<IsBaseOf<T, TObject>::value> = nullptr>
      inline Processor<ChildObject_t>& FindObject(const std::string& name) {
        if (fObject) {
          ChildObject_t* child = dynamic_cast<ChildObject_t*>(fObject->FindObject(name.data()));
          if (!child) {
            std::cout << "[error] object is not found, "
                      << "\"" << name << "\"" << std::endl;
          }
          return child;
        }
        return nullptr;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TObject>::value> = nullptr>
      inline Processor& Draw(const std::string& option = "") {
        if (fObject) {
          fObject->Draw(option.data());
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TNamed>::value> = nullptr>
      inline Processor& SetName(const std::string& name) {
        if (fObject) {
          fObject->SetName(name.data());
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TNamed>::value> = nullptr>
      inline Processor& SetTitle(const std::string& title) {
        if (fObject) {
          fObject->SetTitle(title.data());
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TFile>::value> = nullptr>
      inline Processor& Open(const std::string& filename,
                                   const std::string& option = "") {
        if (TFile* file = fObject) {
          file->Open(filename.data(), option.data());
          if (!file->IsOpen()) {
            std::cout << "[error] file is not opened, "
                      << "\"" << filename << "\"" << ", "
                      << "\"" << option   << "\"" << std::endl;
          }
        }
        return *this;
      }

      template <typename ChildObject_t,
                typename T = void,
                EnableIf<IsBaseOf<T, TFile>::value> = nullptr>
      inline Processor<ChildObject_t>& Get(const std::string& name) {
        if (TFile* file = fObject) {
          ChildObject_t* child = dynamic_cast<ChildObject_t*>(file->Get(name.data()));
          if (!child) {
            std::cout << "[error] object is not found, "
                      << "\"" << name << "\"" << std::endl;
          }
          return child;
        }
        return nullptr;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TH1>::value &&
                         not IsBaseOf<T, TH2>::value> = nullptr>
      inline Processor& ReadFile(const std::string& filename,
                                       const std::string& format = "",
                                       const PH::Func1_t& fx     = nullptr) {
        if (fObject) {
          PH::ReadFile(fObject, filename, format, fx);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TGraph>::value ||
                         IsBaseOf<T, TH2>::value> = nullptr>
      inline Processor& ReadFile(const std::string& filename,
                                       const std::string& format = "",
                                       const PH::Func2_t& fx     = nullptr,
                                       const PH::Func2_t& fy     = nullptr) {
        if (fObject) {
          PH::ReadFile(fObject, filename, format, fx, fy);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TGraph2D>::value> = nullptr>
      inline Processor& ReadFile(const std::string& filename,
                                       const std::string& format = "",
                                       const PH::Func3_t& fx     = nullptr,
                                       const PH::Func3_t& fy     = nullptr,
                                       const PH::Func3_t& fz     = nullptr) {
        if (fObject) {
          PH::ReadFile(fObject, filename, format, fx, fy, fz);
        }
        return *this;
      }
      template <typename T = void,
                EnableIf<IsBaseOf<T, TH1>::value> = nullptr>
      inline Processor& Scale(Double_t scale) {
        if (TH1* hist = fObject) {
          hist->Scale(scale);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TF1>::value> = nullptr>
      inline Processor& SetNpx(Int_t npx) {
        if (TF1* func = fObject) {
          func->SetNpx(npx);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TH1>::value> = nullptr>
      inline Processor& SetRangeX(Double_t xmin, Double_t xmax) {
        if (TH1* hist = fObject) {
          hist->GetXaxis()->SetRangeUser(xmin, xmax);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TGraph>::value> = nullptr>
      inline Processor& SetRangeX(Double_t xmin, Double_t xmax) {
        if (TGraph* graph = fObject) {
          graph->GetXaxis()->SetLimits(xmin, xmax);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TVirtualPad>::value> = nullptr>
      inline Processor& SetCanvasSquare(Double_t margin, SH::CanvasAlign align = SH::CanvasAlign::Center) {
        if (fObject) {
          SH::SetCanvasSquare(fObject, margin, align);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TAttMarker>::value> = nullptr>
      inline Processor& SetMarkerStyle(Color_t color, Style_t style, Size_t size) {
        if (fObject) {
          SH::SetMarkerStyle(fObject, color, style, size);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TAttLine>::value> = nullptr>
      inline Processor& SetLineStyle(Color_t color, Style_t style, Width_t width) {
        if (fObject) {
          SH::SetLineStyle(fObject, color, style, width);
        }
        return *this;
      }

      template <typename T = void,
                EnableIf<IsBaseOf<T, TAttFill>::value> = nullptr>
      inline Processor& SetFillStyle(Color_t color, Style_t style) {
        if (fObject) {
          SH::SetFillStyle(fObject, color, style);
        }
        return *this;
      }

    };

    template <typename Object_t>
    inline Processor<Object_t> From(Object_t* object) {
      return Processor<Object_t>(object);
    }

  }

}

#endif
