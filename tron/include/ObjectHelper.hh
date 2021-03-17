#ifndef Tron_ObjectHelper_hh
#define Tron_ObjectHelper_hh

#include <iostream>
#include <string>
#include "TROOT.h"
#include "TParameter.h"
#include "TFile.h"

namespace Tron {

  namespace ObjectHelper {

    inline std::string FindUnusedObjectName(const std::string& name) {
      std::string objectName = name;
      for (Int_t i = 1; gROOT->FindObject(objectName.data()); ++i) {
        objectName = Form("%s_%d", name.data(), i);
      }
      return objectName;
    }

    template <typename Value_t>
    inline void WriteValue(Value_t value, const std::string& name) {
      if (gFile) {
        TParameter<Value_t> v(name.data(), value);
        v.Write();
      } else {
        std::cout << "[warning] file is not opened" << std::endl;
      }
    }

    template <typename Value_t>
    inline Value_t ReadValue(const std::string& name) {
      Value_t result = Value_t();
      if (gFile) {
        if (auto obj = gFile->Get(name.data())) {
          if (auto value = dynamic_cast<TParameter<Value_t>*>(obj)) {
            result = value->GetVal();
          } else {
            std::cout << "[warning] value cannot cast" << std::endl;
          }
          delete obj;
        } else {
          std::cout << "[warning] value is not found" << std::endl;
        }
      } else {
        std::cout << "[warning] file is not opened" << std::endl;
      }
      return result;
    }

    inline TObject* ReadFrom(const std::string& filename,
                             const std::string& objectname) {
      // Open file
      TFile* file = new TFile(filename.data());
      if (!file->IsOpen()) {
        std::cout << "[error] file is not opened, " << "\"" << filename << "\"" << std::endl;
        return nullptr;
      }
      // Get object
      TObject* object = file->Get(objectname.data());
      // Close file
      file->Clone();
      delete file;
      return object;
    }

    inline std::vector<TObject*> ReadFrom(const std::string& filename,
                                          const std::initializer_list<std::string>& objectnames) {
      std::vector<TObject*> objects(objectnames.size(), nullptr);
      // Open file
      TFile* file = new TFile(filename.data());
      if (!file->IsOpen()) {
        std::cout << "[error] file is not opened, " << "\"" << filename << "\"" << std::endl;
        return objects;
      }
      // Get objects
      Int_t i = 0;
      for (auto&& name : objectnames) {
        objects[i++] = file->Get(name.data());
      }
      // Close file
      file->Clone();
      delete file;
      return objects;
    }

  }

}

#endif
