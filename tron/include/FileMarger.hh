#ifndef Tron_FileMarger_hh
#define Tron_FileMarger_hh

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <functional>
#include <algorithm>
#include "String.hh"

namespace Tron {

  class FileMarger {
  private:
    bool fVerbose = false;
    std::vector<std::string> fFilenames;
    unsigned int fHeaderLines;
    std::vector<std::string> fHeaders;
    std::vector<std::string> fLines;

  public:
    FileMarger() = default;

    void SetVerbose(bool verbose) {
      fVerbose = verbose;
    }

    void SetHeaderLines(int headerLines) {
      fHeaderLines = headerLines;
    }

    void Add(const std::string& filename) {
      fFilenames.push_back(filename);
    }

    void Load() {
      std::string buff;

      // Initialize
      fHeaders.clear();
      fLines.clear();

      // Read from files
      for (auto&& ifilename : fFilenames) {
        if (fVerbose) {
          std::cout << "Read from " << String::Wrap(ifilename) << std::endl;
        }

        // Open input filea
        std::ifstream ifile(ifilename);
        if (!ifile) {
          std::cout << "[warning] input file is not opened" << String::Wrap(ifilename) << std::endl;
          continue;
        }

        // Read/Skip header
        if (fHeaders.empty()) {
          for (unsigned int i = 0; i < fHeaderLines; ++i) {
            if (std::getline(ifile, buff)) {
              fHeaders.push_back(buff);
            } else {
              fHeaders.push_back("");
            }
          }
        } else {
          for (unsigned int i = 0; i < fHeaderLines; ++i) {
            std::getline(ifile, buff);
          }
        }

        // Read lines
        while (std::getline(ifile, buff)) {
          fLines.push_back(buff);
        }
      }

      // Shuffle
      if (fVerbose) {
        std::cout << "Shuffling..." << std::endl;
      }
      std::random_device seed_gen;
      std::mt19937 engine(seed_gen());
      std::shuffle(fLines.begin(), fLines.end(), engine);
    }

    void Output(const std::string& ofilename,
                const std::function<std::string(const std::string&)>& headerPredicate = nullptr,
                const std::function<std::string(const std::string&)>& linePredicate   = nullptr) {
      // Open output file
      std::ofstream ofile(ofilename);
      if (!ofile) {
        std::cout << "[error] output file is not opened, " << String::Wrap(ofilename) << std::endl;
        return;
      }

      Output(ofile, headerPredicate, linePredicate);
    }

    void Output(std::ofstream& ofile,
                const std::function<std::string(const std::string&)>& headerPredicate = nullptr,
                const std::function<std::string(const std::string&)>& linePredicate   = nullptr) {
      // Output to file
      if (fVerbose) {
        std::cout << "Outputing..." << std::endl;
      }
      if (headerPredicate) {
        for (auto&& header : fHeaders) {
          ofile << headerPredicate(header) << std::endl;
        }
      } else {
        for (auto&& header : fHeaders) {
          ofile << header << std::endl;
        }
      }
      if (linePredicate) {
        for (auto&& line : fLines) {
          ofile << linePredicate(line) << std::endl;
        }
      } else {
        for (auto&& line : fLines) {
          ofile << line << std::endl;
        }
      }
    }

  };

}

#endif
