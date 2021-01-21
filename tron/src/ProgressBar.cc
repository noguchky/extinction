#include <string>
#include "ProgressBar.hh"
#include "String.hh"

const double Tron::ProgressBar::kInterval = 0.01;
const double Tron::ProgressBar::kTotal    = 100.0;
const int    Tron::ProgressBar::kWidth    =  50;
const char   Tron::ProgressBar::kDone     = '=';
const char   Tron::ProgressBar::kNotYet   = '-';

void Tron::ProgressBar::Show(double progress) {
  fFormat = ShowType::Real;
  const double ratio = progress / fTotal;
  if (progress == 0.0 || progress == fTotal || ratio - fPrevious > fInterval) {
    fPrevious = ratio;
    const int n1 = ratio * fWidth;
    const int n2 = fWidth - n1;
    const std::string done   = std::string(n1, fDone);
    const std::string notYet = std::string(n2, fNotYet);
    const std::string desc   = String::Format("%5.1lf%% (%3.2lf / %3.2lf)", 100.0 * ratio, progress, fTotal);
    std::cout << "[" << done << notYet << "]   " << desc << "    \r" << std::flush;
  }
}

void Tron::ProgressBar::Show(int progress) {
  fFormat = ShowType::Integer;
  const double ratio = (double)progress / fTotal;
  if (progress == 0 || progress == fTotal || ratio - fPrevious > fInterval) {
    fPrevious = ratio;
    const int n1 = ratio * fWidth;
    const int n2 = fWidth - n1;
    const std::string done   = std::string(n1, fDone);
    const std::string notYet = std::string(n2, fNotYet);
    const std::string desc   = String::Format("%5.1lf%% (%d / %d)", 100.0 * ratio, progress, (int)fTotal);
    std::cout << "[" << done << notYet << "]   " << desc << "    \r" << std::flush;
  }
}

void Tron::ProgressBar::Terminate() {
  if (fPrevious != 1.0 && fFormat != ShowType::None) {
    const double ratio = 1.0;
    fPrevious = ratio;
    const int n1 = ratio * fWidth;
    const int n2 = fWidth - n1;
    const std::string done   = std::string(n1, fDone);
    const std::string notYet = std::string(n2, fNotYet);
    const std::string desc   =
      fFormat == ShowType::Real    ? String::Format("%5.1lf%% (%3.2lf / %3.2lf)", 100.0 * ratio, fTotal, fTotal) :
      fFormat == ShowType::Integer ? String::Format("%5.1lf%% (%d / %d)", 100.0 * ratio, (int)fTotal, (int)fTotal) : "";
    std::cout << "[" << done << notYet << "]   " << desc << "    \r" << std::flush;
  }
  std::cout << std::endl;
}
