#include "Console.hh"
#include <unistd.h>
#include <fcntl.h>

#ifdef __APPLE__
#include <termios.h>
#else
#include <termio.h>
#endif

std::string Tron::Console::Prompt(const std::string& appName, const std::string& funcName) {
  std::stringstream prompt;
  if (funcName.empty()) {
    prompt << appName                    << "> ";
  } else {
    prompt << appName << " " << funcName << "> ";
  }
  return prompt.str();
}

std::string Tron::Console::Prompt(const std::string& appName, const std::string& hostName, const std::string& funcName) {
  std::stringstream prompt;
  if (funcName.empty()) {
    prompt << appName << "@" << hostName                    << "> ";
  } else {
    prompt << appName << "@" << hostName << " " << funcName << "> ";
  }
  return prompt.str();
}

Bool_t Tron::Console::Confirm(const std::string& message, Bool_t defaultAnswer) {
  std::string dmp;
  if (defaultAnswer) {
    std::cout << message << " [Y/n] " << std::flush;
    return !(std::cin >> dmp) || !(dmp == "n" || dmp == "N");
  } else {
    std::cout << message << " [y/N] " << std::flush;
    return  (std::cin >> dmp) &&  (dmp == "y" || dmp == "Y");
  }
}

void Tron::Console::Wait(const std::string& message) {
  std::cout << message << std::flush;
  while (std::getchar() != '\n');
}

bool Tron::Console::KeyboardHit() {
  termios oldt, newt;
  int ch;
  int oldf;

  // Get Current Attribute
  tcgetattr(STDIN_FILENO, &oldt);

  // Set New Attribute
  newt = oldt;
  newt.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &newt);

  // Get Current Control
  oldf = fcntl(STDIN_FILENO, F_GETFL, 0);

  // Set New Control
  fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);

  // Get Chara
  ch = getchar();

  // Set Old Attribute
  tcsetattr(STDIN_FILENO, TCSANOW, &oldt);

  // Set Old Attribute
  fcntl(STDIN_FILENO, F_SETFL, oldf);

  // Check Chara
  if (ch != EOF) {
    ungetc(ch, stdin);
    return 1;
  }

  return 0;}

