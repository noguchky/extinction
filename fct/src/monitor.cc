#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <sys/time.h>
#include <csignal>
#include <cstring>
#include <ctime>
#include <cerrno>
#include <unistd.h>
#include <mutex>
#include <limits>
#include <sys/inotify.h>
#include <pthread.h>

#include "MonitorWindow.hh"

#include "ArgReader.hh"

namespace {
  bool           processing = false;
  bool           exiting    = false;
  std::mutex     filenameMutex;
  std::string    path;
  std::string    filename;
  int            fd, wd;
  Extinction::Fct::MonitorWindow* monitor;

  void TimerTickEventHandler(int) {
    static unsigned long msec_cnt = 0;

    msec_cnt++;
    if (exiting) {
      // Nothing to do

    } else if (gPad) {
      // Check processing
      if (processing) {
        return;
      }
      processing = true;

      // Get Filename
      std::string ifilename;
      {
        std::lock_guard<std::mutex> lock(filenameMutex);
        ifilename = filename;
        filename = "";
      }

      if (ifilename.empty()) {
        // File was not created
      } else {
        try {
          std::cout << "Modify plots for " << ifilename << std::endl;
          monitor->UpdatePlots(path + "/" + ifilename);
        } catch (...) {
          monitor->Terminate();
        }
      }

      processing = false;

    } else {
      std::cerr << "[info] canvas was closed" << std::endl;
      exiting = true;
      monitor->Terminate();
    }
  }

  void* InotiryEventListener(void*) {
    constexpr int BufferSize = 65536;

    char buffer[BufferSize];
    for (int aux = 0; !exiting;) {
      // Read events
      int ret = read(fd, buffer + aux, BufferSize - aux);
      if (ret == -1) {
        if (ret == -EINTR) {
          continue;
        }
        perror("read");
        exit(EXIT_FAILURE);
      }
      ret += aux;

      // Read size check
      if (ret < (int)sizeof(inotify_event)) {
        fprintf(stderr, "short of red bytes\n");
        exit(EXIT_FAILURE);
      }

      if (exiting) {
        break;
      }

      for (int i = 0; i < ret;) {
        inotify_event* inotify_p = (inotify_event*)(buffer + i);

        if (ret < i + (int)offsetof(inotify_event, name)) {
          // Get tail at next read
          aux = ret - i;
          memmove(buffer, buffer + i, aux);
          break;
        }

        int size = sizeof(inotify_event) + inotify_p->len;
        if (ret < i + size) {
          // Get tail at next read
          aux = ret - i;
          memmove(buffer, buffer + i, aux);
          break;
        }

        // File opened for writing was closed (*).
        if (inotify_p->mask & IN_CLOSE_WRITE ||
            inotify_p->mask & IN_MOVED_TO) {
          std::cerr << "[info] file was closed \"" << inotify_p->name << "\"" << std::endl;
          std::lock_guard<std::mutex> lock(filenameMutex);
          filename = inotify_p->name;
        }
        i += size;
      }
    }

    // Remove monitor
    {
      int ret = inotify_rm_watch(fd, wd);
      if (ret == -1) {
        perror("inotify_rm_watch");
        exit(EXIT_FAILURE);
      }
    }

    close(fd);

    return 0;
  }
}

Int_t main(Int_t argc, Char_t** argv) {
  Tron::ArgReader* args = new Tron::ArgReader(argv[0]);
  args->AddArg<std::string>("Directory" ,                    "A rawdata directory");
  args->AddArg<std::string>("Offset"    ,                    "A offset filename format");
  args->AddOpt             ("Help"      , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string directory = args->GetValue("Directory");
  const std::string ffilename = args->GetValue("Offset");
  path = directory;

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize monitor window" << std::endl;
  monitor = new Extinction::Fct::MonitorWindow();

  if (monitor->LoadOffset(ffilename) == 0) {
    std::cerr << "[error] offset file was not found" << std::endl;
    exit(1);
  }

  monitor->InitializeWindow();
  monitor->InitializePlots();
  monitor->DrawPlots();

  std::cout << "--- Set signal handler" << std::endl;
  struct sigaction action;
  memset(&action, 0, sizeof(action));

  action.sa_handler = TimerTickEventHandler;
  action.sa_flags = SA_RESTART;
  sigemptyset(&action.sa_mask);
  if (sigaction(SIGALRM, &action, nullptr) < 0) {
    perror("sigaction error");
    exit(1);
  }

  std::cout << "--- Set intarval timer" << std::endl;
  itimerval timer;
  timer.it_value.tv_sec = 0;
  timer.it_value.tv_usec = 20000; // = 20 ms
  timer.it_interval.tv_sec = 0;
  timer.it_interval.tv_usec = 20000;
  if (setitimer(ITIMER_REAL, &timer, nullptr) < 0) {
    // perror("setitimer error");
    std::cerr << "[error] setitimer error" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Initialize inotify
  fd = inotify_init();
  if (fd == -1) {
    std::cerr << "[error] inotify initialize error" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Add monitor directory
  wd = inotify_add_watch(fd, directory.data(),  IN_ALL_EVENTS);
  if (wd < 0) {
    std::cerr << "[error] directory was not opened, " << directory << std::endl;
    exit(EXIT_FAILURE);
  }

  // Set write monitor
  std::cout << "--- Set inotiry event listener" << std::endl;
  pthread_t pthread;
  if (pthread_create(&pthread, nullptr, &InotiryEventListener, nullptr)) {
    std::cerr << "[error] pthread_create error" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "--- Start to run" << std::endl;
  monitor->Run();
  delete monitor;

  return 0;
}
