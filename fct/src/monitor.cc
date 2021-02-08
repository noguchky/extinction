#include <iostream>
#include <fstream>
#include <iomanip>
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

#include "Linq.hh"
#include "ArgReader.hh"
#include "ConfReader.hh"
#include "ScopeSubstituter.hh"

namespace {
  bool                                idle = true;
  std::mutex                          busyMutex;
  std::mutex                          filenamesMutex;
  std::map<int/*board*/, std::string> filenames;
  Extinction::Fct::MonitorWindow*     monitor;

  void TimerTickEventHandler(int) {
    if (monitor->IsTerminated()) {
      // Nothing to do

    } else if (monitor->IsClosed()) {
      std::cerr << "[info] canvas was closed" << std::endl;
      monitor->Terminate();

    } else if (idle) {
      Tron::ScopeSubstituter<Bool_t> lock(idle, false);

      // Get Filenames
      std::map<int/*board*/, std::string> ifilenames;
      {
        std::lock_guard<std::mutex> lock(filenamesMutex);
        for (auto&& pair : filenames) {
          const int         board    = pair.first;
          const std::string filename = pair.second;
          if (filename.empty()) {
            return;
          }
          ifilenames[board] = filename;
        }
        for (auto&& pair : filenames) {
          pair.second = "";
        }
      }

      try {
        std::cout << "Modify plots for " << ifilenames.begin()->second << " and so on" << std::endl;
        monitor->UpdatePlots(ifilenames);
      } catch (...) {
        monitor->Terminate();
      }
    }
  }

  struct InotifyEventArgs {
    int         fd;
    int         wd;
    int         board;
    std::string path;
  };

  void* InotifyEventListener(void* arg) {
    constexpr int BufferSize = 65536;

    const InotifyEventArgs* eventArg = static_cast<InotifyEventArgs*>(arg);
    const int         fd    = eventArg->fd;
    const int         wd    = eventArg->wd;
    const int         board = eventArg->board;
    const std::string path  = eventArg->path;
    // std::cout << "InotifyEventListener" << std::endl
    //           << "  fd     = " << fd    << std::endl
    //           << "  wd     = " << wd    << std::endl
    //           << "  board  = " << board << std::endl
    //           << "  path   = " << path  << std::endl;

    char buffer[BufferSize];
    for (int aux = 0; !monitor->IsTerminated();) {
      // Read events
      int ret = read(fd, buffer + aux, BufferSize - aux);
      if (ret == -1) {
        if (errno == EINTR) {
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

      if (monitor->IsTerminated()) {
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

        if (inotify_p->mask & IN_CLOSE_WRITE ||
            inotify_p->mask & IN_MOVED_TO) {
          // std::cerr << "[info] file was closed \"" << inotify_p->name << "\"" << std::endl;
          std::lock_guard<std::mutex> lock(filenamesMutex);
          filenames[board] = (path + "/" + inotify_p->name);
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
  args->AddArg<std::string>("ConfFilename"  ,                    "A configure filename");
  args->AddOpt<Int_t>      ("Boards"        , 'b', "boards"    , "Set comma-separated input boards");
  args->AddOpt<std::string>("Filenames"     , 'f', "files"     , "Set comma-separated input filenames (no wait to be written)");
  args->AddOpt<std::string>("ChannelAlign"  , 'a', "align"     , "Set channel align (raw/projection)", "raw");
  args->AddOpt<Int_t>      ("MonitorMode"   , 'm', "mode"      , "Set monitor mode (channel/offset/coincidence)", "coincidence");
  args->AddOpt<Int_t>      ("MonitorChannel", 'c', "channel"   , "Set monitor channel by global channel", "-1");
  args->AddOpt             ("Help"          , 'h', "help"      , "Show usage");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    std::cout << "Global channel map" << std::endl
              << "  " << std::setw(3) << Extinction::ExtinctionDetector::GlobalChannelOffset << "~  extinction detector" << std::endl
              << "  " << std::setw(3) << Extinction::Hodoscope         ::GlobalChannelOffset << "~  (old) hodoscope" << std::endl
              << "  " << std::setw(3) << Extinction::TimingCounter     ::GlobalChannelOffset << "~  timing counter" << std::endl
              << "  " << std::setw(3) << Extinction::BeamlineHodoscope ::GlobalChannelOffset << "~  beamline hodoscope" << std::endl
              << "  " << std::setw(3) << Extinction::MrRf              ::GlobalChannelOffset << "~  MR RF" << std::endl
              << "  " << std::setw(3) << Extinction::MrP3              ::GlobalChannelOffset << "~  MR P3" << std::endl
              << "  " << std::setw(3) << Extinction::MrSync            ::GlobalChannelOffset << "~  MR Sync" << std::endl
      ;
    return 0;
  }

  const std::string confFilename   = args->GetValue("ConfFilename");
  const std::string channelAlign   = args->GetValue("ChannelAlign");
  const std::string monitorMode    = args->GetValue("MonitorMode");
  const Int_t       monitorChannel = args->GetValue<Int_t>("MonitorChannel");
  const Bool_t      isBoardsSet    = args->IsSet("Boards");
  const Bool_t      isFilenamesSet = args->IsSet("Filenames");

  const Bool_t isWatchingMode = !isFilenamesSet;

  const std::vector<std::string> argBoards    = Tron::String::Split(args->GetValue("Boards"   ), ",");
  const std::vector<std::string> argFilenames = Tron::String::Split(args->GetValue("Filenames"), ",");

  std::cout << "--- Load configure" << std::endl;
  Tron::ConfReader* conf = new Tron::ConfReader(confFilename);
  if (!conf->IsOpen()) {
    std::cerr << " [error] config file is not opened, " << confFilename << std::endl;
    return 1;
  }
  conf->ShowContents();

  const Int_t windowWidth  = conf->GetValue<Int_t>("Window.Width" );
  const Int_t windowHeight = conf->GetValue<Int_t>("Window.Height");

  const std::vector<int> boards = isBoardsSet ? Tron::Linq::From(argBoards)
    .Select([](const std::string& board) { return Tron::String::Convert<Int_t>(board); })
    .ToVector() :
    conf->GetValues<Int_t>("Boards");

  const std::map<int, std::string> paths = Tron::Linq::From(boards)
    .ToMap([&](int board) { return board; },
           [&](int board) { return conf->GetValue(Form("Path.%d", board)); });
  const std::map<Int_t, Double_t> timePerTdc = Tron::Linq::From(boards)
    .ToMap([&](int board) { return board; },
           [&](int board) { return conf->GetValue<Double_t>(Form("TimePerTdc.%d", board)) * Extinction::nsec; });
  const std::map<Int_t, Double_t> mrSyncInterval = Tron::Linq::From(boards)
    .ToMap([&](int board) { return board; },
           [&](int board) { return conf->GetValue<Double_t>(Form("MrSyncInterval.%d", board)); });

  if (isWatchingMode) {
    filenames = Tron::Linq::From(boards)
      .ToMap([&](Int_t board) { return board;         },
             [&](Int_t      ) { return std::string(); });
  } else {
    if (argFilenames.size() != boards.size()) {
      std::cerr << "[error] # of filenames is not equal to # of boards" << std::endl;
      return 1;
    }
    filenames = Tron::Linq::From(boards)
      .Join (argFilenames.begin())
      .ToMap([&](std::pair<Int_t, std::string> pair) { return pair.first;  },
             [&](std::pair<Int_t, std::string> pair) { return pair.second; });
  }

  const std::string ffilename = conf->GetValue("Offset");
  if (ffilename.empty()) {
    std::cerr << "[error] offset was not configured" << std::endl;
    return 1;
  }

  Extinction::Fct::ChannelMapWithBoard::Load(conf, boards);

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

  if        (channelAlign == "raw"       ) {
    monitor->SetChannelAlign(Extinction::Fct::ChannelAlign::Raw);
  } else if (channelAlign == "projection") {
    monitor->SetChannelAlign(Extinction::Fct::ChannelAlign::Projection);
  } else {
    std::cerr << "[warning] invalid align, set to raw align" << std::endl;
    monitor->SetChannelAlign(Extinction::Fct::ChannelAlign::Raw);
  }
  if        (monitorMode == "channel") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorMode::Channel, monitorChannel);
  } else if (monitorMode == "offset") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorMode::Offset, monitorChannel);
  } else if (monitorMode == "coincidence") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorMode::Coincidence);
  } else {
    std::cerr << "[warning] invalid monitor mode, set to coincidence" << std::endl;
    monitor->SetMonitorMode(Extinction::Fct::MonitorMode::Coincidence);
  }

  monitor->SetTimePerTdc(timePerTdc);
  monitor->SetMrSyncInterval(mrSyncInterval);
  monitor->InitializeWindow(windowWidth ? windowWidth : 1600, windowHeight ? windowHeight : 1200);
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
    std::cerr << "[error] setitimer error" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::map<pthread_t*, InotifyEventArgs*> threads;
  if (isWatchingMode) {
    std::cout << "--- Register directory to inotify" << std::endl;
    for (auto&& pair : paths) {
      const int         board = pair.first;
      const std::string path  = pair.second;
    
      // Initialize inotify
      int fd = inotify_init();
      if (fd == -1) {
        std::cerr << "[error] inotify initialize error" << std::endl;
        exit(EXIT_FAILURE);
      }

      // Add monitor directory
      int wd = inotify_add_watch(fd, path.data(), IN_ALL_EVENTS);
      if (wd < 0) {
        std::cerr << "[error] directory was not opened, " << path << std::endl;
        exit(EXIT_FAILURE);
      }

      // Create argument
      auto eventArg = new InotifyEventArgs();
      eventArg->fd    = fd;
      eventArg->wd    = wd;
      eventArg->board = board;
      eventArg->path  = path;

      // Set write monitor
      pthread_t* pthread = new pthread_t();
      if (pthread_create(pthread, nullptr, &InotifyEventListener, (void*)eventArg)) {
        std::cerr << "[error] pthread_create error" << std::endl;
        exit(EXIT_FAILURE);
      }

      threads[pthread] = eventArg;
    }
  }

  std::cout << "--- Start to run" << std::endl;
  monitor->Run();
  delete monitor;

  return 0;
}
