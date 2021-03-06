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
#include <dirent.h>
#include <pthread.h>
#include <regex>

#include "MonitorWindow.hh"

#include "Linq.hh"
#include "String.hh"
#include "ArgReader.hh"
#include "ConfReader.hh"
#include "ScopeSubstituter.hh"

namespace {

  auto parser =
    [](const std::string& filename) {
      static const std::regex pattern(R"(fct_id\d\d\d\d_(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d))");

      std::cmatch match;
      if (std::regex_search(filename.data(), match, pattern)) {
        const Int_t year   = Tron::String::Convert<Int_t>(match[1].str());
        const Int_t month  = Tron::String::Convert<Int_t>(match[2].str());
        const Int_t day    = Tron::String::Convert<Int_t>(match[3].str());
        const Int_t hour   = Tron::String::Convert<Int_t>(match[4].str());
        const Int_t minute = Tron::String::Convert<Int_t>(match[5].str());
        const Int_t second = Tron::String::Convert<Int_t>(match[6].str());
        return TDatime(year, month, day, hour, minute, second);
      }

      return TDatime(0U);
    };

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
        monitor->UpdatePlots(ifilenames, parser);
      } catch (...) {
        monitor->Terminate();
      }
    }
  }

  void AddSubDirectory(int fd, std::map<int, std::string>& paths, const std::string& path) {
    DIR*    dir;
    dirent* ent;
    if ((dir = opendir(path.data())) != nullptr) {
      // print all the files and directories within directory
      while ((ent = readdir(dir)) != nullptr) {
        const std::string dirname = ent->d_name;
        if (ent->d_type == DT_DIR && dirname != "." && dirname != "..") {
          const std::string newpath = path + "/" + ent->d_name;
          const int newwd = inotify_add_watch(fd, newpath.data(),  IN_ALL_EVENTS);
          paths[newwd] = newpath;
          std::cout << newpath << std::endl;
          AddSubDirectory(fd, paths, newpath);
        }
      }
      closedir(dir);
    } else {
      // could not open directory
      perror("could not open directory");
    }
  }

  struct InotifyEventArgs {
    int                        board;
    int                        fd;
    std::map<int, std::string> paths;
    // std::vector<int> wds; 
    // std::string      path;
  };

  void* InotifyEventListener(void* arg) {
    constexpr int BufferSize = 65536;

    const InotifyEventArgs*    eventArg = static_cast<InotifyEventArgs*>(arg);
    const int                  fd       = eventArg->fd;
    const int                  board    = eventArg->board;
    std::map<int, std::string> paths    = eventArg->paths;

    bool is_created = false;
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
          if (Tron::String::EndWith(inotify_p->name, ".dat")) {
            std::cerr << "[info] file was closed \"" << inotify_p->name << "\"" << std::endl;
            std::lock_guard<std::mutex> lock(filenamesMutex);
            if (paths.find(inotify_p->wd) != paths.end()) {
              const std::string path = paths[inotify_p->wd];
              filenames[board] = (path + "/" + inotify_p->name);
            } else {
              std::cout << "[warning] unknown wd of inotify" << std::endl;
            }
          }
        }

        if (is_created &&
            inotify_p->mask & IN_ISDIR) {
          const std::string path    = paths[inotify_p->wd];
          const std::string dirname = inotify_p->name;
          const std::string newpath = path + "/" + dirname;
          int newwd = inotify_add_watch(fd, newpath.data(), IN_ALL_EVENTS);
          paths[newwd] = newpath;
        }

        if (inotify_p->mask & IN_CREATE) {
          is_created = true;
        } else {
          is_created = false;
        }
        i += size;
      }
    }

    // Remove monitor
    for (auto& pair : paths) {
      const int wd = pair.first;
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
  args->AddArg<std::string>("ConfFilename"    ,                  "A configure filename");
  args->AddOpt<Int_t>      ("Boards"          , 'b', "boards"  , "Set comma-separated input boards");
  args->AddOpt<std::string>("Filenames"       , 'f', "files"   , "Set comma-separated input filenames (no wait to be written)");
  args->AddOpt<std::string>("ChannelAlign"    , 'a', "align"   , "Set channel align (raw/projection)", "raw");
  args->AddOpt<std::string>("MonitorMode"     , 'm', "mode"    , "Set monitor mode (channel/ext/hod/offset/coincidence)", "coincidence");
  args->AddOpt<Int_t>      ("MonitorChannel"  , 'c', "channel" , "Set monitor channel by global channel", "-1");
  args->AddOpt<std::string>("OffsetFilename"  , 'o', "offset"  , "Set output offset filename", "");
  args->AddOpt<std::string>("IntervalFilename", 'i', "interval", "Set output mr-sync interval filename", "");
  args->AddOpt             ("Sum"             , 's', "sum"     , "Sum coincidence hist");
  args->AddOpt             ("Leakage"         , 'l', "leakage" , "Coincidence only leakage");
  args->AddOpt             ("Help"            , 'h', "help"    , "Show usage");

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

  const std::string confFilename     = args->GetValue("ConfFilename");
  const std::string channelAlign     = args->GetValue("ChannelAlign");
  const std::string monitorMode      = args->GetValue("MonitorMode");
  const Int_t       monitorChannel   = args->GetValue<Int_t>("MonitorChannel");
  const Bool_t      isBoardsSet      = args->IsSet("Boards");
  const Bool_t      isFilenamesSet   = args->IsSet("Filenames");
  const std::string offsetFilename   = args->GetValue("OffsetFilename");
  const std::string intervalFilename = args->GetValue("IntervalFilename");
  const Bool_t      accumulate       = args->IsSet("Sum");
  const Bool_t      leakageMode      = args->IsSet("Leakage");

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

  Extinction::Fct::MonitorWindow::PlotsProfiles profile;
  profile.TimeInSpill   .NbinsX   = conf->GetValue<Double_t>("TimeInSpill.NbinsX");
  profile.TimeInSpill   .Xmin     = conf->GetValue<Double_t>("TimeInSpill.Xmin");
  profile.TimeInSpill   .Xmax     = conf->GetValue<Double_t>("TimeInSpill.Xmax");
  profile.TimeInSpill   .Unit     = conf->GetValue<Double_t>("TimeInSpill.Unit");
  profile.TimeInSync    .Xmin     = conf->GetValue<Double_t>("TimeInSync.Xmin");
  profile.TimeInSync    .Xmax     = conf->GetValue<Double_t>("TimeInSync.Xmax");
  profile.TimeInSync    .BinWidth = conf->GetValue<Double_t>("TimeInSync.BinWidth");
  profile.MrSyncInterval.Mean     = conf->GetValue<Double_t>("MrSyncInterval.Mean");
  profile.MrSyncInterval.Xmin     = conf->GetValue<Double_t>("MrSyncInterval.Xmin");
  profile.MrSyncInterval.Xmax     = conf->GetValue<Double_t>("MrSyncInterval.Xmax");
  profile.TimeDiff      .Xmin     = conf->GetValue<Double_t>("TimeDiff.Xmin");
  profile.TimeDiff      .Xmax     = conf->GetValue<Double_t>("TimeDiff.Xmax");

  const std::vector<int> boards = isBoardsSet ? Tron::Linq::From(argBoards)
    .Select([](const std::string& board) { return Tron::String::Convert<Int_t>(board); })
    .ToVector() :
    conf->GetValues<Int_t>("Boards");

  const std::map<int, std::string> paths = Tron::Linq::From(boards)
    .ToMap([&](int board) { return board; },
           [&](int board) { return conf->GetValue(Form("Path.%d", board)); });

  const Double_t    historyWidth     = conf->GetValue<Double_t   >("HistoryWidth");
  const Double_t    coinTimeWidth    = conf->GetValue<Double_t   >("CoinTimeWidth");
  const Long64_t    spillLimit       = conf->GetValue<Long64_t   >("SpillLimit");
  const std::size_t readBufferSize   = conf->GetValue<std::size_t>("ReadBufferSize");
  const std::size_t readBufferMargin = conf->GetValue<std::size_t>("ReadBufferMargin");

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

  Extinction::Fct::ChannelMapWithBoard::Load(conf, boards);

  std::cout << "--- Initialize style" << std::endl;
  gStyle->SetPalette(1);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetNdivisions(505, "Y");

  std::cout << "--- Initialize monitor window" << std::endl;
  monitor = new Extinction::Fct::MonitorWindow();

  std::map<Int_t, Double_t> timePerTdc;
  do {
    const std::string ifilename = conf->GetValue("TimePerTdc");
    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[warning] time resolution file was not opened, " << ifilename << std::endl;
      break;
    }

    std::string buff;
    Int_t board; Double_t resolution;
    while (std::getline(ifile, buff)) {
      std::stringstream line(buff);
      if (line >> board >> resolution) {
        timePerTdc[board] = resolution * Extinction::nsec;
      }
    }
  } while (false);

  std::map<Int_t, Double_t> mrSyncInterval;
  do {
    const std::string ifilename = conf->GetValue("MrSyncInterval");
    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[warning] MR Sync interval file was not opened, " << ifilename << std::endl;
      break;
    }

    std::string buff;
    Int_t board; Double_t interval;
    while (std::getline(ifile, buff)) {
      std::stringstream line(buff);
      if (line >> board >> interval) {
        mrSyncInterval[board] = interval;
      }
    }
  } while (false);

  Double_t bunchCenters[Extinction::SpillData::kNofBunches] = { 0 };
  Double_t bunchWidths [Extinction::SpillData::kNofBunches] = { 0 };
  do {
    const std::string ifilename = conf->GetValue("BunchProfile");
    std::ifstream ifile(ifilename);
    if (!ifile) {
      std::cout << "[warning] bunch profile file was not opened, " << ifilename << std::endl;
      break;
    }

    std::string buff;
    Int_t bunch; Double_t center, width;
    while (std::getline(ifile, buff)) {
      std::stringstream line(buff);
      if (line >> bunch >> buff >> center >> buff >> width) {
        bunchCenters[bunch] = center * Extinction::nsec;
     // bunchWidths [bunch] = width  * Extinction::nsec;
        bunchWidths [bunch] = 200.0  * Extinction::nsec;
      }
    }
  } while (false);

  {
    const std::string ifilename = conf->GetValue("Offset");
    if (monitor->LoadOffset(ifilename) == 0) {
      std::cerr << "[warning] offset file was not found, " << ifilename << std::endl;
    }
  }

  if        (channelAlign == "raw"       ) {
    monitor->SetChannelAlign(Extinction::Fct::MonitorWindow::ChannelAlign::Raw);
  } else if (channelAlign == "projection") {
    monitor->SetChannelAlign(Extinction::Fct::MonitorWindow::ChannelAlign::Projection);
  } else {
    std::cerr << "[warning] invalid align, set to raw align" << std::endl;
    monitor->SetChannelAlign(Extinction::Fct::MonitorWindow::ChannelAlign::Raw);
  }
  if        (monitorMode == "channel") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorWindow::MonitorMode::Channel, monitorChannel);
  } else if (monitorMode == "ext") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorWindow::MonitorMode::Ext);
  } else if (monitorMode == "hod") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorWindow::MonitorMode::Hod);
  } else if (monitorMode == "offset") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorWindow::MonitorMode::Offset, monitorChannel);
    monitor->SetMrSyncIntervalFilename(intervalFilename);
    monitor->SetCoinDiffsFilename(offsetFilename);
  } else if (monitorMode == "coincidence") {
    monitor->SetMonitorMode(Extinction::Fct::MonitorWindow::MonitorMode::Coincidence);
  } else {
    std::cerr << "[warning] invalid monitor mode, set to coincidence" << std::endl;
    monitor->SetMonitorMode(Extinction::Fct::MonitorWindow::MonitorMode::Coincidence);
  }

  monitor->SetAccumulate(accumulate);
  monitor->SetLeakageMode(leakageMode);
  
  monitor->SetHistoryWidth(historyWidth);
  monitor->SetCoinTimeWidth(coinTimeWidth);
  monitor->SetSpillLimit(spillLimit);
  monitor->SetReadBufferSize(readBufferSize);
  monitor->SetReadBufferMargin(readBufferMargin);

  monitor->SetTimePerTdc(timePerTdc);
  monitor->SetMrSyncInterval(mrSyncInterval);
  monitor->SetBunchCenters(bunchCenters);
  monitor->SetBunchWidths(bunchWidths);
  monitor->InitializeWindow(windowWidth ? windowWidth : 1600, windowHeight ? windowHeight : 1200);
  monitor->InitializePlots(profile);
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
      std::map<int, std::string> paths;
      int wd = inotify_add_watch(fd, path.data(), IN_ALL_EVENTS);
      if (wd < 0) {
        std::cerr << "[error] directory was not opened, " << path << std::endl;
        exit(EXIT_FAILURE);
      }
      paths[wd] = path;
      AddSubDirectory(fd, paths, path.data());

      // Create argument
      auto eventArg = new InotifyEventArgs();
      eventArg->board = board;
      eventArg->fd    = fd;
      // eventArg->wds   = { wd };
      // eventArg->path  = path;
      eventArg->paths = paths;

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
