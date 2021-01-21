#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <limits>
#include <stdlib.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <time.h>
#include <sys/time.h>
#include <signal.h>
#include <errno.h>
#include "ArgReader.hh"

constexpr unsigned int BUFSIZE  = 256000;
constexpr int          THRESIZE = BUFSIZE * 3 / 4;
constexpr int          UNITSIZE = 13;

// int comSignal = 1;
long long evlimit = 0;
void alarm(int) {
  puts("Detect ctrl+c");
  // comSignal = 0;
  evlimit = 0;
}

struct footer_t {
  unsigned int  hoge0;
  unsigned char hoge1;
  unsigned int  magicWord0;
  unsigned int  magicWord1;
  inline bool IsFooter() const {
    static const unsigned int footerWord = htonl(0xAAAAAAAA);
    return magicWord0 == footerWord;
  }
} __attribute__((__packed__));

struct header_t {
  unsigned int  magicWord0;
  unsigned int  magicWord1;
  unsigned int  hoge0;
  unsigned char hoge1;
  inline bool IsHeader() const {
    static const unsigned int headerWord = htonl(0x34567012);
    return magicWord1 == headerWord;
  }
} __attribute__((__packed__));

struct DataType {
  enum {
        Header = -1,
        Data   =  0,
        Footer =  1,
  };
};

inline int GetDataType(unsigned char* rcvdBuffer) {
  if     (((header_t*)rcvdBuffer)->IsHeader()) return DataType::Header;
  else if(((footer_t*)rcvdBuffer)->IsFooter()) return DataType::Footer;
  return DataType::Data;
}

inline const char* filename(const std::string* dir) {
  static char fname[128];
  static long count = 0;
  sprintf(fname, "kc705_%012ld.dat", count++);
  static std::stringstream stream;
  stream.str(*dir);
  stream << fname;
  return stream.str().data();
}

inline const char* commandNull(const std::string* , const std::string*) {
  return nullptr;
}

inline const char* commandMove(const std::string* fname, const std::string* dir) {
  static std::string cmd(256, ' ');
  cmd = "mv ";
  cmd += *fname;
  cmd += " ";
  cmd += *dir;
  cmd += " &";
  return cmd.data();
}

inline const char* commandCompress(const std::string* fname, const std::string* dir) {
  static std::string cmd(256, ' ');
  cmd = "tar -czf ";
  cmd += *dir;
  cmd += "/";
  cmd += *fname;
  cmd += ".tar.gz --remove-files ";
  cmd += *fname;
  cmd += " &";
  return cmd.data();
}

// for Debug
int read(unsigned char* buff, int size) {
  const int n = std::min(size, 2600);
  // const int n = std::min(size, 1326);
  // const int n = std::min(size,  130);
  for (int i = 0; i < n; ++i) {
    buff[i] = i;
  }
  return n;
}

void Usage(){
  //printf("Usage: %s <IP address> <TCP Port #> <# of Event> <File name> \n\n", argv[0]);
  printf("Usage: ./simple-daq <IP address> <TCP Port #>\n");
  printf("     -n  (int)  Set # of events\n");
  printf("     -f  (char) Set output file name\n");
  printf("     -h         Help\n\n");
}

int main(int argc, char** argv){
  const std::string AppName = "daq";

  // Load Arguments
  Tron::ArgReader* args = new Tron::ArgReader(AppName);
  args->AddArg<std::string>("IpAddress"     ,                   "IP address");
  args->AddArg<UInt_t>     ("TcpPort#"      ,                   "TCP port number");
  args->AddOpt<Long64_t>   ("NofEvents"     , 'n', "number"   , "Set # of events in a run" , "1000");
  args->AddOpt<Long64_t>   ("NofSpills"     , 's', "spill"    , "Set # of spills in a file", "100");
  args->AddOpt<Int_t>      ("Pattern"       , 'p', "pattern"  , "Set daq pattern"          , "1");
  args->AddOpt<Int_t>      ("ReadSize"      , 'r', "readsize" , "Set read size [bytes]"    , "2600");
  args->AddOpt<std::string>("Directory"     , 'd', "directory", "Set output directory for raw files" , "");
  args->AddOpt             ("Move"          , 'm', "move"     , "move file to forward directory");
  args->AddOpt             ("Compress"      , 'c', "compress" , "compress and move file to forward directory");
  args->AddOpt<std::string>("ForwardDir"    , 'f', "forward"  , "Set forward directory", ".");
  args->AddOpt             ("Help"          , 'h', "help"     , "Help");

  if (!args->Parse(argc, argv) || args->IsSet("Help") || args->HasUnsetRequired()) {
    args->ShowUsage();
    return 0;
  }

  const std::string  ipAddress   = args->GetValue("IpAddress");
  const unsigned int tcpPort     = args->GetValue<unsigned int>("TcpPort#");
  const long long    nofEvents   = args->GetValue<long long>("NofEvents");
  const long long    nofSpills   = args->GetValue<long long>("NofSpills");
  const int          pattern     = args->GetValue<int>("Pattern");
  const int          readSize    = args->GetValue<int>("ReadSize");
  const std::string  directory   = args->GetValue("Directory");
  const bool         move        = args->IsSet("Move");
  const bool         compress    = args->IsSet("Compress");
  const std::string  forwardDir  = args->GetValue("ForwardDir");
  const char* (*command)(const std::string*,const std::string*) = nullptr;
  
  std::cout << "IP Address: " << ipAddress << std::endl;
  std::cout << "TCP Port: " << tcpPort << std::endl;
  std::cout << "Set # of events: " << nofEvents << std::endl;
  std::cout << "Set # of spills: " << nofSpills << std::endl;
  std::cout << "Set pattern: " << pattern << std::endl;
  std::cout << "Set read size: " << readSize << std::endl;
  if (readSize <= 0 || (unsigned int)readSize > BUFSIZE) {
    std::cout << "[error] invalid read size" << std::endl;
    exit(1);
  }
  std::cout << "Set output directory: \"" << directory << "\"" << std::endl;
  if (move) {
    std::cout << "Set move" << std::endl;
    command = commandMove;
  } else if (compress) {
    std::cout << "Set compress" << std::endl;
    command = commandCompress;
  } else {
    std::cout << "Set no move and no compress" << std::endl;
    command = commandNull;
  }
  std::cout << "Set forward directory: \"" << forwardDir << "\"" << std::endl;
  
  // create a socket
  std::cout << std::endl
            << "Create socket..." << std::endl;
  int sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
  std::cout << "socket created" << std::endl;

  sockaddr_in sitcpAddr;
  sitcpAddr.sin_family      = AF_INET;
  sitcpAddr.sin_port        = htons(tcpPort);
  sitcpAddr.sin_addr.s_addr = inet_addr(ipAddress.data());

  std::cout << "connecting..." << std::endl;
  if (connect(sock, (sockaddr*)&sitcpAddr, sizeof(sitcpAddr)) < 0) {
    puts("Connect() faild");
    close(sock);
    exit(EXIT_FAILURE);
  }
  std::cout << "connection succeded" << std::endl;

  std::cout << "Set signal" << std::endl;
  if (signal(SIGINT, alarm) == SIG_ERR) {
    std::cerr << "error signal" << std::endl;
  }

  std::cout << "Define daq variables" << std::endl;
  unsigned char  rcvdBuffer[BUFSIZE];
  int            filledLength  = 0;
  int            recvLength    = 0;
  int            writtenLength = 0;
  int            writeBegin    = 0;
  int            writeEnd      = 0;
  long long      ievent        = 0;
  long long      ispill        = 0;
  const char*    cmd           = nullptr;
  
  std::cout << "Open output file" << std::endl;
  FILE* fout;
  std::string foutname;
  {
    foutname = filename(&directory);
    std::cout << "File Open: " << foutname << std::endl;
    if ((fout = fopen(foutname.data(), "wb")) == nullptr) {
      std::cerr << "Can't open " << foutname << std::endl; 
      close(sock);
      exit(1);
    }
  }

  std::cout << "Start daq" << std::endl;
  time_t t_start = time(nullptr);

  switch (pattern) {
  case 0:
    {
      // int k;
      // const long long evlimit = nofEvents > 0 ? nofEvents : std::numeric_limits<long long>::max();
      // while (ievent <= evlimit && comSignal) {
      evlimit = nofEvents > 0 ? nofEvents : std::numeric_limits<long long>::max();
      while (ievent <= evlimit) {
        if ((filledLength = read(sock, rcvdBuffer, readSize)) <= 0) {
          std::cerr << " recv() faild" << std::endl;
          close(sock);
          fclose(fout);
          // Compress and move output file
          if ((cmd = command(&foutname, &forwardDir))) {
            system(cmd);
          }
          exit(EXIT_FAILURE);
        }

        // for (k = 0; k < filledLength; ++k) {
        //   printf("%02x ", rcvdBuffer[k]);
        //   if (k % UNITSIZE == UNITSIZE - 1) {
        //     puts("");
        //   }
        // }

        fwrite(rcvdBuffer,  sizeof(char), filledLength, fout);
        if (++ispill % nofSpills == 0) {
          fclose(fout);

          // Compress and move output file
          if ((cmd = command(&foutname, &forwardDir))) {
            system(cmd);
          }

          // Open next file
          foutname = filename(&directory);
          std::cout << "File Open: " << foutname << std::endl;
          if ((fout = fopen(foutname.data(), "wb")) == nullptr) {
            std::cerr << "[error] output file was not opened, " << foutname << std::endl;
            close(sock);
            exit(1);
          }
        }

        ievent += filledLength / UNITSIZE;
      }
    }
    break;

  case 1:
    {
      // int k;
      // const long long evlimit = nofEvents > 0 ? nofEvents : std::numeric_limits<long long>::max();
      // while (ievent <= evlimit && comSignal) {
      evlimit = nofEvents > 0 ? nofEvents : std::numeric_limits<long long>::max();
      const int eventunit = readSize / UNITSIZE;
      while (ievent <= evlimit) {
        for (filledLength = 0; filledLength < readSize;) {
          if ((recvLength = read(sock, rcvdBuffer + filledLength, readSize - filledLength)) <= 0) {
            std::cerr << " recv() faild" << std::endl;
            close(sock);
            fclose(fout);
            // Compress and move output file
            if ((cmd = command(&foutname, &forwardDir))) {
              system(cmd);
            }
            exit(EXIT_FAILURE);
          } else {
            filledLength += recvLength;
          }
        }

        // for (writtenLength = 0; writtenLength < filledLength; writtenLength += UNITSIZE) {
        //   for (k = 0; k < UNITSIZE; ++k) {
        //     printf("%02x ", (rcvdBuffer + writtenLength)[k]);
        //   }
        //   puts("");
        // }

        fwrite(rcvdBuffer,  sizeof(char), filledLength, fout);
        if (++ispill % nofSpills == 0) {
          fclose(fout);

          // Compress and move output file
          if ((cmd = command(&foutname, &forwardDir))) {
            system(cmd);
          }

          // Open next file
          foutname = filename(&directory);
          std::cout << "File Open: " << foutname << std::endl;
          if ((fout = fopen(foutname.data(), "wb")) == nullptr) {
            std::cerr << "[error] output file was not opened, " << foutname << std::endl;
            close(sock);
            exit(1);
          }
        }

        ievent += eventunit;
      }
    }
    break;

  case 2:
    {
      // int k;
      // const long long evlimit = nofEvents > 0 ? nofEvents : std::numeric_limits<long long>::max();
      // while (ievent <= evlimit && comSignal) {
      evlimit = nofEvents > 0 ? nofEvents : std::numeric_limits<long long>::max();
      while (ievent <= evlimit) {

        // std::cout << "Read socket" << std::endl;
        if ((recvLength = read(sock, rcvdBuffer + writtenLength + filledLength, BUFSIZE - (writtenLength + filledLength))) <= 0) {
          std::cerr << " recv() faild" << std::endl;
          close(sock);
          fclose(fout);
          // Compress and move output file
          if ((cmd = command(&foutname, &forwardDir))) {
            system(cmd);
          }
          exit(EXIT_FAILURE);
        } else {
          filledLength += recvLength;
        }

        // std::cout << "// Data processing" << std::endl;
        if (filledLength >= UNITSIZE) {
          for (writeBegin = 0, writeEnd = UNITSIZE; writeEnd <= filledLength;) {
            // for (k = 0; k < UNITSIZE; ++k) {
            //   printf("%02x ", (rcvdBuffer + writtenLength + writeBegin)[k]);
            // }
            // puts("");

            if (GetDataType(rcvdBuffer + writtenLength + writeBegin) == DataType::Footer) {
              if (++ispill % nofSpills == 0) {
                fwrite(rcvdBuffer + writtenLength,  sizeof(char), writeEnd, fout);
                writtenLength += writeEnd;
                filledLength  -= writeEnd;
                writeBegin     = 0;
                writeEnd       = UNITSIZE;

                fclose(fout);

                // Compress and move output file
                if ((cmd = command(&foutname, &forwardDir))) {
                  system(cmd);
                }

                // Open next file
                foutname = filename(&directory);
                std::cout << "File Open: " << foutname << std::endl;
                if ((fout = fopen(foutname.data(), "wb")) == nullptr) {
                  std::cerr << "[error] output file was not opened, " << foutname << std::endl;
                  close(sock);
                  exit(1);
                }
              } else {
                writeBegin  = writeEnd;
                writeEnd   += UNITSIZE;
              }
            } else {
              writeBegin  = writeEnd;
              writeEnd   += UNITSIZE;
            }
            ++ievent;
          }

          if (writeBegin) {
            fwrite(rcvdBuffer + writtenLength,  sizeof(char), writeBegin, fout);
            writtenLength += writeBegin;
            filledLength  -= writeBegin;
          }

          if (filledLength == 0) {
            writtenLength = 0;
          } else if (writtenLength > THRESIZE) {
            std::memcpy(rcvdBuffer, rcvdBuffer + writtenLength, filledLength);
            writtenLength = 0;
          }
        }
      }
    }
    break;

  default:
    std::cerr << "invalid pattern" << std::endl;
    close(sock);
    exit(1);
    break;
  }

  // Finish daq
  time_t t_stop = time(nullptr);
  std::cout << std::endl
            << "Finished" << std::endl;

  // Close socket
  std::cout << "socket closing" << std::endl;
  close(sock);
  std::cout << "socket closed" << std::endl;

  // Close file
  fclose(fout);

  // Compress and move output file
  if ((cmd = command(&foutname, &forwardDir))) {
    system(cmd);
  }

  // Verbose daq time
  const int t_diff = t_stop - t_start;
  std::cout << "DAQ time  " << t_stop << " - " << t_start << " = " << t_diff << " sec" << std::endl;

  // const std::string logname = "../log/run-summary.log";
  // FILE *flog;
  // if ((flog = fopen(logname.data(), "a")) == nullptr){
  //   std::cout << "Can't open " << logname << std::endl;
  // } else {
  //   fprintf(flog, "%s %s %ld -> %ld : %d\n", argv[0], filename.data(), t_start, t_stop, t_diff);
  //   fclose(flog);
  // }
  
  return 0;

}
