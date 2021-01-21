#include <iostream>
#include <fstream>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "usage: " << argv[0] << " [filename]" << std::endl;
  }
  std::ifstream ifile(argv[1], std::ios::binary);

  unsigned int  ibuf = 0;
  unsigned char cbuf = 0;

  while (true) {
    if (ifile.read((char*)&ibuf, sizeof(ibuf))) {
      if (ifile.read((char*)&cbuf, sizeof(cbuf))) {
        printf("Header 0x%02X  Channel  %d  TDC  %d\n", (cbuf & 0xF0), ((ibuf & 0x1F80000) >> 19), (ibuf & 0x7FFFF));
      } else {
        break;
      }
    } else {
      break;
    }
  }

  return 0;
}
