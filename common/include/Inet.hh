#ifndef Extinction_Inet_hh
#define Extinction_Inet_hh

#include <arpa/inet.h>
#include "Rtypes.h"

namespace Extinction {

  inline ULong64_t pntohll(UChar_t* p) {
    return
      (ULong64_t)ntohl(*(UInt_t*) p     ) * (1ULL << 32) +
      (ULong64_t)ntohl(*(UInt_t*)(p + 4));
  }

  inline UInt_t pntohl(UChar_t* p) {
    return ntohl(*(UInt_t*)p);
  }

  inline UShort_t pntohs(UChar_t* p) {
    return ntohs(*(UShort_t*)p);
  }

  inline UChar_t pntohc(UChar_t* p) {
    return *p;
  }

}

#endif
