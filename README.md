## My Environment
- Mac
  - Apple LLVM version 10.0.1
  - cmake version 3.16.2
  - ROOT Version: 6.18/04

- KEK CC
  - gcc (GCC) 8.4.0
  - cmake3 version 3.17.5
  - ROOT Version: 6.22/02

## Build
1. Build tron library
```
$ cd path/to/extinction/tron
```

Make build directory
```
$ mkdir build
$ cd build
```

Cmake
```
$ cmake ..
```
or
```
$ cmake3 ..
```

Build library and install
```
$ make
$ make install
```
libtron is installed at path/to/extinction/tron/lib

2. Build FCT
```
$ cd path/to/extinction/fct
```

The same as tron (mkdir -> cmake -> make -> make install)

The executables are installed at path/to/extinction/fct

4. Build HUL
```
$ cd path/to/extinction/hul
```

The same as tron (mkdir -> cmake -> make -> make install)

The executables are installed at path/to/extinction/hul

3. Build KC705
```
$ cd path/to/extinction/kc705
```

The same as tron (mkdir -> cmake -> make -> make install)

The executables are installed at path/to/extinction/kc705
