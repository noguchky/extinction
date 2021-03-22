#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)

project(fct)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# message(STATUS "CMAKE_EXPORT_COMPILE_COMMANDS = ${CMAKE_EXPORT_COMPILE_COMMANDS}")

if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

#----------------------------------------------------------------------------
# Find CERN/ROOT
#
list(APPEND CMAKE_PREFIX_PATH  $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/cmake)
find_package(ROOT REQUIRED)
# message(STATUS "ROOT_INCLUDE_DIRS = ${ROOT_INCLUDE_DIRS}")
# message(STATUS "ROOT_LIBRARY_DIR  = ${ROOT_LIBRARY_DIR}")
# message(STATUS "ROOT_CXX_FLAGS    = ${ROOT_CXX_FLAGS}")
# message(STATUS "ROOT_LIBRARIES    = ${ROOT_LIBRARIES}")

#----------------------------------------------------------------------------
# Find Tron
#
set(TRON_INCLUDE_DIR       ${CMAKE_SOURCE_DIR}/../tron/include)
set(TRON_LIBRARY_DIR       ${CMAKE_SOURCE_DIR}/../tron/lib    )
set(TRON_CXX_FLAGS         ""                                 )
set(TRON_LIBRARIES         tron                               )

#----------------------------------------------------------------------------
# Locate include directories
#
include_directories(include ../common/include)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${TRON_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Add compile option
#
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -W -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -W -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRON_CXX_FLAGS}")

#----------------------------------------------------------------------------
# Locate sources and headers
#
file(GLOB headers include/* ../common/include/*)

#----------------------------------------------------------------------------
# Add the executable & link libraries
#
link_directories(${ROOT_LIBRARY_DIR})
link_directories(${TRON_LIBRARY_DIR})

# Add the executable
add_executable(genMock src/genMock.cc ${headers})
add_executable(decoder src/decoder.cc ${headers})
add_executable(anaRaw  src/anaRaw.cc  ${headers})
add_executable(anaBeam src/anaBeam.cc ${headers})
add_executable(genHist src/genHist.cc ${headers})
add_executable(watcher src/watcher.cc ${headers})
if(LINUX)
  add_executable(monitor src/monitor.cc ${headers})
endif()

# Add link libraries
target_link_libraries(genMock ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(decoder ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaRaw  ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaBeam ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(genHist ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(watcher ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
if(LINUX)
  target_link_libraries(monitor ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
endif()

#----------------------------------------------------------------------------
# Add install directories
#
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS genMock DESTINATION .)
install(TARGETS decoder DESTINATION .)
install(TARGETS anaRaw  DESTINATION .)
install(TARGETS anaBeam DESTINATION .)
install(TARGETS genHist DESTINATION .)
install(TARGETS watcher DESTINATION .)
if(LINUX)
  install(TARGETS monitor DESTINATION .)
endif()
