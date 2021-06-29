#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)

project(kc705)

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
set(TRON_CXX_FLAGS         ""             )
set(TRON_LIBRARIES         tron           )

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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -W -Wall")
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
add_executable(emcount src/emcount.cc ${headers})
add_executable(decoder src/decoder.cc ${headers})
add_executable(genHist src/genHist.cc ${headers})
add_executable(genCoin src/genCoin.cc ${headers})
add_executable(repHist src/repHist.cc ${headers})
add_executable(repCoin src/repCoin.cc ${headers})
add_executable(getLeak src/getLeak.cc ${headers})

# Add link libraries
target_link_libraries(emcount ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(decoder ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(genHist ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(genCoin ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(repHist ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(repCoin ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(getLeak ${ROOT_LIBRARIES} ${TRON_LIBRARIES})

#----------------------------------------------------------------------------
# Add install directories
#
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS emcount DESTINATION .)
install(TARGETS decoder DESTINATION .)
install(TARGETS genHist DESTINATION .)
install(TARGETS genCoin DESTINATION .)
install(TARGETS repHist DESTINATION .)
install(TARGETS repCoin DESTINATION .)
install(TARGETS getLeak DESTINATION .)
