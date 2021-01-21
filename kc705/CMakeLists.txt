#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)

project(kc705)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# message(STATUS "CMAKE_EXPORT_COMPILE_COMMANDS = ${CMAKE_EXPORT_COMPILE_COMMANDS}")

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
add_executable(daq     src/daq.cc     ${headers})
add_executable(decoder src/decoder.cc ${headers})
add_executable(checker src/checker.cc ${headers})
add_executable(monitor src/monitor.cc ${headers})
add_executable(anaRaw  src/anaRaw.cc  ${headers})
add_executable(anaTree src/anaBeam.cc ${headers})
add_executable(anaBeam src/anaBeam.cc ${headers})
add_executable(anaCoin src/anaCoin.cc ${headers})
add_executable(anaDark src/anaDark.cc ${headers})

# Add link libraries
target_link_libraries(daq     ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(decoder ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(checker ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(monitor ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaRaw  ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaTree ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaBeam ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaCoin ${ROOT_LIBRARIES} ${TRON_LIBRARIES})
target_link_libraries(anaDark ${ROOT_LIBRARIES} ${TRON_LIBRARIES})

#----------------------------------------------------------------------------
# Add install directories
#
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS daq     DESTINATION .)
install(TARGETS decoder DESTINATION .)
install(TARGETS checker DESTINATION .)
install(TARGETS monitor DESTINATION .)
install(TARGETS anaRaw  DESTINATION .)
install(TARGETS anaTree DESTINATION .)
install(TARGETS anaBeam DESTINATION .)
install(TARGETS anaCoin DESTINATION .)
install(TARGETS anaDark DESTINATION .)
