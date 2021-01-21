#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
cmake_policy(SET CMP0042 NEW)

project(tron)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#message(STATUS "CMAKE_EXPORT_COMPILE_COMMANDS = ${CMAKE_EXPORT_COMPILE_COMMANDS}")

#----------------------------------------------------------------------------
# Find CERN/ROOT
#
list(APPEND CMAKE_PREFIX_PATH  $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/cmake)
find_package(ROOT REQUIRED)
message(STATUS "ROOT_INCLUDE_DIRS = ${ROOT_INCLUDE_DIRS}")
message(STATUS "ROOT_LIBRARY_DIR  = ${ROOT_LIBRARY_DIR}")
message(STATUS "ROOT_CXX_FLAGS    = ${ROOT_CXX_FLAGS}")
message(STATUS "ROOT_LIBRARIES    = ${ROOT_LIBRARIES}")

#----------------------------------------------------------------------------
# Add compile option
#
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -W -Wall -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -W -Wall -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

#----------------------------------------------------------------------------
# Locate include directories
#
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(include)

#----------------------------------------------------------------------------
# Locate sources and headers
#
file(GLOB sources     src/*.cc     src/Function/*.cc     src/Random/*.cc)
file(GLOB headers include/*.hh include/Function/*.hh include/Random/*.hh)

#----------------------------------------------------------------------------
# Add the executable
#
link_directories(${ROOT_LIBRARY_DIR})

add_library(tron-object OBJECT ${sources} ${headers})
set_property(TARGET tron-object PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(tron-static STATIC $<TARGET_OBJECTS:tron-object>)
add_library(tron-shared SHARED $<TARGET_OBJECTS:tron-object>)

target_link_libraries(tron-shared ${ROOT_LIBRARIES})

set_target_properties(tron-static PROPERTIES OUTPUT_NAME tron)
set_target_properties(tron-shared PROPERTIES OUTPUT_NAME tron)

#----------------------------------------------------------------------------
# Add install directories
#
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS tron-static DESTINATION lib)
install(TARGETS tron-shared DESTINATION lib)
