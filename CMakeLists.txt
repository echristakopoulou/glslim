cmake_minimum_required(VERSION 2.8)

# project anme
PROJECT(SLIM)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# The version number
SET (Tutorial_VERSION_MAJOR 1)
SET (Tutorial_VERSION_MINOR 1)

# GKlib path
SET(GKLIB_PATH "${PROJECT_SOURCE_DIR}/GKlib" CACHE PATH "path to GKlib")
INCLUDE_DIRECTORIES(${GKLIB_PATH})

# BCLSlib path
set(BCLSLIB_PATH "${PROJECT_SOURCE_DIR}/BCLSlib" CACHE PATH "path to BCLSlib")
INCLUDE_DIRECTORIES(${BCLSLIB_PATH})

# always install
SET(SLIM_INSTALL TRUE)
IF (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local")
   # default installation direcotry
  SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/build")
ENDIF()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -std=c99 -Werror -Wall")

add_definitions(-DNDEBUG=1)
add_definitions(-DNDEBUG2=1)
		  
# messages
MESSAGE(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
MESSAGE(STATUS "GKLIB_PATH = ${GKLIB_PATH}")
MESSAGE(STATUS "BCLSLIB_PATH = ${BCLSLIB_PATH}")

# add subdirectories which contain the necesary code
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(include)
#ADD_SUBDIRECTORY(examples)



