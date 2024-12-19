 #The gurobi part of this file is based on the code supplied by gurobi here https://support.gurobi.com/hc/en-us/articles/360039499751-How-do-I-use-CMake-to-build-Gurobi-C-C-projects-

 # Only use GUROBI_HOME if GUROBI_DIR is not set
 if (NOT GUROBI_DIR)
     message(STATUS "GUROBI_DIR not set, trying to find Gurobi installation directory from environment variable GUROBI_HOME : $ENV{GUROBI_HOME}")
     set(GUROBI_DIR $ENV{GUROBI_HOME} CACHE PATH "Path to Gurobi installation directory" FORCE)
 endif ()
 message(STATUS "GUROBI_DIR = ${GUROBI_DIR}")

 find_path(
         GUROBI_INCLUDE_DIRS
         NAMES gurobi_c.h
         HINTS ${GUROBI_DIR}
         PATH_SUFFIXES include)

 find_library(
         GUROBI_LIBRARY
         NAMES gurobi gurobi100 gurobi110 gurobi120
         HINTS ${GUROBI_DIR}
         PATH_SUFFIXES lib)

 if (MSVC)
     set(MSVC_YEAR "2017")

     if (MT)
         set(M_FLAG "mt")
     else ()
         set(M_FLAG "md")
     endif ()

     find_library(GUROBI_CXX_LIBRARY
             NAMES gurobi_c++${M_FLAG} ${MSVC_YEAR}
             HINTS ${GUROBI_DIR}
             PATH_SUFFIXES lib)
     find_library(GUROBI_CXX_DEBUG_LIBRARY
             NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
             HINTS ${GUROBI_DIR}
             PATH_SUFFIXES lib)
 else ()
     find_library(GUROBI_CXX_LIBRARY
             NAMES gurobi_c++
             HINTS ${GUROBI_DIR}
             PATH_SUFFIXES lib)
     set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})
 endif ()

 include(FindPackageHandleStandardArgs)
 find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY
         GUROBI_INCLUDE_DIRS)