# This file is part of irtools.
# SPDX-Identifier: MIT
#
# MIT License
#   
# Copyright (c) 2024 Philipp Pracht
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
 

#########################################################################################
#########################################################################################
# Add modules to the CMake build 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

## specify module installation directory
#install(
#  DIRECTORY
#  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
#  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
#)

# Options for enabling or disabling features
option(WITH_OpenMP  "Enable OpenMP support"       TRUE)
option(WITH_TESTS   "Enable unit tests"           TRUE)
option(PYTHON_BINDINGS "Compile python bindings"  FALSE)


#########################################################################################
#########################################################################################

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()

#########################################################################################
#########################################################################################

# Compiler settings for GNU and Intel Fortran compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O0 -fbacktrace -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check -ffpe-trap=invalid,zero,overflow")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-g -O2 -r8 -align array64byte -traceback")
  set(bounds "-check all -fpe0")
else()
  message(FATAL_ERROR "Please use an Intel or GNU compiler!")
endif()

# Apply the compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}" PARENT_SCOPE)

#########################################################################################
#########################################################################################

# Populate irtools_metadata.fh
set(version ${PROJECT_VERSION}${SOVERSION})
execute_process(COMMAND git show -s --format=%h
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(git_return)
  set(commit "unknown-commit")
endif()
string(TIMESTAMP date "%a, %d %B %H:%M:%S, %m/%d/%Y")
set(author $ENV{USERNAME})
execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
set(origin ${HOSTNAME})
set(fcid ${CMAKE_Fortran_COMPILER_ID})
set(fcver ${CMAKE_Fortran_COMPILER_VERSION})
set(ccid ${CMAKE_C_COMPILER_ID})
set(ccver ${CMAKE_C_COMPILER_VERSION})
set(bsystem "cmake ${CMAKE_VERSION}")

configure_file(
  "${PROJECT_SOURCE_DIR}/assets/metadata.f90"
  "${PROJECT_BINARY_DIR}/include/irtools_metadata.fh"
  @ONLY
)

#########################################################################################
#########################################################################################
