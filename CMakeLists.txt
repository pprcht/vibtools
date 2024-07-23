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
# 

cmake_minimum_required(VERSION 3.17)
option(INSTALL_MODULES "Install Fortran module files to include directory." FALSE)

# Setup the irtools Project
project(
  irtools
  LANGUAGES "C" "Fortran"
  VERSION  0.1 
  DESCRIPTION "Tool for creating IR spectra in the double harmonic approximation"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Include further configurations
add_subdirectory("config")


###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################
# Check a specific CMake targets

# LAPACK, BLAS, OpenMP (usually via shared dependencies)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
  set(LINK_OpenMP TRUE)
else()
  set(LINK_OpenMP FALSE) 
endif()

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)
add_subdirectory("src")

if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
############################### Static Library ################################
###############################################################################
add_library(
   "lib-${PROJECT_NAME}-static"
   STATIC
   ${srcs}
)

target_link_libraries(
   "lib-${PROJECT_NAME}-static"
   PUBLIC
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   $<$<BOOL:${LINK_OpenMP}>:OpenMP::OpenMP_Fortran>
)


set_target_properties(
   "lib-${PROJECT_NAME}-static"
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
   ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
   POSITION_INDEPENDENT_CODE ON
   OUTPUT_NAME "${PROJECT_NAME}"
 )


target_include_directories(
  "lib-${PROJECT_NAME}-static"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)


###############################################################################
############################### Executables ###################################
###############################################################################
add_executable(
  ${PROJECT_NAME}-exe
  ${prog}
)

target_link_libraries(
   ${PROJECT_NAME}-exe
   PRIVATE
   "lib-${PROJECT_NAME}-static"
)

set_target_properties(
  ${PROJECT_NAME}-exe
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)

target_include_directories(
   ${PROJECT_NAME}-exe 
   PRIVATE 
   ${PROJECT_SOURCE_DIR}/include
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)


###############################################################################
########################### Install Instructions ##############################
###############################################################################

# Build output artifacts
if (WITH_OBJECT)
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
else()
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

##############################################################################
############################# Unit Tests #####################################
##############################################################################

#if (WITH_TESTS)
#  enable_testing()
#  add_subdirectory("test")
#endif()
