# Install script for directory: /homes/4/greve/l/sp1/fsdev.github.local/mri_segment_hypothalamic_subunits

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/homes/4/greve/l/sp1/fsh.github.local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/python/scripts" TYPE FILE FILES "/homes/4/greve/l/sp1/fsdev.github.local/mri_segment_hypothalamic_subunits/mri_segment_hypothalamic_subunits")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  
      message(STATUS "Configuring python wrapper: /homes/4/greve/l/sp1/fsh.github.local/bin/mri_segment_hypothalamic_subunits")
      set(SCRIPTNAME mri_segment_hypothalamic_subunits)
      configure_file(/homes/4/greve/l/sp1/fsdev.github.local/python/wrapper /homes/4/greve/l/sp1/fsh.github.local/bin/mri_segment_hypothalamic_subunits @ONLY)
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/models" TYPE FILE FILES "/homes/4/greve/l/sp1/fsdev.github.local/mri_segment_hypothalamic_subunits/hypothalamic_subunits.h5")
endif()

