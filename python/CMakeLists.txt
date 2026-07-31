project(fspython)

# build the cxx/python bindings for the below modules
add_subdirectories(fsbindings)
## gems in atiic as of 6/2026
# if(NOT INTEGRATE_SAMSEG)
#  add_subdirectories(gems/cxx)
# endif()

## FSPYTHON_VERSION set in top level CMakeLists.txt
# set(FSPYTHON_VERSION "3.8")

set(UV_INSTEAD_OF_PIP_INSTALL OFF)
# set(UV_INSTEAD_OF_PIP_INSTALL ON)

if(UV_INSTEAD_OF_PIP_INSTALL)
   ## 04/2024 allow option for uv package installer (added to fspython base distro)
   # set(PIP_SUBSTRING "uv pip install -v")
   set(PIP_SUBSTRING "uv pip install")
   set(PIP_NODEP_ARG "--no-deps")
   set(PIP_NOVER_ARG "")
   set(PIP_LAST_ARG "--python ${CMAKE_INSTALL_PREFIX}/python")
   set(pip_target_option "--target")
else()
   # set(PIP_SUBSTRING "pip install -v")
   # set(PIP_SUBSTRING "pip install --no-cache-dir -vvv")
   set(PIP_SUBSTRING "pip install")
   set(PIP_NODEP_ARG "--no-dependencies")
   set(PIP_NOVER_ARG "--disable-pip-version-check")
   set(PIP_LAST_ARG "")
   set(pip_target_option "--target")
endif()

if(HOST_OS MATCHES "Ubuntu24")
   if (NOT PYTHON_EXECUTABLE)
      set(PYTHON_EXECUTABLE "${FS_PACKAGES_DIR}/fspython/${FSPYTHON_VERSION}/bin/python")
      message(STATUS "For HOST_OS=${HOST_OS} set PYTHON_EXEUTABLE to ${PYTHON_EXECUTABLE}")
   endif()
endif()

 ## should not be used anymore (use requirements-build.txt)
 # set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt")
 ## should not be used anymore (use requirements-build.txt)
 # set(requirements_extra "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-extra.txt")
 ## cannot be subsumed into requirements-build.txt until INTEGRATE_SAMSEG option is no longer used
 set(requirements_samseg_path "${CMAKE_CURRENT_SOURCE_DIR}")
 set(requirements_samseg_file "requirements-samseg.txt")
 set(requirements_samseg "-r ${requirements_samseg_path}/${requirements_samseg_file}")
 set(pybind11_DIR           "${CMAKE_SOURCE_DIR}/packages/pybind11")

 ## default linux python currently 3.8
 if(NOT FSPYTHON_VERSION)
    set(FSPYTHON_VERSIOM "3.8")
 endif()

 if("${FSPYTHON_VERSION}" MATCHES "3.8")
    ## previous behavior before differentiated by python version (requirements-build.txt and requirements-build-py38.txt started out as same file)
    # set(FSPYTHON_REQFILE "requirements-build.txt")
    if(FSPYTHON_INSTALL_CUDA)
       set(FSPYTHON_REQFILE "requirements-build-py38-cuda.txt")
    else()
       set(FSPYTHON_REQFILE "requirements-build-py38.txt")
    endif() 
 elseif("${FSPYTHON_VERSION}" MATCHES "3.9")
    if(FSPYTHON_INSTALL_CUDA)
       set(FSPYTHON_REQFILE "requirements-build-py39-cuda.txt")
    else()
       set(FSPYTHON_REQFILE "requirements-build-py39.txt")
    endif()
 endif()

 if(FSPYTHON_BUILD_REQ)
    # pip install current snapshot of (hardcoded) package revisions into fspython
    if(NOT APPLE)
       ## python requirments file for all linux OS
       # set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-linux.txt")
       set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/${FSPYTHON_REQFILE}")
       set(requirements_extra "")
    else()
       if(NOT APPLE_ARM64)
          ## darwin_x86_64
          if("${FSPYTHON_VERSION}" MATCHES "3.8")
             set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-build.darwin_x86_64-py38.txt")
          elseif("${FSPYTHON_VERSION}" MATCHES "3.9")
             set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-build.darwin_x86_64-py39.txt")
          endif()
          set(requirements_extra "")
          message(STATUS "For HOST_OS=${HOST_OS} using requirements file ${requirements}")
       else()
          ## darwin_arm64
          if("${FSPYTHON_VERSION}" MATCHES "3.8")
             set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-build.arm_x86_64-py38.txt")
          elseif("${FSPYTHON_VERSION}" MATCHES "3.9")
             set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-build.arm_x86_64-py39.txt")
          endif()
          set(requirements_extra "")
          message(STATUS "For HOST_OS=${HOST_OS} using requirements file ${requirements}")
       endif()
    endif()
 endif()
 if(NOT APPLE_ARM64)
    if(INFANT_MODULE)
       ## provide requirements file for infant module to allow for different package revisions compared to what
       ## freesurfer uses (add cmake option -DFSPYTHON_BUILD_REQ=ON to infant options -DINFANT_MODULE=ON -DMINIMAL=ON)
       # set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-infant.txt")
       # set(requirements_extra "")
    endif()
 endif()
 if(NOT APPLE_ARM64)
    message(STATUS "For HOST_OS=${HOST_OS} fspython pip set to use requirement arguments: ${requirements} ${requirements_extra}")
 endif() 

# Definition of FSPYTHON_INSTALL_PREFIX is in top level CMakeLists.txt

function(requirements_file_install)
  if(FSPYTHON_BUILD_REQ)
      if(NOT APPLE_ARM64)
         set(num ${FSPYTHON_VERSION})
         string(REPLACE "." "" num "${num}")
         if(FSPYTHON_INSTALL_TREE)
            # install the approriate requirement file(s) for the python version per FSPYTHON_VERSION
            install(PROGRAMS requirements-build-py${num}.txt DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/requirements-build-py${num}-cuda.txt")
               install(PROGRAMS requirements-build-py${num}-cuda.txt DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python)
            endif()
         else()
            install(PROGRAMS requirements-build-py${num}.txt DESTINATION python)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/requirements-build-py${num}-cuda.txt")
               install(PROGRAMS requirements-build-py${num}-cuda.txt DESTINATION python)
            endif()
         endif()
      endif()
   endif()
endfunction()

function(fspython_script_install)
   # install the fspython script - FOR NOW install in both trees if split tree install 
   if(APPLE_ARM64)
      # Allow script to use local packages with the system or Apple recommended Anaconda python
      install(PROGRAMS fspython.local DESTINATION bin RENAME fspython)
      if(FSPYTHON_INSTALL_TREE)
         install(PROGRAMS fspython.local DESTINATION ${FSPYTHON_INSTALL_PREFIX}/bin RENAME fspython)
      endif() 
   else()
      install(PROGRAMS fspython DESTINATION bin)
      if(FSPYTHON_INSTALL_TREE)
         install(PROGRAMS fspython DESTINATION ${FSPYTHON_INSTALL_PREFIX}/bin)
      endif()
   endif()
endfunction()

# TO DO  - write a python_3_10_base_install

## FSPYTHON_DIR set in top level ../CMakeLists.txt
function(python_3_9_base_install)
  if(EXISTS ${FSPYTHON_DIR}/bin/python3)
    # specify the python to use
    set(PYTHON_EXECUTABLE ${FSPYTHON_DIR}/bin/python3)
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.9m)
    # install the barebones python distribution
    INSTALL(DIRECTORY ${FSPYTHON_DIR}/bin ${FSPYTHON_DIR}/lib ${FSPYTHON_DIR}/include
      DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python USE_SOURCE_PERMISSIONS
    )
  else()
    message(FATAL_ERROR "Cannot find fspython distribution in FS_PACKAGES_DIR")
  endif()
  message(STATUS "For HOST_OS=${HOST_OS} FSPYTHON_VERSION=${FSPYTHON_VERSION}")
endfunction()

function(python_3_8_base_install)
  if(EXISTS ${FSPYTHON_DIR}/bin/python3)
    # specify the python to use
    set(PYTHON_EXECUTABLE ${FSPYTHON_DIR}/bin/python3)
    set(PYTHON_LIBRARY ${FSPYTHON_DIR}/lib/libpython3.8m)
    # install the barebones python distribution
    INSTALL(DIRECTORY ${FSPYTHON_DIR}/bin ${FSPYTHON_DIR}/lib ${FSPYTHON_DIR}/include
      DESTINATION ${FSPYTHON_INSTALL_PREFIX}/python USE_SOURCE_PERMISSIONS
    )
  else()
    message(FATAL_ERROR "Cannot find fspython distribution in FS_PACKAGES_DIR")
  endif()
  message(STATUS "For HOST_OS=${HOST_OS} FSPYTHON_VERSION=${FSPYTHON_VERSION}")
endfunction()

function(packages_install)
   # run 'pip install' and direct package installs to freesurfer/python/packages
   # pass INTEGRATE_SAMSEG_OPTION=${INTEGRATE_SAMSEG} environment variable for setup.py to act accordingly
   install(CODE "
     message(STATUS \"Installing python packages\")
     message(STATUS \"INTEGRATE_SAMSEG_OPTION=${INTEGRATE_SAMSEG} ${PYTHON_EXECUTABLE} -m ${PIP_SUBSTRING}  ${CMAKE_CURRENT_SOURCE_DIR} ${PIP_NODEP_ARG} --upgrade ${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages ${PIP_LAST_ARG}\")
     execute_process(COMMAND bash -c \"INTEGRATE_SAMSEG_OPTION=${INTEGRATE_SAMSEG} ${PYTHON_EXECUTABLE} -m ${PIP_SUBSTRING} ${CMAKE_CURRENT_SOURCE_DIR} ${PIP_NODEP_ARG} --upgrade ${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages ${PIP_LAST_ARG}\" RESULT_VARIABLE retcode)
     if(NOT \${retcode} STREQUAL 0)
       message(FATAL_ERROR \"Could not install python modules\")
     endif()"
     COMPONENT fspy
   )
endfunction()

# pip can install under python/packages but many/most packages from the requirements files are installed under site-package
function(sitepackages_install)
  install(CODE "
    message(STATUS \"Installing python package dependencies into fspython with cmd:\")
    message(STATUS \"${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} --upgrade ${PIP_NOVER_ARG} ${EXTRA_PACKAGES} ${requirements} ${requirements_extra} ${PKG_TARGET} ${PIP_LAST_ARG}\")
    execute_process(COMMAND bash -c \"${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} --upgrade ${PIP_NOVER_ARG} ${EXTRA_PACKAGES} ${requirements} ${requirements_extra} ${PKG_TARGET} ${PIP_LAST_ARG}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install python dependencies\")
    endif()
  ")
endfunction()

# Expand packages to install from all requirements files directly on to a single pip command line via req.sh.
# This helps pip solve for all dependenices.  Also apparently needed to avoid an issue in docker for linux pip installs.
# Samseg install complicates setup by adding env vars with cmake config on command line; would be nice to eliminate this.

function(sitepackages_install_rev2)
  install(CODE "
    message(STATUS \"Installing python package dependencies into fspython with cmd:\")
    message(STATUS \"FSPYTHON_VERSION=${FSPYTHON_VERSION} ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} SAMSEG_REQ=${requirements_samseg_path}/${requirements_samseg_file} INTEGRATE_SAMSEG=${INTEGRATE_SAMSEG} FSPYTHON_INSTALL_CUDA=${FSPYTHON_INSTALL_CUDA} FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --expand ${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} ${PIP_NODEP_ARG} --upgrade ${PIP_NOVER_ARG} ${EXTRA_PACKAGES} ${requirements} ${requirements_extra} ${PKG_TARGET} ${PIP_LAST_ARG}\")
    execute_process(COMMAND bash -c \"FSPYTHON_VERSION=${FSPYTHON_VERSION} ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} SAMSEG_REQ=${requirements_samseg_path}/${requirements_samseg_file} INTEGRATE_SAMSEG=${INTEGRATE_SAMSEG} FSPYTHON_INSTALL_CUDA=${FSPYTHON_INSTALL_CUDA} FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --expand ${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} ${PIP_NODEP_ARG} --upgrade ${PIP_NOVER_ARG} ${EXTRA_PACKAGES} ${requirements} ${requirements_extra} ${PKG_TARGET} ${PIP_LAST_ARG}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install python dependencies\")
    endif()
  ")
endfunction()

function(clean_and_optional_link)
  install(CODE "
    execute_process(COMMAND bash -c \"find ${FSPYTHON_INSTALL_PREFIX}/python/bin ! -name 'python3*' ! -name 'ds' -type f -exec rm -f {} +\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not remove pip-installed scripts\")
    endif()

    # Create the symbolic link after pip install
    if(EXISTS ${FSPYTHON_INSTALL_PREFIX}/python/bin/ds)
       execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink \"${FSPYTHON_INSTALL_PREFIX}/python/bin/ds\" \"${FSPYTHON_INSTALL_PREFIX}/bin/ds\")
    else()
       message(STATUS \"Cannot create soft link because ${FSPYTHON_INSTALL_PREFIX}/python/bin/ds was not installed by pip install deepsurfer\")
    endif()
  ")
endfunction()

function(fspython_patch)
   # Use to workaround a python module not installing, e.g., problem with setup.py.
   # Enable with -DPATCH_FSPYTHON=ON as arg to cmake command to copy install from previous build.
   # The path to the previous build must be available/mounted in the build container.
   set(FSPYTHON_PREVIOUS_BUILD_BASE "/space/freesurfer/build/garnet/sandboxes/utility/freesurfer/dev.20221207")
   install(CODE "
     message(STATUS \"patching fspython distribution\")
     execute_process(COMMAND bash -c \" ( cd ${FSPYTHON_PREVIOUS_BUILD_BASE}/python/lib/python${FSPYTHON_VERSION}/site-packages && tar cpf - surfa* xxhash* ) | ( cd ${FSPYTHON_INSTALL_PREFIX}/python/lib/python${FSPYTHON_VERSION}/site-packages && tar xpf - ) \" RESULT_VARIABLE patch)
     if(NOT \${patch} STREQUAL 0)
       message(FATAL_ERROR \"Could not patch fspython\")
     endif()"
   )
endfunction()

function(darwin_arm64_py_install)
   install(CODE "
     message(STATUS \"Installing darwin_arm64 miniforge3 python ${FSPYTHON_VERSION} for fspython\")
     message(STATUS \"( cd $HOME/space_freesurfer && bash -x fs_install_darwin_arm64_conda ${FSPYTHON_VERSION} ${FSPYTHON_INSTALL_PREFIX} )\")
     execute_process(COMMAND bash -c \" ( cd $HOME/space_freesurfer && bash -x fs_install_darwin_arm64_conda ${FSPYTHON_VERSION} ${FSPYTHON_INSTALL_PREFIX} ) \" RESULT_VARIABLE retcode)
     if(NOT \${retcode} STREQUAL 0)
       message(FATAL_ERROR \"Could not install darwin_arm64 fspython ${FSPYTHON_VERSION}\")
     endif()"
   )
endfunction()

# export FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} && 

function(prune_cuda)
  if(NOT APPLE)
     install(CODE "
       message(STATUS \"FSPYTHON_VERSION=${FSPYTHON_VERSION} FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --torchcpu\")
       execute_process(COMMAND bash -c \"FSPYTHON_VERSION=${FSPYTHON_VERSION} FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} FSPYTHON_INSTALL_PREFIX=${FSPYTHON_INSTALL_PREFIX} bash ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --torchcpu\" RESULT_VARIABLE retcode)
       if(NOT \${retcode} STREQUAL 0)
         message(FATAL_ERROR \"Could not replace torch module\")
       endif()
       message(STATUS \"FSPYTHON_VERSION=${FSPYTHON_VERSION} FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --uninstall\")
       execute_process(COMMAND bash -c \"FSPYTHON_VERSION=${FSPYTHON_VERSION} FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} FSPYTHON_INSTALL_PREFIX=${FSPYTHON_INSTALL_PREFIX} bash ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --uninstall\" RESULT_VARIABLE retcode)
       if(NOT \${retcode} STREQUAL 0)
         message(FATAL_ERROR \"Could not uninstall python modules\")
       endif()
     ")
   else()
      message(STATUS \"Nothing to prune for Cuda in fspython on MacOS\")
   endif()
endfunction()

function(integrate_samseg)
  # setup STANDALONE_SAMSEG_PATH and pip install flags
  # set(PIP_INSTALL_FLAGS "${PIP_NODEP_ARG} ${PIP_NOVER_ARG} --upgrade")
  set(PIP_INSTALL_FLAGS "${PIP_NOVER_ARG} --upgrade")
  ## global now as may use with req.sh - can install samseg at same time with all other python pkgs
  # set(requirements_samseg "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-samseg.txt")
  set(STANDALONE_SAMSEG_PATH ${requirements_samseg})

  # set samseg package install location
  if(NOT "${USER_SAMSEG_PATH}" STREQUAL "")
    set(SAMSEG_INSTALL_DIR "${FSPYTHON_INSTALL_PREFIX}/python/lib/python${FSPYTHON_VERSION}/site-packages")
    set(STANDALONE_SAMSEG_PATH "${USER_SAMSEG_PATH}")
    # for editable install to work, the .pth needs to be python/lib/pythonX.Y/site-packages
    # editable samseg install doesn't work in the following build environments:
    # 1. if DISTRIBUTE_FSPYTHON=OFF,
    #       module samseg cannot be found
    # 2. if DISTRIBUTE_FSPYTHON=ON and INSTALL_PYTHON_DEPENDENCIES=OFF,
    #       samseg dependencies modules cannot be found
    set(PIP_INSTALL_FLAGS "${PIP_INSTALL_FLAGS} --editable")
  endif()
 
  ## make global now to use with req.sh
  # set(pybind11_DIR           "${CMAKE_SOURCE_DIR}/packages/pybind11")

  if(APPLE)
    set(APPLE_BUILD_SETTING "APPLE=ON APPLE_ARM64=${APPLE_ARM64}")
  endif()
  set(INTEGRATE_SAMSEG_COMMAND "ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} ${APPLE_BUILD_SETTING} ${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} ${PIP_INSTALL_FLAGS} ${STANDALONE_SAMSEG_PATH} ${PKG_TARGET} ${PIP_LAST_ARG}")

  # pip install samseg via separate command
  install(CODE "
    message(STATUS \" for HOST_OS=${HOST_OS} \")
    message(STATUS \" ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} \")
    message(STATUS \" PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} PYTHON_INSTALL=${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 PKG_TARGET=${PKG_TARGET} \")

    # all dependencies should be already satisfied at this point
    if(NOT ${USER_SAMSEG_PATH} STREQUAL \"\")
      message(STATUS \" removing previous install ${SAMSEG_INSTALL_DIR}/*samseg* ... \")
      execute_process(COMMAND bash -c \"rm -rf ${SAMSEG_INSTALL_DIR}/*samseg*\")
    endif()
    message(STATUS \" installing standalone samseg ... \")
    message(STATUS \" ${INTEGRATE_SAMSEG_COMMAND} \")
    execute_process(COMMAND bash -c \"${INTEGRATE_SAMSEG_COMMAND}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install Standalone Samseg\")
    endif()"
  )
endfunction()

function(rm_fspython_packages_bin)
   # there should be no python/packages/bin subdriectory after a package install
   install(CODE "
     message(STATUS \"nuke python/packages/bin\")
     execute_process(COMMAND bash -c \" ( cd ${FSPYTHON_INSTALL_PREFIX}/python/packages && rm -rf bin) \" RESULT_VARIABLE patch)
     if(NOT \${patch} STREQUAL 0)
       message(FATAL_ERROR \"Could not rm python/pakcages/bin\")
     endif()"
   )
endfunction()

## CASE 1
# This is the default for all linux and darwin_x86_64 builds where all python packages listed
# in the requirements file should be installed under ./site-packages.  These packages should
# not be duplicated under python/packages.
if(DISTRIBUTE_FSPYTHON AND INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS "For HOST_OS=${HOST_OS} CASE 1: DISTRIBUTE_FSPYTHON ON AND INSTALL_PYTHON_DEPENDENCIES ON")
   # install base python distribution
   if("${FSPYTHON_VERSION}" MATCHES "3.8")
      python_3_8_base_install()
   elseif("${FSPYTHON_VERSION}" MATCHES "3.9")
      python_3_9_base_install()
   else()
      message(FATAL_ERROR "Don\'t know how to build with FSPYTHON_VERSION=${FSPYTHON_VERSION}.")
   endif()
   # install C bindings under python/packages
   packages_install()
   # install packages with pip using requirements file under python/lib/python<rev>/site-packages
   if(NOT APPLE)
      sitepackages_install_rev2()
   else()
      sitepackages_install()
   endif()
   clean_and_optional_link()
   # add fspython wrapper
   fspython_script_install()
   # add copy of requirmeents file under python/lib/python<rev>/site-packages
   requirements_file_install()
   # remove any NVIDIA related files we should not distribute
   if(NOT FSPYTHON_INSTALL_CUDA)
      prune_cuda()
   endif()
   if(INTEGRATE_SAMSEG)
      integrate_samseg()
      # rm_fspython_packages_bin()
   endif()    
   if(PATCH_FSPYTHON)
      fspython_patch()
   endif()
endif()

## CASE 2
# useful for basic debugging w/o anything pre-installed under site-packages
if(DISTRIBUTE_FSPYTHON AND NOT INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS "For HOST_OS=${HOST_OS} CASE 2: DISTRIBUTE_FSPYTHON ON AND INSTALL_PYTHON_DEPENDENCIES OFF")
   # install base python distribution
   if("${FSPYTHON_VERSION}" STREQUAL "3.8")
      python_3_8_base_install()
   elseif("${FSPYTHON_VERSION}" STREQUAL "3.9")
      python_3_9_base_install()
   else()
      message(FATAL_ERROR "Don't know how to build with FSPYTHON_VERSION=${FSPYTHON_VERSION}.")
   endif()
   # install C bindings under python/packages
   packages_install()
   #
   # NO PACKAGES INSTALLED VIA PIP AND REQUIREMENTS FILES UNDER python/lib/python<rev>/site-packages
   # 
   if(INTEGRATE_SAMSEG)
      integrate_samseg()
      # rm_fspython_packages_bin()
   endif()
   # add fspython wrapper
   fspython_script_install()
endif()

## CASE 3
# use what is normally installed for under site-packages with an external/3rd party, e.g., anaconda, python distribution
if(NOT DISTRIBUTE_FSPYTHON AND INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS For "HOST_OS=${HOST_OS} CASE 3: DISTRIBUTE_FSPYTHON OFF AND INSTALL_PYTHON_DEPENDENCIES ON")
   #
   # NO BASE PYTHON DISTRIBUTION INSTALLED
   # 
   set(PKG_TARGET "${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages")
   message(STATUS For "HOST_OS=${HOST_OS} PKG_TARGET=${PKG_TARGET}")
   # install C bindings under python/packages
   packages_install()
   # install packages with pip using requirements file under python/lib/python<rev>/site-packages
   if(NOT APPLE)
      sitepackages_install_rev2()
   else()
      sitepackages_install()
   endif()
   clean_and_optional_link()
   # add fspython wrapper
   fspython_script_install()
   # add copy of requirmeents file under python/lib/python<rev>/site-packages
   requirements_file_install()
   # remove any NVIDIA related files we should not distribute
   if(NOT FSPYTHON_INSTALL_CUDA)
      prune_cuda()
   endif()
   if(INTEGRATE_SAMSEG)
      integrate_samseg()
      # rm_fspython_packages_bin()
   endif()
   if(PATCH_FSPYTHON)
      fspython_patch()
   endif()
endif()

## CASE 4
# do not install a base python distribution, but instead create a soft link from FSPYTHON_INSTALL_PREFIX/python/bin/python3 to the systemn python.
# do not install anything under site packages using pip and requirements files.
if(NOT DISTRIBUTE_FSPYTHON AND NOT INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS "For HOST_OS=${HOST_OS} CASE 4: DISTRIBUTE_FSPYTHON OFF AND INSTALL_PYTHON_DEPENDENCIES OFF")
   #
   # NO BASE PYTHON DISTRIBUTION INSTALLED
   #
   # install C bindings under python/packages
   packages_install()
   #   
   # NO PACKAGES INSTALLED VIA PIP AND REQUIREMENTS FILES UNDER python/lib/python<rev>/site-packages
   #
   if(NOT APPLE_ARM64)
      if(INTEGRATE_SAMSEG)
         set(PKG_TARGET "${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages")
         message(STATUS For "HOST_OS=${HOST_OS} PKG_TARGET=${PKG_TARGET}")
         integrate_samseg()
         # rm_fspython_packages_bin()
      endif()
   else()
      # However for darwin_arm64 MacOS (-DAPPLE_ARM64) replace soft link to system python under ./python/bin with a
      # link to the miniforge3 python and do a normal install via pip and requirements files under site-packages, i.e.,
      # THIS IS CASE 1 ABOVE ONLY FOR darwin_arm64
      #
      # install base python distribution
      # install packages with pip using requirements file under python/lib/python<rev>/site-packages
      message(STATUS "For HOST_OS=${HOST_OS} CASE 4: will instead install darwin_arm64 miniforge3 python")
      darwin_arm64_py_install()
      if(INTEGRATE_SAMSEG)
         integrate_samseg()
      endif()
   endif()
   # add fspython wrapper
   fspython_script_install()
endif()

# WIP: Custom packages to modify existing install - install over pre-existing package with --ignore-installed and rebuild (--no-cache-dir)
# surfa
set(SURFA_LOCAL_PATH "/absolute/path/to/local/surfa/code")
# ?? ifdef on ENV setting for SURFA_LOCAL_PATH, else use default PKG_TARGET setting for site-packages ??
# set(PKG_TARGET "${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages")
#  message(STATUS For "HOST_OS=${HOST_OS} PKG_TARGET=${PKG_TARGET}")

add_custom_target(setup_surfa)

file(GLOB_RECURSE surfa_dists ${FSPYTHON_INSTALL_PREFIX}/python/packages/surfa-*.dist-info)

foreach(dist ${surfa_dists})
  file(REMOVE_RECURSE ${dist})
endforeach()

# message(STATUS "Installing Surfa from local path: ${SURFA_LOCAL_PATH}")
add_custom_command(
  TARGET setup_surfa 
  COMMAND 
    ${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} --ignore-installed --no-cache-dir --upgrade -e ${SURFA_LOCAL_PATH} ${PKG_TARGET} ${PIP_LAST_ARG}
  WORKING_DIRECTORY 
  ${CMAKE_CURRENT_SOURCE_DIR}
)

