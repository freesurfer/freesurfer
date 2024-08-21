project(fspython)

# build the cxx/python bindings for the below modules
add_subdirectories(fsbindings)
if(NOT INTEGRATE_SAMSEG)
  add_subdirectories(gems/cxx)
endif()

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
   set(PIP_SUBSTRING "pip install")
   set(PIP_NODEP_ARG "--no-dependencies")
   set(PIP_NOVER_ARG "--disable-pip-version-check")
   set(PIP_LAST_ARG "")
   set(pip_target_option "--target")
endif()

function(requirements_file_install)
  if(FSPYTHON_BUILD_REQ)
      if(NOT APPLE_ARM64)
         if(FSPYTHON_INSTALL_TREE)
            install(PROGRAMS requirements-build.txt DESTINATION ../fspython/python)
         else()
            install(PROGRAMS requirements-build.txt DESTINATION python)
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
         install(PROGRAMS fspython.local DESTINATION ../fspython/bin RENAME fspython)
      endif() 
   else()
      install(PROGRAMS fspython DESTINATION bin)
      if(FSPYTHON_INSTALL_TREE)
         install(PROGRAMS fspython DESTINATION ../fspython/bin)
      endif()
   endif()
endfunction()

if(FSPYTHON_INSTALL_TREE)
   set(FSPYTHON_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/\.\./fspython" CACHE PATH "fspython install path" FORCE)
   # install(DIRECTORY DESTINATION ${FSPYTHON_INSTALL_PREFIX}) 
else()
   set(FSPYTHON_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH "fspython install path" FORCE)
endif()

# TO DO  - write a python_3_10_base_install

function(python_3_8_base_install)
  set(FSPYTHON_DIR ${FS_PACKAGES_DIR}/fspython/3.8)
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
  # default is least version specific or sparse requirements to allow for latest pkg revisions
  set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt")
  set(requirements_extra "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-extra.txt")
  if(FSPYTHON_BUILD_REQ)
     # pip install current snapshot of (hardcoded) package revisions into fspython
     if(NOT APPLE)
        ## all linux OS
        # set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-linux.txt")
        set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-build.txt")
        set(requirements_extra "")
     else()
        if(NOT APPLE_ARM64)
           ## darwin_x86_64
           # set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-darwin-intel.txt")
  	 message(STATUS "For HOST_OS=${HOST_OS} using default requirements files")
        else()
           ## darwin_arm64
  	 message(STATUS "For HOST_OS=${HOST_OS} - no requirements files currently in use with mininconda darwin_arm64.")
        endif()
     endif()
  endif()
  if(NOT APPLE_ARM64)
     if(INFANT_MODULE)
        # provide requirements file for infant module to allow for different package revisions compared to what
        # freesurfer uses (add cmake option -DFSPYTHON_BUILD_REQ=ON to infant options -DINFANT_MODULE=ON -DMINIMAL=ON)
        set(requirements "-r ${CMAKE_CURRENT_SOURCE_DIR}/requirements-infant.txt")
        set(requirements_extra "")
     endif()
  endif()
  if(NOT APPLE_ARM64)
     message(STATUS "For HOST_OS=${HOST_OS} fspython pip set to use requirement arguments: ${requirements} ${requirements_extra}")
  endif() 

  install(CODE "
    message(STATUS \"Installing python package dependencies into fspython with cmd:\")
    message(STATUS \"${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} --upgrade ${PIP_NOVER_ARG} ${EXTRA_PACKAGES} ${requirements} ${requirements_extra} ${PKG_TARGET} ${PIP_LAST_ARG}\")

    execute_process(COMMAND bash -c \"${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} --upgrade ${PIP_NOVER_ARG} ${EXTRA_PACKAGES} ${requirements} ${requirements_extra} ${PKG_TARGET} ${PIP_LAST_ARG}\" RESULT_VARIABLE retcode)

    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install python dependencies\")
    endif()

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
     execute_process(COMMAND bash -c \" ( cd ${FSPYTHON_PREVIOUS_BUILD_BASE}/python/lib/python3.8/site-packages && tar cpf - surfa* xxhash* ) | ( cd ${FSPYTHON_INSTALL_PREFIX}/python/lib/python3.8/site-packages && tar xpf - ) \" RESULT_VARIABLE patch)
     if(NOT \${patch} STREQUAL 0)
       message(FATAL_ERROR \"Could not patch fspython\")
     endif()"
   )
endfunction()

# export FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} && 

function(prune_cuda)
  if(NOT APPLE)
     install(CODE "
       message(STATUS \"FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --torchcpu\")
       execute_process(COMMAND bash -c \"FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --torchcpu\" RESULT_VARIABLE retcode)
       if(NOT \${retcode} STREQUAL 0)
         message(FATAL_ERROR \"Could not replace torch module\")
       endif()
       message(STATUS \"FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --uninstall\")
       execute_process(COMMAND bash -c \"FS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR}/req.sh --uninstall\" RESULT_VARIABLE retcode)
       if(NOT \${retcode} STREQUAL 0)
         message(FATAL_ERROR \"Could not uninstall python modules\")
       endif()
     ")
   else()
      message(STATUS \"Nothing to prune for Cuda in fspython on MacOS\")
   endif()
endfunction()

function(integrate_samseg)
  # pip install samseg
  set(PIP_INSTALL_FLAGS "${PIP_NODEP_ARG} ${PIP_NOVER_ARG} --upgrade")
  set(STANDALONE_SAMSEG_PATH "${CMAKE_CURRENT_BINARY_DIR}/git-samseg.standalone")
  if(NOT "${USER_SAMSEG_PATH}" STREQUAL "")
    set(STANDALONE_SAMSEG_PATH "${USER_SAMSEG_PATH}")
    # editable samseg install doesn't work in the following build environments:
    # 1. if DISTRIBUTE_FSPYTHON=OFF,
    #       module samseg cannot be found
    # 2. if DISTRIBUTE_FSPYTHON=ON and INSTALL_PYTHON_DEPENDENCIES=OFF,
    #       samseg dependencies modules cannot be found
    set(PIP_INSTALL_FLAGS "${PIP_INSTALL_FLAGS} --editable")
  endif()
  set(STANDALONE_SAMSEG_URL  "https://github.com/freesurfer/samseg.git")
  
  set(pybind11_DIR           "${CMAKE_SOURCE_DIR}/packages/pybind11")
        
  install(CODE "
    message(STATUS \" for HOST_OS=${HOST_OS} \")
    message(STATUS \" ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} \")
    message(STATUS \" PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} PYTHON_INSTALL=${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 PKG_TARGET=${PKG_TARGET} \")

    if(\"${USER_SAMSEG_PATH}\" STREQUAL \"\")
      message(STATUS \" git clone ${STANDALONE_SAMSEG_URL} ${STANDALONE_SAMSEG_PATH} \")
      execute_process(COMMAND bash -c \"git clone --quiet ${STANDALONE_SAMSEG_URL} ${STANDALONE_SAMSEG_PATH}\" result_variable retcode)
      if(NOT ${retcode} STREQUAL 0)
        message(FATAL_ERROR \"could not git clone standalone samseg\")
      endif()
    endif()

    # no check/install samseg dependencies, all dependencies will be handled in requirement files
    message(STATUS \" installing standalone samseg from ${STANDALONE_SAMSEG_PATH} \") 
    message(STATUS \" ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} CMAKE_C_COMPILER=${CMAKE_C_COMPILER} CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} APPLE_ARM64=${APPLE_ARM64} ${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} ${PIP_INSTALL_FLAGS} ${STANDALONE_SAMSEG_PATH} ${PKG_TARGET} ${PIP_LAST_ARG}\")
    execute_process(COMMAND bash -c \"ITK_DIR=${ITK_DIR} pybind11_DIR=${pybind11_DIR} CMAKE_C_COMPILER=${CMAKE_C_COMPILER} CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} APPLE_ARM64=${APPLE_ARM64} ${FSPYTHON_INSTALL_PREFIX}/python/bin/python3 -m ${PIP_SUBSTRING} ${PIP_INSTALL_FLAGS} ${STANDALONE_SAMSEG_PATH} ${PKG_TARGET} ${PIP_LAST_ARG}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not install Standalone Samseg\")
    endif()"
  )
endfunction()

# used as current default for intel linux and mac builds
if(DISTRIBUTE_FSPYTHON AND INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS "For HOST_OS=${HOST_OS} DISTRIBUTE_FSPYTHON ON AND INSTALL_PYTHON_DEPENDENCIES ON")
   # install base python distribution
   # install C bindings, etc. under python/packages
   # install packages under python/lib/python<rev>/site-packages, add fspython wrapper
   # install a copy of the requirements-build.txt file if the build was set to use it for the revs of pkgs to install
   python_3_8_base_install()
   packages_install()
   sitepackages_install()
   if(NOT FSPYTHON_INSTALL_CUDA)
      prune_cuda()
   endif()
   fspython_script_install()
   requirements_file_install()
   if(PATCH_FSPYTHON)
      fspython_patch()
   endif()
endif()

# useful for basic debugging
if(DISTRIBUTE_FSPYTHON AND NOT INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS "For HOST_OS=${HOST_OS} DISTRIBUTE_FSPYTHON ON AND INSTALL_PYTHON_DEPENDENCIES OFF")
   # install base python distribution
   # install C bindings, etc. under python/packages
   # NO pre-installed site-packages, add fspython wrapper
   python_3_8_base_install()
   packages_install()
   fspython_script_install()
endif()

# used for linking what's installed under site-packages to an external/3rd party, e.g., anaconda, python distribution
# FSPYTHON_INSTALL_PREFIX/python/bin/python3 soft link not used
if(NOT DISTRIBUTE_FSPYTHON AND INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS For "HOST_OS=${HOST_OS} DISTRIBUTE_FSPYTHON OFF AND INSTALL_PYTHON_DEPENDENCIES ON")
   # Use an existing 3rd party python instead of base python distribution
   # install C bindings, etc. under python/packages
   # install packages under python/lib/python<rev>/site-packages, add fspython wrapper
   # install a copy of the requirements-build.txt file if the build was set to use it for the revs of pkgs to install
   set(PKG_TARGET "${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages")
   message(STATUS For "HOST_OS=${HOST_OS} PKG_TARGET=${PKG_TARGET}")
   packages_install()
   sitepackages_install()
   if(NOT FSPYTHON_INSTALL_CUDA)
      prune_cuda()
   endif()
   fspython_script_install()
   requirements_file_install()
   if(PATCH_FSPYTHON)
      fspython_patch()
   endif()
endif()

# used for silicon Mac, -DAPPLE_ARM64, with download and install of external/3rd party miniconda arm64 python
# FSPYTHON_INSTALL_PREFIX/python/bin/python3 soft link changed to point to anaconda python
if(NOT DISTRIBUTE_FSPYTHON AND NOT INSTALL_PYTHON_DEPENDENCIES)
   message(STATUS "For HOST_OS=${HOST_OS} DISTRIBUTE_FSPYTHON OFF AND INSTALL_PYTHON_DEPENDENCIES OFF")
   # Use an existing 3rd party python instead of base python distribution
   # install C bindings, etc. under python/packages
   # NO pre-installed site-packages, add fspython wrapper
   # set PKG_TARGET, this is needed for integrate_samseg()
   set(PKG_TARGET "${pip_target_option}=${FSPYTHON_INSTALL_PREFIX}/python/packages")
   message(STATUS For "HOST_OS=${HOST_OS} PKG_TARGET=${PKG_TARGET}")
   packages_install()
   fspython_script_install()
endif()


if(INTEGRATE_SAMSEG)
  integrate_samseg()
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

