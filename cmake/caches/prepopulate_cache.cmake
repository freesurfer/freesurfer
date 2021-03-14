set(CMAKE_INSTALL_PREFIX
    "/Users/aboualiaa/Desktop/.temp/install/fs/"
    CACHE STRING "install prefix"
    )

set(CMAKE_MESSAGE_LOG_LEVEL
    "WARNING"
    CACHE STRING ""
    )

set(CMAKE_BUILD_TYPE
    "Release"
    CACHE STRING "Default: Release"
    )

set(CMAKE_GENERATOR
    "Ninja Multi-Config"
    CACHE STRING "Ninja, Unix Makefiles, Xcode, Ninja Multi-Config"
    )

string(TIMESTAMP TODAY "%Y%m%d")

set(BUILD_STAMP
    "freesurfer-local-build-${TODAY}"
    CACHE STRING "Distribution build stamp"
    )

set(CMAKE_EXPORT_COMPILE_COMMANDS
    ON
    CACHE BOOL "create a json database of compile commands for tooling"
    )

set(CMAKE_CXX_STANDARD
    20
    CACHE STRING "20;17;14;11"
    )

set(CMAKE_POSITION_INDEPENDENT_CODE
    ON
    CACHE BOOL ""
    )

set(FS_BUILD_ATTIC
    ON
    CACHE BOOL "Build deprecated tools from the attic"
    )

set(FS_BUILD_DNG
    ON
    CACHE BOOL "Build Doug's testing tools"
    )

set(FS_BUILD_GUIS
    ON
    CACHE BOOL "Build GUI tools"
    )

set(FS_BUILD_DOCS
    ON
    CACHE BOOL "Build Documentation"
    )

set(FS_BUILD_TESTING
    ON
    CACHE BOOL "build test targets"
    )

set(FS_PACKAGES_DIR
    "/Users/aboualiaa/Downloads/packages"
    CACHE STRING ""
    )

set(CMAKE_INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/lib"
    CACHE STRING "for shared Libraries"
    )

# TODO: implement
set(FS_OSX_BUILD_UNIVERSAL_BINARIES
    OFF
    CACHE BOOL "sets CMAKE_OSX_ARCHITECTURES=arm64;x86_64"
    )

# Unfortunately, the python version used to run pybind c-libraries must be equivalent to
# the version used to build the libraries. The easiest and least intrusive way of making freesurfer
# python scripts run out-of-the-box (and to help guarantee reproducibility) requires
# distributing a minimal, custom python installation called fspython. This fspython package
# will get installed to "freesurfer/python", but an external python can be used instead when the
# FS_DISTRIBUTE_FSPYTHON option is turned off. Turning this off will create a freesurfer/bin/fspythonlink
# symlink during install that points to the external python executable located by pybind
set(FS_DISTRIBUTE_FSPYTHON
    OFF
    CACHE BOOL "Include the fspython distribution in the installation"
    )

set(FS_COVERAGE_STYLE
    "all"
    CACHE STRING ""
    )

set(FS_DOCS_FORMAT
    "all"
    CACHE STRING ""
    )

set(FS_DOCS_GENERATOR
    "all"
    CACHE STRING ""
    )

set(FS_EXP_BUILD_CONFIGURATIONS
    "Release;Debug;RelWithDebInfo;MinSizeRel;asan;ubsan;tsan;msan;sssan;cfisan;coverage;profile;lsan;thinlto;fulllto"
    CACHE STRING ""
    )

set(CMAKE_CONFIGURATION_TYPES ${FS_EXP_BUILD_CONFIGURATIONS})

set(FS_ENABLE_LTO
    "Full"
    CACHE STRING "off thin full"
    )

set(FS_INFANT_MODULE
    ON
    CACHE BOOL "Include infant recon-all"
    )

set(FS_QATOOLS_MODULE
    OFF
    CACHE BOOL "Include quality assurance tools"
    )

set(FS_FREEVIEW_LINEPROF
    OFF
    CACHE BOOL "Enable Lineprof"
    )

set(FS_INSTALL_PYTHON_DEPS
    ON
    CACHE BOOL "python deps"
    )

set(FS_INTEGRATION_TESTING
    ON
    CACHE BOOL "Copy files for integration tests"
    )

set(FS_MINIMAL_BUILD
    OFF
    CACHE BOOL "Only build core components"
    )

set(FS_REPO_ENVIRONMENT
    "develop"
    CACHE STRING ""
    )

set(FS_SUPRESS_WARNINGS
    ON
    CACHE BOOL "Suppress some selected warnings"
    )

set(FS_USE_CCACHE
    ON
    CACHE BOOL "Use ccache (if present) to reduce build times"
    )

set(FS_USE_GCC
    ""
    CACHE STRING ""
    )

set(FS_USE_LLVM
    ""
    CACHE STRING ""
    )

set(FS_VERSION
    "$ENV{USER}-local"
    CACHE STRING "Distribution version"
    )

if(NOT APPLE)

  set(FS_GEMS_BUILD_CUDA
      OFF
      CACHE BOOL "Build CUDA stuff"
      )

  set(FS_TKTOOLS_MODULE
      OFF
      CACHE BOOL "Install old Linux TK GUIs"
      )

endif()

set(FS_BUILD_MATLAB
    ON
    CACHE BOOL "build matlab"
    )

set(FS_GEMS_BUILD_EXECUTABLES
    ON
    CACHE BOOL "Build command line executables"
    )

set(FS_GEMS_BUILD_GUI
    ON
    CACHE BOOL "build gems guis"
    )

set(FS_GEMS_BUILD_MATLAB
    ON
    CACHE BOOL "build matlab wrappers"
    )

set(FS_GEMS_BUILD_PYTHON
    ON
    CACHE BOOL ""
    )

set(FS_GEMS_BUILD_SHARED_LIBS
    ON
    CACHE BOOL "Build GEMS with shared libraries"
    )

set(FS_GEMS_BUILD_TESTING
    ON
    CACHE BOOL "Build tests"
    )

set(FS_GEMS_MAKE_SPARSE_INITIAL_MESHES
    OFF
    CACHE BOOL "Make sparse initial meshes"
    )

set(FS_USE_ARRAYFIRE
    ON
    CACHE BOOL "Use Arrayfire for parallelism"
    )

set(FS_USE_HPX
    ON
    CACHE BOOL "Use Arrayfire for parallelism"
    )

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${FS_EXP_BUILD_CONFIGURATIONS})
set_property(CACHE CMAKE_GENERATOR PROPERTY STRINGS "Ninja;Unix Makefiles;Xcode;Ninja Multi-Config")
set_property(CACHE CMAKE_CXX_STANDARD PROPERTY STRINGS "20;17;14;11")
set_property(CACHE FS_REPO_ENVIRONMENT PROPERTY STRINGS "develop;travis;azure;cirrus;deploy;martinos")
set_property(CACHE FS_DOCS_FORMAT PROPERTY STRINGS "html;md;yaml;man;dash;all")
set_property(CACHE FS_DOCS_GENERATOR PROPERTY STRINGS "doxygen;sphinx;clang-doc;all")
set_property(CACHE FS_ENABLE_LTO PROPERTY STRINGS "Off;Thin;Full")
set_property(CACHE FS_COVERAGE_STYLE PROPERTY STRINGS "gcov;llvm;all")

#set(CMAKE_CONFIGURATION_TYPES ${FS_EXP_BUILD_CONFIGURATIONS} CACHE STRING "")
