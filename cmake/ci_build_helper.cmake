# ci_build_helper.cmake
#
# Injected via -DCMAKE_PROJECT_freesurfer_INCLUDE so the sphere-registration
# tools (mris_register, mris_convert, mris_curvature, ...) build out-of-the-box
# on stock Ubuntu and macOS runners WITHOUT editing any FreeSurfer source.
#
# It addresses two issues seen when building against system packages with a
# modern toolchain:
#
#  1. System VTK (>=9) advertises an X11::X11 imported target that FreeSurfer
#     never explicitly find_package(X11)'s, so configuration aborts. We resolve
#     X11 here and provide the imported target if CMake's module did not.
#
#  2. Newer GCC/clang promote warnings (e.g. -Warray-compare in gca.cpp) to
#     errors under FreeSurfer's -Werror. add_compile_options() is applied AFTER
#     CMAKE_CXX_FLAGS on the command line, so -Wno-error wins and pre-existing
#     files compile. This changes no source and does not affect our own code.

find_package(X11)
if(X11_FOUND AND NOT TARGET X11::X11)
  add_library(X11::X11 INTERFACE IMPORTED)
  set_target_properties(X11::X11 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${X11_X11_INCLUDE_PATH}"
    INTERFACE_LINK_LIBRARIES "${X11_X11_LIB}")
endif()

# -Wno-error relaxes FreeSurfer's -Werror generally. The remaining entries
# downgrade warnings that Apple/LLVM clang promotes to errors BY DEFAULT
# (int-conversion, implicit-function-declaration, implicit-int,
# incompatible-function-pointer-types), which the vendored C packages
# (e.g. packages/nifti) trip on macOS. They are scoped to C compiles, and the
# clang-only spelling is scoped to clang, so the Linux gcc build is unaffected.
add_compile_options(-Wno-error)
add_compile_options($<$<COMPILE_LANGUAGE:C>:-Wno-error=int-conversion>)
add_compile_options($<$<COMPILE_LANGUAGE:C>:-Wno-error=implicit-function-declaration>)
add_compile_options($<$<COMPILE_LANGUAGE:C>:-Wno-error=implicit-int>)
add_compile_options(
  $<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:AppleClang,Clang>>:-Wno-error=incompatible-function-pointer-types>)

# Homebrew's ITK (and similar) record a full C SDK include dir
# (.../SDKs/MacOSX*.sdk/usr/include) in ITK_INCLUDE_DIRS. FreeSurfer adds those
# via include_directories(SYSTEM ${ITK_INCLUDE_DIRS}); when the active toolchain
# is a *different* SDK (e.g. Xcode vs CommandLineTools), that C usr/include
# shadows libc++'s <limits.h>/<math.h> and the C++ build fails. Override
# include_directories to drop any such SDK C include (the real ITK-5.x include
# is kept). The original command remains available as _include_directories.
macro(include_directories)
  set(_fs_inc_args)
  foreach(_a ${ARGV})
    if(_a MATCHES "\\.sdk/usr/include$" OR _a STREQUAL "/usr/include")
      message(STATUS "ci_build_helper: dropping C sysroot include: ${_a}")
    else()
      list(APPEND _fs_inc_args "${_a}")
    endif()
  endforeach()
  _include_directories(${_fs_inc_args})
endmacro()


# Minimal build: only mris_register (and the sibling surface tools used by the
# sphere->sphere.reg->std-mesh pipeline) are wanted. FreeSurfer's root adds
# hundreds of unrelated tool subdirectories (mri_robust_register, mri_probedicom,
# Fortran tools, GUIs, ...) whose missing deps break CMake's *generate* step even
# though we never build them. Override add_subdirectory so that, AT THE ROOT
# ONLY, just the directories mris_register needs are configured. Subdirectories
# added deeper (e.g. packages/* from packages/CMakeLists.txt) are untouched.
macro(add_subdirectory _fs_dir)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(_fs_keep packages utils mris_register mris_convert mris_curvature mris_make_template)
    list(FIND _fs_keep "${_fs_dir}" _fs_idx)
    if(_fs_idx GREATER -1)
      _add_subdirectory(${ARGV})
    else()
      message(STATUS "ci_build_helper: skipping root subdir not needed for mris_register: ${_fs_dir}")
    endif()
  else()
    _add_subdirectory(${ARGV})
  endif()
endmacro()
