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

# -Wno-error relaxes FreeSurfer's -Werror; the -Wno-error=<name> entries
# additionally downgrade warnings that Apple/LLVM clang promotes to errors BY
# DEFAULT (int-conversion, implicit-function-declaration, implicit-int), which
# the vendored C packages (e.g. packages/nifti) still trip. All are accepted
# (and harmless) by the Linux gcc toolchain too.
add_compile_options(
  -Wno-error
  -Wno-error=int-conversion
  -Wno-error=implicit-function-declaration
  -Wno-error=implicit-int
  -Wno-error=incompatible-function-pointer-types)
