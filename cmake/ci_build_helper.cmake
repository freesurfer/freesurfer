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
