project(mris_convert)

include_directories(${FS_INCLUDE_DIRS})

# to set additional debug cxxflags:
#   export MRIS_CONVERT_DEBUG_CXXFLAG="-D__MRISURFOVERLAY_DEBUG"
#   touch CMakeLists.txt to trigger re-configuration
if(DEFINED ENV{MRIS_CONVERT_DEBUG_CXXFLAG})
  message(WARNING "additional mris_convert cxxflags: $ENV{MRIS_CONVERT_DEBUG_CXXFLAG}")
  add_compile_options($ENV{MRIS_CONVERT_DEBUG_CXXFLAG})
endif()

add_executable(mris_convert mris_convert.cpp)
add_help(mris_convert mris_convert.help.xml)
target_link_libraries(mris_convert utils)

add_test_script(NAME mris_convert_test SCRIPT test.sh)

install(TARGETS mris_convert DESTINATION bin)
