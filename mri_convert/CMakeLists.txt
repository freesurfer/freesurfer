project(mri_convert)

# ATH: temporarily turning off this warning
add_compile_options(-Wno-self-assign)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom ${CMAKE_SOURCE_DIR}/packages/dcm2niix)

add_executable(mri_convert mri_convert.cpp)
add_help(mri_convert mri_convert.help.xml)
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-dead_strip -Wl,-stack_size -Wl,0x1000000")
endif()
target_link_libraries(mri_convert utils)

install(TARGETS mri_convert DESTINATION bin)

add_test_script(NAME mri_convert_test SCRIPT test.sh)

add_executable(mri_make_uchar mri_make_uchar.cpp)
target_link_libraries(mri_make_uchar utils)
install(TARGETS mri_make_uchar DESTINATION bin)
