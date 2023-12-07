project(mri_ctab_fix)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ctab_fix mri_ctab_fix.cpp)
target_link_libraries(mri_ctab_fix utils)

install(TARGETS mri_ctab_fix DESTINATION bin)
