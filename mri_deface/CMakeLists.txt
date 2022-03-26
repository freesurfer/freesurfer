project(mri_deface)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_deface mri_deface.cpp)
target_link_libraries(mri_deface utils)

add_executable(mri_defacer mri_defacer.cpp)
target_link_libraries(mri_defacer utils)

install(TARGETS mri_deface mri_defacer DESTINATION bin)
