project(mri_voldiff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_voldiff mri_voldiff.cpp)
target_link_libraries(mri_voldiff utils)

install(TARGETS mri_voldiff DESTINATION bin)
