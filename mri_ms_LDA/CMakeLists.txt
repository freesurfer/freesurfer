project(mri_ms_LDA)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_ms_LDA mri_ms_LDA.cpp)
target_link_libraries(mri_ms_LDA utils)

install(TARGETS mri_ms_LDA DESTINATION bin)
