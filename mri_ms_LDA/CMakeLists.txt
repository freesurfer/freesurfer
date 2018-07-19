project(mri_ms_LDA)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_ms_LDA mri_ms_LDA.c)
target_link_libraries(mri_ms_LDA utils)

install(TARGETS mri_ms_LDA DESTINATION bin)
