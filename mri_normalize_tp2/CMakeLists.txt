project(mri_normalize_tp2)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_normalize_tp2 mri_normalize_tp2.c)
target_link_libraries(mri_normalize_tp2 utils)

install(TARGETS mri_normalize_tp2 DESTINATION bin)
