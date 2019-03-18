project(mri_normalize_tp2)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_normalize_tp2 mri_normalize_tp2.cpp)
target_link_libraries(mri_normalize_tp2 utils)

install(TARGETS mri_normalize_tp2 DESTINATION bin)
