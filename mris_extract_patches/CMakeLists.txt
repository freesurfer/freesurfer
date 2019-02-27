project(mris_extract_patches)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_extract_patches mris_extract_patches.cpp)
target_link_libraries(mris_extract_patches utils)

install(TARGETS mris_extract_patches DESTINATION bin)
