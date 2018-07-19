project(mris_extract_patches)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_extract_patches mris_extract_patches.c)
target_link_libraries(mris_extract_patches utils)

install(TARGETS mris_extract_patches DESTINATION bin)
