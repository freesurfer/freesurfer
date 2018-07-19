project(mris_merge_parcellations)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_merge_parcellations mris_merge_parcellations.c)
target_link_libraries(mris_merge_parcellations utils)

install(TARGETS mris_merge_parcellations DESTINATION bin)
