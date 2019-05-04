project(mris_merge_parcellations)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_merge_parcellations mris_merge_parcellations.cpp)
target_link_libraries(mris_merge_parcellations utils)

install(TARGETS mris_merge_parcellations DESTINATION bin)
