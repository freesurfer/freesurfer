project(mris_average_parcellation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_average_parcellation mris_average_parcellation.cpp)
target_link_libraries(mris_average_parcellation utils)

install(TARGETS mris_average_parcellation DESTINATION bin)
