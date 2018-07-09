project(mris_average_parcellation)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_average_parcellation mris_average_parcellation.c)
target_link_libraries(mris_average_parcellation utils)

install(TARGETS mris_average_parcellation DESTINATION bin)
