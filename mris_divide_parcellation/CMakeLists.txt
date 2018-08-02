project(mris_divide_parcellation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_divide_parcellation mris_divide_parcellation.c)
target_link_libraries(mris_divide_parcellation utils)

install(TARGETS mris_divide_parcellation DESTINATION bin)
