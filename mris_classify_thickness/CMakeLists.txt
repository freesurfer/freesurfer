project(mris_classify_thickness)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_classify_thickness mris_classify_thickness.c)
target_link_libraries(mris_classify_thickness utils)

install(TARGETS mris_classify_thickness DESTINATION bin)
