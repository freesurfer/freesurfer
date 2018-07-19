project(mri_otl)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(list_otl_labels list_otl_labels.c)
target_link_libraries(list_otl_labels utils)
install(TARGETS list_otl_labels DESTINATION bin)
