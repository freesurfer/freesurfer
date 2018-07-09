project(mris_translate_annotation)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_translate_annotation mris_translate_annotation.c)
target_link_libraries(mris_translate_annotation utils)

install(TARGETS mris_translate_annotation DESTINATION bin)
