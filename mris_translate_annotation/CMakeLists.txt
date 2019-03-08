project(mris_translate_annotation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_translate_annotation mris_translate_annotation.cpp)
target_link_libraries(mris_translate_annotation utils)

install(TARGETS mris_translate_annotation DESTINATION bin)
