project(mris_glm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_glm mris_glm.c)
target_link_libraries(mris_glm utils)
install(TARGETS mris_glm DESTINATION bin)
