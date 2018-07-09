project(mris_glm)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_glm mris_glm.c)
target_link_libraries(mris_glm utils)

install(TARGETS mris_glm DESTINATION bin)
