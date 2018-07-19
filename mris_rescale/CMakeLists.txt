project(mris_rescale)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_rescale mris_rescale.c)
target_link_libraries(mris_rescale utils)

install(TARGETS mris_rescale DESTINATION bin)
