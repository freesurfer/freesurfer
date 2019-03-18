project(mris_rescale)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_rescale mris_rescale.cpp)
target_link_libraries(mris_rescale utils)

install(TARGETS mris_rescale DESTINATION bin)
