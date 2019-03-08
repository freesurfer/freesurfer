project(mris_entropy)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_entropy mris_entropy.cpp)
target_link_libraries(mris_entropy utils)

install(TARGETS mris_entropy DESTINATION bin)
