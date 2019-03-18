project(mris_ca_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_ca_train mris_ca_train.cpp)
target_link_libraries(mris_ca_train utils)

install(TARGETS mris_ca_train DESTINATION bin)
