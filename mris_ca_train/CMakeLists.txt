project(mris_ca_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_ca_train mris_ca_train.c)
target_link_libraries(mris_ca_train utils)

install(TARGETS mris_ca_train DESTINATION bin)
