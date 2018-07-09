project(mris_rf_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_rf_train mris_rf_train.c)
target_link_libraries(mris_rf_train utils)

install(TARGETS mris_rf_train DESTINATION bin)
