project(mris_rf_label)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_rf_label mris_rf_label.c)
target_link_libraries(mris_rf_label utils)

install(TARGETS mris_rf_label DESTINATION bin)
