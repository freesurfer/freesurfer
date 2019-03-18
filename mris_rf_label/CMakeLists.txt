project(mris_rf_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_rf_label mris_rf_label.cpp)
target_link_libraries(mris_rf_label utils)

install(TARGETS mris_rf_label DESTINATION bin)
