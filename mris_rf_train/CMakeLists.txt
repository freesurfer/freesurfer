project(mris_rf_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_rf_train mris_rf_train.cpp)
target_link_libraries(mris_rf_train utils)

install(TARGETS mris_rf_train DESTINATION bin)
