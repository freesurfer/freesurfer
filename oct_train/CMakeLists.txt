project(oct_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(oct_train oct_train.cpp)
target_link_libraries(oct_train utils)
install(TARGETS oct_train DESTINATION bin)

add_executable(oct_rf_train oct_rf_train.cpp)
target_link_libraries(oct_rf_train utils)
install(TARGETS oct_rf_train DESTINATION bin)
