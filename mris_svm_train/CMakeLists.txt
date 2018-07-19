project(mris_svm_train)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_svm_train mris_svm_train.c)
target_link_libraries(mris_svm_train utils svm)

install(TARGETS mris_svm_train DESTINATION bin)
