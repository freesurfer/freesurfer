project(mris_svm_train)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_svm_train mris_svm_train.cpp)
target_link_libraries(mris_svm_train utils svm)

install(TARGETS mris_svm_train DESTINATION bin)
