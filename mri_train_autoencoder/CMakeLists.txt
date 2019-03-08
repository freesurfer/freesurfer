project(mri_train_autoencoder)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_train_autoencoder mri_train_autoencoder.cpp)
target_link_libraries(mri_train_autoencoder utils)
install(TARGETS mri_train_autoencoder DESTINATION bin)

add_executable(mri_apply_autoencoder mri_apply_autoencoder.cpp)
target_link_libraries(mri_apply_autoencoder utils)
install(TARGETS mri_apply_autoencoder DESTINATION bin)
