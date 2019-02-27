project(mri_relabel_nonwm_hypos)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_relabel_nonwm_hypos mri_relabel_nonwm_hypos.cpp)
target_link_libraries(mri_relabel_nonwm_hypos utils)

install(TARGETS mri_relabel_nonwm_hypos DESTINATION bin)
