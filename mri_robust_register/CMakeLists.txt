project(mri_robust_register)

if(APPLE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else()
  if(NOT LAPACK_LIBRARIES OR NOT BLAS_LIBRARIES OR NOT GFORTRAN_LIBRARIES)
    message(WARNING "NOT building robust register tools due to missing fortran libs")
    return()
  endif()
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -z muldefs")
endif()

include_directories(${FS_INCLUDE_DIRS})

set(FORTRAN_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${GFORTRAN_LIBRARIES} ${QUADMATH_LIBRARIES})

# mri_robust_register
add_executable(
  mri_robust_register
  mri_robust_register.cpp
  JointHisto.cpp
  CostFunctions.cpp
  MyMatrix.cpp
  MyMRI.cpp
  Quaternion.cpp
  Registration.cpp
  RegRobust.cpp
  RegPowell.cpp
)
target_link_libraries(mri_robust_register utils ${FORTRAN_LIBS})
add_help(mri_robust_register mri_robust_register.help.xml)
install(TARGETS mri_robust_register DESTINATION bin)

# mri_robust_template
add_executable(
  mri_robust_template
  mri_robust_template.cpp
  Registration.cpp
  RegRobust.cpp
  CostFunctions.cpp
  MyMatrix.cpp
  MyMRI.cpp
  Quaternion.cpp
  MultiRegistration.cpp
)
target_link_libraries(mri_robust_template utils ${FORTRAN_LIBS})
add_help(mri_robust_template mri_robust_template.help.xml)
add_test_script(NAME mri_robust_template_test SCRIPT test.sh)
install(TARGETS mri_robust_template DESTINATION bin)

# lta_diff
add_executable(
  lta_diff
  lta_diff.cpp
  Registration.cpp
  CostFunctions.cpp
  MyMatrix.cpp
  MyMRI.cpp
  Quaternion.cpp
)
target_link_libraries(lta_diff utils ${FORTRAN_LIBS})
add_help(lta_diff lta_diff.help.xml)
install(TARGETS lta_diff DESTINATION bin)

# mri_create_tests
add_executable(
  mri_create_tests
  mri_create_tests.cpp
  CostFunctions.cpp
  MyMatrix.cpp
  MyMRI.cpp
  Quaternion.cpp
)
target_link_libraries(mri_create_tests utils ${FORTRAN_LIBS})
install(TARGETS mri_create_tests DESTINATION bin)

# mri_gradient_info
add_executable(mri_gradient_info mri_gradient_info.cpp)
target_link_libraries(mri_gradient_info utils ${FORTRAN_LIBS})
install(TARGETS mri_gradient_info DESTINATION bin)
