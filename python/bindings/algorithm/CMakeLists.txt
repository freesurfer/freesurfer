pybind11_add_module(
  algorithm_python
  python.cpp
  cmf3d_ml.cpp
  fib.cpp
  perform_front_propagation_3d.cpp
)

set_target_properties(algorithm_python PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${FS_PYTHON_DIR}/algorithm)
