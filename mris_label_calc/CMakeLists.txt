project(mris_label_calc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_label_calc mris_label_calc.cpp)
target_link_libraries(mris_label_calc utils)

install(TARGETS mris_label_calc DESTINATION bin)
