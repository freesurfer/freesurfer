project(oct_register_mosaic)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(oct_register_mosaic oct_register_mosaic.c)
target_link_libraries(oct_register_mosaic utils)

install(TARGETS oct_register_mosaic DESTINATION bin)
