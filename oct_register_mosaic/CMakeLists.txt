project(oct_register_mosaic)

include_directories(${FS_INCLUDE_DIRS})

add_executable(oct_register_mosaic oct_register_mosaic.cpp)
target_link_libraries(oct_register_mosaic utils)

install(TARGETS oct_register_mosaic DESTINATION bin)
