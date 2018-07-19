# NetCDF Find Module

if(NOT NetCDF_DIR)
  set(NetCDF_DIR ${FS_PACKAGES_DIR}/netcdf/3.6.0-p1)
endif()

find_path(NetCDF_INCLUDE_DIR PATHS ${NetCDF_DIR} NAMES netcdf.h PATH_SUFFIXES include)
find_library(NetCDF_LIBRARIES PATHS ${NetCDF_DIR} NAMES libnetcdf.a PATH_SUFFIXES lib)
find_package_handle_standard_args(NetCDF DEFAULT_MSG NetCDF_INCLUDE_DIR NetCDF_LIBRARIES)
