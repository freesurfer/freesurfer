# try to find the MNI library

find_library(MNI_LIBRARY NAMES minc volume_io netcdf
                         PATH_SUFFIXES lib)
