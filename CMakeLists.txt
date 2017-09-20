project(freesurfer)
cmake_minimum_required(VERSION 2.6)

SET(FS_CORE utils rgb dicom hipsstubs log image unix utilscpp svm-lib)
SET(FS_LIBS ${FS_AUX_LIBS} ${FS_CORE})
