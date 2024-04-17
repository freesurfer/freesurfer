TEMPLATE = subdirs
CONFIG += static
CONFIG -= shared
SUBDIRS = retrospective_correction connected_components \
  dissection_photo \
  fiducials_calibration \
  fiducials_correction \
  mask_extraction
