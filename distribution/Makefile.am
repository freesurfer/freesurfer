#
# Makefile.am
#

if SMALL_DIST_INSTALL
SUBDIRS=average lib docs fsafd subjects
else
SUBDIRS=average lib docs \
	fsafd subjects sessions diffusion fsl-extra trctrain
endif

foodir=$(prefix)
foo_DATA= \
	tkmeditParcColorsCMA \
	Simple_surface_labels2005.txt	\
	Simple_surface_labels2009.txt	\
	FreeSurferColorLUT.txt \
	WMParcStatsLUT.txt \
	ASegStatsLUT.txt \
	DefectLUT.txt \
	SegmentNoLUT.txt

foo_SCRIPTS= \
	FreeSurferEnv.csh \
	SetUpFreeSurfer.csh \
	sources.csh \
	FreeSurferEnv.sh \
	SetUpFreeSurfer.sh \
	sources.sh

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=.license .xdebug_tk* error.log surfer.log csurf.log
include $(top_srcdir)/Makefile.extra
