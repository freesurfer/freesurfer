#
# Makefile.am
#

SUBDIRS=average lib docs \
	fsafd subjects sessions diffusion fsl-extra trctrain

foodir=$(prefix)
foo_DATA= \
	tkmeditParcColorsCMA \
	Simple_surface_labels2005.txt	\
	Simple_surface_labels2009.txt	\
	FreeSurferColorLUT.txt \
	WMParcStatsLUT.txt \
	numaris4-protocols.unpackcfg \
	ASegStatsLUT.txt \
	DefectLUT.txt \
	SegmentNoLUT.txt

foo_SCRIPTS= \
	FreeSurferEnv.csh \
	SetUpFreeSurfer.csh \
	FreeSurferEnv.sh \
	SetUpFreeSurfer.sh

clean-local:
	rm -f SetUpFreeSurfer.csh SetUpFreeSurfer.sh

EXTRA_DIST=$(foo_DATA) $(foo_SCRIPTS) \
	SetUpFreeSurfer-SAMPLE.csh SetUpFreeSurfer-SAMPLE.sh

SetUpFreeSurfer.csh: SetUpFreeSurfer-SAMPLE.csh
	cat $< > $@
	chmod a+x $@

SetUpFreeSurfer.sh: SetUpFreeSurfer-SAMPLE.sh
	cat $< > $@
	chmod a+x $@

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=.license .xdebug_tk* error.log surfer.log csurf.log
include $(top_srcdir)/Makefile.extra
