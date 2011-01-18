#
# Makefile.am
#

SUBDIRS=

foodir=$(prefix)/fsafd
foo_DATA=\
	TalairachingCovariance.adf \
	TalairachingMean.adf \
	TalairachingProbas.adf

# not yet used:
#	lh.CorpusCallosumCutDice.adf \
#	lh.CorticalLabelingPercentArea_txt.adf \
#	lh.SurfRegistrationResidualVar.adf \
#	PonsCutDice.adf \
#	rh.CorpusCallosumCutDice.adf \
#	rh.CorticalLabelingPercentArea_txt.adf \
#	RibbonPlacementDice.adf \
#	SubcorticalLabelingPercentVol.adf \
#	WMsegmentationDice.adf \
#	lh.EV90.mgh \
#	lh.ThicknessMatrix.mgh

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
