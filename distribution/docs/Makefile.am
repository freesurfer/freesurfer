#
# Makefile.am
#
# licence files should be copied over to the install directory
foodir=$(prefix)/docs
foo_DATA=license.freesurfer \
	citation.minc \
	copywrite.minc \
	license.minc \
	registration.notes \
	license.mesa \
	license.fsl \
	license.gnu \
	license.tcltk \
	license.html_lib \
	license.tix \
	license.afni

# big pdf and .ps files are only installed when release is specified
.PHONY: release

release: articles.html \
	recon1.pdf \
	recon2.pdf \
	morphing.pdf \
	cortical-thickness.pdf \
	topology-fixing.pdf \
	labeling.pdf \
	parcellation.pdf
	$(INSTALL_DATA) articles.html \
	recon1.pdf \
	recon2.pdf \
	morphing.pdf \
	cortical-thickness.pdf \
	topology-fixing.pdf \
	labeling.pdf \
	parcellation.pdf \
	$(foodir)

uninstall-release:
	cd $(prefix)/docs; rm -f articles.html \
	recon1.pdf \
	recon2.pdf \
	morphing.pdf \
	cortical-thickness.pdf \
	topology-fixing.pdf \
	labeling.pdf \
	parcellation.pdf

