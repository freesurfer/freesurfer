#
# Makefile.am
#
# licence files should be copied over to the install directory
foodir=$(prefix)/docs
foo_DATA=citation.minc license.freesurfer  license.mesa  copywrite.minc license.fsl license.minc \
  registration.notes license.gnu license.tcltk license.html_lib license.tix license.afni

# big pdf and .ps files are only installed when release is specified
.PHONY: release

release: recon1.ps recon2.ps FreeSurferManual.pdf 
	$(INSTALL_DATA) recon1.ps recon2.ps FreeSurferManual.pdf $(foodir)

uninstall-release:
	cd $(prefix)/docs; rm -f recon1.ps recon2.ps FreesurferManual.pdf
