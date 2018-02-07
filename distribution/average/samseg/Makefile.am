#
# Makefile.am
#

foodir=$(prefix)/average/samseg

foo_tarballs=20Subjects_smoothing2_down2_smoothingForAffine2.tar.gz

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	for tarball in $(foo_tarballs); do \
		gunzip -c $(top_srcdir)/distribution/average/samseg/$$tarball | tar xvf - -C $(foodir);\
	done

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
