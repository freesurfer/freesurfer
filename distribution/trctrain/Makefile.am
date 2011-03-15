#
# Makefile.am
#

foodir=$(prefix)/trctrain

install: trctraindata.tar.gz
	rm -Rf $(DESTDIR)$(prefix)/trctrain
	rm -Rf $(top_builddir)/distribution/trctrain/trc???
	rm -Rf $(top_builddir)/distribution/trctrain/README
	rm -Rf $(top_builddir)/distribution/trctrain/trainlist.txt
	gunzip -c $(top_srcdir)/distribution/trctrain/trctraindata.tar.gz | tar xvf -
	mkdir -p $(DESTDIR)$(prefix)/trctrain
	mv -v $(top_builddir)/distribution/trctrain/trc??? \
		$(DESTDIR)$(prefix)/trctrain
	mv -v $(top_builddir)/distribution/trctrain/README \
		$(DESTDIR)$(prefix)/trctrain
	mv -v $(top_builddir)/distribution/trctrain/trainlist.txt \
		$(DESTDIR)$(prefix)/trctrain
	chmod -R 775 $(DESTDIR)$(prefix)/trctrain

uninstall-hook:
	rm -Rf $(DESTDIR)$(prefix)/trctrain

EXTRA_DIST=trctraindata.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=

include $(top_srcdir)/Makefile.extra
