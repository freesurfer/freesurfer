#
# Makefile.am
#
foodir=$(prefix)/trctrain

all:
	cd $(top_srcdir)/distribution/trctrain
	rm -Rf README trainlist.txt trc???
	gunzip -c trctraindata.tar.gz | tar xvf -

install:
	rm -Rf $(DESTDIR)$(prefix)/trctrain
	cp -R $(top_srcdir)/distribution/trctrain $(DESTDIR)$(prefix)
	chmod -R 775 $(DESTDIR)$(prefix)/trctrain

uninstall-hook:
	rm -Rf $(DESTDIR)$(prefix)/trctrain

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=

include $(top_srcdir)/Makefile.extra
