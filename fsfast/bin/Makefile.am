#
# Makefile.am for bin
#

install-exec-local:
	$(mkinstalldirs) $(DESTDIR)$(prefix)/fsfast/bin; \
	d="$(DESTDIR)$(prefix)/fsfast/bin/"; \
	for p in *; do \
	val=`file $$p | grep script`;\
	if test "$$val" != "" ; then \
	  echo " $(INSTALL_SCRIPT) $$p $$d$$p"; \
	    $(INSTALL_SCRIPT) $$p $$d$$p; \
	fi; \
	done

uninstall-local:
	$(mkinstalldirs) $(DESTDIR)$(prefix)/fsfast/bin; \
	d="$(DESTDIR)$(prefix)/fsfast/bin/"; \
	for p in *; do \
	val=`file $$p | grep script`;\
	if test "$$val" != "" ; then \
	  echo " rm -f $$d$$p"; \
	  rm -f $$d$$p; \
	fi; \
	done
 
#release:
#	$(MAKE) install
#	@list=$(BIN_REMOVE_LIST); for p in $$list ; do \
#	  d="$(DESTDIR)$(prefix)/fsfast/bin/"; \
#	  if test -x $$d$$p; then \
#	  echo "rm -f $$d$$p"; \
#	  rm -f $$d$$p; \
#	  fi; \
#	done

#uninstall-release:
#	$(MAKE) uninstall

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
