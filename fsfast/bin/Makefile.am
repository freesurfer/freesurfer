#
# Makefile.am for bin
#
# no file list is given. Therefore I have to check whether a file is a script or not
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
 
include exclude_files.sh

release:
	$(MAKE) install
	@list=$(BIN_REMOVE_LIST); for p in $$list ; do \
	  d="$(DESTDIR)$(prefix)/fsfast/bin/"; \
	  if test -x $$d$$p; then \
	  echo "rm -f $$d$$p"; \
	  rm -f $$d$$p; \
	  fi; \
	done

uninstall-release:
	$(MAKE) uninstall
