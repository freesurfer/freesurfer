#
# Makefile.am for fsfast/toolbox
#
# here even though file list is not given, but .m is enough to discriminate files
#
foodir=$(prefix)/fsfast/toolbox
foo_DATA= $(wildcard *.m)

include exclude_files.sh

release:
	$(MAKE) install
	@list=$(REMOVE_LIST); for p in $$list ; do \
	  d="$(DESTDIR)$(prefix)/fsfast/toolbox/"; \
	  if test -r $$d$$p; then \
	  echo "rm -f $$d$$p"; \
	  rm -f $$d$$p; \
	  fi; \
	done

uninstall-release:
	$(MAKE) uninstall


