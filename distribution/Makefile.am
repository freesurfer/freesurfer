#
# Makefile.am
#
SUBDIRS=. average lib docs

foodir=$(prefix)
foo_DATA=surface_labels.txt tkmeditColorsCMA tkmeditParcColorsCMA Simple_surface_labels2002.txt FreeSurferEnv.csh SetUpFreeSurfer.csh


SetUpFreeSurfer.csh: SetUpFreeSurfer-SAMPLE.csh
	sed < $< -e 's,REPLACE_WITH_PREFIX,$(prefix),' > $@
	chmod a+x $@

.PHONY: release

# nothing here goto docs directory
release:
	$(MAKE) install
	$(MAKE) release -C docs

uninstall-release:
	$(MAKE) uninstall 
	$(MAKE) uninstall-release -C docs
