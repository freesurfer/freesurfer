#
# Makefile.am
#
SUBDIRS=. average lib docs grad_unwarp_tables fsafd subjects sessions

foodir=$(prefix)
foo_DATA=surface_labels.txt \
	tkmeditColorsCMA \
	tkmeditParcColorsCMA \
	Simple_surface_labels2002.txt \
        Simple_surface_labels2005.txt	\
	FreeSurferEnv.csh \
	SetUpFreeSurfer.csh \
	FreeSurferEnv.sh \
	SetUpFreeSurfer.sh \
	jeans_labels.txt \
	FreeSurferColorLUT.txt \
	numaris4-protocols.unpackcfg

SetUpFreeSurfer.csh: SetUpFreeSurfer-SAMPLE.csh
	sed < $< -e 's,REPLACE_WITH_PREFIX,$(prefix),' > $@
	chmod a+x $@

SetUpFreeSurfer.sh: SetUpFreeSurfer-SAMPLE.sh
	sed < $< -e 's,REPLACE_WITH_PREFIX,$(prefix),' > $@
	chmod a+x $@

.PHONY: release

# nothing here goto docs directory
release:
	$(MAKE) install
	$(MAKE) release -C docs
	$(MAKE) release -C grad_unwarp_tables

uninstall-release:
	$(MAKE) uninstall 
	$(MAKE) uninstall-release -C docs
	$(MAKE) uninstall-release -C grad_unwarp_tables
