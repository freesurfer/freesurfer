#
# Makefile.am
#
SUBDIRS=. average lib docs

foodir=$(prefix)
foo_DATA=surface_labels.txt tkmeditColorsCMA tkmeditParcColorsCMA Simple_surface_labels2002.txt

.PHONY: release

# nothing here goto docs directory
release:
	$(MAKE) release -C docs
