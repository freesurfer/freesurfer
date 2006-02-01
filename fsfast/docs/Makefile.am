#
# Makefile.am for docs
#

foodir=$(prefix)/fsfast/docs
# cannot specify *, since Makefile etc will be copied
TEXFILES=$(wildcard *.tex)
foo_DATA = $(TEXFILES:%.tex=%.ps)

# latex FILE.tex
# dvips FILE.dvi -o FILE.ps
# rm FILE.dvi

.dvi.ps:
	@dvips $< -o $@ >& /dev/null
	@rm -f $<

.tex.dvi:
	@latex $< >& /dev/null

clean:
	rm -f *.aux *.log *.ps

#release:
#	$(MAKE) install

#uninstall-release:
#	$(MAKE) uninstall 

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
