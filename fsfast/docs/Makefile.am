#
# Makefile.am for docs
#

# Note: the pdf files should be committed to CVS, so as not to require a
# latex installation on all build systems.

foodir=$(prefix)/fsfast/docs
# cannot specify *, since Makefile etc will be copied
TEXFILES=$(wildcard *.tex)
foo_DATA = $(TEXFILES:%.tex=%.pdf)

.ps.pdf:
	@ps2pdf $< $@
	@rm -f $<

.dvi.ps:
	@dvips $< -o $@ >& /dev/null
	@rm -f $<

.tex.dvi:
	@latex $< >& /dev/null

clean:
	rm -f *.aux *.log *.ps

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
