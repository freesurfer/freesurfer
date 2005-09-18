#
# Makefile.am for docs
#
foodir=$(prefix)/fsfast/docs
# I cannot specify *, since Makefile etc will be copied
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
distclean:
	rm -f *.aux *.log *.ps

release:
	$(MAKE) install

uninstall-release:
	$(MAKE) uninstall 

