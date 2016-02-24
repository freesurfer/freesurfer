#
# Makefile.am for docs
#

# Note: the pdf files should be committed to CVS, so as not to require a
# latex installation on all build systems.

# NJS: notice the .pdf files are not removed during clean or distclean.
# this is done so as not to play havoc with cvs updates since the pdf
# files are committed.  so this means that if the tex files change, 
# someone must manually delete the pdfs then run make.   

foodir=$(prefix)/fsfast/docs
# cannot specify *, since Makefile etc will be copied
TEXFILES=\
	bert-functional.howto.tex\
	bfile.tex\
	condfit.tex\
	fmri-analysis-theory.tex\
	func2roi-sess.tex\
	inorm.tex\
	INSTALL.tex\
	intergroup-sess.tex\
	isxavg-fe.tex\
	isxavg-re.tex\
	mc-afni.tex\
	mkcontrast.tex\
	mkmosaic16.tex\
	mkmosaic.tex\
	mtp-repack.tex\
	optseq.tex\
	overview.tex\
	paint.tex\
	paradigm.tex\
	selfreq.tex\
	selxavg.tex\
	seqefficiency.tex\
	stxgrinder.tex\
	univar-analysis.tex\
	vol2roi.tex\
	yakview.tex

foo_datafiles = $(TEXFILES:%.tex=%.pdf)

clean-local:
	rm -f *.aux *.log *.ps

distclean-local:
	rm -f Makefile

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
