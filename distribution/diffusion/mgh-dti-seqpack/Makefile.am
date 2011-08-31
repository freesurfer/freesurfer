#
# Makefile.am
#

foodir=$(prefix)/diffusion/mgh-dti-seqpack
foo_DATA=README \
	gradient_mgh_dti6.gdt   gradient_mgh_dti12.gdt \
	gradient_mgh_dti13.gdt  gradient_mgh_dti57.gdt \
	gradient_mgh_dti14.gdt  gradient_mgh_dti58.gdt \
	gradient_mgh_dti15.gdt  gradient_mgh_dti59.gdt \
	gradient_mgh_dti16.gdt  gradient_mgh_dti60.gdt \
	gradient_mgh_dti17.gdt  gradient_mgh_dti61.gdt \
	gradient_mgh_dti18.gdt  gradient_mgh_dti62.gdt \
	gradient_mgh_dti19.gdt  gradient_mgh_dti63.gdt \
	gradient_mgh_dti20.gdt  gradient_mgh_dti64.gdt \
	gradient_mgh_dti21.gdt  gradient_mgh_dti65.gdt \
	gradient_mgh_dti22.gdt  gradient_mgh_dti66.gdt \
	gradient_mgh_dti23.gdt  gradient_mgh_dti67.gdt \
	gradient_mgh_dti24.gdt  gradient_mgh_dti68.gdt \
	gradient_mgh_dti25.gdt  gradient_mgh_dti69.gdt \
	gradient_mgh_dti26.gdt  gradient_mgh_dti70.gdt \
	gradient_mgh_dti27.gdt  gradient_mgh_dti71.gdt \
	gradient_mgh_dti28.gdt  gradient_mgh_dti72.gdt \
	gradient_mgh_dti29.gdt  gradient_mgh_dti73.gdt \
	gradient_mgh_dti30.gdt  gradient_mgh_dti74.gdt \
	gradient_mgh_dti31.gdt  gradient_mgh_dti75.gdt \
	gradient_mgh_dti32.gdt  gradient_mgh_dti76.gdt \
	gradient_mgh_dti33.gdt  gradient_mgh_dti77.gdt \
	gradient_mgh_dti34.gdt  gradient_mgh_dti78.gdt \
	gradient_mgh_dti35.gdt  gradient_mgh_dti79.gdt \
	gradient_mgh_dti36.gdt  gradient_mgh_dti80.gdt \
	gradient_mgh_dti37.gdt  gradient_mgh_dti81.gdt \
	gradient_mgh_dti38.gdt  gradient_mgh_dti82.gdt \
	gradient_mgh_dti39.gdt  gradient_mgh_dti83.gdt \
	gradient_mgh_dti40.gdt  gradient_mgh_dti84.gdt \
	gradient_mgh_dti41.gdt  gradient_mgh_dti85.gdt \
	gradient_mgh_dti42.gdt  gradient_mgh_dti86.gdt \
	gradient_mgh_dti43.gdt  gradient_mgh_dti87.gdt \
	gradient_mgh_dti44.gdt  gradient_mgh_dti88.gdt \
	gradient_mgh_dti45.gdt  gradient_mgh_dti89.gdt \
	gradient_mgh_dti46.gdt  gradient_mgh_dti90.gdt \
	gradient_mgh_dti47.gdt  gradient_mgh_dti91.gdt \
	gradient_mgh_dti48.gdt  gradient_mgh_dti92.gdt \
	gradient_mgh_dti49.gdt  gradient_mgh_dti93.gdt \
	gradient_mgh_dti50.gdt  gradient_mgh_dti94.gdt \
	gradient_mgh_dti51.gdt  gradient_mgh_dti95.gdt \
	gradient_mgh_dti52.gdt  gradient_mgh_dti96.gdt \
	gradient_mgh_dti53.gdt  gradient_mgh_dti97.gdt \
	gradient_mgh_dti54.gdt  gradient_mgh_dti98.gdt \
	gradient_mgh_dti55.gdt  gradient_mgh_dti99.gdt \
	gradient_mgh_dti56.gdt \
	gradient_mgh_dti06.gdt gradient_mgh_dti240.gdt

EXTRA_DIST=$(foo_DATA)

install: 
	rm   -Rf $(DESTDIR)$(prefix)/diffusion/mgh-dti-seqpack/
	mkdir -p $(DESTDIR)$(prefix)/diffusion/mgh-dti-seqpack/
	cd $(top_srcdir)/distribution/diffusion/mgh-dti-seqpack/; cp $(foo_DATA) \
$(DESTDIR)$(prefix)/diffusion/mgh-dti-seqpack/

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
