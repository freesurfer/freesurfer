##
## Makefile.am
##

# Fortran flags
AM_FFLAGS=-ffixed-line-length-132 -fno-second-underscore

# LIB_G2C_A is defined to libg2c.a (static lib) if configure find it
AM_LDFLAGS=$(LIB_G2C_A)

LIBS=

foodir=$(prefix)/average
foo_datafiles= \
	711-2B_as_mni_average_305.4dfp.hdr \
	711-2B_as_mni_average_305_mask.4dfp.ifh \
	711-2C_as_mni_average_305.4dfp.img \
	711-2B_as_mni_average_305.4dfp.ifh \
	711-2B_as_mni_average_305_mask.4dfp.img \
	711-2C_as_mni_average_305.4dfp.img.rec \
	711-2B_as_mni_average_305.4dfp.img \
	711-2B_as_mni_average_305_mask.4dfp.img.rec \
	711-2B_as_mni_average_305.4dfp.img.rec \
	711-2C_as_mni_average_305.4dfp.hdr \
	711-2B_as_mni_average_305_mask.4dfp.hdr \
	711-2C_as_mni_average_305.4dfp.ifh \
	711-2C_as_mni_average_305.4dfp.mat \
	mni_average_305.4dfp.hdr \
	mni_average_305.4dfp.ifh \
	mni_average_305.4dfp.img \
	mni_average_305.4dfp.img.rec \
	3T18yoSchwartzReactN32_as_orig.4dfp.hdr \
	3T18yoSchwartzReactN32_as_orig.4dfp.img \
	3T18yoSchwartzReactN32_as_orig.4dfp.ifh \
	3T18yoSchwartzReactN32_as_orig.4dfp.img.rec \
	3T18yoSchwartzReactN32_as_orig.4dfp.mat \
	3T18yoSchwartzReactN32_as_orig.lst \
	RLB700_atlas_as_orig.4dfp.hdr \
	RLB700_atlas_as_orig.4dfp.ifh \
	RLB700_atlas_as_orig.4dfp.img \
	RLB700_atlas_as_orig.4dfp.img.rec \
	SVIP_Adult_Comp_N24_as_orig.4dfp.hdr \
	SVIP_Adult_Comp_N24_as_orig.4dfp.ifh \
	SVIP_Adult_Comp_N24_as_orig.4dfp.img \
	SVIP_Adult_Comp_N24_as_orig.4dfp.img.rec \
	SVIP_Adult_Control_N12_as_orig.4dfp.hdr \
	SVIP_Adult_Control_N12_as_orig.4dfp.ifh \
	SVIP_Adult_Control_N12_as_orig.4dfp.img \
	SVIP_Adult_Control_N12_as_orig.4dfp.img.rec \
	SVIP_Adult_Patient_N12_as_orig.4dfp.hdr \
	SVIP_Adult_Patient_N12_as_orig.4dfp.ifh \
	SVIP_Adult_Patient_N12_as_orig.4dfp.img \
	SVIP_Adult_Patient_N12_as_orig.4dfp.img.rec \
	SVIP_Child_Comp_N24_as_orig.4dfp.hdr \
	SVIP_Child_Comp_N24_as_orig.4dfp.ifh \
	SVIP_Child_Comp_N24_as_orig.4dfp.img \
	SVIP_Child_Comp_N24_as_orig.4dfp.img.rec \
	SVIP_Full_Comp_N48_as_orig.4dfp.hdr \
	SVIP_Full_Comp_N48_as_orig.4dfp.ifh \
	SVIP_Full_Comp_N48_as_orig.4dfp.img \
	SVIP_Full_Comp_N48_as_orig.4dfp.img.rec \
	SVIP_Young_Control_N12_as_orig.4dfp.hdr \
	SVIP_Young_Control_N12_as_orig.4dfp.ifh \
	SVIP_Young_Control_N12_as_orig.4dfp.img \
	SVIP_Young_Control_N12_as_orig.4dfp.img.rec \
	SVIP_Young_Patient_N12_as_orig.4dfp.hdr \
	SVIP_Young_Patient_N12_as_orig.4dfp.ifh \
	SVIP_Young_Patient_N12_as_orig.4dfp.img \
	SVIP_Young_Patient_N12_as_orig.4dfp.img.rec

xmldir=$(prefix)/docs/xml
xml_DATA=talairach_avi.help.xml
htmldir=$(prefix)/docs/html
html_DATA=talairach_avi.help.xml.html

dist_bin_SCRIPTS = talairach_avi mpr2mni305 brec brec.awk \
	t4img_4dfp tal_QC_AZS

bin_PROGRAMS= \
	analyzeto4dfp \
	4dfptoanalyze \
	flip_4dfp \
	gauss_4dfp \
	zero_lt_4dfp \
	imgreg_4dfp \
	t4imgs_4dfp \
	compute_vox2vox \
	ifh2hdr

TRX_SRCS=rec.c Getifh.c endianio.c

analyzeto4dfp_SOURCES=analyzeto4dfp.c $(TRX_SRCS) cflip.c
analyzeto4dfp_LDADD=$(LIBS_ICC) $(LIB_G2C_A)

4dfptoanalyze_SOURCES=4dfptoanalyze.c $(TRX_SRCS) cflip.c Inithdr.c
4dfptoanalyze_LDADD=$(LIBS_ICC) $(LIB_G2C_A)

flip_4dfp_SOURCES=flip_4dfp.c $(TRX_SRCS) cflip.c
flip_4dfp_LDADD=$(LIBS_ICC) $(LIB_G2C_A)

RMS_SRCS=cgauss3d.c fftsol.f imgpad.f polfit.f matopr.f eigen.f param6opr.f

ACT_SRCS=conc.c

gauss_4dfp_SOURCES=gauss_4dfp.c $(TRX_SRCS) $(RMS_SRCS) $(ACT_SRCS) cgauss3dd.c
gauss_4dfp_LDADD=$(LIBS_ICC) $(FCLIBS) $(LIB_G2C_A)

zero_lt_4dfp_SOURCES=zero_lt_4dfp.c $(TRX_SRCS)
zero_lt_4dfp_LDADD=$(LIBS_ICC) $(LIB_G2C_A)

LIN_SRCS=t4_sub.f t4scale.c

IMG_SRCS=fimgreg.f imgvalm.f imgvalx.f ffind.f

imgreg_4dfp_SOURCES=imgreg_4dfp.c $(TRX_SRCS) $(RMS_SRCS) \
	$(LIN_SRCS) $(IMG_SRCS) cflip.c
imgreg_4dfp_CFLAGS=-DHAVE_GFORTRAN
imgreg_4dfp_LDADD=$(LIBS_ICC) $(FCLIBS) $(LIB_G2C_A)

T4IMGS_SRCS=ft4imgo.f to_711-2B.f ft4imgn.f set_rnan.f

T4A_SRCS=ft4ixyz.f

t4imgs_4dfp_SOURCES=t4imgs_4dfp.c $(TRX_SRCS) $(RMS_SRCS) $(LIN_SRCS) \
	$(IMG_SRCS) $(T4IMGS_SRCS) $(T4A_SRCS) cflip.c spline3dvgh.f
t4imgs_4dfp_CFLAGS=-DHAVE_GFORTRAN
t4imgs_4dfp_LDADD=$(LIBS_ICC) $(FCLIBS) $(LIB_G2C_A)

compute_vox2vox_SOURCES=compute_vox2vox.c $(TRX_SRCS) $(RMS_SRCS) $(LIN_SRCS) \
	$(IMG_SRCS) cflip.c set_rnan.f ft4imgo.f 
compute_vox2vox_CFLAGS=-DHAVE_GFORTRAN
compute_vox2vox_LDADD=$(LIBS_ICC) $(FCLIBS) $(LIB_G2C_A)

ifh2hdr_SOURCES=ifh2hdr.c Inithdr.c $(TRX_SRCS)
ifh2hdr_LDADD=$(LIBS_ICC) $(LIB_G2C_A)

CFILES=\
	4dfptoanalyze.c \
	analyzeto4dfp.c \
	cflip.c \
	cgauss3d.c \
	cgauss3dd.c \
	compute_vox2vox.c \
	conc.c \
	endianio.c \
	flip_4dfp.c \
	gauss_4dfp.c \
	Getifh.c \
	ifh2hdr.c \
	imgreg_4dfp.c \
	Inithdr.c \
	rec.c \
	t4imgs_4dfp.c \
	t4scale.c \
	zero_lt_4dfp.c

HFILES=\
	ANALYZE.h \
	conc.h \
	endianio.h \
	Getifh.h \
	ifh.h \
	librms.h \
	rec.h

FFILES=\
	eigen.f \
	ffind.f \
	fftsol.f \
	fimgreg.f \
	ft4imgn.f \
	ft4imgo.f \
	ft4ixyz.f \
	imgpad.f \
	imgvalm.f \
	imgvalx.f \
	matopr.f \
	param6opr.f \
	polfit.f \
	set_rnan.f \
	spline3dvgh.f \
	t4_sub.f \
	to_711-2B.f

# trick to get test data and scripts into the build dir (during make distcheck)
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/distribution/average/mni305.cor.mgz"; \
	then \
		cp -v $(top_srcdir)/distribution/average/mni305.cor.mgz \
					$(top_builddir)/distribution/average/; fi;
	if ! test -f "$(top_builddir)/scripts/avi2talxfm"; then \
		cp -v $(top_srcdir)/scripts/avi2talxfm \
					$(top_builddir)/scripts/; fi;
	if ! test -f "$(top_builddir)/talairach_avi/nu.mgz"; then \
		cp -v $(top_srcdir)/talairach_avi/nu.mgz \
					$(top_builddir)/talairach_avi/; fi;
	if ! test -f "$(top_builddir)/talairach_avi/reference-tal.xfm"; then \
		cp -v $(top_srcdir)/talairach_avi/reference-tal.xfm \
					$(top_builddir)/talairach_avi/; fi;
	if ! test -f "$(top_builddir)/talairach_avi/711-2C_as_mni_average_305.4dfp.img"; then \
		cp -v $(top_srcdir)/talairach_avi/711-2* \
					$(top_builddir)/talairach_avi/; fi;
	if ! test -f "$(top_builddir)/talairach_avi/3T18yoSchwartzReactN32_as_orig.4dfp.img"; then \
		cp -v $(top_srcdir)/talairach_avi/*_orig* \
					$(top_builddir)/talairach_avi/; fi;
	if ! test -f "$(top_builddir)/talairach_avi/talairach_avi"; then \
		list=( $(dist_bin_SCRIPTS) ); \
		for f in $${list[@]}; do \
			cp -v $(top_srcdir)/talairach_avi/$$f $(top_builddir)/talairach_avi/; \
		done; fi;

check_PROGRAMS=foo

TESTS=test_libg2c test_talairach_avi $(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES) 

EXTRA_DIST=$(dist_bin_SCRIPTS) talairach_avi.help.xml.h \
	$(CFILES) \
	$(HFILES) \
	$(FFILES) \
	test_libg2c test_talairach_avi nu.mgz reference-tal.xfm $(xml_DATA)

listfiles:
	@echo $(EXTRA_DIST) Makefile.am

tarball:
	tar cvf talairach_avi.tar $(EXTRA_DIST) Makefile.am
	gzip talairach_avi.tar
	ls -alF talairach_avi.tar.gz

# if Avi Snyder needs to supply new source files, then a directory
# called 'src' should be created (presumably from a tarball that he
# supplies), containing his source.  then these targets can be used
# to show differences, and copy these files from src

diff_cfiles_from_src:
	list=( $(CFILES) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec diff $$p {} \;; \
	done;

diff_hfiles_from_src:
	list=( $(HFILES) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec diff $$p {} \;; \
	done;

diff_ffiles_from_src:
	list=( $(FFILES) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec diff $$p {} \;; \
	done;

copy_cfiles_from_src:
	list=( $(CFILES) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec cp -iv {} . \;; \
	done;

copy_hfiles_from_src:
	list=( $(HFILES) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec cp -iv {} . \;; \
	done;

copy_ffiles_from_src:
	list=( $(FFILES) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec cp -iv {} . \;; \
	done;

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=talairach_avi.help.xml.h talairach_avi.help.xml.html
