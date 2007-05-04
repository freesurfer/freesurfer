##
## Makefile.am
##

# Fortran flags
AM_FFLAGS=-ffixed-line-length-132 -fno-second-underscore

# declare LDFLAGS and LIBS as empty, to override unnecessary default settings
LDFLAGS=
LIBS=

foodir=$(prefix)/average
foo_DATA= \
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
	mni_average_305.4dfp.hdr \
	mni_average_305.4dfp.ifh \
	mni_average_305.4dfp.img \
	mni_average_305.4dfp.img.rec

dist_bin_SCRIPTS = talairach_avi mpr2mni305 brec brec.awk t4img_4dfp

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

4dfptoanalyze_SOURCES=4dfptoanalyze.c $(TRX_SRCS) cflip.c Inithdr.c

flip_4dfp_SOURCES=flip_4dfp.c $(TRX_SRCS) cflip.c

RMS_SRCS=cgauss3d.c fftsol.f imgpad.f polfit.f matopr.f eigen.f param6opr.f

ACT_SRCS=conc.c

gauss_4dfp_SOURCES=gauss_4dfp.c $(TRX_SRCS) $(RMS_SRCS) $(ACT_SRCS) cgauss3dd.c

zero_lt_4dfp_SOURCES=zero_lt_4dfp.c $(TRX_SRCS)

LIN_SRCS=t4_sub.f t4scale.c

IMG_SRCS=fimgreg.f imgvalm.f imgvalx.f ffind.f

imgreg_4dfp_SOURCES=imgreg_4dfp.c $(TRX_SRCS) $(RMS_SRCS) \
	$(LIN_SRCS) $(IMG_SRCS) cflip.c

T4IMGS_SRCS=ft4imgo.f to_711-2B.f ft4imgn.f set_rnan.f

T4A_SRCS=ft4ixyz.f

t4imgs_4dfp_SOURCES=t4imgs_4dfp.c $(TRX_SRCS) $(RMS_SRCS) $(LIN_SRCS) \
	$(IMG_SRCS) $(T4IMGS_SRCS) $(T4A_SRCS) cflip.c spline3dvgh.f

compute_vox2vox_SOURCES=compute_vox2vox.c $(TRX_SRCS) $(RMS_SRCS) $(LIN_SRCS) \
	$(IMG_SRCS) cflip.c set_rnan.f ft4imgo.f 

ifh2hdr_SOURCES=ifh2hdr.c Inithdr.c $(TRX_SRCS)

TESTS=

EXTRA_DIST=$(dist_bin_SCRIPTS) \
	$(wildcard *.h) $(wildcard *.c) $(wildcard *.f) \
	$(foo_DATA)

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
	list=( $(wildcard *.c) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec diff $$p {} \;; \
	done;

diff_hfiles_from_src:
	list=( $(wildcard *.h) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec diff $$p {} \;; \
	done;

diff_ffiles_from_src:
	list=( $(wildcard *.f) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec diff $$p {} \;; \
	done;

copy_cfiles_from_src:
	list=( $(wildcard *.c) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec cp -iv {} . \;; \
	done;

copy_hfiles_from_src:
	list=( $(wildcard *.h) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec cp -iv {} . \;; \
	done;

copy_ffiles_from_src:
	list=( $(wildcard *.f) ); for p in $${list[@]}; do \
		find src -follow -type f -name $$p -exec cp -iv {} . \;; \
	done;

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
