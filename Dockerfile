# Docker file for creating a container with freesurfer 7.2.0 release (07/2021). Note that installs which do
# not use the linux packaging tools (to install via an .rpm or .deb file) are not officially supported.
#
# Users are expected to fill out the registration form at https://surfer.nmr.mgh.harvard.edu/registration.html
# in order to download a license.txt file received thru the email address provided in the form.  Subsequently
# you can copy the license.txt file to live under $FREESURFER_HOME as .license (see commented out COPY command below).
#
# Alternately, you can place license.txt anywhere inside the container, e.g., under $HOME, and set
# the FS_LICENSE environment variable to point to it (see commented out COPY and ENV commands below).
# Then your license file will always be found even if you install/remove other freesurfer distributions.

FROM centos:7

# shell settings
WORKDIR /root

# install utils
RUN yum -y update
RUN yum -y install bc libgomp perl tar tcsh wget vim-common
RUN yum -y install mesa-libGL libXext libSM libXrender libXmu
RUN yum clean all

# install fs
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.2.0/freesurfer-linux-centos7_x86_64-7.2.0.tar.gz -O fs.tar.gz && \
    tar --no-same-owner -xzvf fs.tar.gz && \
    mv freesurfer /usr/local && \
    rm fs.tar.gz

# setup fs env
ENV OS Linux
ENV PATH /usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV FREESURFER_HOME /usr/local/freesurfer
# COPY license.txt $FREESURFER_HOME/.license
## alternately
# COPY license.txt $HOME/license.txt
# ENV FS_LICENSE=$HOME/license.txt
ENV FREESURFER /usr/local/freesurfer
ENV SUBJECTS_DIR /usr/local/freesurfer/subjects
ENV LOCAL_DIR /usr/local/freesurfer/local
ENV FSFAST_HOME /usr/local/freesurfer/fsfast
ENV FMRI_ANALYSIS_DIR /usr/local/freesurfer/fsfast
ENV FUNCTIONALS_DIR /usr/local/freesurfer/sessions

# set default fs options
ENV FS_OVERRIDE 0
ENV FIX_VERTEX_AREA ""
ENV FSF_OUTPUT_FORMAT nii.gz

# mni env requirements
ENV MINC_BIN_DIR /usr/local/freesurfer/mni/bin
ENV MINC_LIB_DIR /usr/local/freesurfer/mni/lib
ENV MNI_DIR /usr/local/freesurfer/mni
ENV MNI_DATAPATH /usr/local/freesurfer/mni/data
ENV MNI_PERL5LIB /usr/local/freesurfer/mni/share/perl5
ENV PERL5LIB /usr/local/freesurfer/mni/share/perl5
