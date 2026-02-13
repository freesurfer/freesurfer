# Docker file for creating a container with freesurfer 8.1.0 release (07/2025). Note that installs which do
# not use the linux packaging tools (to install via an .rpm or .deb file) are not officially supported.
#
# Users are expected to fill out the registration form at https://surfer.nmr.mgh.harvard.edu/registration.html
# in order to download a license.txt file received thru the email address provided in the form.  Subsequently
# you can copy the license.txt file to live under $FREESURFER_HOME as .license (see commented out COPY command below).
#
# Alternately, you can place license.txt anywhere inside the container, e.g., under $HOME, and set
# the FS_LICENSE environment variable to point to it (see commented out COPY and ENV commands below).
# Then your license file will always be found even if you install/remove other freesurfer distributions.

FROM ubuntu:20.04

# shell settings
WORKDIR /root

# install utils
RUN apt update
RUN apt install -y bc binutils libgomp1 perl psmisc sudo tar tcsh unzip uuid-dev vim-common libjpeg62-dev wget

ENV DEBIAN_FRONTEND=noninteractive

# install fs
ARG FS_version=8.1.0
RUN wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/${FS_version}/freesurfer_ubuntu20-${FS_version}_amd64.deb -O freesurfer.deb
RUN echo "deb http://security.ubuntu.com/ubuntu focal-security main universe restricted multiverse" > /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu focal main universe restricted multiverse" >> /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu focal-updates main universe restricted multiverse" >> /etc/apt/sources.list && \
    echo "deb http://archive.ubuntu.com/ubuntu focal-backports main universe restricted multiverse" >> /etc/apt/sources.list && \
    apt-get update

RUN dpkg -i freesurfer.deb || true && \
    apt-get install -f -y && \
    rm freesurfer.deb
	
# setup fs env
ENV OS Linux
ENV PATH /usr/local/freesurfer/${FS_version}/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV FREESURFER_HOME /usr/local/freesurfer/${FS_version}
# COPY license.txt $FREESURFER_HOME/.license
## alternately
# COPY license.txt $HOME/license.txt
# ENV FS_LICENSE=$HOME/license.txt
ENV FREESURFER /usr/local/freesurfer/${FS_version}
ENV SUBJECTS_DIR /usr/local/freesurfer/${FS_version}
ENV LOCAL_DIR /usr/local/freesurfer/${FS_version}/local
ENV FSFAST_HOME /usr/local/freesurfer/${FS_version}/fsfast
ENV FMRI_ANALYSIS_DIR /usr/local/freesurfer/${FS_version}/fsfast
ENV FUNCTIONALS_DIR /usr/local/freesurfer/${FS_version}/sessions
ENV FS_ALLOW_DEEP 1

# set default fs options
ENV FS_OVERRIDE 0
ENV FIX_VERTEX_AREA ""
ENV FSF_OUTPUT_FORMAT nii.gz

# mni env requirements
ENV MINC_BIN_DIR /usr/local/freesurfer/${FS_version}/mni/bin
ENV MINC_LIB_DIR /usr/local/freesurfer/${FS_version}/mni/lib
ENV MNI_DIR /usr/local/freesurfer/${FS_version}/mni
ENV MNI_DATAPATH /usr/local/freesurfer/${FS_version}/mni/data
ENV MNI_PERL5LIB /usr/local/freesurfer/${FS_version}/mni/share/perl5
ENV PERL5LIB /usr/local/freesurfer/${FS_version}/mni/share/perl5
