# Define base image. Set HOME to avoid Matplotlib warning about non-writable
# MPLCONFIGDIR on Neurite import when running as non-root user.
FROM tensorflow/tensorflow:2.14.0-gpu AS base
ENV FREESURFER_HOME=/freesurfer
ENV PYTHONUSERBASE="$FREESURFER_HOME/env"
ENV PATH="$FREESURFER_HOME:$PATH"
ENV HOME=/tmp


# Intermediate build stage. Install Python packages to user base for easy COPY.
FROM base AS copy

COPY --chmod=0775 mri_synthmorph $FREESURFER_HOME/
COPY --chmod=0664 synthmorph.*.h5 $FREESURFER_HOME/models/

RUN apt-get update && apt-get install -y --no-install-recommends git
RUN python3 -m pip install -U pip
RUN python3 -m pip install --user \
    git+https://github.com/adalca/pystrum.git@ba35d4b357f54e5ed577cbd413076a07ef810a21 \
    git+https://github.com/adalca/neurite.git@d931ac70b12914a80f7e95c486d716d510c1be2c \
    git+https://github.com/freesurfer/surfa.git@3c9cf585d3bde229ef9a8b170ae12a96707f3987 \
    git+https://github.com/voxelmorph/voxelmorph.git@ca3d47a2c254aae9a0c0e1b30c24c324c211ebc8


# Exclude Git and caches from final image to save space. Copy only once to
# avoid unnecessary container layers. Set working directory to /mnt for data
# exchange with the host without having to specify the full path.
FROM base
COPY --from=copy $FREESURFER_HOME $FREESURFER_HOME
WORKDIR /mnt
ENTRYPOINT ["mri_synthmorph"]
