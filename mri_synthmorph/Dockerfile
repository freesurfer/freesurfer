# Define base image. Set HOME to avoid Matplotlib warning about non-writable
# MPLCONFIGDIR on Neurite import when running as non-root user.
FROM tensorflow/tensorflow:2.15.0.post1-gpu AS base
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
    'numpy<2.0' \
    git+https://github.com/adalca/pystrum.git@ba35d4b357f54e5ed577cbd413076a07ef810a21 \
    git+https://github.com/adalca/neurite.git@9167ad8ad4ef5cb9b22ee8c67d3c42c38eea4bda \
    git+https://github.com/freesurfer/surfa.git@ded5f1d3d90e223050ab7792ac4760c3242e43c7 \
    git+https://github.com/voxelmorph/voxelmorph.git@d3b0417ae6f662d8d90a277125062f233c71296c

RUN python3 -V >python.txt
RUN python3 -m pip freeze >requirements.txt
RUN mri_synthmorph -h >help.txt


# Export Python requirements for reference. Build artifacts will only exist in
# in the target stage `export`.
FROM scratch AS export
COPY --from=copy *.txt /


# Exclude Git and caches from final image to save space. Copy only once to
# avoid unnecessary container layers. Set working directory to /mnt for data
# exchange with the host without having to specify the full path.
FROM base
COPY --from=copy $FREESURFER_HOME $FREESURFER_HOME
WORKDIR /mnt
ENTRYPOINT ["mri_synthmorph"]
