# Define base image in a single place.
FROM tensorflow/tensorflow:2.11.0-gpu AS base
ENV FREESURFER_HOME=/freesurfer


# Install SynthMorph and Python packages.
FROM base AS copy

COPY mri_synthmorph $FREESURFER_HOME/
COPY synthmorph*.h5 $FREESURFER_HOME/models/

RUN apt-get update && apt-get install -y --no-install-recommends git
RUN python3 -m pip install -U pip
RUN PYTHONUSERBASE=/packages python3 -m pip install --user \
    git+https://github.com/adalca/pystrum.git@6422922af2607bf740a78526c6be83c388b25c1c \
    git+https://github.com/adalca/neurite.git@152368b260889f7c5b637a5d2d0cc4cc2ea20de3 \
    git+https://github.com/voxelmorph/voxelmorph.git@f55ee9e4118a87f457accbff6726ea1d92475128


# Final image without Git and caches.
FROM base
ENV PATH="$FREESURFER_HOME:$PATH"

COPY --from=copy $FREESURFER_HOME $FREESURFER_HOME
COPY --from=copy \
    /packages/lib/python3.8/site-packages/ \
    /usr/local/lib/python3.8/dist-packages

ENTRYPOINT ["mri_synthmorph"]
