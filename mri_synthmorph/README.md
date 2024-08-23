# SynthMorph

This guide explains how to build SynthMorph container images.
It assumes you execute commands in the `mri_synthmorph` directory.
For general information about SynthMorph, visit [synthmorph.io](https://synthmorph.io).


## Managing weight files with git-annex

Weight files are large and therefore managed with `git-annex`.
Instructions with examples are available elsewhere:

* https://surfer.nmr.mgh.harvard.edu/fswiki/GitAnnex
* https://git-annex.branchable.com/walkthrough


## Building SynthMorph images with Docker

FreeSurfer automatically ships the most recent `mri_synthmorph` and weight files.
Building a standalone container image requires fetching and unlocking the model files with `git-annex`, replacing the symbolic links with the actual files.

```sh
git fetch datasrc
git annex get .
git annex unlock synthmorph.*.h5
```

Build a new image with the appropriate version tag:

```sh
tag=X
docker build -t freesurfer/synthmorph:$tag .
```


## Testing the local Docker image

Update the version reference in the wrapper script and run it to test the local image with Docker.

```sh
sed -i "s/^\(version = \).*/\1$tag/" synthmorph
./synthmorph -h
```


## Testing with Apptainer

Testing the image with Apptainer (Singularity) before making it public requires conversion.
If your home directory has a low quota, set up a cache elsewhere:

```sh
d=$(mktemp -d)
export APPTAINER_CACHEDIR="$d"
export APPTAINER_TMPDIR="$d"
```

On the machine running Docker, convert the image with:

```sh
apptainer build -f synthmorph_$tag.sif docker-daemon://freesurfer/synthmorph:$tag
```

If you want to test the image on another machine, save it first.
After transfer to the target machine, build a SIF file as a non-root user using the fakeroot feature.
This relies on namespace mappings set up in /etc/subuid and /etc/subgid (likely by Help).

```sh
docker save synthmorph:$tag | gzip >synthmorph_$tag.tar.gz
apptainer build -f synthmorph_$tag.sif docker-archive://synthmorph_$tag.tar.gz
```

Finally, run the image.

```sh
apptainer run --nv -e -B /autofs synthmorph_$tag.sif
```


## Pushing to the Docker Hub

Push the new image to the Docker Hub to make it public.
Update the default "latest" tag, so that `docker pull freesurfer/synthmorph` without a tag will fetch the most recent image.

```sh
docker push freesurfer/synthmorph:$tag
docker tag freesurfer/synthmorph:$tag freesurfer/synthmorph:latest
docker push freesurfer/synthmorph:latest
```


## Exporting Python requirements

Export build artifacts for users who wish to create a custom Python environment.

```sh
docker build --target export --output env .
```


## Final steps

Lock the annexed weight files again to prevent modification.

```sh
git annex lock synthmorph.*.h5
```
