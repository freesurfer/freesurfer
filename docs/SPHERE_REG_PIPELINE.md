# Parallelized sphere-registration → standardized-mesh pipeline

This branch parallelizes the slow stages of going from a subject's spherical
surface to a standardized icosahedral mesh:

```
?h.sphere ──mris_register (FreeSurfer)──► ?h.sphere.reg ──MapIcosahedron -ld 125 (AFNI/SUMA)──► std.125 mesh
```

* **FreeSurfer** `mris_register` — the per-iteration correlation/curvature
  gradient terms now run with OpenMP, and the matching SSE terms use the
  deterministic (reproducible) parallel reduction. Output is unchanged.
* **AFNI/SUMA** `MapIcosahedron` — the per-node mapping core (`SUMA_MapSurface`)
  now runs with OpenMP and honors a new `-threads` option. Output is unchanged
  (bitwise identical to the serial result).

There are two ways to use it: **download prebuilt binaries**, or **build from
source with one command**. No FreeSurfer/AFNI build experience is required.

---

## 1. Download prebuilt binaries (no build)

CI builds **macOS Intel (x86_64)** and **macOS Apple Silicon (arm64)** binaries
for both tools on every push, and attaches them to a GitHub Release on tags.

1. Open the repo's **Releases** page (or the **Actions** tab → latest
   "sphere-reg binaries" run → **Artifacts**).
2. Download the archive for your Mac:
   * `mris_register-macos-arm64.tgz` / `mris_register-macos-intel-x86_64.tgz`
   * the matching `MapIcosahedron-macos-*` archive from the **afni** repo.
3. Unpack and run via the bundled `run.sh` (it sets the library paths for you):

   ```bash
   tar xzf mris_register-macos-arm64.tgz
   ./run.sh mris_register --help
   ```

> macOS Gatekeeper: the first time, you may need to allow the binaries with
> `xattr -dr com.apple.quarantine .` inside the unpacked folder.

---

## 2. Build from source with one command

Clone this branch and run the builder. It installs dependencies, configures,
compiles, and packages — no manual steps.

### FreeSurfer tools (mris_register, mris_convert, mris_curvature) — Linux

```bash
git clone -b claude/sphere-reg-pipeline-optimization-xtpgi7 <freesurfer-repo-url>
cd freesurfer
./build_sphere_reg_pipeline.sh          # add --jobs N to set parallelism
# -> dist-freesurfer/bin/{mris_register,mris_convert,mris_curvature}
# -> dist-freesurfer/run.sh
```

The FreeSurfer build script targets **Linux**. On macOS, use FreeSurfer's
official `mris_register` and only build AFNI's `MapIcosahedron` (below) locally.

### AFNI/SUMA MapIcosahedron

```bash
git clone -b claude/sphere-reg-pipeline-optimization-xtpgi7 <afni-repo-url>
cd afni/src
./scripts_install/build_mapicosahedron.sh
# -> dist-afni/bin/{MapIcosahedron,quickspec,ConvertSurface}
# -> dist-afni/run.sh
```

Supported hosts: `build_sphere_reg_pipeline.sh` (FreeSurfer) → Ubuntu/Debian.
`build_mapicosahedron.sh` (AFNI) → Ubuntu/Debian and macOS Intel + Apple Silicon.
Pass `--no-deps` to skip dependency installation if you already have them.

---

## 3. Run the pipeline

`mris_register` is threaded via `-threads N` (or the `OMP_NUM_THREADS`
environment variable); `recon-all -threads N` passes this through automatically.

```bash
# (a) subject sphere -> registered sphere, using a FreeSurfer atlas .tif
./dist-freesurfer/run.sh mris_register -threads 8 \
    $SUBJECTS_DIR/SUBJ/surf/lh.sphere \
    $FREESURFER_HOME/average/lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif \
    $SUBJECTS_DIR/SUBJ/surf/lh.sphere.reg

# (b) convert surfaces to SUMA-readable ASCII and make a spec
./dist-freesurfer/run.sh mris_convert lh.sphere.reg lh.sphere.reg.asc
./dist-freesurfer/run.sh mris_convert lh.white     lh.white.asc
./dist-afni/run.sh quickspec \
    -tsnad FS white      lh.white.asc      Y lh.white.asc \
    -tsnad FS sphere.reg lh.sphere.reg.asc N lh.white.asc \
    -spec SUBJ.spec

# (c) registered sphere -> standardized -ld 125 mesh (156,252 nodes)
./dist-afni/run.sh MapIcosahedron -spec SUBJ.spec -ld 125 -threads 8 \
    -morph lh.sphere.reg.asc -prefix std.125.
# -> std.125.lh.white.asc  (the standardized white surface) + node mapping
```

`-ld 125` gives `2 + 10·125² = 156,252` nodes; choose `-ld 141` (≈200k) etc. as
needed. Use all cores by passing `-threads $(nproc)` (Linux) or
`-threads $(sysctl -n hw.ncpu)` (macOS).

---

## 4. Correctness

* `mris_register` recovers a known rotation exactly and the threaded SSE term is
  bitwise thread-count-invariant (deterministic reduction).
* `MapIcosahedron` output is bitwise identical between serial and threaded runs;
  threading only changes speed (measured ~2.3–2.6× on 4 cores for the mapping
  core).

The two helper scripts (`build_sphere_reg_pipeline.sh`,
`scripts_install/build_mapicosahedron.sh`) and the CI workflows under
`.github/workflows/` contain the exact, reproducible recipe.
