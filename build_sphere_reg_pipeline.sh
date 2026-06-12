#!/usr/bin/env bash
#
# build_sphere_reg_pipeline.sh
# ---------------------------------------------------------------------------
# One-shot builder for the parallelized FreeSurfer surface-registration tools
# (mris_register, mris_make_template, mris_convert, mris_curvature).
#
# It installs the required dependencies, configures, and compiles WITHOUT any
# manual editing. Run it from a checkout of this repository/branch.
#
#   Linux  (Ubuntu/Debian) : uses apt
#   macOS  (Intel & Apple Silicon) : uses Homebrew
#
# Usage:
#     ./build_sphere_reg_pipeline.sh [--jobs N] [--builddir DIR] [--no-deps]
#
# Result:
#     ./dist-freesurfer/bin/{mris_register,mris_convert,mris_curvature,...}
#     ./dist-freesurfer/run.sh   (sets up library paths and runs a tool)
#
# Example after building:
#     ./dist-freesurfer/run.sh mris_register -threads 8 \
#         lh.sphere lh.atlas.tif lh.sphere.reg
# ---------------------------------------------------------------------------
set -euo pipefail

JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
BUILD_DIR="_build_sphere_reg"
INSTALL_DEPS=1
SRC_DIR="$(cd "$(dirname "$0")" && pwd)"

while [ $# -gt 0 ]; do
  case "$1" in
    --jobs) JOBS="$2"; shift 2;;
    --builddir) BUILD_DIR="$2"; shift 2;;
    --no-deps) INSTALL_DEPS=0; shift;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) echo "unknown option: $1"; exit 1;;
  esac
done

OS="$(uname -s)"; ARCH="$(uname -m)"
echo ">>> Building FreeSurfer sphere-reg tools on ${OS}/${ARCH} with ${JOBS} jobs"

# ---------------------------------------------------------------------------
# 1. Dependencies
# ---------------------------------------------------------------------------
ITK_DIR=""
if [ "$INSTALL_DEPS" = "1" ]; then
  if [ "$OS" = "Linux" ]; then
    echo ">>> Installing Linux dependencies (apt)"
    export DEBIAN_FRONTEND=noninteractive
    SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO="sudo"
    $SUDO apt-get update -q
    $SUDO apt-get install -y -q \
        build-essential cmake git xxd \
        libinsighttoolkit5-dev libgsl-dev zlib1g-dev libtiff-dev \
        libx11-dev libxext-dev libxmu-dev libxt-dev libexpat1-dev \
        python3
  elif [ "$OS" = "Darwin" ]; then
    echo ">>> Installing macOS dependencies (brew)"
    command -v brew >/dev/null || { echo "Homebrew required: https://brew.sh"; exit 1; }
    brew update >/dev/null || true
    # gcc provides gfortran, which FreeSurfer's CMake enables as a language
    brew install cmake itk gsl libomp jpeg-turbo libtiff expat gcc || true
    ITK_DIR="$(echo "$(brew --prefix)"/lib/cmake/ITK-* | tr ' ' '\n' | head -1)"
  fi
fi

# Locate ITK if not set above
if [ -z "$ITK_DIR" ]; then
  for d in /usr/lib/*/cmake/ITK-* /usr/local/lib/cmake/ITK-* "$(brew --prefix 2>/dev/null)"/lib/cmake/ITK-*; do
    [ -d "$d" ] && ITK_DIR="$d" && break
  done
fi
echo ">>> Using ITK_DIR=${ITK_DIR:-<auto>}"

# ---------------------------------------------------------------------------
# 2. Configure (no source edits; helper resolves X11 + relaxes legacy -Werror)
# ---------------------------------------------------------------------------
EXTRA=""
if [ "$OS" = "Darwin" ]; then
  # AppleClang OpenMP via Homebrew libomp
  OMP_PREFIX="$(brew --prefix libomp 2>/dev/null || true)"
  [ -n "$OMP_PREFIX" ] && EXTRA="-DOpenMP_ROOT=$OMP_PREFIX"
  # FreeSurfer's CMake enables the Fortran language; point it at Homebrew gfortran
  GFORTRAN="$(ls "$(brew --prefix gcc 2>/dev/null)"/bin/gfortran* 2>/dev/null | head -1 || true)"
  [ -z "$GFORTRAN" ] && GFORTRAN="$(command -v gfortran || true)"
  [ -n "$GFORTRAN" ] && EXTRA="$EXTRA -DCMAKE_Fortran_COMPILER=$GFORTRAN"
fi

echo ">>> Configuring (cmake)"
cmake -S "$SRC_DIR" -B "$BUILD_DIR" \
  -DCMAKE_PROJECT_freesurfer_INCLUDE="$SRC_DIR/cmake/ci_build_helper.cmake" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DFS_PACKAGES_DIR=/usr \
  ${ITK_DIR:+-DITK_DIR="$ITK_DIR"} \
  -DBUILD_GUIS=OFF -DMINIMAL=ON \
  -DINSTALL_PYTHON_DEPENDENCIES=OFF -DDISTRIBUTE_FSPYTHON=OFF -DFSPYTHON_BUILD_REQ=OFF \
  -DBUILD_FORTRAN=OFF -DINTEGRATE_SAMSEG=OFF -DBUILD_ATTIC=OFF -DQATOOLS_MODULE=OFF \
  -DPYTHON_EXECUTABLE="$(command -v python3)" \
  $EXTRA

# ---------------------------------------------------------------------------
# 3. Build the tools the sphere->sphere.reg->std-mesh pipeline needs
# ---------------------------------------------------------------------------
echo ">>> Building tools"
cmake --build "$BUILD_DIR" -j"$JOBS" \
  --target mris_register mris_convert mris_curvature
# mris_make_template is excluded from MINIMAL; build it directly if present
if [ -d "$BUILD_DIR/mris_make_template" ]; then
  cmake --build "$BUILD_DIR" -j"$JOBS" --target mris_make_template || true
fi

# ---------------------------------------------------------------------------
# 4. Package
# ---------------------------------------------------------------------------
DIST="dist-freesurfer"
rm -rf "$DIST"; mkdir -p "$DIST/bin" "$DIST/lib"
for t in mris_register mris_convert mris_curvature mris_make_template; do
  f="$(find "$BUILD_DIR" -name "$t" -type f -perm -u+x 2>/dev/null | head -1)"
  [ -n "$f" ] && cp "$f" "$DIST/bin/" && echo "    packaged $t"
done
# bundle non-system shared libs
for b in "$DIST"/bin/*; do
  if [ "$OS" = "Linux" ]; then
    ldd "$b" 2>/dev/null | awk '/=>/{print $3}' \
      | grep -iE "ITK|libgomp|libtiff|libexpat|libgsl|libz|libnetcdf" || true
  else
    otool -L "$b" 2>/dev/null | awk 'NR>1{print $1}' \
      | grep -iE "ITK|libomp|libtiff|libexpat|libgsl" || true
  fi
done | sort -u | while read -r so; do [ -f "$so" ] && cp -n "$so" "$DIST/lib/" 2>/dev/null || true; done

cat > "$DIST/run.sh" <<'RUN'
#!/usr/bin/env bash
HERE="$(cd "$(dirname "$0")" && pwd)"
case "$(uname -s)" in
  Linux)  export LD_LIBRARY_PATH="$HERE/lib:${LD_LIBRARY_PATH:-}";;
  Darwin) export DYLD_LIBRARY_PATH="$HERE/lib:${DYLD_LIBRARY_PATH:-}";;
esac
export SURFER_SIDEDOOR=1   # lets surface tools run without a recon license
exec "$HERE/bin/$@"
RUN
chmod +x "$DIST/run.sh"

echo ""
echo ">>> DONE. Binaries in $DIST/bin :"
ls -1 "$DIST/bin"
echo ">>> Run with, e.g.:"
echo "    $DIST/run.sh mris_register -threads $JOBS lh.sphere lh.atlas.tif lh.sphere.reg"
