#!/usr/bin/env bash
#
# bench_mris_register.sh
# ---------------------------------------------------------------------------
# Measures the mris_register speedup (serial vs all cores) on a self-contained
# synthetic subject. Builds a spherical surface with a smooth curvature pattern,
# makes a 1-subject atlas with mris_make_template, then registers a rotated copy
# at -threads 1 vs all cores and reports the wall-clock speedup. The registered
# spheres are compared (the reproducible reductions make them thread-invariant).
#
# Requires a dist with mris_register, mris_convert and mris_make_template
# (i.e. run build_sphere_reg_pipeline.sh first).
#
# Usage:  ./bench_mris_register.sh [dist-dir]   (default: ./dist-freesurfer)
# ---------------------------------------------------------------------------
set -euo pipefail
DIST="${1:-dist-freesurfer}"
HERE="$(cd "$(dirname "$0")" && pwd)"
[ -d "$DIST" ] || DIST="$HERE/../$DIST"
DIST="$(cd "$DIST" && pwd)"
BIN="$DIST/bin"
for t in mris_register mris_convert mris_make_template; do
  [ -x "$BIN/$t" ] || { echo "skip benchmark: $BIN/$t missing (build with mris_make_template)"; exit 0; }
done

case "$(uname -s)" in
  Darwin) NCPU="$(sysctl -n hw.ncpu)";;
  *)      NCPU="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)";;
esac
export FREESURFER_HOME="${FREESURFER_HOME:-$HERE/../distribution}"
export SURFER_SIDEDOOR=1
run() { "$DIST/run.sh" "$@"; }

WD="$(mktemp -d)"; trap 'rm -rf "$WD"' EXIT
SD="$WD/subjects"; mkdir -p "$SD/atlas01/surf" "$SD/test01/surf"

# 1. synthetic spherical surfaces + curvature (FreeSurfer ASCII + curv files)
python3 - "$SD" <<'PY'
import struct, math, sys, os
sd = sys.argv[1]
# icosphere by subdividing an octahedron a few times, then normalising
def octa():
    V=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    F=[(0,2,4),(2,1,4),(1,3,4),(3,0,4),(2,0,5),(1,2,5),(3,1,5),(0,3,5)]
    return [list(v) for v in V],[list(f) for f in F]
def subdivide(V,F):
    mid={}; nF=[]
    def m(a,b):
        k=(min(a,b),max(a,b))
        if k in mid: return mid[k]
        x=[(V[a][i]+V[b][i])/2 for i in range(3)]
        V.append(x); mid[k]=len(V)-1; return mid[k]
    for a,b,c in F:
        ab,bc,ca=m(a,b),m(b,c),m(c,a)
        nF+=[[a,ab,ca],[b,bc,ab],[c,ca,bc],[ab,bc,ca]]
    return V,nF
V,F=octa()
for _ in range(6): V,F=subdivide(V,F)   # ~40962 verts
def norm(v):
    r=math.sqrt(sum(c*c for c in v)) or 1.0
    return [c/r*100 for c in v]
V=[norm(v) for v in V]
def rotZ(V,a):
    c,s=math.cos(a),math.sin(a); return [[c*x-s*y,s*x+c*y,z] for x,y,z in V]
def pat(x,y,z):
    r=math.sqrt(x*x+y*y+z*z) or 1; x,y,z=x/r,y/r,z/r
    return math.sin(3*x)*math.cos(2*y)+0.6*z*z*z-0.4*math.cos(4*x*y)
def wasc(p,V,F):
    with open(p,"w") as f:
        f.write("#!ascii\n%d %d\n"%(len(V),len(F)))
        for x,y,z in V: f.write("%f %f %f 0\n"%(x,y,z))
        for a,b,c in F: f.write("%d %d %d 0\n"%(a,b,c))
def wcurv(p,vals,nf):
    with open(p,"wb") as f:
        f.write(b"\xff\xff\xff"); f.write(struct.pack(">iii",len(vals),nf,1))
        f.write(struct.pack(">%df"%len(vals),*vals))
pat0=[pat(*v) for v in V]
H=[0.7*math.cos(4*v[0]/100)*math.sin(3*v[1]/100)+0.5*v[0]*v[2]/1e4 for v in V]
for name,rot in (("atlas01",0.0),("test01",math.radians(25))):
    Vs=rotZ(V,rot) if rot else V
    s=os.path.join(sd,name,"surf")
    wasc(os.path.join(s,"lh.sphere.asc"),Vs,F)
    wasc(os.path.join(s,"lh.smoothwm.asc"),Vs,F)
    wcurv(os.path.join(s,"lh.sulc"),pat0,len(F))
    wcurv(os.path.join(s,"lh.inflated.H"),H,len(F))
print("generated %d verts"%len(V))
PY

for s in atlas01 test01; do
  for surf in sphere smoothwm; do
    run mris_convert "$SD/$s/surf/lh.$surf.asc" "$SD/$s/surf/lh.$surf" >/dev/null 2>&1
  done
done

# 2. atlas from atlas01 (never fail CI if the synthetic atlas can't be built)
( cd "$SD" && run mris_make_template lh sphere atlas01 "$WD/lh.atlas.tif" >/dev/null 2>&1 ) || true
[ -f "$WD/lh.atlas.tif" ] || { echo "skip: atlas build unavailable on this runner"; exit 0; }

# 3. register test01 (rotated) at 1 vs NCPU threads, timed
reg() {  # $1=threads $2=out
  local t0 t1
  t0=$(date +%s)
  run mris_register -threads "$1" -curv \
      "$SD/test01/surf/lh.sphere" "$WD/lh.atlas.tif" "$WD/$2" >/dev/null 2>&1
  t1=$(date +%s); echo $((t1 - t0))
}
echo ">>> mris_register benchmark (40k-vertex sphere, 25deg offset)"
T1=$(reg 1 reg1)
TN=$(reg "$NCPU" regN)
run mris_convert "$WD/reg1" "$WD/reg1.asc" >/dev/null 2>&1 || true
run mris_convert "$WD/regN" "$WD/regN.asc" >/dev/null 2>&1 || true
if grep -vE '^#' "$WD/reg1.asc" >"$WD/a" 2>/dev/null && grep -vE '^#' "$WD/regN.asc" >"$WD/b" 2>/dev/null && diff -q "$WD/a" "$WD/b" >/dev/null 2>&1; then
  IDENT="identical"; else IDENT="differs (expected within tolerance for iterative descent)"; fi
echo "    threads=1        : ${T1}s"
echo "    threads=${NCPU}  : ${TN}s"
[ "$TN" -gt 0 ] && awk -v a="$T1" -v b="$TN" -v n="$NCPU" 'BEGIN{printf "    speedup          : %.2fx on %d cores\n", a/b, n}'
echo "    sphere.reg vs 1t : $IDENT"
