#!/bin/bash
# ============================================================================
#  build_jewel.sh [TARBALL] [BUILD_DIR]
#
#  Builds JEWEL 2.4.0 with the centrality patch (jewel-2.4.0-centrality.patch)
#  exactly as the production binaries jewel-2.4.0-gcc14-cent/jewel-2.4.0-{vac,simple}:
#    gfortran 14.2 from the LCG_107 x86_64-el9-gcc14-opt view, -O2, linked to
#    LHAPDF 6 of the same view (rpath set, no LD_LIBRARY_PATH needed at run time).
#
#  TARBALL    jewel-2.4.0.tar.gz (default: download from hepforge)
#  BUILD_DIR  output directory (default: ./jewel-2.4.0-gcc14-cent)
#
#  Do NOT use gcc 16 (LCG_110): with -O2 it miscompiles JEWEL and medium runs
#  hang after a few events (NaN gluon momenta, endless loop in MAKECASCADE /
#  DOINSTATESCAT). gcc 14 -O2 (this script) or gcc 16 -O0 run normally.
# ============================================================================
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
TARBALL=${1:-}
BUILD_DIR=${2:-$PWD/jewel-2.4.0-gcc14-cent}
URL=https://www.hepforge.org/archive/jewel/jewel-2.4.0.tar.gz
SHA256=8038e76d54df03d380cea58e5c0a5a6aee582b19e6d62e7753b99a34dedbe417   # tarball used in production
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc14-opt

source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_107 x86_64-el9-gcc14-opt >/dev/null 2>&1 || true
gfortran --version | head -1 | grep -q " 14\." || { echo "need gfortran 14 (LCG_107 gcc14 view)"; exit 1; }

mkdir -p "$BUILD_DIR"
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
if [ -z "$TARBALL" ]; then TARBALL=$WORK/jewel-2.4.0.tar.gz; curl -sSL -o "$TARBALL" "$URL"; fi
echo "$SHA256  $TARBALL" | sha256sum -c - || echo "WARNING: tarball checksum differs from the production tarball"

# unpack and apply the centrality patch (writes the sampled centrality into the
# impact-parameter slot of the HepMC H record; see README_generation.md)
tar xzf "$TARBALL" -C "$WORK"
patch -p1 -d "$WORK" < "$HERE/jewel-2.4.0-centrality.patch"
cp "$WORK"/jewel-2.4.0/*.f "$BUILD_DIR"/

cd "$BUILD_DIR"
for f in jewel-2.4.0 medium-simple medium-vac pythia6425mod-lhapdf6 meix; do
  gfortran -O2 -c $f.f -o $f.o
done
LINK="-L$LCG/lib -Wl,-rpath,$LCG/lib -lLHAPDF -lstdc++"
gfortran -o jewel-2.4.0-simple jewel-2.4.0.o medium-simple.o pythia6425mod-lhapdf6.o meix.o $LINK
gfortran -o jewel-2.4.0-vac    jewel-2.4.0.o medium-vac.o    pythia6425mod-lhapdf6.o meix.o $LINK
echo "built $BUILD_DIR/jewel-2.4.0-{vac,simple}"
