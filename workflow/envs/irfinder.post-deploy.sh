#!/usr/bin/env bash
set -euo pipefail

VERSION="2.0.1"

cd "$CONDA_PREFIX"

wget -O "IRFinder-${VERSION}.tar.gz" "https://github.com/RitchieLabIGH/IRFinder/archive/refs/tags/v${VERSION}.tar.gz"
tar -xzf "IRFinder-${VERSION}.tar.gz"

cd "IRFinder-${VERSION}/src/irfinder"
"${CXX:-g++}" -O3 -march=x86-64-v2 -include cstdint src/IRFinder2.cpp src/*/*.cpp \
    -I"$CONDA_PREFIX/include" -L"$CONDA_PREFIX/lib" -Wl,-rpath,"$CONDA_PREFIX/lib" \
    -lboost_iostreams -o "$CONDA_PREFIX/IRFinder-${VERSION}/bin/util/irfinder"
"${OBJCOPY:-objcopy}" --remove-section .note.gnu.property "$CONDA_PREFIX/IRFinder-${VERSION}/bin/util/irfinder"
cd "$CONDA_PREFIX"

chmod -R a+x "IRFinder-${VERSION}/bin"

ln -sf "$CONDA_PREFIX/IRFinder-${VERSION}/bin/IRFinder" "$CONDA_PREFIX/bin/IRFinder"

"$CONDA_PREFIX/bin/IRFinder" -v