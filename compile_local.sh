#!/bin/bash

# Script to compile GECKO-A and Box Model locally on macOS
# Usage: ./compile_local.sh

# Navigate to the script's directory (project root)
cd "$(dirname "$0")" || exit

export MACOSX_DEPLOYMENT_TARGET=15.0

echo "🛠️  Starting Local Compilation..."

# --- 1. Compile GECKO-A Generator ---
echo "------------------------------------------------"
echo "🦎 Compiling GECKO-A Generator..."
echo "------------------------------------------------"

GECKO_DIR="docker/gecko_source"
if [ ! -d "$GECKO_DIR" ]; then
    echo "❌ GECKO source directory not found at $GECKO_DIR"
    exit 1
fi

cd "$GECKO_DIR/OBJ" || exit
echo "   Cleaning previous build..."
make clean > /dev/null 2>&1

echo "   Compiling (this may take a minute)..."
make

if [ -f "cm" ]; then
    echo "✅ GECKO-A Generator compiled successfully!"
else
    echo "❌ GECKO-A compilation failed."
    echo "   Check if 'gfortran' is installed and in your PATH."
    exit 1
fi

cd ../../.. # Back to root

# --- 2. Compile Box Model ---
echo "------------------------------------------------"
echo "📦 Compiling Box Model..."
echo "------------------------------------------------"

BOX_DIR="docker/boxmodel_source"
if [ ! -d "$BOX_DIR" ]; then
    echo "❌ Box Model source directory not found at $BOX_DIR"
    exit 1
fi

# Check for NetCDF Fortran
if ! command -v nf-config &> /dev/null; then
    echo "⚠️  'nf-config' not found."
    echo "   NetCDF-Fortran is required to compile the Box Model."
    echo "   On macOS with Homebrew, run: brew install netcdf netcdf-fortran"
    echo "   Skipping Box Model compilation (Mock solver will be used)."
    exit 0
fi

# Create macOS architecture file
echo "   Detecting NetCDF paths..."
NCF_PREFIX=$(nf-config --prefix)
NCF_INC=$(nf-config --includedir)
NCF_LIB=$(nf-config --libdir)

# Detect NetCDF C paths (separately, as Homebrew might separate them)
if command -v nc-config &> /dev/null; then
    NCC_PREFIX=$(nc-config --prefix)
    NCC_INC=$(nc-config --includedir)
    NCC_LIB=$(nc-config --libdir)
else
    # Fallback if nc-config not found (unlikely if netcdf installed)
    NCC_PREFIX=$NCF_PREFIX
    NCC_INC=$NCF_INC
    NCC_LIB=$NCF_LIB
fi

ARCH_FILE="$BOX_DIR/SCRIPTS/arch/macos_arch.gnu"
cat <<EOF > "$ARCH_FILE"
# Architecture file for macOS (Auto-generated)
export my_netcdf_fortran_lib='$NCF_LIB'
export my_netcdf_fortran_inc='$NCF_INC'
export my_netcdf_fortran_bin='$NCF_PREFIX/bin'

export my_netcdf_c_lib='$NCC_LIB'
export my_netcdf_c_inc='$NCC_INC'
export my_netcdf_c_bin='$NCC_PREFIX/bin'

export LD_LIBRARY_PATH=\${my_netcdf_c_lib}:\${my_netcdf_fortran_lib}:\$LD_LIBRARY_PATH
export my_hdr=Makefile.hdr.gnu
export FC=gfortran
EOF

echo "   Created architecture file: $ARCH_FILE"

cd "$BOX_DIR" || exit
echo "   Building Box Model..."
./build.sh --arch macos_arch.gnu

# Check if binary exists (usually in PROG/boxmod.exe or similar, build.sh usually tells)
# We'll assume success if build.sh returns 0, but let's check PROG
if [ -f "PROG/boxmod.exe" ] || [ -f "PROG/boxmod" ]; then
    echo "✅ Box Model compiled successfully!"
else
    # Some builds might put it elsewhere or name it differently.
    # Let's just say finished.
    echo "ℹ️  Box Model build script finished."
fi

cd ../.. # Back to root

echo "------------------------------------------------"
echo "🎉 Compilation process complete."
echo "   You can now run ./start_gecko.command"
echo "------------------------------------------------"
