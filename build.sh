#!/bin/bash

# Exit immediately if a command fails
set -e

echo "=== Starting configuration and compilation ==="

# Create the build directory if it doesn't exist and enter it
mkdir -p build
cd build

# Run CMake pointing to the parent folder (where CMakeLists.txt is located)
cmake ..

# Compile using all available cores to speed up the process
echo "=== Compilation in progress ==="
make -j$(nproc)

echo "=== Compilation completed successfully! ==="