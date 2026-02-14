#!/bin/bash

# Exit immediately if a command fails
set -e

# Set default values if not provided from command line
NUM_PROCS=${1:-2}
PARAM_FILE=${2:-"../src/case_2.prm"}

# ---------------------------------------------------------
# CLEAN OUTPUT FOLDER
# ---------------------------------------------------------
echo "=== Cleaning 'output' folder ==="
# Create the folder if it doesn't exist, otherwise empty its contents
mkdir -p output
rm -rf output/*
echo "Cleanup completed."
echo "---------------------------------------------------------"

# ---------------------------------------------------------
# START SIMULATION
# ---------------------------------------------------------
echo "=== Starting simulation ==="
echo "  MPI processes: $NUM_PROCS"
echo "  Parameter file: $PARAM_FILE"
echo "---------------------------------------------------------"

# Enter build folder and execute
cd build
mpirun -np $NUM_PROCS --allow-run-as-root ./navier_stokes "$PARAM_FILE"

echo "=== Simulation completed successfully! ==="