#!/bin/bash

# Usage: ./gen_flamegraph.sh <python_program> [output_name] [sample_frequency]
# Example: ./gen_flamegraph.sh my_script.py my_profile 99

set -e  # Exit on any error

FLAMEGRAPH_DIR="$(dirname "$0")/FlameGraph"

# Check if FlameGraph directory exists
if [ ! -d "$FLAMEGRAPH_DIR" ]; then
    echo "Error: FlameGraph submodule not found. Please run:"
    echo "git submodule update --init --recursive"
    exit 1
fi
# Parse arguments
PYTHON_PROGRAM="$1"
OUTPUT_NAME="${2:-flamegraph}"
SAMPLE_FREQ="${3:-99}"

# Validate input
if [ -z "$PYTHON_PROGRAM" ]; then
    echo "Usage: $0 <python_program> [output_name] [sample_frequency]"
    echo "  python_program: Path to the Python script to profile"
    echo "  output_name: Name for output files (default: flamegraph)"
    echo "  sample_frequency: Sampling frequency in Hz (default: 99)"
    exit 1
fi

if [ ! -f "$PYTHON_PROGRAM" ]; then
    echo "Error: Python program '$PYTHON_PROGRAM' not found"
    exit 1
fi

# Check if running as root or with sudo privileges
if [ "$EUID" -ne 0 ]; then
    echo "Error: This script requires root privileges for perf record"
    echo "Please run with sudo: sudo $0 $@"
    exit 1
fi

echo "Profiling Python program: $PYTHON_PROGRAM"
echo "Sample frequency: ${SAMPLE_FREQ} Hz"
echo "Output name: $OUTPUT_NAME"

# Cleanup function
cleanup() {
    echo "Cleaning up temporary files..."
    rm -f perf.data perf.data.old
}
trap cleanup EXIT

# Record performance data
echo "Recording performance data..."
perf record -g -F "$SAMPLE_FREQ" --call-graph dwarf python3 "$PYTHON_PROGRAM"

# Generate flame graph
echo "Generating flame graph..."
cd "$FLAMEGRAPH_DIR"
perf script -i ../perf.data | ./stackcollapse-perf.pl > "${OUTPUT_NAME}.perf-folded"
./flamegraph.pl "${OUTPUT_NAME}.perf-folded" > "../../${OUTPUT_NAME}.svg"
cd -

echo "Flame graph generated: ${OUTPUT_NAME}.svg"
echo "Folded data saved: FlameGraph/${OUTPUT_NAME}.perf-folded"