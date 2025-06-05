#!/bin/bash

# Usage: ./gen_flamegraph.sh <python_program> [python_args...] [--output <output_name>] [--freq <sample_frequency>]
# Example: ./gen_flamegraph.sh my_script.py arg1 arg2 --output my_profile --freq 99

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
shift

# Default values
OUTPUT_NAME="flamegraph"
SAMPLE_FREQ="99"
PYTHON_ARGS=()

# Parse remaining arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output)
            OUTPUT_NAME="$2"
            shift 2
            ;;
        --freq)
            SAMPLE_FREQ="$2"
            shift 2
            ;;
        *)
            PYTHON_ARGS+=("$1")
            shift
            ;;
    esac
done

# Validate input
if [ -z "$PYTHON_PROGRAM" ]; then
    echo "Usage: $0 <python_program> [python_args...] [--output <output_name>] [--freq <sample_frequency>]"
    echo "  python_program: Path to the Python script to profile"
    echo "  python_args: Arguments to pass to the Python script"
    echo "  --output: Name for output files (default: flamegraph)"
    echo "  --freq: Sampling frequency in Hz (default: 99)"
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
perf record -g -F "$SAMPLE_FREQ" --call-graph dwarf python3.12 "$PYTHON_PROGRAM" "${PYTHON_ARGS[@]}"

# Generate flame graph
echo "Generating flame graph..."
cd "$FLAMEGRAPH_DIR"
perf script -i ../perf.data | ./stackcollapse-perf.pl > "${OUTPUT_NAME}.perf-folded"

# Generate interactive SVG with explicit options
./flamegraph.pl "${OUTPUT_NAME}.perf-folded" > "../${OUTPUT_NAME}.svg"

cd -
