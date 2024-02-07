#!/bin/bash

# Define the primality test name
PRIM_TEST_NAME=("aks")

# Define the list of files for deterministic tests
DET_TEST_FILES="fermat_little_theorem.txt"

# Output file
output_file="det_speed.txt"

# Run deterministic tests
for test in "${PRIM_TEST_NAME[@]}"; do
    echo "Running deterministic test: $test"
    echo "Test: $test" >> $output_file
    cargo run --bin $test $DET_TEST_FILES >> $output_file
    echo "---------------------------------" >> $output_file
done
