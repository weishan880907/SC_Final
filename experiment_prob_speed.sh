#!/bin/bash

# Define the primality test binary name
PRIM_TEST_NAME=("fermat_little_theorem")

# Define the list of files for deterministic tests
UPPER_BOUND=24

# Output file
output_file="prob_speed.txt"

# Run probabilistic tests
for test in "${PRIM_TEST_NAME[@]}"; do
    echo "Running probabilistic test: $test"
    echo "Test: $test" >> $output_file
    cargo run --bin $test $UPPER_BOUND >> $output_file
    echo "---------------------------------" >> $output_file
done
