#!/bin/bash

# Define the primality test binary name
PRIM_TEST_NAME=("fermat_little_theorem")

# Define the number of upper bound for probabilistic tests
UPPER_BOUND=24

# Run probabilistic tests
for test in "${PRIM_TEST_NAME[@]}"; do
    echo "Running probabilistic test: $test"
    cargo run --bin $test $UPPER_BOUND 
done

# Define the primality test name
DET_TEST_NAME="brute_force"

# Run deterministic tests
for test in "${PRIM_TEST_NAME[@]}"; do
    echo "Running deterministic test: $test"
    cargo run --bin $DET_TEST_NAME $test.txt 
done