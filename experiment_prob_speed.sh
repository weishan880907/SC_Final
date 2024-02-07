#!/bin/bash

# Define the primality test binary name
PRIM_TEST_NAME=("fermat_little_theorem" "miller_rabin" "strong_lucas" "baillie_psw")

# Define the number of upper bound for probabilistic tests
UPPER_BOUNDS=(10 100 1000 10000 100000 1000000)

# Output file
output_file="prob_speed.txt"

# Number of iterations
NUM_ITERATIONS=5

# Run probabilistic tests
for upper in "${UPPER_BOUNDS[@]}"; do
    for iteration in $(seq 1 $NUM_ITERATIONS); do
        for test in "${PRIM_TEST_NAME[@]}"; do
            echo "Running test: $test, with upper bound of $upper in iteration $iteration."
            # echo "---------------------------------" >> $output_file
            echo "Iteration: $iteration" >> $output_file
            echo "Upper Bound: $upper" >> $output_file
            echo "Test: $test" >> $output_file
            cargo run --bin $test $upper >> $output_file
            # echo "---------------------------------" >> $output_file
        done
        echo "---------------------------------" >> $output_file
    done
    echo "---------------------------------" >> $output_file
done
