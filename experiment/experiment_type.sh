#!/bin/bash

# Define the primality test binary name
PRIM_TEST_NAME=("brute_force" "brute_force_u64")

# Define the file of random number for type tests
RANDOM_NUM="input/type_random_input.txt"

# Output file
output_file="output/type.txt"

# Number of iterations
NUM_ITERATIONS=5

# Run probabilistic tests
for upper in "${RANDOM_NUM[@]}"; do
    for iteration in $(seq 1 $NUM_ITERATIONS); do
        for test in "${PRIM_TEST_NAME[@]}"; do
            echo "Running test: $test, with upper bound of $upper in iteration $iteration."
            # echo "---------------------------------" >> $output_file
            echo "Iteration: $iteration" >> $output_file
            echo "Upper Bound: $upper" >> $output_file
            echo "Test: $test" >> $output_file
            cargo run --bin $test $upper >> $output_file
            echo "---------------------------------" >> $output_file
        done
        echo "---------------------------------" >> $output_file
    done
    echo "---------------------------------" >> $output_file
done
