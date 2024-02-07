#!/bin/bash

# Define the primality test binary name
PRIM_TEST_NAME=("aks" "brute_force" "sieve_of_eratosthenes" )

# # Define the file name for deterministic tests (speed)
FILENAME_S=("input/det_test_inputs_1.txt" "input/det_test_inputs_2.txt" "input/det_test_inputs_3.txt" "input/det_test_inputs_4.txt" "input/det_test_inputs_5.txt" "input/det_test_inputs_6.txt" "input/det_test_inputs_7.txt")

# Output file for deterministic tests (speed)
output_file_S="output/det_speed.txt"

# Define the file name for deterministic tests (accuracy)
FILENAME_A=("output/baillie_psw.txt" "output/fermat_little_theorem.txt" "output/miller_rabin.txt" "output/strong_lucas.txt")

# Output file for deterministic tests (speed)
output_file_A="output/det_acc.txt"

# Number of iterations
NUM_ITERATIONS=5

# Run deterministic tests
for file in "${FILENAME_S[@]}"; do
    for iteration in $(seq 1 $NUM_ITERATIONS); do
        for test in "${PRIM_TEST_NAME[@]}"; do
            echo "Running test: $test, with file name of $file in iteration $iteration."
            echo "Iteration: $iteration" >> $output_file_S
            echo "File name: $file" >> $output_file_S
            echo "Test: $test" >> $output_file_S
            cargo run --bin $test $file >> $output_file_S
            echo "---------------------------------" >> $output_file_S
        done
        echo "---------------------------------" >> $output_file_S
    done
    echo "---------------------------------" >> $output_file_S
done
