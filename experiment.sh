#!/bin/bash

# List of test names
test_names=("brute_force" "fermat_little_theorem")

# Test numbers list
test_numbers="10 23 37"

# Output file
output_file="test_results.txt"

# Loop over each test name
for test_name in "${test_names[@]}"; do
    echo "Running tests for $test_name"

    # Save the test name to the file
    echo "Test: $test_name" >> $output_file
    
    # Run the test and append the output to the file
    cargo run --bin $test_name $test_numbers >> $output_file

    echo "---------------------------------" >> $output_file
done

echo "Results saved to $output_file"
