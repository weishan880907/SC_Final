#!/bin/bash

# Test names list
test_names=("brute_force" "fermat_little_theorem" "sieve_of_eratosthenes" "aks")

# Test numbers list
test_numbers=(3	17	577	9721	8145299	6759414029		940096417244711
                 38474065577182139287		340681920367515012250940028947
                 14477964987562292665377353639322377415787332354497        
                   8086121675349298378918633775461057603945233990664964430350610220444038837310829557805357343925593697
                        100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000118000000080101811009000118101080000000811000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)

# Output file
output_file="test_results.txt"


for num in "${test_numbers[@]}"; do
    for test_name in "${test_names[@]}"; do
        echo "Running tests for $test_name"
        # Save the test name to the file
        echo "Test: $test_name" >> $output_file
        # Run the test and append the output to the file
        cargo run --bin $test_name $num >> $output_file
        echo "---------------------------------" >> $output_file
    done
    echo "Number Test End: $num" >> $output_file
    echo "---------------------------------" >> $output_file
done

echo "Results saved to $output_file"
