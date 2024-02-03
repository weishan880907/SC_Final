// 499263 Wei-Shan Chang

use std::{env, time, fs};
use num_bigint::BigUint;
use num_traits::{One, Zero};

fn is_prime(n: &BigUint) -> bool {
    if n <= &BigUint::one() {
        return false;
    }

    let mut i = BigUint::from(2u32);
    let sqrt_n = n.sqrt();

    while &i <= &sqrt_n {
        if n % &i == BigUint::zero() {
            return false;
        }
        i += BigUint::one();
    }

    true
}

fn find_prime_numbers(input: Vec<BigUint>) -> Vec<BigUint> {
    input.into_iter().filter(|n| is_prime(n)).collect()
}

fn main() {
    // Get the input from the command line as Vec<BigUint>.
    let args: Vec<String> = env::args().collect();

    // Extract the input file name from command-line arguments.
    let input_filename = &args[1];

    // Read the input from the specified file.
    let input_content = fs::read_to_string(input_filename).expect("Failed to read file");
    
    // Parse the input content into a Vec<BigUint>.
    let input: Vec<BigUint> = input_content
        .trim()
        .split_whitespace()
        .map(|num| num.parse().unwrap())
        .collect();

    // Find the prime numbers in the input vector.
    let start_time = time::Instant::now();
    let primes: Vec<BigUint> = find_prime_numbers(input.clone());
    let elapsed_time = start_time.elapsed();

    let pseudo_prime_length = primes.len() - input.len();
    // Output the prime numbers and the time taken.
    println!("The count of pseudo primes from the previous test: {}", pseudo_prime_length);
    println!("Time: {:?}", elapsed_time);
}
