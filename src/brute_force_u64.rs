// 499263 Wei-Shan Chang

use std::{env, fs, time};

fn is_prime(n: u64) -> bool {
    if n <= 1 {
        return false;
    }

    let sqrt_n = (n as f64).sqrt() as u64;

    for i in 2..=sqrt_n {
        if n % i == 0 {
            return false;
        }
    }

    true
}

fn find_prime_numbers(input: Vec<u64>) -> Vec<u64> {
    input.into_iter().filter(|&n| is_prime(n)).collect()
}

fn main() {
    // Get the input from the command line as Vec<u64>.
    let args: Vec<String> = env::args().collect();

    // Extract the input file name from command-line arguments.
    let input_filename = &args[1];

    // Read the input from the specified file.
    let input_content = fs::read_to_string(input_filename).expect("Failed to read file");

    // Parse the input content into a Vec<u64>.
    let input: Vec<u64> = input_content
        .trim()
        .split_whitespace()
        .map(|num| num.parse().unwrap())
        .collect();

    // Find the prime numbers in the input vector.
    let start_time = time::Instant::now();
    let primes: Vec<u64> = find_prime_numbers(input.clone());
    let elapsed_time = start_time.elapsed();

    let pseudo_prime_length = input.len() - primes.len();
    // Output the prime numbers and the time taken.
    println!("The count of pseudo primes from the previous test: {}", pseudo_prime_length);
    println!("Time: {:?}", elapsed_time);
}
