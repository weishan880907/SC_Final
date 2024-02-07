// 499263 Wei-Shan Chang

use std::{env, time, fs};
use num_bigint::BigUint;
use num_traits::{One, Zero, ToPrimitive};

fn sieve_of_eratosthenes(limit: &BigUint) -> Vec<bool> {
    let limit_usize: usize = limit.to_usize().unwrap_or(usize::MAX);
    let mut is_prime = vec![true; limit_usize + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    for i in 2..=limit_usize {
        if is_prime[i] {
            let mut j = i * i;
            while j <= limit_usize {
                is_prime[j] = false;
                j += i;
            }
        }
    }
    is_prime
}

fn find_primes_sieve(numbers: Vec<BigUint>) -> Vec<BigUint> {
    let max_number = numbers.iter().fold(BigUint::zero(), |acc, x| acc.max(x.clone()));
    let sieve = sieve_of_eratosthenes(&max_number);

    numbers
        .into_iter()
        .filter(|n| is_prime_sieve(n, &sieve))
        .collect()
}
fn is_prime_sieve(n: &BigUint, sieve: &Vec<bool>) -> bool {
    if n <= &BigUint::one() {
        return false;
    }

    let n_usize: usize = n.to_usize().unwrap_or(usize::MAX);
    if n_usize >= sieve.len() {
        return false;
    }

    sieve[n_usize]
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
    let primes: Vec<BigUint> = find_primes_sieve(input.clone());
    let elapsed_time = start_time.elapsed();

    let pseudo_prime_length = input.len() - primes.len();
    // Output the prime numbers and the time taken.
    println!("The count of pseudo primes from the previous test: {}", pseudo_prime_length);
    println!("Time: {:?}", elapsed_time);
}
