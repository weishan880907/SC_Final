// 499263 Wei-Shan Chang

use std::{env, time};
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
    let input: Vec<BigUint> = args
        .iter()
        .skip(1) // Skip the program name
        .map(|arg| arg.parse().unwrap())
        .collect();

    // Find the prime numbers in the input vector.
    let start_time = time::Instant::now();
    let primes: Vec<BigUint> = find_prime_numbers(input.clone());
    let elapsed_time = start_time.elapsed();

    // Output the prime numbers.
    println!("Prime numbers: {:?}", primes);
    println!("Time: {:?}", elapsed_time);
}
