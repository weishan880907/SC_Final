// 499263 Wei-Shan Chang

use std::{env, time};
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
    let input: Vec<BigUint> = args
        .iter()
        .skip(1) // Skip the program name
        .map(|arg| arg.parse().unwrap())
        .collect();

    // Find the prime numbers in the input vector.
    let start_time = time::Instant::now();
    let primes = find_primes_sieve(input.clone());
    let elapsed_time = start_time.elapsed();

    // Output the prime numbers and the time taken.
    println!("Prime numbers: {:?}", primes);
    println!("Time: {:?}", elapsed_time);
}
