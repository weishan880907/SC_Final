// 499263 Wei-Shan Chang

use std::{env, time};
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive};

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
    // Get the input from the command line.
    let args: Vec<String> = env::args().collect();

    // Parse the input as BigUint and set default value as 23.
    let num_to_test: BigUint = args.get(1)
        .and_then(|arg| arg.parse().ok())
        .unwrap_or_else(|| {
            eprintln!("Invalid or missing input. Using default value: 23");
            BigUint::from(23u32)
        });

    // Generate a sieve of Eratosthenes up to the given number
    let sieve = sieve_of_eratosthenes(&num_to_test);

    // Measure the time taken to check if the number is prime using the sieve
    let start_time = time::Instant::now();
    if is_prime_sieve(&num_to_test, &sieve) {
        println!("{} is prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
    let elapsed_time = start_time.elapsed();

    println!("Time taken: {:?}", elapsed_time);
}
