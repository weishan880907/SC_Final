// 499263 Wei-Shan Chang
// TODO: add the restriction of input with Big Int.

use std::env;

fn sieve_of_eratosthenes(limit: usize) -> Vec<bool> {
    let mut is_prime = vec![true; limit + 1];
    is_prime[0] = false;
    is_prime[1] = false;

    for i in 2..=limit {
        if is_prime[i] {
            for j in (i * i..=limit).step_by(i) {
                is_prime[j] = false;
            }
        }
    }

    is_prime
}

fn is_prime_sieve(n: usize, sieve: &Vec<bool>) -> bool {
    if n <= 1 || n >= sieve.len() {
        return false;
    }

    sieve[n]
}

fn main() {
    // Get the input from the command line.
    let args: Vec<String> = env::args().collect();

    // Parse the input as u64 and set default value as 23.
    let num_to_test: usize = args.get(1)
        .and_then(|arg| arg.parse().ok())
        .unwrap_or_else(|| {
            eprintln!("Invalid or missing input. Using default value: 23");
            23
        });

    // Generate a sieve of Eratosthenes up to the given number
    let sieve = sieve_of_eratosthenes(num_to_test);

    // Check if the number is prime using the sieve
    if is_prime_sieve(num_to_test, &sieve) {
        println!("{} is prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
}
