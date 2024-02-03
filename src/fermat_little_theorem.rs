// 499263 Wei-Shan Chang

use std::{env, time, fs};
use num_bigint::BigUint;
use num_traits::{One, Zero, ToPrimitive};
use rand::Rng;

fn mod_exp(base: &BigUint, exponent: &BigUint, modulus: &BigUint) -> BigUint {
    if modulus.is_one() {
        return BigUint::zero();
    }

    let mut result = BigUint::one();
    let mut base_power = base.clone() % modulus;

    let mut exp = exponent.clone();
    while exp > BigUint::zero() {
        if exp.clone() % BigUint::from(2u32) == BigUint::one() {
            result = (result * &base_power) % modulus.clone();
        }
        base_power = (base_power.clone() * base_power) % modulus.clone();
        exp = exp >> 1;
    }

    result
}

fn gcd(a: &BigUint, b: &BigUint) -> BigUint {
    if a.is_zero() {
        b.clone()
    } else {
        gcd(&(b % a), a)
    }
}

fn is_prime_fermat(n: &BigUint, k: usize) -> bool {
    if n.is_one() {
        return false;
    }

    let mut rng = rand::thread_rng();

    for _ in 0..k {
        let a_u64: u64 = rng.gen_range(1..n.to_u64().unwrap_or(u64::MAX));
        let a = BigUint::from(a_u64);
        // Check gcd(a, n) != 1
        if gcd(&a, n) != BigUint::one() {
            return false;
        }

        // Check a^(n-1) ≡ 1 (mod n)
        if mod_exp(&a, &(n - BigUint::one()), n) != BigUint::one() {
            return false;
        }
    }
    true
}

fn main() {
    // Get the input from the command line as Vec<BigUint>.
    let args: Vec<String> = env::args().collect();

    // Extract the upper bound for primality testing from command-line arguments.
    let upper_bound: BigUint = args[1].parse().unwrap();

    // Check if the upper bound is less than 2.
    if upper_bound < BigUint::from(2u32) {
        println!("Upper bound should be at least 2.");
        return;
    }
    // Create a vector from 1 to the upper bound using a while loop.
    let mut current_num = BigUint::one();
    let mut input_vector = Vec::new();
    while current_num <= upper_bound.clone() {
        input_vector.push(current_num.clone());
        current_num += BigUint::one();
    }

    // Choose a value of k, which determines the number of tests.
    let k = 2;

    // Find the probably prime numbers in the input vector.
    let start_time: time::Instant = time::Instant::now();
    let primes: Vec<BigUint> = input_vector
        .into_iter()
        .filter(|num| is_prime_fermat(num, k))
        .collect();
    let elapsed_time = start_time.elapsed();

    // Output the probably prime numbers and the time taken.
    println!("Time: {:?}", elapsed_time);

    // Split the string by "/" and collect it into a vector of parts.
    let parts: Vec<&str> = args[0].split('/').collect(); // target/debug/fermat_little_theorem
    let output_filename = format!("{}.txt", parts[parts.len() - 1]);
    let primes_str: String = primes.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(" ");
    fs::write(&output_filename, primes_str).expect("Failed to write to file");
}