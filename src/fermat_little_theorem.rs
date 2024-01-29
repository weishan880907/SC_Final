// 499263 Wei-Shan Chang

use std::{env, time};
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
        let a_u64: u64 = rng.gen_range(2..n.to_u64().unwrap_or(u64::MAX));
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
    // Get the input from the command line.
    let args: Vec<String> = env::args().collect();

    // Parse the input as Vec<BigUint> and set the default value as [23u32].
    let input: Vec<BigUint> = args
        .iter()
        .skip(1) // Skip the program name
        .map(|arg| arg.parse().unwrap())
        .collect();

    // Choose a value of k, which determines the number of tests.
    let k = 5;

    // Measure the time taken to check if each number in the input vector is probably prime.
    let start_time = time::Instant::now();
    let primes: Vec<BigUint> = input
        .into_iter()
        .filter(|num| is_prime_fermat(num, k))
        .collect();
    let elapsed_time = start_time.elapsed();

    // Output the prime numbers and the time taken.
    println!("Prime numbers: {:?}", primes);
    println!("Time: {:?}", elapsed_time);
}
