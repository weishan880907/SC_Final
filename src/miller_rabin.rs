// 501900 Ximing Zhang

use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, Zero};
use num_bigint::RandBigInt;
use std::{env, time, fs};

// Utility function to do modular exponentiation.
// It returns (x^y) % p
fn power(x: &BigInt, y: &BigInt, p: &BigInt) -> BigInt {
    let mut result = One::one();
    let mut base = x % p;
    let mut exponent = y.clone();
    //let two = 2.to_bigint().unwrap();

    while exponent > Zero::zero() {
        if &exponent & &One::one() == One::one() {
            result = (result * &base) % p;
        }
        exponent >>= 1;
        base = (&base * &base) % p;
    }

    result
}

// This function is called for all k trials. It returns
// false if n is composite and returns true if n is
// probably prime.
fn miller_test(d: BigInt, n: &BigInt) -> bool {
    let mut rng=rand::thread_rng();
    let two = 2.to_bigint().unwrap();

    let n_minus_3 = n - 1.to_bigint().unwrap();
    let a = rng.gen_bigint_range(&two, &n_minus_3);

    let mut x = power(&a, &d, n);

    if x == One::one() || x == (n - 1) {
        return true;
    }

    let mut d = d;
    while d < (n - 1) {
        x = (&x * &x) % n;
        d <<= 1;

        if x == One::one() {
            return false;
        }
        if x == (n - 1) {
            return true;
        }
    }

    false
}

// It returns false if n is composite and returns true if n
// is probably prime. k is an input parameter that determines
// accuracy level. Higher value of k indicates more accuracy.
fn is_prime(n: &BigInt, k: i64) -> bool {
    if n <= &One::one() || n == &4.to_bigint().unwrap() {
        return false;
    }
    if n <= &3.to_bigint().unwrap() {
        return true;
    }

    let mut d:BigInt = n - 1;
    while &d & &One::one() == Zero::zero() {
        d >>= 1;
    }

    for _ in 0..k {
        if !miller_test(d.clone(), n) {
            return false;
        }
    }

    true
}

// Driver program
fn main() {
    // Get the input from the command line as Vec<BigInt>.
    let args: Vec<String> = env::args().collect();

    // Extract the upper bound for primality testing from command-line arguments.
    let upper_bound: BigInt = args[1].parse().unwrap();

    // Check if the upper bound is less than 2.
    if upper_bound < BigInt::from(2u32) {
        println!("Upper bound should be at least 2.");
        return;
    }
    // Create a vector from 1 to the upper bound using a while loop.
    let mut current_num = BigInt::one();
    let mut input_vector = Vec::new();
    while current_num <= upper_bound.clone() {
        input_vector.push(current_num.clone());
        current_num += BigInt::one();
    }

    // Choose a value of k, which determines the number of tests.
    let k = 2;

    // Find the probably prime numbers in the input vector.
    let start_time: time::Instant = time::Instant::now();
    let primes: Vec<BigInt> = input_vector
        .into_iter()
        .filter(|num| is_prime(num, k))
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