// 499263 Wei-Shan Chang

use std::{env, time, fs};
use num_bigint::BigUint;
use num_traits::{One, Zero, ToPrimitive};

fn binomial(n: &BigUint, k: &BigUint) -> BigUint {
    let mut result = BigUint::one();
    let mut num = n.clone();
    let mut den = BigUint::one();

    let mut i = BigUint::zero();
    while i < k.clone() {
        result *= &num;
        result /= &den;
        num -= BigUint::one();
        den += BigUint::one();
        i += BigUint::one();
    }

    result
}

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

fn phi(n: &BigUint) -> BigUint {
    let mut count = BigUint::zero();
    let mut i = BigUint::one();

    while i < *n {
        if gcd(&i, n) == BigUint::one() {
            count += BigUint::one();
        }
        i += BigUint::one();
    }

    count
}

fn aks_prime_test(numbers: Vec<BigUint>) -> Vec<BigUint> {
    numbers
        .into_iter()
        .filter(|num| {
            let n = num.to_usize().unwrap_or(0);
            if n <= 1 {
                return false;
            }

            // Check if there exists a non-trivial divisor
            if (2..n).any(|i| gcd(&BigUint::from(i), num) > BigUint::one()) {
                return false;
            }

            // Check AKS primality condition for a range of 'a' values
            let phi_n = phi(&num);
            let mut a = BigUint::one();
            let upper_bound = phi_n.clone() * BigUint::from(
                num.bits() / 64 + if num.bits() % 64 == 0 { 0 } else { 1 }
            );
                        
            while &a <= &upper_bound {
                if gcd(&a, &num) > BigUint::one() {
                    continue;
                }
            
                let bin_coeff = binomial(&num, &a);
                let a_exp_n_mod_n = mod_exp(&a, &(num.clone() - BigUint::one()), &num);
                let bin_coeff_mod_n = bin_coeff % num;
            
                // Check AKS primality condition
                if bin_coeff_mod_n != a_exp_n_mod_n {
                    return true;
                }
            
                a += BigUint::one();
            }

            false
        })
        .collect()
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
    let primes: Vec<BigUint> = aks_prime_test(input.clone());
    let elapsed_time = start_time.elapsed();

    let pseudo_prime_length = primes.len() - input.len();
    // Output the prime numbers and the time taken.
    println!("The number of pseudo primes: {}", pseudo_prime_length);
    println!("Time: {:?}", elapsed_time);
}
