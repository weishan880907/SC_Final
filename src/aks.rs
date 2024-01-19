// 499263 Wei-Shan Chang
// TODO: add the restriction of input with Big Int.

use std::env;

// Function to calculate binomial coefficient
fn binomial(n: usize, k: usize) -> usize {
    let mut result = 1;
    let mut num = n;
    let mut den = 1;

    for _ in 0..k {
        result *= num;
        result /= den;
        num -= 1;
        den += 1;
    }

    result
}

// Function to calculate the greatest common divisor (GCD) using Euclid's algorithm
fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

// Function for modular exponentiation
fn power_mod(base: usize, exponent: usize, modulus: usize) -> usize {
    let mut result = 1;

    for _ in 0..exponent {
        result = (result * base) % modulus;
    }

    result
}

// Euler's totient function
fn phi(n: usize) -> usize {
    (1..n).filter(|&i| gcd(i, n) == 1).count()
}

// AKS primality testing function
fn is_aks_prime(n: usize) -> bool {
    if n <= 1 {
        return false;
    }

    // Check if there exists a non-trivial divisor
    let r = (2..n).find(|&i| gcd(i, n) > 1);
    if let Some(_r) = r {
        return false;
    }

    // Check AKS primality condition for a range of 'a' values
    for a in 1..(phi(n) * (n as f64).log2() as usize) {
        if gcd(a, n) > 1 {
            continue;
        }

        let bin_coeff = binomial(n, a);
        let a_exp_n_mod_n = power_mod(a, n, n);
        let bin_coeff_mod_n = bin_coeff % n;

        // Check AKS primality condition
        if bin_coeff_mod_n != a_exp_n_mod_n {
            return true;
        }
    }

    false
}

// Main function
fn main() {
    // Get the input from the command line.
    let args: Vec<String> = env::args().collect();

    // Parse the input as u64 and set default value as 23.
    let num_to_test: u64 = args.get(1)
        .and_then(|arg| arg.parse().ok())
        .unwrap_or_else(|| {
            eprintln!("Invalid or missing input. Using default value: 23");
            23
        });

    // Perform AKS primality test and print the result
    if is_aks_prime(num_to_test as usize) {
        println!("{} is prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
}
