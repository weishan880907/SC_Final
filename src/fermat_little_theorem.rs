// 499263 Wei-Shan Chang
// TODO: add the restriction of input with Big Int.

use std::env;

fn mod_exp(base: u64, exponent: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }

    let mut result = 1;

    // This operation ensures that base remains within the range of 0 to modulus - 1.
    let base = base % modulus;

    for _ in 0..exponent {
        result = (result * base) % modulus;
    }

    result
}

fn is_prime_fermat(n: u64, k: u64) -> bool {
    if n <= 1 {
        return false;
    }

    for _ in 0..k {
        let a = rand::random::<u64>() % (n - 1) + 1; // Randomly choose 'a' between 1 and n-1
        if mod_exp(a, n - 1, n) != 1 {
            return false;
        }
    }

    true
}

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

    // Choose a value of k, which determines the number of tests.
    let k = 5;

    // Check if the number is prime by using Fermat's Little Theorem.
    if is_prime_fermat(num_to_test, k) {
        println!("{} is probably prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
}
