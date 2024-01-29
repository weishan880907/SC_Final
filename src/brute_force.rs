use std::env;
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

    // Measure the time taken to check if the number is prime.
    let start_time = std::time::Instant::now();
    if is_prime(&num_to_test) {
        println!("{} is prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
    let elapsed_time = start_time.elapsed();

    println!("Time taken: {:?}", elapsed_time);
}
