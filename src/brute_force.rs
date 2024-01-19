// 499263 Wei-Shan Chang
use std::env;

fn is_prime(n: u64) -> bool {
    if n <= 1 {
        return false;
    }

    let sqrt_n = (n as f64).sqrt() as u64;

    for i in 2..=sqrt_n {
        if n % i == 0 {
            return false;
        }
    }

    true
}

fn main() {
    // Get the input from the environment variable
    let input_str = env::var("NUM_TO_TEST").unwrap_or_else(|_| String::from("23"));

    // Parse the input as u64
    let num_to_test: u64 = input_str.parse().unwrap_or_else(|_| {
        eprintln!("Invalid input. Using default value: 23");
        23
    });

    // Check if the number is prime
    if is_prime(num_to_test) {
        println!("{} is prime.", num_to_test);
    } else {
        // In this simple example, we're considering any non-prime as composite.
        // You might want to use a more sophisticated algorithm for larger inputs.
        println!("{} is composite.", num_to_test);
    }
}
