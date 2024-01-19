// 499263 Wei-Shan Chang
// TODO: add the restriction of input with Big Int.

use std::env;

fn is_prime(n: u64) -> bool {
    if n <= 1 {
        return false;
    }

    for i in 2..= n-1 {
        if n % i == 0 {
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

    // Check if the number is prime. There will not be probably prime here.
    if is_prime(num_to_test) {
        println!("{} is prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
}
