use std::{env, time};
use num_bigint::BigUint;
use num_traits::{One, Zero, ToPrimitive};
use rand::Rng;

fn mod_exp(base: &BigUint, exponent: &BigUint, modulus: &BigUint) -> BigUint {
    let mod_exp_time = time::Instant::now();
    if modulus.is_one() {
        return BigUint::zero();
    }

    let mut result = BigUint::one();
    let base = base % modulus;

    let mut current_exp = BigUint::zero();
    while current_exp < *exponent {
        result = (result * &base) % modulus;
        current_exp += BigUint::one();
    }
    println!("Mod Exp Time: {:?}", mod_exp_time.elapsed());

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
    let is = time::Instant::now();
    if n.is_one() {
        return false;
    }

    let mut rng = rand::thread_rng();

    for _ in 0..k {
        let a_u64: u64 = rng.gen_range(2..n.to_u64().unwrap_or(u64::MAX));
        let a = BigUint::from(a_u64);
        println!("{}, {}", a_u64,a );
        // Check gcd(a, n) != 1
        if gcd(&a, n) != BigUint::one() {
            return false;
        }

        // Check a^(n-1) ≡ 1 (mod n)
        if mod_exp(&a, &(n - BigUint::one()), n) != BigUint::one() {
            return false;
        }
    }
    println!("is fermat Time: {:?}", is.elapsed());
    true
    
}

fn main() {
    // Get the input from the command line.
    let args: Vec<String> = env::args().collect();

    // Parse the input as BigUint and set the default value as 23.
    let num_to_test: BigUint = args.get(1)
        .and_then(|arg| arg.parse().ok())
        .unwrap_or_else(|| {
            eprintln!("Invalid or missing input. Using default value: 23");
            BigUint::from(23u32)
        });

    // Choose a value of k, which determines the number of tests.
    let k = 5;

    // Measure the time taken to check if the number is probably prime.
    let start_time = time::Instant::now();
    if is_prime_fermat(&num_to_test, k) {
        println!("{} is probably prime.", num_to_test);
    } else {
        println!("{} is composite.", num_to_test);
    }
    let elapsed_time = start_time.elapsed();

    println!("Time taken: {:?}", elapsed_time);
}
