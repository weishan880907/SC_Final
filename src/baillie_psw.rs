// 403581 Yunzhou Lu

use num_bigint::{ToBigInt, BigInt};
use num_traits::{One, Zero};
use std::time::Instant;
use std::env;
use std::fs::File;
use std::io::Write;

/* This program decides whether a number passes the Baillie-PSW test.
   It is the combination of strong lucas test and a Miller-Rabin test with base 2.*/
fn main() {
    let start_time = Instant::now();

    let args :Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: cargo run --bin baillie_psw <upper bound>");
        std::process::exit(1);
    }

    let mut file = match File::create("baillie_psw.txt") {
        Ok(file) => file,
        Err(err) => {
            panic!("Failed to create file: {}", err);
        }
    };

    let upper_bound: &str = &args[1];
    let upper_bound = BigInt::parse_bytes(upper_bound.as_bytes(), 10).unwrap();

    let mut r = BigInt::from(1);
    // let mut c = 0;

    while &r <= &upper_bound {
        r += 1;
        if baillie(&r) {
            // c+=1;
            if let Err(err) = write!(&mut file, "{} ", r) {
                panic!("Failed to write to file: {}", err);
            }
        }
        
    }

    // println!("{c}");
    //print the total time
    let duration = start_time.elapsed();
    println!("Time : {:?}", duration);
}

fn baillie(r: &BigInt) -> bool {
    if is_small_prime(r) {
        return true
    } else {
        if !is_perfect_square(&r) && not_devisable_by_small_prime(r) 
        && strong_lucas_tests(r) && is_prime_miller(r) {
            return true
        } else {
            return false
        }
    }
}

fn miller_test_base_two(d: BigInt, n: &BigInt) -> bool {
    let mut x = power(&BigInt::from(2), &d, n);

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

fn is_prime_miller(n: &BigInt) -> bool {
    if n == &One::one() || n == &4.to_bigint().unwrap() {
        return false;
    }
    let mut d:BigInt = n - 1;
    while &d & &One::one() == Zero::zero() {
        d >>= 1;
    }

    if !miller_test_base_two(d.clone(), n) {
        return false;
    }
    true
}

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

fn strong_lucas_tests(num: &BigInt) -> bool {
/*  This function perform a strong Lucas probable prime test with following steps:
    1. decomposite n + 1 = 2 ^ s * d
    2. choose D, set P, Q
    3. calculate the d-th and (d+1)-th element V_d, V_{d+1}
    4. calculate U_d, decide whether U_d ≡ 0 (mod n)
        if true, n is probably prime
    5. calculate V_{d*2^r}, decide whether V_{d*2^r} ≡ 0 (mod n), for 0 <= r < s
        if true, n is probably prime */

    let iterations = 2;  //number of tests with different P, Q, D
    
    let (s, d) = decomposition(&(num+1));  //1.
    let mut delta = 1;

    'outer: for _ in 0..iterations {
        delta = search_for_d(num, delta);  //2.
        let (p,q) = suitabale_p_q_to_d(&delta);
        let (mut u,v, mut b_d) =lucas_chain(&d, &p,& q, num); //3.

        let big_u_d = (2 * &v - &p * &u) % num;   //4.
        if big_u_d == BigInt::from(0) {
            continue 'outer;
        }

        for _ in 0..s {                     //5.
            if u == BigInt::from(0) {
                continue 'outer;
            }
            u = &u * &u - &b_d * 2;
            b_d = &b_d * &b_d;
            (u, b_d) = ((u % num + num) % num, b_d % num);
        }
        return false
    }
    true
}

fn jacobi_symbol(odd_int: &BigInt, int: &BigInt) -> i32 {
/*  Given positive odd integer odd_int, and integer int, this algorithm returns the
    Jacobi symbol (int / odd_int), which for odd_int prime is also the Legendre symbol.

    See Definition 2.3.2, Algorithm 2.3.5 from book: Richard E. Crandall, Carl Pomerance
    Prime Numbers: A Computational Perspective. 2. ed.*/

    if odd_int % 2 == BigInt::from(0) {
        panic!("Expected the lower argument of Jacobi symbol being odd, recieved even");
    }

    let mut numerator = ((int % odd_int) + odd_int) % odd_int;
    let mut denominator = odd_int.clone();
    let mut t: i32 = 1;

    while numerator != BigInt::from(0) {
        while &numerator % 2 == BigInt::from(0) {

            numerator /= 2;

            if &denominator % 8 == BigInt::from(3) 
            || &denominator % 8 == BigInt::from(5) {  // denominator mod 8 ∈ {3,5}
                t *= -1; 
            }
        }

        (numerator, denominator) = (denominator, numerator);  //swap variables

        if (&numerator % 4 == BigInt::from(3) || &numerator % 4 == BigInt::from(-1))
        && &denominator % 4 == BigInt::from(3) { //numerator ≡ denominator ≡ 3 (mod 4)
            t *= -1;
        }
        numerator = numerator % &denominator;

        if denominator == BigInt::from(1) {
            return t;
        }
    }
    0
}

fn search_for_d(num: &BigInt, start: i32) -> i32 {
/*  Given n, one technique for choosing D is to use trial and error to find the
    first D in the sequence 5,9,13,17...(for start=1) such that (D / n) = −1*/
    let mut test_d = start;
    loop {
        test_d += 4;

        if jacobi_symbol(num, &BigInt::from(test_d)) == -1 
        && !is_perfect_square(&BigInt::from(test_d)) {
            return test_d
        }
    }
}

fn suitabale_p_q_to_d(d: &i32) -> (i32, i32) {
/*  Once we have D, we set P = least odd number exceeding sqrt(d) and Q = (P^2 − D) / 4.
    It is a good idea to check that n has no prime factors in common with P or Q. 
    This method of choosing D, P, and Q was suggested by Baillie and Wagstaff.  */
    let sqrt_d = (*d as f32).sqrt();
    let mut big_p = sqrt_d.floor() as i32;
    if big_p % 2 == 0 {
        big_p += 1;
    } else {
        big_p += 2;
    }
    (big_p, (big_p * big_p - d) / 4)
}

fn decomposition(num: &BigInt) -> (u64, BigInt) {
/*  This function decomposite for input num, num = 2 ^ s * d, and return (s, d)  */
    let mut d = num.clone();
    let mut s = 0;

    while &d % 2 == BigInt::from(0) {
        d /= 2;
        s += 1;
    }

    (s, d)
}

fn lucas_chain(computing_position: &BigInt, a: &i32, b: &i32, num: &BigInt) -> (BigInt, BigInt, BigInt) {
/*  For a sequence x_0, x_1, . . . with the computation rules
      1. x_{2j} = x_j^2 - 2* b^j,
      2. x_{2j+1} = x_j * x_{j+1} - a * b^j,
      3. x_0 = 2, x_1 = a,
    this algorithm computes the pair (x_n, x_{n+1}) for a given positive integer n.

    We have n in binary as (n_0, n_1, . . . , n_{B−1}),
        for n_j == 0, we travel from (x_j,x_{j+1}) to (x_{2j},x_{2j+1}),        (1)
        for n_j == 1, we travel from (x_j,x_{j+1}) to (x_{2j+1},x_{2j+2}),      (2)
    for example m is 97. We travel from 0,1 to 97,98 as follows:
            0,1 → 1,2 → 3,4 → 6,7 → 12,13 → 24,25 → 48,49 → 97,98 
    while 97 = (1     1     0     0       0       0       1)_2
    
    See Section 3.6.3, Algorithm 3.6.7 from the book mentioned above*/

    let mut u: BigInt = 2.to_bigint().unwrap();  //initial x_0
    let mut v: BigInt = a.to_bigint().unwrap();  //initial x_1
    let mut b_j: BigInt = 1.to_bigint().unwrap();  //b^0

    let position_binary = format!("{:b}", computing_position); //compute the binary form for position

    for (index, digit) in position_binary.chars().enumerate() {
        match digit {
            '0' => {                        //match (1)
                v = &u * &v - a * &b_j;
                u = &u * &u - &b_j * 2;
                b_j = b_j.pow(2);
            }
            '1' => {                        //match (2)
                u = &u * &v - a * &b_j;
                v = &v * &v - &b_j * b * 2;
                b_j = b_j.pow(2);
                b_j = b_j * b;         
            }
            _ => {                          //I hope nothing shity would happen here
                panic!("Failed at the {} step of computing the lucas chain", index) 
            }
        }
        (u, v, b_j) = ((u % num + num) % num, (v % num + num) % num, (b_j % num + num) % num);
    }
    (u,v,b_j)
}

fn is_perfect_square(num: &BigInt) -> bool {
/*  This function decides efficiently, whether num is square of some number
    using Newton's method */
    let mut x0 = num.clone();

    while {
        let x1 = (&x0 + num / &x0) / 2;
        if &x0 - &x1 < BigInt::from(1) {
            false
        } else {
            x0 = x1;
            true
        }
    } {}

    &x0 * &x0 == *num
}

fn is_small_prime(num: &BigInt) -> bool {
    for i in SMALL_PRIME_LIST{
        if num == &BigInt::from(i){
            return true
        }
    }
    false
}

fn not_devisable_by_small_prime(num: &BigInt) -> bool {
/*  This function decides whether num is not devided by some small prime up to 150 */

    for i in SMALL_PRIME_LIST {
        if num % i == BigInt::from(0) {
            return false
        }
    }
    true
}

const SMALL_PRIME_LIST: [i32; 25] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 
47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];