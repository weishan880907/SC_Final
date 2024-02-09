// Yunzhou 403581 Feb. 2023

/* This project yields to prove if a given number is prime using ECPP method.
   Unfortunately, due to some time limitation, this rust implementation
   does not work yet. I will modify it in the future and make the program 
   run successfully*/

use num::integer::Roots;
use num_bigint::{ToBigInt,  BigInt, RandBigInt};
use num_traits::{Signed, One, Zero};
use rand::Error;
use std::collections::{HashSet,HashMap};
use std::f64::consts::{PI,E};
use num::complex::{Complex, ComplexFloat};
use std::time::Instant;
use poly::FieldElement;
use std::{env,fs};

fn main() {
    let start_time = Instant::now();

    let args: Vec<String> = env::args().collect();

    // Extract the input file name from command-line arguments.
    let input_filename = &args[1];

    // Read the input from the specified file.
    let input_content = fs::read_to_string(input_filename).expect("Failed to read file");
    
    // Parse the input content into a Vec<BigUint>.
    let input: Vec<BigInt> = input_content
        .trim()
        .split_whitespace()
        .map(|num| num.parse().unwrap())
        .collect();

    let mut count = 0; //count the number of composite numbers in input flie
    for n in input {
        if !ecpp_prime_test(&n) {
            count += 1;
        }
    }
    println!("The count of pseudo primes from the previous test: {}", count);
    
    let elapsed_time = start_time.elapsed();
    println!("Time : {:?}", elapsed_time);
}

fn ecpp_prime_test(n: &BigInt) -> bool {
    //run ecpp test for n
    let mut n = n.clone();
    let (m, q, _, _, _, a, b) = choose_discriminant(&n);
    if &n % 2 == BigInt::zero() {
        return true;
    }

    loop {
        match decide_q(&a, &b, &n, &m, &q) {
            Some(x) => {
                n = x;
                for i in SMALL_PRIME_LIST {
                    if n == BigInt::from(i){
                        return true;
                    }
                }
            }
            None => {
                return false;
            }
        }
    }
}

//list of small primes
const SMALL_PRIME_LIST: [i32; 95] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 
233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499];


fn choose_discriminant(n: &BigInt) -> (BigInt, BigInt, i64, BigInt, BigInt, BigInt, BigInt) {
    /* Select a fundamental discriminant D by increasing value of h(D) for which (D/n) = 1 and
       for which we are successful in finding a solution u^2 + |D|v^2 = 4n via modified_cornacchia_smith,
       yielding possible curve orders m 
       
       See Algorithm 7.6.3 from book: Richard E. Crandall, Carl Pomerance
       Prime Numbers: A Computational Perspective. 2. ed.*/

    let list = [-3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163];
/*     let list = [-3, -4, -7, -8, -11, -19, -43, -67, -163,
    -15, -20, -24, -35, -40, -51, -52, -88, -91, -115, -123, -148, -187, -232, -235, -267, -403, -427,
    -23, -31, -59, -83, -107, -139, -211, -283, -307, -331, -379, -499, -547, -643, -883, -907,
    -39, -55, -56, -68, -84, -120, -132, -136, -155, -168, -184, -195, -203, -219, -228, -259, -280, -291, -292,
    -312, -323, -328, -340, -355, -372, -388, -408, -435, -483, -520, -532, -555, -568, -595, -627, -667, -708,
     -715, -723, -760, -763, -772, -795, -955, -1003, -1012, -1027, -1227, -1243, -1387, -1411, -1435, -1507, -1555]; */
    for i in list {
        if let Some((u, v)) = modified_cornacchia_smith(n, &BigInt::from(i)) {
            if i == -4 {
                for m in [n+1+&u,n+1-&u,n+1+2*&v,n+1-2*&v] {
                    let result = find_a_b_q_for_m(n, &m, i);
                    match result {
                        Some((a, b, q)) => return (m, q, i, u, v, a, b),
                        None => continue
                    }
                }
            } else if i == -3 {
                for m in [n+1+&u,n+1-&u,n+1+(&u+3*&v)/2,n+1-(&u+3*&v)/2,n+1+(&u-3*&v)/2,n+1-(&u-3*&v)/2] {
                    let result = find_a_b_q_for_m(n, &m, i);
                    match result {
                        Some((a, b, q)) => return (m, q, i, u, v, a, b),
                        None => continue
                    }
                }
            } else {
                for m in [n+1+&u,n+1-&u] {
                    let result = find_a_b_q_for_m(n, &m, i);
                    match result {
                        Some((a, b, q)) => return (m, q, i, u, v, a, b),
                        None => continue
                    }
                }
            }
        } else {
            continue;
        }
    }
    panic!("Failed to find a discriminant D");
}

fn kq_factor(m: &BigInt, n: &BigInt) -> Option<(BigInt, BigInt)> {
    /* This function factorizes m as m = kq, where k > 1 and q is a probable prime > n  */
    let mut k = BigInt::from(1);
    let bound = integer_part_of_square_root(n) 
        + 2 * integer_part_of_square_root(&integer_part_of_square_root(n)) + 1;
    let mut q = m.clone();
    for i in SMALL_PRIME_LIST {
        if i != 2 && baillie(&q) == true {
            return Some((k, q))
        }
        if q == BigInt::from(i) && &BigInt::from(i) > &bound {
            return Some((k, q))
        }
        if &q <= &bound {
            return None
        } else {
            while &q % i == BigInt::from(0) {
                q = q / i;
                k = k * i;
            }
        }
    }
    
    if k == BigInt::from(1) && m > n {
        return None
    }

    if baillie(&q) == false {
        return None
    }

    Some((k, q))
}

fn modified_cornacchia_smith(n: &BigInt, big_d: &BigInt) -> Option<(BigInt, BigInt)> {
    /* Input n prime, bid_b negative integer. Represent 4n as x^2 + |D|y^2 (modified 
       Cornacchia– Smith)) Given a prime n and −4n < D < 0 with D ≡ 0, 1 (mod 4), this 
       algorithm either reports that no solution exists, or returns a solution (x, y).

       See Algorithm 2.3.13 from book: Richard E. Crandall, Carl Pomerance
       Prime Numbers: A Computational Perspective. 2. ed.*/
    if n * -4 >= *big_d  && *big_d >= BigInt::from(0) {
        panic!("-4n < D < 0 not true")
    } else if (big_d % 4 + 4) % 4 < BigInt::from(0) || (big_d % 4 + 4) % 4 > BigInt::from(1) {
        panic!("D ≡ 0, 1 (mod 4) not true")
    }

    if n == &BigInt::from(2) {
        if is_perfect_square(&(big_d + 8)) {
            return Some((square_root(&(big_d + 8)),BigInt::from(1)))
        } else {
            return None
        }
    }

    if jacobi_symbol(n, big_d) < 1 {
        return None;
    }

    let mut x_0 = match square_root_mod_prime(n, big_d) {
        Ok(x) => x,
        Err(err_msg) => {
            panic!("Failed at modified_cornacchia_smith(): {}",err_msg);
        }
    };

    if &x_0 % 2 != big_d % 2 && &x_0 % 2 != (big_d % 2) * (-1) {
        x_0 = n - x_0;
    }

    let (mut a,mut b) = (n * 2, x_0);
    let c = integer_part_of_square_root(&(4 * n));

    while b > c {
        let a_mod_b = a % &b;
        (a, b) = (b, a_mod_b);
    }

    let t = n * 4 - &b *&b;
    
    if (n * 4 - &b *&b) % (big_d *(-1)) != BigInt::from(0)
    || !is_perfect_square(&((n * 4 - &b *&b) / (big_d *(-1)))) {
        return None
    }

    Some((b, square_root(&(t / (big_d.abs())))))
}

fn jacobi_symbol(odd_int: &BigInt, int: &BigInt) -> i64 {
    /*  Given positive odd integer odd_int, and integer int, this algorithm returns the
        Jacobi symbol (int / odd_int), which for odd_int prime is also the Legendre symbol.
    
        See Definition 2.3.2, Algorithm 2.3.5 from book: Richard E. Crandall, Carl Pomerance
        Prime Numbers: A Computational Perspective. 2. ed.*/
    
    if odd_int % 2 == BigInt::from(0) {
        panic!("Expected the lower argument of Jacobi symbol being odd, recieved even");
    }

    let mut numerator = ((int % odd_int) + odd_int) % odd_int;
    let mut denominator = odd_int.clone();
    let mut t: i64 = 1;

    while numerator != BigInt::from(0) {
        while &numerator % 2 == BigInt::from(0) {
            numerator /= 2;

            if &denominator % 8 == BigInt::from(3) 
            || &denominator % 8 == BigInt::from(5) {  // denominator mod 8 ∈ {3,5}
                t *= -1; 
            }
        }

        (numerator, denominator) = (denominator, numerator);  //swap variables

        if (&numerator % 4 + 4) % 4 == BigInt::from(3)
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

fn square_root(num: &BigInt) -> BigInt {
    /*  This function computes efficiently the square root of some number
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
    
    x0
}

fn big_b(n: &BigInt) -> u64 {
    /* This function compute for nonnegative integers N the number of bits
       in the binary representation of N , except that B(0) = 0.  */

    if n == &BigInt::from(0) {
        return 0
    } else {
        let binary_string = n.to_str_radix(2);
        return binary_string.len() as u64
    }
}

fn integer_part_of_square_root(num: &BigInt) -> BigInt {
    /*  This function computes efficiently the integer part of square root
        of some number using Newton's method */

    let big_b = big_b(num);
    let power = match big_b % 2 {
        1 => big_b / 2 + 1,
        0 => big_b / 2,
        _ => panic!("Failed with computing integer part of square root")
    };
    let mut x = a_power_b(2, power);
    loop {
        let y = (&x + num / &x) / 2;
        if y >= x {
            return x
        }
        x = y;
    }
}

fn square_root_mod_prime(odd_prime: &BigInt, a: &BigInt) -> Result<BigInt, Error> {
    /* Given an odd prime and an integer a with (a / odd_prime) = 1, this algorithm
       returns a solution x to x^2 ≡ a (mod odd_prime)
       (Notice a here needs not to be positive) 

       See Algorithm 2.3.8 from book: Richard E. Crandall, Carl Pomerance
       Prime Numbers: A Computational Perspective. 2. ed.*/
    if jacobi_symbol(odd_prime, a) != 1 {
        return Err(rand::Error::new("No result"));
    }

    let a = (a % odd_prime + odd_prime) % odd_prime;
    if odd_prime % 8 == BigInt::from(3) || odd_prime % 8 == BigInt::from(7) {
        let b = (odd_prime + 1) / 4;
        let y = a_power_b_mod_c(&a, &b, odd_prime);
        if &y * &y % odd_prime == a {
            return Ok(y)
        }
    }

    if odd_prime % 8 == BigInt::from(5) {
        let b = (odd_prime + 3) / 8;
        let mut y = a_power_b_mod_c(&a, &b, odd_prime);
        let c = (&y * &y) % odd_prime;
        if c != a {
            let b = (odd_prime - 1) / 4;
            let y_y = a_power_b_mod_c(&BigInt::from(2), &b, odd_prime);
            y = (y * y_y) % odd_prime;
        }
        if &y * &y % odd_prime == a {
            return Ok(y)
        }
    }

    if odd_prime % 8 == BigInt::from(1) {
        let mut d = BigInt::from(2);
        while jacobi_symbol(odd_prime, &d) != -1 {
            d = rand::thread_rng()
            .gen_bigint_range(&BigInt::from(2), &BigInt::from(odd_prime - 1));
        }
        let (s,t) = decomposition(&(odd_prime - 1));
        let big_a = a_power_b_mod_c(&a, &t, odd_prime);
        let big_d = a_power_b_mod_c(&d, &t, odd_prime);

        let mut big_d_power_two_i = big_d.clone();
        let mut big_a_d_m = big_a.clone();
        let mut b = (odd_prime - 1) / &t;
        let mut m = BigInt::zero();
        let mut i = BigInt::one();

        for _ in 0..s {
            b = b / 2;
            let res = a_power_b_mod_c(&big_a_d_m, &b, odd_prime);
            if res == BigInt::from(-1) || res == odd_prime - 1 {
                big_a_d_m= (big_a_d_m * &big_d_power_two_i) % odd_prime;
                m += &i;
            }
            big_d_power_two_i = (&big_d_power_two_i * &big_d_power_two_i) % odd_prime;
            i = i * 2;
        }

        let y =(a_power_b_mod_c(&a, &((t + BigInt::one()) / 2), odd_prime) * a_power_b_mod_c(&big_d, &(m/2), odd_prime)) % odd_prime;
        if &y * &y % odd_prime == a {
            return Ok(y)
        }
    }

    return Err(rand::Error::new("Input not odd by computing square root mod prime"));
}

fn decomposition(num: &BigInt) -> (u64, BigInt) {
    /*  This function decomposite for input num, num = 2 ^ s * d, d odd, and return (s, d)  */

    let mut d: BigInt = num.clone();
    let mut s = 0;

    loop {
        let e = d.clone();
        d /= 2;
        if &d * 2 != e {
            return (s, e)
        }
        s += 1;
    }
}

fn a_power_b_mod_c(a: &BigInt, b: &BigInt, c: &BigInt) -> BigInt {
    //computes a^b mod c

    let mut x = a.clone();
    let mut y = BigInt::from(1);
    let mut p = b.clone();
    while p > BigInt::from(0) {
        if &p % 2 != BigInt::from(0) {
            y *= &x;
            y = y % c;
        }
        x = &x * &x;
        x = x % c;
        p = p / 2;
    }
    y
}

fn a_power_b(a: u64, b: u64) -> BigInt {
    let mut x = BigInt::from(a);
    let mut y = BigInt::from(1);
    let mut p = b;
    while p > 0 {
        if &p % 2 != 0 {
            y *= &x;
        }
        x = &x * &x;
        p = p / 2;
    }
    y
}

fn hilbert(big_d: i64) -> (i64, Vec<i64>, HashSet<(i64, i64, i64)>) {
    /* Given a (negative) fundamental discriminant D, this algorithm returns any desired
       combination of the class number h(D), the Hilbert class polynomial T ∈ Z[X]
       (whose degree is h(D)), and the set of reduced forms (a, b, c) of discriminant D
       (whose cardinality is h(D)). 
       
       See Algorithm 7.5.8 from book: Richard E. Crandall, Carl Pomerance
       Prime Numbers: A Computational Perspective. 2. ed.*/

    let reduced_form = reduced_form(big_d);
    let mut a_inverse_sum = 0.0;
    for tuple in &reduced_form {
        let reciprocal = 1.0 / tuple.0 as f64;
        a_inverse_sum += reciprocal;
    }

    let precision = (PI * ((big_d.abs() as f64).sqrt()) * a_inverse_sum / 10.0_f64.log(E)).round() as i64 + 10;
    let mut t = vec![Complex::new(1.0, 0.0)];
    let mut b = ((big_d % 2) + 2) % 2;
    let r = (big_d.abs() / 3).sqrt();
    let mut h = 0;
    let mut red = HashSet::new();
    while b <= r {
        let m = (b * b - big_d) / 4;
        for a in 1..m.sqrt() + 1 {
            if m % a != 0 {
                continue; 
            }
            let c = m / a;
            if b > a {
                continue;
            }

            let tau = Complex::new(-b as f64 / (2.0 * a as f64), (big_d as f64).abs().sqrt() / (2.0 * a as f64));
            /* let numerator = dedekind_eta((Complex::new(0.0, 4.0 * PI ) * tau).exp(), precision);
            let denominator = dedekind_eta((Complex::new(0.0, 2.0 * PI ) * tau).exp(), precision); */
            let numerator = dedekind_eta(Complex::new(2.0, 0.0)* tau, precision);
            let denominator = dedekind_eta(tau, precision);
            let ratio = numerator / denominator;
            let f = ratio.powi(24);


            let j = (256.0 * f + Complex::new(1.0, 0.0)).powi(3) / f;

            if b == a || c == a|| b == 0 {
                t = polynomial_mul(&t, &vec![-j,Complex::new(1.0, 0.0)]);
                h += 1;
                red.insert((a,b,c));
            } else {
                t = polynomial_mul(&t, &vec![Complex::new(j.norm().powf(2.0),0.0), Complex::new(-2.0 * j.re,0.0), Complex::new(1.0, 0.0)]);
                h += 2;
                red.insert((a,b,c));
                red.insert((a,-b,c));
            }
        }
        b += 2;
    }
    /* if red != reduced_form {
        panic!("Reduced form inconsistent");
    } */
    let mut t_x = Vec::new();
    for p in t {
        t_x.push((p.re + 0.5).floor() as i64);
    }
    (h,t_x, red)
}

fn reduced_form(d: i64) -> HashSet<(i64, i64, i64)> {
    //Given discriminant D compute its reduced forms. Used to calculate preicion

    let mut red = HashSet::new();
    let mut b = ((d % 2)+2) % 2;
    let r = (d.abs() / 3).sqrt();

    while b <= r {
        let m = (b * b - d) / 4;
        let m_sqrt = (m as f64).sqrt() as i64;
        
        for a in 1..=m_sqrt {
            if m % a != 0 {
                continue;
            }
            
            let c = m / a;
            
            if b > a {
                continue;
            }
            
            if b == a || c == a || b == 0 {
                red.insert((a, b, c));
            } else {
                red.insert((a, b, c));
                red.insert((a, -b, c));
            }
        }
        
        b += 2;
    }
    red
}

fn dedekind_eta(tau: Complex<f64>, precision: i64) -> Complex<f64> {
    /* Implementation of dedekind's eta function.
       This implementation follows the idea in NZMATH's implementation */

    let mut tau = tau;
    let x = Complex::new(0.0, 2.0 * PI / 24.0).exp();
    let mut outer = Complex::new(1.0, 0.0);
    let mut absolute = 0.0;

    while absolute <= 1.0 - 0.1_f64.powi(5) {
        let real_tau = tau.re.floor() as i64;
        if real_tau != 0 {
            let real_tau_complex = Complex::new(real_tau as f64, 0.0);
            tau -= real_tau_complex;
            outer *= x.powf(real_tau as f64);
        }
        absolute = tau.norm();
        if absolute > 1.0 - 0.1_f64.powi(5) {
            break;
        }
        let mut ro = (-Complex::new(1.0, 0.0) / tau * Complex::new(0.0, 1.0)).sqrt();
        if ro.re < 0.0 {
            ro = -ro;
        }
        outer *= ro;
        tau = (-outer.re + outer.im * Complex::new(0.0, 1.0)) / absolute;
    }

    let q1 = (PI / 12.0 * tau * Complex::new(0.0, 1.0)).exp();
    let q = q1.powf(24.0);
    let mut sum = Complex::new(1.0, 0.0);
    let mut qs = Complex::new(1.0, 0.0);
    let mut qn = Complex::new(1.0, 0.0);
    let bound = f64::powf(10.0, -(precision + 2) as f64);

    while qs.norm() > bound {
        let t = -q * qn * qn * qs;
        qn *= q;
        qs = qn * t;
        sum += t + qs;
    }

    outer * q1 * sum
}

fn polynomial_mul(p1: &Vec<Complex<f64>>, p2: &Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    if p1.is_empty() || p2.is_empty() {
        panic!("Polynomial array empty.");
    }
    let mut m = vec![Complex::new(0.0, 0.0); p1.len() + p2.len() - 1];
    for i in 0..p1.len() {
        for j in 0..p2.len() {
            m[i+j] += p1[i] * p2[j];
        }
    }
    m
}

fn gen_quardratic_residue(p: &BigInt) -> BigInt {
    loop {
        let g = rand::thread_rng()
        .gen_bigint_range(&BigInt::from(2), &(p - 1)); 
        if jacobi_symbol(p, &g) != -1 {
            continue;
        }
        if p % 3 != BigInt::from(1) {
            return g
        }
        let cube = a_power_b_mod_c(&g, &((p-1)/3), p);
        if (&cube + p) % p == BigInt::from(1) {
            continue;
        }
        return g
    }    
}

fn cm_generating_curves(n: &BigInt, big_d: i64) -> Vec<(BigInt, BigInt)> {
    /* CM method for generating curves and orders: We assume a list of fundamental
       discriminants {Dj < 0 : j = 1, 2, 3, . . .} ordered, say, by increasing class
       number h(D), and within the same class number by increasing |D|. We are given
       a prime p > 3. The algorithm reports (optionally) possible curve orders or
       (also optionally) curve parameters for CM curves associated with the various Dj. 
       
       See Algorithm 7.5.9 from book: Richard E. Crandall, Carl Pomerance
       Prime Numbers: A Computational Perspective. 2. ed.*/

    let g = gen_quardratic_residue(n);

    if big_d == -3 {
        let a_b: Vec<(BigInt, BigInt)> = (0..6)
        .map(|i| (BigInt::from(0), (-a_power_b_mod_c(&g, &BigInt::from(i), n) + n) % n))
        .collect();
        return a_b
    } else if big_d == -4 {
        let a_b: Vec<(BigInt, BigInt)> = (0..4)
        .map(|i| ((-a_power_b_mod_c(&g, &BigInt::from(i), n) + n) % n, BigInt::from(0)))
        .collect();
        return a_b
    } else {
        let (_,t,_) = hilbert(big_d);
        let mut s = HashMap::new();
        for i in 0..t.len() {
            s.insert(BigInt::from(i), FieldElement::new((t[i]+n) % n, Some(n.clone())));
        }
        let j = root_of_poly(&s, n).pop().unwrap();
        let c = j.clone() / (j.clone() - FieldElement::new(BigInt::from(1728) % n, Some(n.clone())));
        let r: BigInt = ((-3 * &c.num) % n + n) % n;
        let s: BigInt = ((2 * &c.num) % n + n) % n;
        let a_b = vec![(r.clone(), s.clone()),((r * &g * &g) % n,(s * &g * &g * &g) % n)];
        return a_b
    }
}

fn find_a_b_q_for_m(n: &BigInt, m: &BigInt, big_d: i64) -> Option<(BigInt, BigInt, BigInt)> {
    let a_b = cm_generating_curves(n, big_d);
    for i in a_b {
        let ellip = EllipticCurve {
            a: i.0.clone(),
            b: i.1.clone(),
            field_p: n.clone(),
        };
        if test_a_b_order_m(&ellip, &m) == false {
            continue
        }
        let factor = kq_factor(&m, n);
        if let Some((_, q)) = factor {
            return Some((i.0, i.1, q));
        }
    }
    None
}

fn test_a_b_order_m(ellip: &EllipticCurve, m: &BigInt) -> bool {
    for _ in 0..2 {
        let random_point = ellip.random_point();
        if let PointForm::Infinity = ellip.scalar_mul(&random_point, m) {
            continue
        } else {
            return false
        }
    }
    true
}

fn decide_q(a: &BigInt, b: &BigInt, n: &BigInt, m: &BigInt, q: &BigInt) -> Option<BigInt> {
    let ellip = EllipticCurve{
        a: a.clone(),
        b: b.clone(),
        field_p: n.clone(),
    };
    let point = ellip.random_point();
    let u = ellip.scalar_mul(&point, &(m / q));
    if let PointForm::Infinity = u {
        return None
    }
    let v = ellip.scalar_mul(&u, q);
    if let PointForm::Infinity = v {
        return Some(q.clone())
    } else {
        return None
    }
}

//compute roots of polynomial over finite field
fn root_of_poly(g: &HashMap<BigInt, FieldElement>, p: &BigInt) -> Vec<FieldElement> {
    let mut r = Vec::new();
    let x = FieldElement {
        num: BigInt::one(),
        prime: p.clone()
    };
    let y = FieldElement {
        num: p - BigInt::one(),
        prime: p.clone()
    };
    
    let f_x = HashMap::from([(BigInt::one(), x)]);
    let mut g_new = poly_x_power_b_mod_poly_c(&f_x, &(p-2), &g);
    g_new = mul_f_g(&f_x, &g_new);
    match g_new.get(&BigInt::zero()) {
        Some(value) => g_new.insert(BigInt::zero(), value.clone()+y),
        None => g_new.insert(BigInt::zero(), FieldElement::new(p-1, Some(p.clone()))),
    };

    let mut g = polynomials_gcd(&mul_f_g(&f_x, &g_new), &g);

    while let None = g.get(&BigInt::zero()) {
        r.push(FieldElement::new(BigInt::zero(), Some(p.clone())));
        let mut new_g = HashMap::new();
        for i in g {
            new_g.insert(i.0 - 1, i.1);
        }
        g = new_g;
    }
    let a = roots_g(&g, p);
    r.extend(a);
    r
}

//saparate roots to subproblems
fn roots_g(g: &HashMap<BigInt, FieldElement>, p: &BigInt) -> Vec<FieldElement> {
    if max_key(&g) == BigInt::from(0) {
        return Vec::new();
    }
    if max_key(&g) == BigInt::from(1) {
        return root_deg_one(&g);
    }
    if max_key(&g) == BigInt::from(2) {
        return roots_deg_two(&g, &p);
    }
    let mut res: Vec<FieldElement> = Vec::new();
    let mut h = g.clone();
    while max_key(&h) == BigInt::one() || &h == g {
        let a = rand::thread_rng()
        .gen_bigint_range(&BigInt::from(2), p);
        let x_plus_a = 
        HashMap::from([(BigInt::one(), FieldElement {num: BigInt::one(),prime: p.clone()}),
        (BigInt::zero(), FieldElement {num: a,prime: p.clone()})]);
        let mut x_plus_a = poly_x_power_b_mod_poly_c(&x_plus_a, &((p-1)/2), &g);
        match x_plus_a.get(&BigInt::zero()) {
            Some(value) => x_plus_a.insert(BigInt::zero(), value.clone()+FieldElement::new(p-1, Some(p.clone()))),
            None => x_plus_a.insert(BigInt::zero(), FieldElement::new(p-1, Some(p.clone()))),
        };
        if x_plus_a == HashMap::from([(BigInt::zero(), FieldElement {num: BigInt::zero(),prime: p.clone()})]) {
            continue
        }
        h = polynomials_gcd(&x_plus_a, &g);
    }
    res.extend(roots_g(&h, p));
    res.extend(roots_g(&div_f_g(&g, &h), p));
    res
}

//root for degree 2
fn roots_deg_two(g: &HashMap<BigInt, FieldElement>, p: &BigInt) -> Vec<FieldElement> {
    let mut res = Vec::new();
    if let Some(two_value) = g.get(&BigInt::from(2)) {
        if let Some(one_value) = g.get(&BigInt::one()) {
            if let Some(zero_value) = g.get(&BigInt::zero()) {
                let delta = one_value.clone() * one_value.clone() - FieldElement::new(BigInt::from(4), Some(p.clone())) * zero_value.clone() * two_value.clone();
                if let Ok(ro) = square_root_mod_prime(p, &delta.num) {
                    res.push((FieldElement::new(&ro % p, Some(p.clone()))- one_value.clone())/two_value.clone());
                    res.push((FieldElement::new((p-ro) % p, Some(p.clone()))- one_value.clone())/two_value.clone());
                }
            } else {
                res.push(two_value.clone()/one_value.clone());
            }
        } else {
            if let Some(zero_value) = g.get(&BigInt::zero()) {
                let l = FieldElement::new((p-zero_value.num.clone()) % p, Some(p.clone())) / two_value.clone();
                if let Ok(ro) = square_root_mod_prime(p, &l.num) {
                    res.push(FieldElement::new(&ro % p, Some(p.clone())));
                    res.push(FieldElement::new((p-ro) % p, Some(p.clone())));
                }
            }
        }
    }
    res
}

//root for degree 1
fn root_deg_one(g: &HashMap<BigInt, FieldElement>) -> Vec<FieldElement> {
    let mut res = Vec::new();
    if let Some(one_value) = g.get(&BigInt::one()) {
        if let Some(zero_value) = g.get(&BigInt::zero()) {
            res.push(one_value.clone() / zero_value.clone());
        } 
    }
    res
}

//searches the gcd for two polynomials
fn polynomials_gcd(f: &HashMap<BigInt, FieldElement>, g: &HashMap<BigInt, FieldElement>)
    -> HashMap<BigInt, FieldElement> {
   

    if f.is_empty() && g.is_empty() {
        panic!("Cannot compute gcd for zero polynomials");
    }

    let (mut u,mut v) = match max_key(f) >= max_key(g) {
        true  => (f.clone(), g.clone()),
        false => (g.clone(), f.clone())
    };

    loop {
        if u.is_empty() && !v.is_empty() {
            return v
        }
        if !u.is_empty() && v.is_empty() {
            return u
        }

        let (u_deg, v_deg) = (max_key(&u), max_key(&v));
        let a = u.get(&u_deg).unwrap().clone() / v.get(&v_deg).unwrap().clone();
        for i in &v {
            match u.get(&(&u_deg - &v_deg + i.0)) {
                Some(value) => {
                    let u_value = value.clone() - a.clone() * i.1.clone();
                    u.insert(&u_deg - &v_deg + i.0, u_value);
                }
                None =>  {
                    let zero = FieldElement {
                        num: BigInt::zero(),
                        prime: i.1.prime.clone()
                    };
                    let u_value = zero - a.clone() * i.1.clone();
                    u.insert(&u_deg - &v_deg + i.0, u_value);
                }
            }
        }

        for i in u.clone() {
            if i.1.num == BigInt::zero() {
                u.remove(&i.0);
            }
        }

        if max_key(&u) < max_key(&v) || u.is_empty() {
            let u_new = v;
            v = u;
            u = u_new;
        }
    }
}

//get highest coefficient
fn max_key (f: &HashMap<BigInt, FieldElement>) -> BigInt {
    if f == &HashMap::new() {
        return BigInt::from(-1);
    }
    if let Some(max_key) = f.iter().max_by_key(|&(key, _value)| key) {
        return max_key.0.clone();
    } else {
        return BigInt::zero();
    }
}

//compute f mod g
fn f_mod_g(f: &HashMap<BigInt, FieldElement>, g: &HashMap<BigInt, FieldElement>) -> HashMap<BigInt, FieldElement> {
    let (mut u, v) = (f.clone(), g.clone());
    while max_key(&u) >= max_key(&v) {
        let (u_deg, v_deg) = (max_key(&u), max_key(&v));
        let a = u.get(&u_deg).unwrap().clone() / v.get(&v_deg).unwrap().clone();
        for i in &v {
            match u.get(&(&u_deg - &v_deg + i.0)) {
                Some(value) => {
                    let u_value = value.clone() - a.clone() * i.1.clone();
                    u.insert(&u_deg - &v_deg + i.0, u_value);
                }
                None =>  {
                    let zero = FieldElement {
                        num: BigInt::zero(),
                        prime: i.1.prime.clone()
                    };
                    let u_value = zero - a.clone() * i.1.clone();
                    u.insert(&u_deg - &v_deg + i.0, u_value);
                }
            }
        }

        for i in u.clone() {
            if i.1.num == BigInt::zero() {
                u.remove(&i.0);
            }
        }
    }
    u
}

//compute f^b mod g
fn poly_x_power_b_mod_poly_c(f: &HashMap<BigInt, FieldElement>, b: &BigInt, g: &HashMap<BigInt, FieldElement>)
-> HashMap<BigInt, FieldElement> {
    let mut y = f.clone();

    let scalar_binary = format!("{:b}", b);

    for (i, digit) in scalar_binary.chars().enumerate() {
        if i == 0 {
            y = f_mod_g(&y, &g);
        } else {
            match digit {
                '0' => {                        
                    y = mul_f_g(&y, &y);
                    y = f_mod_g(&y, &g);
                }
                '1' => {                        
                    y = mul_f_g(&y, &y);
                    y = mul_f_g(&y, f);
                    y = f_mod_g(&y, &g);   
                }
                _ => {
                    panic!("Failed multiply")
                }
            }
        }       
    }

    y
} 

//f * g
fn mul_f_g(f: &HashMap<BigInt, FieldElement>, g: &HashMap<BigInt, FieldElement>) -> HashMap<BigInt, FieldElement> {
    if f == &HashMap::new() || g == &HashMap::new() {
        return HashMap::new();
    }
    let mut f_g: HashMap<BigInt, FieldElement> = HashMap::new();
    for i in f {
        for j in g {
            match f_g.get(&(i.0 + j.0)) {
                Some(value) => {
                    let new_value = value.clone() + i.1.clone() * j.1.clone();
                    f_g.insert(i.0.clone() + j.0.clone(), new_value);
                }
                None => {
                    f_g.insert(i.0.clone() + j.0.clone(), i.1.clone() * j.1.clone());
                }
            }
        }
    }
    f_g
}

//f / g
fn div_f_g(f: &HashMap<BigInt, FieldElement>, g: &HashMap<BigInt, FieldElement>) -> HashMap<BigInt, FieldElement> {
    let (mut u, v) = (f.clone(), g.clone());
    let mut res = HashMap::new();
    while max_key(&u) >= max_key(&v) {
        let (u_deg, v_deg) = (max_key(&u), max_key(&v));
        let a = u.get(&u_deg).unwrap().clone() / v.get(&v_deg).unwrap().clone();
        for i in &v {
            match u.get(&(&u_deg - &v_deg + i.0)) {
                Some(value) => {
                    let u_value = value.clone() - a.clone() * i.1.clone();
                    u.insert(&u_deg - &v_deg + i.0, u_value);
                }
                None =>  {
                    let zero = FieldElement {
                        num: BigInt::zero(),
                        prime: i.1.prime.clone()
                    };
                    let u_value = zero - a.clone() * i.1.clone();
                    u.insert(&u_deg - &v_deg + i.0, u_value);
                }
            }
            
        }

        res.insert(BigInt::from(u_deg- v_deg), a);

        for i in u.clone() {
            if i.1.num == BigInt::zero() {
                u.remove(&i.0);
            }
        }
    }
    res
}



#[derive(Clone)]
#[derive(Debug)]
enum PointForm {
    Infinity,
    Point(BigInt, BigInt),
}

#[warn(dead_code)]
struct EllipticCurve {
    a: BigInt,
    b: BigInt,
    field_p: BigInt,
}

impl EllipticCurve {
    fn random_point(&self) -> PointForm {
        loop {
            let x = rand::thread_rng()
            .gen_bigint_range(&BigInt::from(0), &(&self.field_p));
            match self.y(&x) {
                Some(y) => return PointForm::Point(x, y),
                None => continue
            }
        }
    }

    fn y(&self, x: &BigInt) -> Option<BigInt>{
        /* Given x compute y such that (x,y) on E(Fp) */
        let big_q = (x * x * x + &self.a * x + &self.b) % &self.field_p;
        if jacobi_symbol(&self.field_p, &big_q) != -1 {
            let y = match square_root_mod_prime(&self.field_p, &big_q) {
                Ok(some) => some,
                Err(_) => return None
            };
            if &y * &y % &self.field_p != big_q {
                return None
            }
            return Some(y)
        }
        None
    }

    fn add(&self, p: &PointForm, q: &PointForm) -> PointForm {
        let mut m = BigInt::from(0);
        match (p, q) {
            (PointForm::Infinity, _) => q.clone(),
            (_, PointForm::Infinity) => p.clone(),
            (PointForm::Point(x1, y1), PointForm::Point(x2, y2)) => {
                if x1 == x2 && y1 != y2 {
                    return PointForm::Infinity
                } else if x1 == x2 {
                    let inv = a_power_b_mod_c(&(2 * y1), &(&self.field_p-2), &self.field_p);
                    m = m + (((3 * x1 *x1 + &self.a) % &self.field_p) * inv) % &self.field_p;
                } else {
                    let inv = a_power_b_mod_c(&(x2 - x1), &(&self.field_p-2), &self.field_p);
                    m = m + (((y2 - y1) % &self.field_p) * inv) % &self.field_p;
                }
                let x3 = (((&m * &m - x1 - x2) % &self.field_p) + &self.field_p) % &self.field_p;
                PointForm::Point(x3.clone(), ((m * (x1 - x3) - y1) % &self.field_p + &self.field_p) % &self.field_p)
            }
        }
    }

    /* fn sub(&self, p: &PointForm, q: &PointForm) -> PointForm {
        match (p, q) {
            (PointForm::Infinity, PointForm::Infinity) => PointForm::Infinity,
            (PointForm::Infinity, PointForm::Point(x2, y2)) => PointForm::Point(x2.clone(), -y2),
            (_, PointForm::Infinity) => p.clone(),
            (_, PointForm::Point(x2, y2)) => {
                self.add(p, &PointForm::Point(x2.clone(), -y2))
            }
        }
    } */

    fn scalar_mul(&self, p: &PointForm, scalar: &BigInt) -> PointForm {
        let mut res = PointForm::Infinity;
        if scalar == &BigInt::from(0) {
            return PointForm::Infinity
        } else if scalar == &BigInt::from(1) {
            return p.clone()
        } else {
            let scalar_binary = format!("{:b}", scalar);

            for (_, digit) in scalar_binary.chars().enumerate() {
                match digit {
                    '0' => {                        
                        res = self.add(&res, &res);
                    }
                    '1' => {                        
                        res = self.add(&res, &res);
                        res = self.add(&res, p);    
                    }
                    _ => {
                        panic!("Failed multiply")
                    }
                }
                
            }
        }
        res
    }
}

fn baillie(r: &BigInt) -> bool {
    /* Baillie primality test */
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

fn strong_lucas_tests(num: &BigInt) -> bool {
/*  This function perform a strong Lucas probable prime test with following steps:
    1. decomposite n + 1 = 2 ^ s * d
    2. choose D, set P, Q
    3. calculate the d-th and (d+1)-th element V_d, V_{d+1}
    4. calculate U_d, decide whether U_d ≡ 0 (mod n)
        if true, n is probably prime
    5. calculate V_{d*2^r}, decide whether V_{d*2^r} ≡ 0 (mod n), for 0 <= r < s
        if true, n is probably prime */

    let iterations = 5;  //number of tests with different P, Q, D
    
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
                b_j = &b_j * &b_j;
            }
            '1' => {                        //match (2)
                u = &u * &v - a * &b_j;
                v = &v * &v - &b_j * b * 2;
                b_j = &b_j * &b_j;
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

fn miller_test_base_two(d: BigInt, n: &BigInt) -> bool {
    let mut x = a_power_b_mod_c(&BigInt::from(2), &d, n);

    if x == BigInt::from(1) || x == (n - 1) {
        return true;
    }

    let mut d = d;
    while d < (n - 1) {
        x = (&x * &x) % n;
        d <<= 1;

        if x == BigInt::from(1) {
            return false;
        }
        if x == (n - 1) {
            return true;
        }
    }

    false
}

fn is_prime_miller(n: &BigInt) -> bool {
    if n == &BigInt::from(1) || n == &4.to_bigint().unwrap() {
        return false;
    }
    let mut d:BigInt = n - 1;
    while &d & &BigInt::from(1) == BigInt::from(0) {
        d >>= 1;
    }

    if !miller_test_base_two(d.clone(), n) {
        return false;
    }
    true
}


mod poly {
    use std::ops::{Add, Div, Mul, Sub};

    use num_bigint::BigInt;
    use num_traits::{Num, One, Zero};

    pub const P: &str = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";

    #[derive(Debug, Clone)]
    pub struct FieldElement {
        pub num: BigInt,
        pub prime: BigInt,
    }

    impl FieldElement {
        pub fn new(num: BigInt, prime: Option<BigInt>) -> Self {
            let prime = if prime.is_none() {
                BigInt::from_str_radix(P, 16).unwrap()
            } else {
                prime.unwrap()
            };

            if num >= prime {
                panic!(
                    "Num {} not in field range 0 to {}",
                    num,
                    prime - BigInt::one()
                );
            }
            Self { num, prime }
        }

        // pub fn zero(prime: BigInt) -> Self {
        //     Self {
        //         num: BigInt::zero(),
        //         prime,
        //     }
        // }

        // pub fn get_prime(&self) -> BigInt {
        //     self.prime.clone()
        // }

        // pub fn get_number(&self) -> BigInt {
        //     self.num.clone()
        // }

        pub fn to_the_power_of(&self, exponent: BigInt) -> Self {
            let exp = exponent % (self.prime.clone() - BigInt::one());
            let new_num = Self::mod_pow(&self.num, exp, &self.prime);
            Self {
                num: new_num,
                prime: self.prime.clone(),
            }
        }

        pub fn mod_pow(base: &BigInt, mut exp: BigInt, modulus: &BigInt) -> BigInt {
            if modulus == &BigInt::one() {
                return BigInt::zero();
            }

            let mut result = BigInt::one();
            let mut base = base % modulus;
            while exp > BigInt::zero() {
                if exp.clone() % (BigInt::one() + BigInt::one()) == BigInt::one() {
                    result = result * base.clone() % modulus;
                }
                exp = exp / (BigInt::one() + BigInt::one());
                base = base.clone() * base.clone() % modulus;
            }

            result
        }

        // pub fn ne(&self, other: &FieldElement) -> bool {
        //     self.num != other.num || self.prime != other.prime
        // }

        // pub fn pow(&self, exp: u32) -> Self {
        //     let num = self.modulo(&self.num.pow(exp));
        //     Self {
        //         num,
        //         prime: self.prime.clone(),
        //     }
        // }

        fn modulo(&self, b: &BigInt) -> BigInt {
            let result = b % self.prime.clone();
            if result < BigInt::zero() {
                result + self.prime.clone()
            } else {
                result
            }
        }

        // pub fn sqrt(&self) -> Self {
        //     let p = BigInt::from_str_radix(P, 16).unwrap();
        //     self.to_the_power_of((p + BigInt::one()) / (BigInt::from_u8(4).unwrap()))
        // }
    }

    impl PartialEq for FieldElement {
        fn eq(&self, other: &FieldElement) -> bool {
            self.num == other.num && self.prime == other.prime
        }
    }

    impl Eq for FieldElement {}

    impl Add for FieldElement {
        type Output = Self;

        fn add(self, rhs: Self) -> Self::Output {
            if self.prime != rhs.prime {
                panic!("cannot add two numbers in different Fields");
            }

            let num = self.modulo(&(self.num.clone() + rhs.num));
            Self {
                num,
                prime: self.prime.clone(),
            }
        }
    }

    impl Sub for FieldElement {
        type Output = Self;

        fn sub(self, rhs: Self) -> Self::Output {
            if self.prime != rhs.prime {
                panic!("cannot subtract two numbers in different Fields");
            }

            let difference = BigInt::from(self.num.clone()) - BigInt::from(rhs.num.clone());
            let big_prime = BigInt::from(self.prime.clone());
            let remainder = difference % big_prime.clone();
            if remainder < BigInt::zero() {
                let new_num = remainder + big_prime;
                Self {
                    num: new_num.try_into().unwrap(),
                    prime: self.prime.clone(),
                }
            } else {
                Self {
                    num: remainder.try_into().unwrap(),
                    prime: self.prime.clone(),
                }
            }
        }
    }

    impl Mul for FieldElement {
        type Output = Self;

        fn mul(self, rhs: Self) -> Self::Output {
            if self.prime != rhs.prime {
                panic!("cannot multiply two numbers in different Fields");
            }

            let num = self.modulo(&(self.num.clone() * rhs.num));
            Self {
                num,
                prime: self.prime.clone(),
            }
        }
    }

    impl Div for FieldElement {
        type Output = Self;

        fn div(self, rhs: Self) -> Self::Output {
            if self.prime != rhs.prime {
                panic!("cannot divide two numbers in different Fields");
            }

            // use Fermat's little theorem
            // self.num.pow(p-1) % p == 1
            // this means:
            // 1/n == pow(n, p-2, p) in Python
            let exp = rhs.prime.clone() - (BigInt::one() + BigInt::one());
            let num_pow = rhs.to_the_power_of(exp);
            let result = self.num.clone() * num_pow.num;
            Self {
                num: result % self.prime.clone(),
                prime: self.prime.clone(),
            }
        }
    }
}