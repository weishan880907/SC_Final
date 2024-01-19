// We will find all primes in the range 1 to 120
const N: usize = 121;
static mut IS_PRIME: [i32; N] = [1; N];

fn sieve() {
    unsafe {
        // We cross out all composites from 2 to sqrt(N)
        let mut i = 2;
        // This will loop from 2 to int(sqrt(x))
        while i * i <= N {
            // If we already crossed out this number, then continue
            if IS_PRIME[i] == 0 {
                i += 1;
                continue;
            }
            let mut j = 2 * i;
            while j < N {
                // Cross out this as it is composite
                IS_PRIME[j] = 0;
                // j is incremented by i, because we want to cover all multiples of i
                j += i;
            }
            i += 1;
        }
    }
}

fn main() {
    unsafe {
        // is_prime[i] = 1 means that i is prime and is_prime[i] = 0 means that i is composite
        // Initially, we say all of them are prime
        for i in IS_PRIME.iter_mut() {
            *i = 1;
        }
        // We know 0 and 1 are composites
        IS_PRIME[0] = 0;
        IS_PRIME[1] = 0;
        sieve();
        // Print all the primes in between 1 and 121
        for i in 1..N {
            if IS_PRIME[i] == 1 {
                print!("{} ", i);
            }
        }
        // Output: 2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113
    }
}
