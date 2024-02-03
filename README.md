# SC_Final
Exploring Prime Number Testing Algorithms: A Comparative Analysis

## Introduction

This study explores and compares prime number testing algorithms, focusing on their implementation in the Rust programming language. By assessing computational complexity, accuracy, and practical implications, we aim to provide insights into the strengths and weaknesses of these algorithms in the context of Rust. This research contributes to the understanding of algorithmic efficiency in Rust, assisting developers in making informed choices for diverse applications.

## Table of Contents

- [SC_Final](#SC_Final)
    - [Introduction](#introduction)
    - [Table of Contents](#table-of-contents)
    - [Built with](#built-with)
- [Installation Guide](#installation-guide)
    - [Workspace](#workspace)
- [Usage Instructions](#usage-instructions)
    - [Primality Test Tool](#primality-test-tool)
        - [Overview](#overview)
        - [How to Run Each Test](#how-to-run-each-test)
        - [Example](#example)
        - [Output Interpretation](#output-interpretation)
        - [How to Loop Over the Test](#how-to-loop-over-the-test)
  

## Built with
* [Rust](https://www.rust-lang.org) - A systems programming language that combines performance and safety, empowering developers to write reliable and efficient software.

# Installation Guide

## Workspace
1. To clone this repository, execute the following command in your terminal:
    ```sh
    git clone https://github.com/weishan880907/SC_Final.git
    ```

2. Initialize the repository with the following commands:
    ```sh
    git init
    git remote add SC_Final https://github.com/weishan880907/SC_Final.git
    ```

3. After cloning the repository, your project workspace should have the following folder structure:
    ```
    SC_Final 
    ├── src
        ├── aks.rs
        ├── brute_force.rs
        ├── fermat_little_theorem.rs
        ├── sieve_of_eratosthenes.rs
    ├── Cargo.toml
    ├── README.md
    ├── experiment.sh
    ```

4. Switch to your personal development branch using the following command:
    ```sh
    git checkout -b <your_name>-dev
    git push -u <remote_name> <your_name>-dev
    ```

# Usage Instructions

## Primality Test Tool

### Overview
This tool provides functionality for testing the primality of given BigInt numbers using different primality tests. The updated version now supports testing a vector of BigInt numbers and provides results in the form of prime numbers along with the execution time.

### How to Run Each Test
#### Probabilistic Test

To execute a probabilistic primality test for a vector of BigInt numbers, use the following command format:

```sh
cargo run --bin <probabilistic_primality_test_name> <num_upperbound>
```

#### Example
For instance, to apply the Fermat's Little Theorem primality test to check numbers up to 24, you can use the following command:
```sh
cargo run --bin fermat_little_theorem 24
```

Replace <primality_test_name> with the specific probabilistic primality test you want to run, and provide an upper bound <num_upperbound> to test for primality.

#### Deterministic Test

To execute a deterministic primality test for a vector of BigInt numbers, use the following command format:
```sh
cargo run --bin <primality_test_name> <probabilistic_primality_test_name>.txt

```

### Output Interpretation
#### Probabilistic Test
The output includes a file named <Prob_test_name>.txt, containing probable prime numbers, and the corresponding execution time is displayed in the terminal.

#### Deterministic Test

The output indicates the count of pseudo primes obtained from the previous test and the corresponding execution time.
### How to Run the Experiments
We have designed three experiments to assess different aspects of our primality tests. To execute these experiments, follow the steps below:

#### Experiment 1: Deterministic Test Speed

Run the following command to evaluate the speed of the deterministic test:
```sh
bash experiment_det_speed.sh
```

#### Experiment 2: Probabilistic Test Speed

Execute the following command to measure the speed of the probabilistic test:
```sh
bash experiment_prob_speed.sh
```

#### Experiment 3: Correctness Test

For this experiment, we first conduct the probabilistic test to identify probable primes. Then, we use these probable primes as input for the deterministic test to check for pseudo primes.

To perform the correctness test, use the following command:
```sh
bash experiment_correctness.sh
```
These commands will execute the respective experiments and provide insights into the speed and correctness of the implemented primality tests. Adjustments can be made to the scripts or experiments as needed based on specific requirements.