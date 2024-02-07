# SC_Final
Exploring Prime Number Testing Algorithms: A Comparative Analysis

## Authors
- Wei-Shan Chang
- Ximing Zhang
- Yunzhou Lu

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
        - [How to Run the Experiments](#how-to-run-the-experiments)
  

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
        ├── ...
    ├── experiment
        ├── experiment_det_speed.sh
        ├── experiment_prob_speed.sh
        ├── experiment_type.sh
    ├── Cargo.toml
    ├── README.md
    
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

#### Experiment 1: Efficiency of Probabilistic Primality Testing
This experiment evaluates the efficiency of probabilistic primality testing algorithms, including Fermat's Little Theorem, Lucas Primality Test, Baillie–PSW Primality Test, and Miller-Rabin Primality Test. Tests are conducted with a 600-second time limit and an upper bound of 10^6, providing insights into the reliability and computational efficiency of these algorithms. Pseudo primes identified in the output are examined for potential limitations and areas for improvement.

Execute the following command to measure the speed of the probabilistic test under the <SC_Final> directory:
```sh
bash experiment/experiment_prob_speed.sh
```

#### Experiment 2: Efficiency of Deterministic Primality Testing
The study explores the efficiency of deterministic primality testing algorithms, encompassing the Adleman–Pomerance–Rumely Primality Test, Brute-Force method, Sieve of Eratosthenes, and Elliptic Curve Primality Testing. With a 600-second time limit and a list of random number chosed from 1 to 10^14, the tests aim to assess reliability and computational efficiency under specific constraints.


Run the following command to evaluate the speed of the deterministic test under the <SC_Final> directory:
```sh
bash experiment/experiment_det_speed.sh
```


#### Experiment 3: Accuracy Assessment and Enhancing Accuracy of Probabilistic Testing
To address potential inaccuracies in probabilistic primality tests, this experiment combines probabilistic tests with the fastest deterministic test. The goal is to assess accuracy and understand the effectiveness of deterministic tests in validating results, aiming to enhance the reliability of primality testing.

This test is a combination of the previous two tests.

#### Experiment 4: Comparative Efficiency Analysis of Different Types - BigInt vs. u64
Initially utilizing BigInt for handling large numbers, this experiment compares the efficiency of tests using BigInt and u64. Employing brute force methods, the study highlights performance differences between these data types in prime testing scenarios.

Execute the following command to measure the type test under the <SC_Final> directory:
```sh
bash experiment/experiment_type.sh
```
#### Program Running Environment
The Rust program, developed with Rust version 1.73.0, relies on key dependencies such as rand 0.8 and num-bigint 0.4 with the "rand" feature. Tested on a MacBook Pro running macOS with Darwin Kernel Version 22.3.0, built on Jan 30, 2023, the system supports the ARM64 architecture. Managed with rustup, the program's Cargo.toml file ensures dependency reproducibility and considers specific configurations for the target operating system.

These commands will execute the respective experiments and provide insights into the speed and correctness of the implemented primality tests. Adjustments can be made to the scripts or experiments as needed based on specific requirements.