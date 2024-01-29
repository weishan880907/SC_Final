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
- [Usage Intructions](#usage-instructions)
    - [Performing a Primality Test](#performing-a-primality-test)
  

## Built with
* [Rust](https://www.rust-lang.org) - A systems programming language that combines performance and safety, empowering developers to write reliable and efficient software.

# Installation Guide

## Workspace
1. To clone this repository, execute the following command in your terminal:
    ```
    git clone https://github.com/weishan880907/SC_Final.git
    ```


2. Initialize the repository with the following commands:
    ```
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
    ```
    
4. Switch to your personal development branch using the following command:
    ```
    git checkout -b <your_name>-dev
    git push -u <remote_name> <your_name>-dev
    ```

# Usage instruction

# Primality Test Tool

## Usage Instructions

### Overview
This tool provides functionality for testing the primality of given BigInt numbers using different primality tests. The updated version now supports testing a vector of BigInt numbers and provides results in the form of prime numbers along with the execution time.

### How to Run

To execute a primality test for a vector of BigInt numbers, use the following command format:

```sh
cargo run --bin <primality_test_name> <number1> <number2> ...
```

Replace <primality_test_name> with the specific primality test you want to run, and provide a list of BigInt numbers (<number1>, <number2>, ...) to test for primality.

### Example
For example, to test whether 10, 23, and 37 are prime or probably prime using the brute-force primality test, run the following command:

```sh
cargo run --bin brute_force 10 23 37
```
### Output Interpretation
The output of the test will indicate the prime numbers found in the provided list along with the execution time.
```
Prime numbers: [23, 37]
Time: 36.583µs
```

Adjust the command and interpret the output accordingly based on your specific primality test requirements.

Feel free to explore different primality tests by replacing <primality_test_name> with the desired test name.