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
    - [How to run a primality test]


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
    ```
    
4. Switch to your personal development branch using the following command:
    
    ```
    git checkout -b <your_name>-dev
    git push -u <remote_name> <your_name-dev
    ```


# Usage instruction

## Performing a Primality Test
To execute a primality test, use the following command format: 

    cargo run --bin <primality_test_name> <number_to_test>
    
Replace <primality_test_name> with the specific primality test you want to run, and <number_to_test> with the numerical value you wish to test for primality.

If no argument is provided, the default input is set to 23.

The output of the test will indicate whether the specified number is determined to be a prime, composite, or probable prime.

For example, the output should resemble:

1. <number_to_test> is prime. (if the number is prime)

2. <number_to_test> is composite. (if the number is composite)

3. <number_to_test> is probably prime. (if the number is a probably prime)

Adjust the command and interpret the output accordingly based on your specific primality test requirements.
