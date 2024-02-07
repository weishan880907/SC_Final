import random

# Set the seed for reproducibility
random.seed(42)

# Define the range
start_range = 1
end_range = 100000000000000

# Number of random numbers to generate
num_random_numbers = 10

# Generate random numbers
random_numbers = [random.randrange(start_range, end_range) for _ in range(num_random_numbers)]

# Write random numbers to a file
output_file_path = f"random_numbers_{end_range}.txt"
with open(output_file_path, "w") as output_file:
    for number in random_numbers:
        output_file.write(str(number) + "\n")

print(f"Random numbers written to {output_file_path}")
