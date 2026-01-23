import random

def generate_random_genome(filename, size_bp=100000):
    print(f"Generating random genome of size {size_bp} bp to {filename}...")
    bases = ['A', 'C', 'G', 'T']
    
    with open(filename, 'w') as f:
        f.write(">chr1 random_sequence\n")
        # Write in lines of 80 chars
        chunk_size = 80
        for i in range(0, size_bp, chunk_size):
            line = "".join(random.choices(bases, k=min(chunk_size, size_bp - i)))
            f.write(line + "\n")
            
    print("Done.")

if __name__ == "__main__":
    generate_random_genome("test_genome.fasta")
