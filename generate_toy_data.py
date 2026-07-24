import random
import os
import gzip

os.makedirs("reference", exist_ok=True)
os.makedirs("data", exist_ok=True)

# 1. Generate a synthetic 10kb reference genome
random.seed(42)
ref_seq = "".join(random.choices("ACGT", k=10000))

with open("reference/genome.fa", "w") as f:
    f.write(">chr1\n")
    # Wrap text to 80 characters per line
    f.write("\n".join(ref_seq[i:i+80] for i in range(0, len(ref_seq), 80)) + "\n")

# 2. Generate simulated paired-end reads (100bp) with artificial coverage
print("Generating reads...")
with gzip.open("data/ToySample_R1.fastq.gz", "wt") as f1, gzip.open("data/ToySample_R2.fastq.gz", "wt") as f2:
    for i in range(20000):
        # Pick a random start site, leaving room for a 300bp insert size
        start = random.randint(0, 9500)
        
        # Extract read 1 and read 2 sequences
        read1 = ref_seq[start : start+100]
        # Reverse complement logic is skipped here for pure toy speed; 
        # BWA will still map it, just with different flags.
        read2 = ref_seq[start+200 : start+300] 
        
        # Introduce a fake SNP in a few reads to trigger the variant callers
        if i % 2 == 0:
            read1 = read1[:50] + random.choice([c for c in "ACGT" if c != read1[50]]) + read1[51:]
        
        # Write FASTQ format with dummy high quality scores ('I')
        f1.write(f"@read_{i}/1\n{read1}\n+\n{'I'*100}\n")
        f2.write(f"@read_{i}/2\n{read2}\n+\n{'I'*100}\n")

print("Toy data generated successfully.")
