import subprocess
import pandas as pd

# 1. Configuration for ABF1
REF = "S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa"
VCF = "abf1_region.vcf.gz"
FASTA_CHR = "ref|NC_001143|" 
START, END = 458951, 459030
REGION = f"{FASTA_CHR}:{START}-{END}"
ENVIRONMENTS = ["YPD", "SD"]

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def main():
    strains = subprocess.check_output(f"bcftools query -l {VCF}", shell=True).decode().splitlines()
    print(f"Extracting ABF1 sequences for {len(strains)} strains...")

    all_data = []
    for i, strain in enumerate(strains):
        cmd = f"samtools faidx {REF} \"{REGION}\" | bcftools consensus -s {strain} {VCF}"
        try:
            fasta_output = subprocess.check_output(cmd, shell=True).decode()
            raw_seq = "".join(fasta_output.splitlines()[1:])
            # ABF1 is on the negative strand -> Reverse Complement
            final_seq = reverse_complement(raw_seq)[:80]

            for env in ENVIRONMENTS:
                all_data.append({
                    "Strain": strain, "Gene": "ABF1", 
                    "Environment": env, "Sequence_80bp": final_seq
                })
        except:
            continue

    pd.DataFrame(all_data).to_csv("WT_ABF1_Sequences.csv", index=False)
    print("Finished! Created 'WT_ABF1_Sequences.csv'.")

if __name__ == "__main__":
    main()
